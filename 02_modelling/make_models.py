#!/usr/bin/env python

import argparse
import collections
import configparser
import copy
import itertools
import os
import pathlib
import shutil
import subprocess
import time

from Bio import SeqIO
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import PDBParser, NeighborSearch, PDBIO


protein_letters_3to1 = {
    k.upper(): v for k, v in protein_letters_3to1.items()
}

def read_seqfile(fpath):
    seqs = {}
    with open(fpath) as handle:
        for r in SeqIO.parse(handle, 'fasta'):
            isolate_id = r.description.split("|")[0]
            seqs[isolate_id] = str(r.seq)

    print(f'Read {len(seqs)} sequences from {fpath}')
    return seqs


def read_pdb(fpath):
    fpath = pathlib.Path(fpath)

    parser = PDBParser(QUIET=1)
    return parser.get_structure(fpath.name, str(fpath))


def map_seq_to_pdbchain(seq, pdb, chain):
    """Crude but works for us here."""
    aaset = set(protein_letters_3to1)

    pdbseq = [
        (protein_letters_3to1.get(aa.resname, "X"), aa.id[1])
        for aa in pdb[0][chain] if aa.resname in aaset
    ]

    # We assume seq is at least the same size as pdbseq. Might be longer
    # because of gaps.
    mapping = {}

    for idx, aa in enumerate(seq):
        if aa == "-":
            continue
        pdbname, pdbnum = pdbseq.pop(0)
        assert pdbname == aa, \
            f"Sequence mismatch at position {idx}: {aa} is not {pdbname}"
        mapping[idx] = pdbnum

    return mapping


def get_glycans(structure):
    """Returns a dictionary of the glycans structures in the PDB."""
    # Used to remove glycans from model when there is no
    # ASN due to mutation.

    aaset = set(protein_letters_3to1)

    atoms = [a for a in structure.get_atoms() if a.element != "H"]
    ns = NeighborSearch(atoms)

    # First pass. Find ASNs bound to sugars (roots).
    roots = set()
    nags = [r for r in structure.get_residues() if r.resname not in aaset]
    for nag in nags:
        c1_nag = nag['C1']
        for n in ns.search(c1_nag.coord, 2.5, 'A'):
            if n.parent.resname == 'ASN' and n.name == 'ND2':
                roots.add(n.parent)
                break

    # Second pass. Find branched structures.
    branches = {r: set() for r in roots}

    for r in roots:
        pool = set()  # pool of neighbors to visit

        # seed pool with immediate neighbors
        for atom in r.get_unpacked_list():
            neighbors = ns.search(atom.coord, 2.0, "R")
            for n in neighbors:
                is_aa = n.resname in aaset
                is_root = n is r
                is_same_chain = n.parent is r.parent
                have_visited = n in branches[r]

                if not is_aa and not is_root and is_same_chain and not have_visited:
                    pool.add(n)

        # Search neighbors "recursively"
        while True:
            for p in list(pool):
                branches[r].add(p)
                for atom in p.get_unpacked_list():
                    neighbors = ns.search(atom.coord, 2.0, "R")
                    for n in neighbors:
                        is_aa = n.resname in aaset
                        is_root = n is r
                        is_same_chain = n.parent is r.parent
                        have_visited = n in branches[r]

                        if not is_aa and not is_root and is_same_chain and not have_visited:
                            pool.add(n)
                pool.remove(p)

            if not pool:
                break
    # debug
    # for r, g in branches.items():
    #     print(
    #         f"{r.id[1]} {r.parent.id}",
    #         "=>",
    #         ",".join(sorted(str(gg.id[1]) for gg in g))
    #         )

    # Change branches to index by resi id/chain tuple
    return {
        (r.id[1], r.parent.id): g for r, g in branches.items()
    }


ap = argparse.ArgumentParser()
ap.add_argument('config', type=pathlib.Path)
args = ap.parse_args()

# Read config file
config = configparser.ConfigParser()
config.read(args.config)

# Read in MODELLER Python scripts
modloop, modbb, fixatom = None, None, None
moddir = pathlib.Path(config['modeller']['datadir'])
for pyf in moddir.rglob('*.py'):
    if 'loops' in str(pyf):
        modloop = pyf
    elif 'bb' in str(pyf):
        modbb = pyf
    elif 'fixatom' in str(pyf):
        fixatom = pyf

assert all([modloop, modbb, fixatom])

# Read template structure
template = read_pdb(config['modeller']['template'])
pdb_glycans = get_glycans(template)

# Read alignment files
# First sequence is template
seqs_ace2 = read_seqfile(config['sequences']['ace2_seqfile'])
seqs_vrbd = read_seqfile(config['sequences']['vrbd_seqfile'])

ref_ace2 = list(seqs_ace2)[0]
ref_vrbd = list(seqs_vrbd)[0]
assert "6M0J" in ref_ace2, f"Wrong reference: {ref_ace2}"
assert "6M0J" in ref_vrbd, f"Wrong reference: {ref_vrbd}"

# Hack to change ACE2 reference name here to match species
nname = "Homo_sapiens"
seqs_ace2[nname] = seqs_ace2[ref_ace2]
del seqs_ace2[ref_ace2]
ref_ace2 = nname

# Map ref sequence numbering to template. We can be lazy and not do
# an alignment since we use the template sequences as references in
# the alignments.
ace2_mapping = map_seq_to_pdbchain(seqs_ace2[ref_ace2], template, "A")
vrbd_mapping = map_seq_to_pdbchain(seqs_vrbd[ref_vrbd], template, "E")

# If firstonly is False, model all species, otherwise just first
do_all = bool(config['modeller']['firstonly'])
if do_all:
    paired_seqs = sorted(itertools.product(seqs_ace2, seqs_vrbd))
else:
    paired_seqs = sorted(itertools.product([ref_ace2], seqs_vrbd))
print(f'Will generate {len(paired_seqs)} models')

# Read template alignment
alnfpath = pathlib.Path(config['modeller']['datadir']) / 'ali.fasta'
with alnfpath.open() as handle:
    alnfile = handle.read()

# Create directories if necessary
curdir = pathlib.Path('.').resolve()
modeldir = pathlib.Path(config['modeller']['modeldir'])
modeldir.mkdir(exist_ok=True)

# Prepare and launch modelling
writer = PDBIO()

ace2_refseq = seqs_ace2[ref_ace2]
vrbd_refseq = seqs_vrbd[ref_vrbd]
for ace2, vrbd in paired_seqs:

    title = f"{ace2}-{vrbd}"
    subdir = modeldir / title

    if subdir.exists():
        done_flag = subdir / 'DONE'
        if done_flag.exists():
            print(f"Skipping {subdir.name}", flush=True)
            continue

        print(f"Emptying old contents for {subdir.name}", flush=True)

        # Remove all contents and re-run
        for f in subdir.rglob('.*'):
            f.unlink()
            pass
    else:
        subdir.mkdir()

    #
    # Prepare template and alignment
    #
    ace2_glycosites = []
    vrbd_glycosites = []
    n_ace2_glycans = 0
    n_vrbd_glycans = 0

    # Compare model and template to see if we have mutations
    # on sugar-linked ASNs.
    for idx, (aa_seq, aa_ref) in enumerate(zip(ace2_refseq, seqs_ace2[ace2])):
        if aa_seq == aa_ref == "N":
            pdbnum = ace2_mapping.get(idx)
            if pdbnum is not None:
                rid = (pdbnum, "A")
                glycans = pdb_glycans.get(rid)
                if glycans:
                    ace2_glycosites.append(rid)
                    n_ace2_glycans += len(glycans)

    for idx, (aa_seq, aa_ref) in enumerate(zip(vrbd_refseq, seqs_vrbd[vrbd])):
        if aa_seq == aa_ref == "N":
            pdbnum = vrbd_mapping.get(idx)
            if pdbnum is not None:
                rid = (pdbnum, "E")
                glycans = pdb_glycans.get(rid)
                if glycans:
                    vrbd_glycosites.append(rid)
                    n_vrbd_glycans += len(glycans)

    # Replace aligned template seqs
    a = alnfile.replace('ACE2_TEMPLATE', ace2_refseq + "." * n_ace2_glycans)
    a = a.replace('VRBD_TEMPLATE', vrbd_refseq + "." * n_vrbd_glycans)
    # Replace model sequences
    a = a.replace('ACE2_SEQ', seqs_ace2[ace2]  + "." * n_ace2_glycans)
    a = a.replace('VRBD_SEQ', seqs_vrbd[vrbd] + "." * n_vrbd_glycans)
    # Replace title
    a = a.replace('SEQNAME', title)

    fpath = subdir / "ali.pir"
    with fpath.open("w") as handle:
        handle.write(a)
    
    # Now edit template to remove glycans if necessary.
    template_copy = copy.deepcopy(template)
    model_glycosites = set(ace2_glycosites + vrbd_glycosites)
    n_gly_removed = 0
    for root, glycans in pdb_glycans.items():
        if root in model_glycosites:
            continue  # leave it be

        # For debug only.
        # print(f"Will remove glycans attached to: {root}")
        for g in glycans:
            chain = g.parent.id
            pdbchain = template_copy[0][chain]
            pdbchain.detach_child(g.id)
        n_gly_removed += 1
    
    # Save template PDB
    writer.set_structure(template_copy)
    writer.save(str(subdir / "template.pdb"))

    #
    # Finally, copy MODELLER scripts.
    #
    # Do we need loop modelling (have insertions in the model?)
    have_ace2_inserts = ace2_refseq.count("-") < seqs_ace2[ace2].count("-")
    have_vrbd_inserts = vrbd_refseq.count("-") < seqs_vrbd[vrbd].count("-")
    if have_ace2_inserts or have_vrbd_inserts:
        shutil.copy(modloop, subdir / 'cmd_modeller.py')
    else:
        shutil.copy(modbb, subdir / 'cmd_modeller.py')
    shutil.copy(fixatom, subdir / fixatom.name)  # copy fixatoms.py


    # Launch MODELLER
    os.chdir(subdir)
    print(f'[{time.ctime()}] Running {subdir.name}', flush=True)
    print(f"\tACE2 glycans: {n_ace2_glycans}")
    print(f"\tvRBD glycans: {n_vrbd_glycans}")
    print(f"\tTotal N-glycosites removed: {n_gly_removed}")
    with open('modeller.log', 'w') as logfile:
        with open('modeller.err', 'w') as logerr:
            p = subprocess.Popen(
                'python cmd_modeller.py',
                shell=True,
                stdout=logfile,
                stderr=logerr,
                close_fds=True
            )
            p.communicate()

    os.chdir(curdir)
    
print("Done")
