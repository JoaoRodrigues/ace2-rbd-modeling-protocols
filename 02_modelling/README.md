# Modelling of ACE2:RBD complexes

## Template Structure
The template structure is a refined model based on the crystal structure
6m0j and including full-length glycans, made available by John Chodera
and co-workers at https://github.com/choderalab/rbd-ace2-contact-analysis.

To minimize run time for our refinements, we cut the C-terminal portion of
the template (keep region 19-615) and any associated glycans. We also removed
the Cl and Zn ions associated with ACE2, since these are not at the interface
with the viral RBD and as such, would be restrained during refinement.

Sugars were renamed to match HADDOCK topologies.

## Protocol

The modelling will proceeed in an all-vs-all manner, using the settings defined
in the `modeling.cfg` file. In a first step, we generate all alignments for
MODELLER, to allow visual inspection if necessary. Then, we use MODELLER to
generate models for each ACE2:RBD ortholog/variant pair.

### Sequences and Template/Model Alignments
We used the aligned RBD sequences in  `00_gisaid-processing`) and the aligned
ACE2 sequences from `01_ace2-processing`. No re-alignment is done here.

### Modelling
For each ACE2:RBD ortholog/variant per, we will generate 10 models. To keep the models
as true to the template as possible (and since MODELLER doesn't handle glycans), we
restrain the entire structure except for residues that differ between the sequence and
the template (to relax/avoid clashes). We also delete any glycans if the corresponding
ASN residue was mutated (e.g. plenty in Mouse ACE2) before creating the models.