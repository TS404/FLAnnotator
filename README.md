# FLAnnotator
A set of functions for annotation of Fasciclin-like Arabinogalactan sequences

### Background

Fasciclin-like arabinogalactan-proteins (FLAs) are cell-wall associated glycoproteins found across the pant kingdom. They contain both globular fasciclin domains (Pfam PF02469) as well as long disordered regions. Portions of these disordered regions contain 'glycomotifs' which direct the hydroxylation and glycosylation of the prolines within those motifs.

The high sequence diversity of the fasciclin domains, and the _extremely_ high sequence diversity of the disordered regions means that traiditional bioinformtaic methods have to be supplemented with techniqes specialised to sequences with such properties. 

### Functions

- Split FLAs into:
  - fasciclin domains
  - arabinogalactan regions 
  - non-arabinogalactan regions
- Numericise sequence data
  - biophysical properties of residues in an alignment of fasciclin domain sequences
  - kmer fequencies of regions outside of fasciclin domains (non-alignable)
- Transform and project multidimensional numericised data down into a simplified sequence space and identify clusters
  - by PCA & bayesian clustering
  - by UMAP & HDBscan
- Visualise reulsts
  - sequence space projections
  - domain diagrams
