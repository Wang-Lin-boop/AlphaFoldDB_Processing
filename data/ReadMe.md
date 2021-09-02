Disorder Protein Sequences
---
A fasta file contains all the disorder regions in AlphaFold Human Protein Structure Database. 

This fasta file are produce in two steps:
1. The models were segment into different fargments by pLDDT of residues within 1-50, but the last 2 residues per fargments near to next ordered fargment were discard.
2. Writing the sequences of all disorder fargments into `Human_Disorder_Protein_Sequences.fasta`. 
 

