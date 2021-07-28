# AlphaFoldDB_Processing
This script can clean up the disordered region in the AlphaFold DB structures to adapt it to the molecular docking, extract the sequence of the disordered region to make it available for sequence motif searching.

Introduction
===

Traditional methods explore potential Protein-Protein Inteaction and Protein-Ligand Interaction by molecular docking. The molecular docking is highly dependent on the protein structure, the [Alphafold DB](https://alphafold.ebi.ac.uk/) gives us a large amount of high quality structures. Howerver, many low confidence structures of disordered regions make it difficult for these structures to be used for molecular docking. 

Here, we provide a scirpt (Domain_Parser-BioStructures.jl) to extract ordered structures from the structure of alphafold DB for molecular docking.

<div align="center">
<img src=https://user-images.githubusercontent.com/58931275/127297567-36ab0c1e-f97b-44fa-a0a3-ee2f9c715699.png width=60% />
</div>
    
Many disordered regions play critical roles in protein-protein interactions, and we also provide an additional script (Disorder_Reader.jl) to extract the sequences of disordered regions from alphafold DB for motif searching.

Dependencies
===

1. [Julia](https://julialang.org/downloads/)
2. [BioStructures](https://github.com/BioJulia/BioStructures.jl) library in julia. 

Usage
===

To read the disorder region from AF2 structures.

    julia -p 8 Disorder_Reader.jl

![image](https://user-images.githubusercontent.com/58931275/127270117-22e8ce78-e542-4161-8066-2ec7d3b1c9ca.png)


To extract the ordered regions from AF2 structures.

    julia -p 8 Domain_Parser-BioStructures.jl

![image](https://user-images.githubusercontent.com/58931275/127270091-068d684e-620b-4b60-bcd5-4ed11247faae.png)



