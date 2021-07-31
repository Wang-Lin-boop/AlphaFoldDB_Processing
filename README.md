# AlphaFoldDB_Processing
This script can clean up the disordered region in the AlphaFold DB structures to adapt it to the molecular docking, extract the sequence of the disordered region to make it available for sequence motif searching.

Introduction
===

Traditional methods explore potential Protein-Protein Inteaction and Protein-Ligand Interaction by molecular docking. The molecular docking is highly dependent on the protein structure, the [Alphafold DB](https://alphafold.ebi.ac.uk/) gives us a large amount of high quality structures. Howerver, many low confidence structures of disordered regions let it difficult for these structures to be used for molecular docking. 

Here, we provide a scirpt (Domain_Parser-BioStructures.jl) to extract ordered structures from the structure of alphafold DB for molecular docking.

<div align="center" style="border-top:1px solid black;border-botton:1px">
<img src=https://user-images.githubusercontent.com/58931275/127298352-f0da558d-9486-4f1c-aace-a4beb45b9f10.png width=60% />
<p> </p>
</div>


Many disordered regions play critical roles in protein-protein interactions, and we also provide an additional script (Disorder_Reader.jl) to extract the sequences of disordered regions from alphafold DB for motif searching.

Dependencies
===

1. [Julia](https://julialang.org/downloads/)
2. [BioStructures](https://github.com/BioJulia/BioStructures.jl) library in julia. 

Installtion
====
1. Install [Julia](https://julialang.org/downloads/) on your computer, Linux, macOS or Windows.
2. run Julia, click julia.exe on windows, run ```julia``` on Linux or macOS.
3. typing ```]add BioStructures``` in RPEL and pressing enter to add [BioStructures](https://github.com/BioJulia/BioStructures.jl) library in julia. 

Usage
===

Linux:
---

To read the disorder region from AF2 structures.

    julia -p 8 Disorder_Reader.jl

![image](https://user-images.githubusercontent.com/58931275/127270117-22e8ce78-e542-4161-8066-2ec7d3b1c9ca.png)


To extract the ordered regions from AF2 structures.

    julia -p 8 Domain_Parser-BioStructures.jl

![image](https://user-images.githubusercontent.com/58931275/127270091-068d684e-620b-4b60-bcd5-4ed11247faae.png)

Windows:
---

Drag scripts to your working directory containing your AF2 DB folder， then clicking scripts and select julia.exe to run:

![image](https://user-images.githubusercontent.com/58931275/127741637-b34e6686-ea58-44c8-bfb0-99d0091f8ac7.png)

Typing your AF2_DB foldername in the cmd and running it. 
