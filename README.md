# AlphaFold_Processing
This script can clean up the disordered region in the AlphaFold DB structures to adapt it to the molecular docking, extract the sequence of the disordered region to make it available for sequence motif searching.

Usage
===

To read the disorder region from AF2 structures.

    julia -p 8 Disorder_Reader.jl

To extract the ordered regions from AF2 structures.

    julia -p 8 Domain_Parser-BioStructures.jl


Dependencies
===

1. [Julia](https://julialang.org/downloads/)
2. [BioStructures](https://github.com/BioJulia/BioStructures.jl) library in julia. 

