# Data and Code for 8-mer TcdA peptides
This repository contains the data and codes associated with the paper-Design of 8-mer peptides that Block Clostridoides difficile Toxin A in Intestinal Cells, submitted to Communications Biology.


# Requirement and Installation
The source codes are in /code/. Compiler such as gfortran or ifort is required

# Getting Started
/code/-This directory contains the following files and directories:
1. main.f90-source code
2. /lib/-This directory contains the forcefield and rotamer parameters of all the 20 natural amino acids
3. comp.pdb-An example of the input PDB file needed to start the design
4. input.txt-Input file to start design
5. pdbfiles-This is an empty directory where the pdbfiles will be stored of the new peptides that are designed

# Data for paper
/SA1 to SA7/-This directory contains the following files:
1. SA1-SA7: Each folder contains the topology file and final PDB file from a 100 ns simulation for runs 1 to 3

# Reproducing figures
/Figures/-This directory contains the raw data for all plots in the Main Text and Supplemnetary Material. It also contains the matlab code for image processing  in the bead assay


