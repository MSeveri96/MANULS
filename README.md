# Optimal External Electric Field

## Disclaimer 
Please cite **J. Chem. Theory Comput. 2022, 18, 2, 935–952** and **xx--xx--xx** if you use any part of the following code in your work!

##
This simple repository contains the script used to find the Optimal External Electric Field (OEEF) to make the Wittig reaction barrierless, as described in xx-xx-xx. This code also calculates the Optimal Bond Breaking Point (OBBP). The OBBP is a special point of the potential energy surface. If the system is in the OBBP and is perturbed with the OEFF the reaction becomes barrierless.

This first version of the code works with 2-D scans of the potential energy surface.

## Preparation of the input

The data need to be grouped in a simple text file which contains the following:

1. First column: parameter that describe the x-axis (in the example the P-C distance)
2. Second column: parameter that describe the y-axis (in the example the C-O distance)
3. Third column: energy of printed by the quantum chemical program
4. Fourth, fifth and sixth columns: x, y and z component of the dipole moment vector
