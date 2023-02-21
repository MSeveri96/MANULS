# Optimal External Electric Field

A program to calculate the external electric field that makes a chemical reaction barrierless.



This first version of the code works with 1-D and  2-D scans of the potential energy surface.

## Prerequisites
The script has been tested with:
1. Python 3.8
2. Numpy 1.19.5
3. Scipy 1.9.3
4. Matplotlib 3.3.4 (optional)

Usually these packages can be installed using the 'pip'

## Preparation of the input
( These are some general guidelines, in the follwing sections there is a detailed discussion of the example)

The data need to be grouped in a simple text file which contains the following:

1. First column: parameter that describe the x-axis (in the example the P-C distance for each point)
2. Second column: parameter that describe the y-axis (in the example the C-O distance for each point)
3. Third column: energy of printed by the quantum chemical program
4. Fourth, fifth and sixth columns: x, y and z component of the dipole moment vector. Since the code computes the second derivatives of the dipole it is very important that the origin of the dipole moment vector is the same for each point of the scan.

## Preparation of the script

The preparation of the script is minimal: the user just need to modify the step size and the number of points of the grid. The step size must be constant in the calculation.

## Running the script

Now you're ready to find the OEEF!

## Detailed discussion of the example

### Scanning the PES

