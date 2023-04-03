#!/bin/bash

grep -A3 'raw cartesian tensor' $1.out |sed '/raw/d'|sed '/--/d'> $1_polarizability.txt

grep 'Total Dipole' $1.out > tmp.txt

cut -d " " -f 10- tmp.txt > dipole.txt

paste $1.relaxscanact.dat dipole.txt > $1_grid.txt

rm tmp.txt dipole.txt

