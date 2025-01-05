# MANULS (Make chemicAl reactioNs spontaneoUs via optimaL fieldS)


<div align="center">
<img src="./MANULS_logo.png" alt="MANULS_logo" width="450">
</div>

Version 1.0.1 is out!

A program to calculate the smallest electric field that makes a chemical reaction barrierless. In other words: with MANULS it is possible to compute the external electric field that enables the perfect catalysis of a chemical reaction.

This first version of the code works with 1-D and  2-D scans of the potential energy surface. A script is provided to extract the relevant data from ORCA relaxed scan calculations.

## Capabilities
MANULS computes the smallest electric field that removes the energy barrier between the reactants and the transition state of a chemical reaction. At the moment the code work with potential energy curves and 2-D potential energy surfaces. 

The following figures show one prototypical application of MANULS. The code starts from a potential energy curve computed without any field effecs, in vacuo (on the left). In this case the blue curve describes a trans-to-cis isomerization, that is the reactants are on the right of the blue curve. The application of MANULS on the blue curve gives an external electric field that transforms the blue curve into the orange one. In other words the external electric field computed by MANULS is able to remove the reaction barrier and catalyze the reaction. 


| Without Electric Field                | With Electric Field                          |
| ----------------------------------- | ----------------------------------- |
| ![no_field](https://github.com/MSeveri96/MANULS/blob/main/original_pes.png) | ![field](https://github.com/MSeveri96/MANULS/blob/main/perturbed_pes_polar.png) |


## Jump-Start Guide 

### Documentation

You can find the documentation in PDF format, along with a worked example, [here](https://github.com/MSeveri96/MANULS/blob/main/MANULS_documentation.pdf) in this repository 

### Prerequisites

1. Python 3
2. Numpy
3. Matplotlib (optional)

Usually the `numpy` and `matplotlib` libraries can be installed using the `pip` or `pip3` (for python3) package. 

For instance: `pip install numpy matplotlib`

### Download and Usage

MANULS can be dowloaded using the `git clone` command or using the green button at the top of this page. After pressing the button you can choose "download ZIP", unzip the folder and run the code.

To run the the code you need to:
1. Have a scan of the potential energy surface, along with the dipole moment vector and the polarizability for each point of the scan.
2. Prepare the data for the code (see the [documentation](https://github.com/MSeveri96/MANULS/blob/main/MANULS_documentation.pdf)). We provide a script to automatically extract the data from a ORCA calculation.
3. Unzip MANULS in a directory of your choice.
4. Modify the file MANULS_2D.py or MANULS_1D.py files putting the correct location of the files containing your data and the scan.
5. Run the code with `python3 MANULS_2D.py` or `python3 MANULS_1D.py`. If you prefer you can run the code in your favorite IDE as well.

### Output
 
 The code produces two files: `external_electric_fields.txt` and `coordinates_bbps.txt`.
 
`external_electric_fields.txt` contains, at its top, the coordinates of the optimal bond-breaking-point and the optimal external electric field, the optimal external electric field is the field with the smallest amplitude capable of removing the energy barrier of a chemical reaction.
 
`coordinates_bbps.txt` is mainly for the interested user and for debugging purposes. More information is available in the [documentation](https://github.com/MSeveri96/MANULS/blob/main/MANULS_documentation.pdf).

We provide the `MANULS_plot_2D.py` and `MANULS_plot_1D.py` scripts to plot the results.
 
 

## Citation
The main publication related to MANULS is "An algorithm to find the optimal oriented external electrostatic field for annihilating a reaction barrier in a polarizable molecular system", J. Chem. Phys. 159, 114112 (2023). 

The user is encouraged to cite this paper whenever results are obtained with the MANULS program.  J. Chem. Phys. 159, 114112 (2023) is also the reference from which the description of the model is taken. Other works regarding the optimal bond-breaking point and the optimal external electric field are
1. Josep Maria Bofill, Wolfgang Quapp, Guillermo Albareda, Ibério de P. R. Moreira, and Jordi Ribas-Ariño
Journal of Chemical Theory and Computation 2022 18 (2), 935-952
DOI: 10.1021/acs.jctc.1c00943 
2. Bofill, J.M., Quapp, W., Albareda, G. et al. A catastrophe theory-based model for optimal control of chemical reactions by means of oriented electric fields. Theor Chem Acc 142, 22 (2023). https://doi.org/10.1007/s00214-023-02959-0
3. Josep Maria Bofill, Marco Severi, Wolfgang Quapp, Jordi Ribas-Ariño, Ibério de P. R. Moreira, Guillermo Albareda. Optimal oriented external electric elds to trigger a barrierless oxaphosphetane ring opening step of the Wittig reaction. Chem. Eur. J. 2024, e202400173. https://doi.org/10.1002/chem.202400173


## Questions?

If something is unclear or for any kind of question you are encouraged to write to the maintainer of this repository (marco.severi6@unibo.it) or to post in the [Discussions section](https://github.com/MSeveri96/MANULS/discussions). Every comment, question or suggestion is more than welcome!

If you spot a bug in the code (thanks!) you can post in the [Issues section](https://github.com/MSeveri96/MANULS/issues) or write to the maintainer.
