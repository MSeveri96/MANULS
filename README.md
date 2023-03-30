# MANULS (Make chemicAl reactioNs spontaoUs via optimaL fieldS)


<div align="center">
<img src="./MANULS_logo_small.png" alt="MANULS_logo" width="1500">
</div>


A program to calculate the smallest electric field to make a chemical reaction barrierless.

This first version of the code works with 1-D and  2-D scans of the potential energy surface. A script is provided to extract the relevant data from ORCA relaxed scan calculations.

## Jump-Start Guide 

### Documentation

You can find the documentation in PDF format, along with a worked example, [here](https://github.com/MSeveri96/MANULS/MANULS_documentaton.pdf) in this repository 

### Prerequisites

1. Python 3
2. Numpy
3. Matplotlib (optional)

Usually the `numpy` and `matplotlib` libraries can be installed using the `pip` or `pip3` (for python3) package. 

For instance: `pip install numpy matplotlib`

### Download and Usage

XXX can be dowloaded using the `git clone` command or using the green button at the top of this page. After pressing the button you can choose "download ZIP", unzip the folder and run the code.

To run the the code you need to:
1. Prepare the data for the code (see the [documentation](https://github.com/MSeveri96/Optimal-external-electric-field/blob/main/Documentation/README.md)). We provide a script to automatically extract the data from a ORCA calculation.
2. Unzip in a directory of your choice.
3. Modify the file xxx.py and yyy.py files putting the correct location of the files containing your data.
4. Run the code with `python3 xxxx.py` or `python xxx.py`. If you prefer you can run the code in your favorite IDE as well.

### Output
 
 The code produces two files: `optimal_field_and_optimal_bbp.txt` and `fields_and_bbps.txt`.
 
 `optimal_field_and_optimal_bbp.txt` contains, at its top, the coordinates of the optimal bond-breaking-point and the optimal external electric field, along with other data to control the quality of the result. It also contains four "sub-optimal" fields and "sub-optimal" bond breaking points. These "sub-optimal" values are not optimal but close to the optimality. 
 
 `fields_and_bbps.txt` is mainly for the interested user and for debugging purposes. More information is available in the [documentation](https://github.com/MSeveri96/Optimal-external-electric-field/blob/main/Documentation/README.md).
 
 

## Citation

## Questions?

If something is unclear or for any kind of question you are encouraged to write to the maintainer of this repository (marco.severi6@unibo.it) or to post in the [Discussions section](https://github.com/MSeveri96/PALLAS/discussions). Every comment, question or suggestion is more than welcome!

If you spot a bug in the code (thanks!) you can post in the [Issues section](https://github.com/MSeveri96/PALLAS/issues) or write to the maintainer.
