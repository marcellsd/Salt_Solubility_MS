## Introduction

This is a Python program to calculate the salt solubility in a mixture of solvents (water + MEG).

## Requirements

You will need Python 3.7.0 or higher version installed on your OS with the following dependecies:
* Numpy
* Scipy
* Pandas
* Matplotlib
* tkinter
* pip
* pyswarm

 To install these dependencies on Ubuntu 20.04:

```
sudo apt-get install python3-numpy python3-scipy python3-matplotlib python3-pandas python3-tk python3-pip
```
And so:

```
pip install --upgrade pyswarm
```

The input data file must be in a .xlsx format, following the examples of Data folder for single or multi temperatures. 

## Usage

Clone this repository, cd into the folder directory:

```
cd /SALT_SOLUBILITY_MS
```

then run:

```
python3 multi_isothermal.py
```
for data with 3 or more temperatures.
Or:

```
python3 single_isothermal.py
```

for data with just one temperature.

The outputs are pdf file with the graphics and a xlsx file with the calculated values.

## Comments

Feel free to let me know if you have any questions or comments about this code.

## Credits

