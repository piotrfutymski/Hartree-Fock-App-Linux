# Hartree Fock App

## Table of contents
* [General info](#general-info)
* [Setup](#setup)
* [Input file](#input-file)
* [H2O orbital](#H2O-orbital)

## General info

Calculates energies and drows plot of molecular orbitals for small molecules (up to 3-4 nucleons) form hydrogen to neon.
	
## Setup

First you need to download and install armadillo.
You can download this from: http://arma.sourceforge.net/download.html

Next install it:
1. You will need cmake

```
$ sudo apt-get install cmake
```
2. If not present already, install LAPACK and BLAS:
```
$ sudo apt-get install liblapack-dev
$ sudo apt-get install libblas-dev
```
3. Open terminal into the directory that was created by unpacking the Armadillo archive
Type:
```
$ cmake .
$ make
$ sudo make install
```

After this you can download repository:

```
$ git clone https://github.com/piotrfutymski/Hartree-Fock-Attempt-2.git
$ cd Hartree-Fock-Attempt-2
$ cmake .
$ make
```
To run script you will need matplotlib and tornado:

```
$ pip install matplotlib
$ pip install tornado
```

Now you can run bash script with example input like this:

```
$ chmod +x HF.sh
$ ./HF.sh examples/H2O.txt out.txt
```
Output will be find in out.txt and molecular orbitals grqphs in pic folder.

## Input file

For example:
```
HF MO_GRAPHS DENSITY_GRAPH 
*xyz 0 1
Li 	0.0000 	0.0000 	0.0000
H 	0.0000 	0.0000 	1.5949

*
*name LiH *
*plain yz *
```
You can change/add nucleons and positions. Choose from: [ H He Li Be B C N O F Ne ].
All positions are in angstroms.


## H2O orbital

https://github.com/piotrfutymski/Hartree-Fock-App-Linux/blob/master/pic/H2O_2.png
