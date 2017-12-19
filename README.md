![code coverage badge here](img/coverage.svg)

# CHE477 1D Numerical Schroedinger Solver

*Author: Clyde Overby*


## Installation

Install using 
```sh
$ pip install setup.py
```
## Usage

After installation, calling
```sh
$ schro
```
Will initate the script, upon which it will as for an input file and a destination file.

## Overview

This is a 1D Schrodinger equation solver for Fourier and Legendre Polynomial basis sets using the variational method.  Input files are given in a CSV as :
Input file contains index, potential energy(V0), kinetic energy constant (c), basis set cooefficient number (bs-coefficients), basis set choice (basis-set), and domain.

```sh
###Index:
(integer)
###V0:
(integer or float)
###constant:
(integer or float)
###bs_coefficients:
(integer)
###basis set:
(char) either l for Legendre or f for Fourier
###domain:
(float, float) (-1, 1) for nearly all cases
```

Output file is in the form of index: [coefficients] where the coefficients are those that satisfy the Schrodinger equation for the basis set chosen.

(c) 2017