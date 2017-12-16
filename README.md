![code coverage badge here](img/coverage.svg)

CHE477 Langevin Project
====

*Clyde Overby*

Overview
======

This is a Langevin dynamics project for CHE477.  Inputs are given as two files: one that includes all of the particles in the format of index, position, potential energy, and force and a second file that includes the temperature, damping coefficient, time step size, and total integration time.

Output file is given as lines containing index, time, position, and velocity separated by tabs.

To use, install setup.py and enter "langevin" into the command line while in this directory.  You will be prompted for a position/energy file (specified above), a parameter file (specified above), particle mass, particle initial velocity, output file name, and how many timesteps (counting backwards from final) for each particle you want included in the output.  

Example input files (positionfile1.txt and paramfile1.txt) as well as a sample output (testout1.txt) are provided.  testout1.txt was generated with mass = 1, velocity =1, and output lines = 5.

(c) 2017