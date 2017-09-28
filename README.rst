
=======================
Poiseuille Flow Example
=======================

This example demostrates Poiseuille flow in a tube. The Navier-Stokes equations are solved in a multi-block mesh of the tube. The velocity uses triquadradtic-Lagrange elements and the pressure trilinear-Lagrange elements. Residual base stabilisation is used. A CellML model is used to specify the time-varying Poiseuille flow input velocities. 

Command Line arguments
======================

Up to four command line arguments can be specified. They are (in order):
* number of square elements in the multi-block mesh
* number of arm elements in the multi-block mesh
* number of length elements in the multi-block mesh
* the Reynolds number
* the maximum input flow velocity
* the start time
* the stop time
* the time step


