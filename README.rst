
===============
Poiseuille Flow
===============

This example demostrates Poiseuille flow in a tube. The Navier-Stokes equations are solved in a multi-block mesh of the tube. The velocity uses triquadradtic-Lagrange elements and the pressure trilinear-Lagrange elements. Residual base stabilisation is used. A CellML model is used to specify the time-varying Poiseuille flow input velocities.

Running the example
===================

Python version::

  source /path/to/opencmisslibs/install/virtaul_environments/oclibs_venv_pyXY_release/bin/activate
  cd ./navierstokes_poiseuilleflow/src/python
  python src/python/poiseuille.py

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

Verifying the example
=====================

Results can be visualised by running `visualise.cmgui <./src/python/visualisePoiseuille.cmgui>`_ with the `Cmgui visualiser <http://physiomeproject.org/software/opencmiss/cmgui/download>`_.

.. figure:: docs/images/PoiseuilleFlow_0.png
   :align: center
   :width: 30%

   **Figure:** Mesh structure for the poiseuille flow example.

The expected results from this example are available in `expected_results <./src/python/expected_results>`_ folder.

Prerequisites
=============

There are no additional input files required for this example as it is self-contained.

License
=======

License applicable to this example is described in `LICENSE <./LICENSE>`_.
