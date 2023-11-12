# PyRosetta source code for PhaC1437 engineering

This source code was implemented based on PyRosetta in a Python 2.7 environment, and was confirmed to work in Ubuntu 12.04 version.
A complex structure in which phenyllactyl-CoA is bound to the active site of PhaC is used as input. The mutation site is determined within a predetermined range near the active site of PhaC (see lines 15-16 of the source code).
PyRosetta can be installed after obtaining a license and is free for academic users.


## Installation instructions for PyRosetta in Linux
GNU/Linux and Mac OS X
1. Obtain a Rosetta license from to receive a username and password.

2. Download the appropriate version of PyRosetta from the links above.

3. Unpack the downloaded file to the location of your choice to create the PyRosetta directory.

4. (From a terminal/console window, you can unpack the archive using the command: tar -vjxf PyRosetta-<version>.tar.bz2. Please note, there is no special install procedure required; after unpacking, PyRosetta is ready to use. So unpack it to the location from where you want to execute it.)

5. From within the new PyRosetta directory, type source SetPyRosettaEnvironment.sh into the command line to set up the PyRosetta library file paths.

6. Start Python.

7. In Python, you should be able to import the PyRosetta library with the command import rosetta; rosetta.init().
   
        (If this step does not produce a complaint or error, your installation has been successful.)

For detailed instructions, you can check the installation and usage instructions through the link below.
        https://www.pyrosetta.org/downloads/legacy-pyrosetta3-download

