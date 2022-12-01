# NIM-2D-NSE-POLAR-CAVITY
The solution of the time-dependent, incompressible, Navier-Stokes equations in curvilinear coordinates using nodal integral method. <br />
<br />
The code is written by Ibrahim Jarrah. If you have any questions, feel free to reach out at [ibrahim.jarrah92@gmail.com](mailto:ibrahim.jarrah92@gmail.com).
# Test case: Lid-Driven Flow in a Polar Cavity
The details of the case can be found [here](https://www.ideals.illinois.edu/items/124600) in page 118.<br />
# Silo Library
[Silo](https://github.com/LLNL/Silo) library is required to run the code. You can istall it and provide the path to silo in the Makefile, or use the provided script [install.sh](3rd_part/install.sh).
# Installation
```
git clone https://github.com/ijarrah2/NIM-2D-NSE-POLAR-CAVITY.git
cd NIM-2D-NSE-POLAR-CAVITY/3rd_party
./install.sh
cd ..
```
# Running the case
```
./run.sh
```
# Running the case with different parameters
To change the test case setup, modify the [user_parameters.f90](user_parameters.f90) file.
