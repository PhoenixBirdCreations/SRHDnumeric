# Approximate Riemann solver for SRHD with a realistic cold equation of state (EoS)

This folder contains a stable implementation in C language of a numerical Riemann solver for Special Relativistic Hydrodynamics (SRHD) with third order of accuracy.
The code was developed by Marina Berbel (me) in 2023 for personal use in her PhD thesis: "On nonconvex Special Relativistic Hydrodynamics".

Copyright 2023, Marina Berbel, All rights reserved

## The numerical flux is Harten-Lax-van Leer (HLL) [REF], adapted to nonconvex dynamics
## The EoS models coded are tabulated, piecewise polytropic (PP) [Read et al. Phys. Rev. D, 2009] and thermodynamically adaptive slope piecewise polytropic (T-ASPP) [Berbel, Serna. Phys. Rev. D, 2023]

This folder contains:
### Libraries
Functions implemented for the solver, separated by topic for readability.

* auxFun: auxiliar functions for basic math, initialization of quantities and output of solution
* EOS_cold: all functions related to the EoS 
* highOrder: functions to reconstruct variables at high order 
* itmethod: also known as con2prim. Recovery of primitive variables from conserved ones. Uses fix point iterations.
* HLLrel: implementation of HLL for nonconvex fluxes at high order, along a Runge-Kutta 3 for the high order in time
* reconstructions: several interpolation methods for the high order spatial reconstruction
* rel_euler: eigenvalues, eigenvectors and fluxes of SRHD equations in 1D

### Main execution script
Reads the par files and start the integrations. Meassures execution time and writes output to a file.

### Technicalities 
**Allowupwind**: this variable takes values 0 (false) or 1 (true). This is not used for HLL flux, but it is inherited from the MFF implementation in order to share the initial conditions file.


**Method**: this variable takes values 1 or 2.

Value 1 selects piecewise hyperbolic method (phm) [Marquina. J. Sci. Comput. 1994] 3rd order reconstruction. This is the one recommended.

Value 2 selects esentially non oscillatory (eno) [Harten, Osher. J. Numer. Anal. 1987] 3rd order reconstruction.


### Use of the code 
Download all files and put them in the same directory.

There is a Makefile distributed with the source code, simply change the name of the C compiler, if you use a different one, and type make.
It needs no external libraries besides the standard math, stdlib, string, time and stdio.

The executable is named Rsolver_MFF and runs without any parameters. The details of the EoS are specified in a par file, eos.par.
The initial conditions and the solver details are specified in an ic.par file.
Examples of both type of par files are included in the folder.

Once the integration is finished, the program creates a txt file in the path specified with the solution. It is stored in columns of the format
1-position &nbsp; &nbsp; &nbsp; 2-density &nbsp; &nbsp; &nbsp; 3-velocity &nbsp; &nbsp; &nbsp;   4-preassure &nbsp; &nbsp; &nbsp; 5-fund_derivative &nbsp; &nbsp; &nbsp; 6-nonlinearityfactor

in scientific notation with 10 decimal places. Any other quantity of the problem can be reconstructed from these.

The use of a tabulated EoS expect a file with the following columns:
1-rho(g/cm3) &nbsp; &nbsp;2-Fundamental derivative &nbsp; &nbsp;3-adiabatic exponent &nbsp; &nbsp;4-Pressure (g/cm3) 
&nbsp; &nbsp;5-internal energy &nbsp; &nbsp;6-relativistic sound speed squared &nbsp; &nbsp;7-classic sound speed squared &nbsp; &nbsp;8-Grüneissen coefficient 
