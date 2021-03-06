# IMEX_SfloW2D-LAHARS

IMEX_SfloW2D-LAHARS is a FORTRAN90 code designed to model mudflows (lahars) over digital elevation models (DEMs) of natural terrain, with the total friction slope given by the sum of yield slope, viscous slope and turbulent slope (FLO-2D, O'Brien et al., 1993). A simple erosion and deposition model based on Fagents and Baloga, 2006 is implemented in the model. The rate of erosion is assumed to be proportional to the velocity, thickness and density (which depends on the solid sediment fraction) of the flow, through an emirical coefficient of proportionality (provided as input parameter). The deposition rate, on the other hand, is proportional to the fraction of sediments in the the flow, through the sediments settling velocity (provided as input parameter).

The system is described by an hyperbolic system of partial differential equations with relaxation and source terms. It is possible to select a simpler rheology in order to mimic the system of equations described in Kurganov and Petrova, 2007. IMEX_SfloW2D-LAHARS is written in Fortran90, which has the advantage to have a simple syntax and to be well known into the scientific community.

The model is discretized in time with an explicit-implicit Runge-Kutta method where the hyperbolic part is solved explicitly and the other terms (friction) are treated implicitly to allow for larger time steps and to enforce the stopping condition. The finite volume solver for the hyperbolic part of the system is based on the Kurganov and Petrova 2007 semi-discrete central scheme and it is not tied on the specific eigen-structure of the model. The implicit part is solved with a Newton-Raphson method where the elements of the Jacobian of the nonlinear system are evaluated numerically with a complex step derivative technique. This automatic procedure allows for easy changes of the friction term.

The code can deal with different scenarios, but its first aim is to treat lahars over topographies described as digital elevation models (DEMs) in the ESRI ascii format. The output files can be handled very well with gnuplot, in particular with the package it is provided a small script that create a video from the output data saved at different times. Moreover it is possible to save the solution as as ESRI ascii files, suitable for GIS softwares.

Several examples are provided with the code to test it with different scenarios (1D/2D, enabling or not erosion and deposition, ...).

### Authors and Contributors

Mattia de' Michieli Vitturi (@demichie)

Giacomo Lari

### Installation and execution

Check first if you have the LAPACK library installed on your system.

Download the IMEX_SfloW package and create the executable with the following commands from a terminal:

>./configure
>
>make
>
>make install

This will create the executable and copy it in the bin folder. You can test the executable copying it in the EXAMPLES folder and running it.

### Documentation

A wiki page describing the model is available at:

[https://github.com/demichie/LAHARS/wiki](https://github.com/demichie/LAHARS/wiki) 

Doxygen generated documentation of the code can be found at:

[http://demichie.github.io/LAHARS/html/](http://demichie.github.io/LAHARS/html/) 

### Acknowledgments

The development of IMEX-SfloW2D-LAHARS has been partially funded by Istituto Nazionale di Geofisica e Vulcanologia and Italian DPC (INGV - DPC B2).

### References


