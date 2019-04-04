This example simulate a 1D lahar over a variable slope topography. 
The rheology model adopted is the same of FLO-2D [1], where the total friction slope is written as the sum of the yield slope, the viscous slope and the turbulent dispersive slope components.

A Python script is provided to create the input file for this example. 
Please provide three arguments:

1) Number of cells

2) Variables to reconstruct: phys or cons

3) Order of the RK scheme

4) Initial solid volume fraction

Usage example of the script:

>> ./create_example5.py 400 phys 2 0.5

Once the input file (IMEX_SfloW2D.inp) is created create a simbolic link of the executable 
in this folder:

>> ln -s ../../../bin/IMEX_SfloW2D .

Finally, launch the solver:

>> ./IMEX_SfloW2D

A Python script to plot the results is provided. With this script you can choose the output and the variable to plot (h,hB,B,u,v)
Usage example:

>> ./plot_p_2d.py exampleLahar_0020.p_2d B hB

A Python script to create an animation (mp4) of the simulation is also provided. This script plot the topography and animate the flow over it. The interval between frames in milliseconds has to be given as input.
Usage example:

>> ./plot_animated.py exampleLahar 100


REFERENCES

[1] O'brien, J. S., P. Y. Julien, and W. T. Fullerton. "Two-dimensional water flood and mudflow simulation." Journal of hydraulic engineering 119.2 (1993): 244-261.
