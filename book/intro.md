# RRBAXI

This project is a recreation of a project I worked on while attempting to complete a PhD degree in Aerospace Engineering during my first assignment as a research scientist at the USAF Aerospace Research Laboratiry at Wright-Patterson Air Force Base in Dayton, Ohio. That project was supposed to create a program that would solve 3-dimensional flow over a hypersonic vehicle similar to the NASA Shuttle. Sadly, I never quite got this finished, despite working on it for five years. (NASA managed to get the idea running, but it took them another eight years to do. My Air Force career forced me to abandon this pursuit.)

The code I did get running solver the flow over a simple ogive cylinder at zero angle of attack. The resulting flow is axisymmetric. The governing equations are the full compressible Navier-Stokes equations. The solution technique uses a spacial marching scheme, ignoring conventional time dependent terms. This scheme is called Parabolized Navier-Stokes (PNS). 

After movinginto my second assignment, teaching Computer Science at the Air Force Institute of Technology, the USAF Geaduate School, also at WPAFB, my code was used by several students pursuing Master's degrees in Aerospace Engineering.

The code was originally written in Fortran in 1972 and ran on a CDC 6600 mainframe. SOlutions normally used about 10 minutes of processing time. In 2003, as part of a second Master's degree program at Texas State, I re-wrote the code in Python and discovered solutions took only 10 seconds on my laptop. COmputer technology had seriously advanced in those 30 years! 

My intent in the current effort is to get the original code running once again, document the code more completely, and move it to a more modern version using parallel processing on current computer systems, including a cluster of Raspberry-Pi processors. 

This document was written using *Jupyter Book* - see [Jupyter Book](https://jupyterbook.org) for more information.

```{tableofcontents}
```
