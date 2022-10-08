RBAXI - 1976 Parabolized Navier-Stokes Solver
#############################################
;Author: Roie R. Black
:Email: roie.black@gmail.com

This project is recreating original code developed by the author while
conducting research in the Aerospace Research Laboratory at Wright-Patterson
AFB in Dayton, Ohio. The code solves the axisymmetric hypersonic flow over an
ogive-cylinder using the parabolized form of the governing equations. The
finite-difference solution uses MacCormack's explicit scheme. 

The code was originally run on a CDC 6600 mainframe, and was even run on an HP
desktop programmable calculator. Solutions took around 10 minutes of time on
the CDC 6600, and untold hours on the calculator. 

In 2003, this code was converted to Python as part of a course in
object-oriented design that was part of my Master's program in Computer Science
at Texas State - San Marcos. That code ran in about 10 seconds on my laptop!

This project returns to the original Fortran code, this time using **gfortran**
on a MacBook Pro laptop, The code is also being tested on a Windows 11 laptop.
