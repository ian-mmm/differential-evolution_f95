# differential-evolution_f95
Fortran 90/95 implementation of a basic differential evolution optimization algorithm.

Files:
DE_main.f95       -main Fortran file to test DE algorithm
DE_mod.f95        -module Fortran file, contains DE subroutine and Griewank Function as an example objective function

For the Griewank Function, this DE tends to lock up fairly quickly. I have found two successful strategies: (1) increase the grid size sunstantially to avoid grid lock, without ad hoc pertubations (e.g. NP=1040, T=1000 for nop=4 and -/+1000 starting grid); (2) turn on the ad hoc pertubations and increase the number of generations substantially (e.g. NP=50, T=10000 with pertubations for nop=4 and -/+1000 starting grid). The provided DE_main.f95 file has the later input values.

Griewank Function:
http://mathworld.wolfram.com/GriewankFunction.html
https://www.sfu.ca/~ssurjano/griewank.html

Tested on GCC 6.1, stable results for -O -O2 -O3
