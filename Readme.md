# clustermcml

Implementation of MCML in MPI and OpenCL.



## Build

Currently only Windows.
Get MSVC, MSMPI and OpenCL SDK.
Set paths in Makefile.
Run nmake.



## Run

mpiexec clustermcml.exe kernel.cl -Werror

"-Werror" is given to OpenCL compiler.
Without this option you won't see warnings because the program prints only errors (currently).



## Changing code

Currently there are two other kernel functions in kernel.cl besides "mcml" for testing purposes:
- "simpson" integrates the function "simpson_f" with the simpson method (integrate from 0 to 1 to approximate PI)
- "mcpi" approximates PI with a monte carlo method

Replace "runMCML" in the Makefile against "runSimpson" or "runMonteCarloPi" to run these kernels.
This works by simply linking different "runCLKernel" functions into the main program.



## References

[Original MCML](https://omlc.org/software/mc/)

[CUDA MCML](http://www.atomic.physics.lu.se/biophotonics/research/monte-carlo-simulations/gpu-monte-carlo/)



## Code TODOs

- consider using -cl-mad-enable


## Project TODOs

- build 64 bit

- port to MPICH (Linux)

- support more compilers