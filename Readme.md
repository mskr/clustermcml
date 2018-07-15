# clustermcml

Implementation of MCML in MPI and OpenCL.



## Build

Currently only Windows.
1. Get MSVC, MSMPI and OpenCL SDK.
2. Set paths in Makefile.
3. Run nmake.



## Run

mpiexec clustermcml.exe kernel.cl -Werror

"-Werror" is given to OpenCL compiler.
Without this option you won't see warnings because the program prints only errors (currently).
For multiple OpenCL compiler options separate them by spaces and wrap the whole string in "".
[A list of all options is found in the spec](https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/clBuildProgram.html#notes).



## Changing code

MPI and OpenCL boilerplate setup is done in "main.cpp".
Real work is done in the "run\*.cpp" files and "kernel.cl".

Currently there are two other kernel functions in kernel.cl besides "mcml" for testing purposes:
- "simpson" integrates the function "simpson_f" with the simpson method (integrate from 0 to 1 to approximate PI)
- "mcpi" approximates PI with a monte carlo method

Replace "runMCML" in the Makefile against "runSimpson" or "runMonteCarloPi" to run these kernels.
This works by simply linking different "getCLKernelName" and "runCLKernel" functions into the main program.
As you can see in the "run\*.cpp" files, the latter function sets the kernel arguments,
places it in a command queue, waits for it to finish and accumulates the results from all threads and processes.



## References

[Original MCML](https://omlc.org/software/mc/)

[CUDA MCML](http://www.atomic.physics.lu.se/biophotonics/research/monte-carlo-simulations/gpu-monte-carlo/)



## Code TODOs

- use all available devices: multiple command queues for multiple GPUs, multiple contexts when adding CPUs

- consider using -cl-mad-enable, native_log() and other kernel optimizations

- only compile kernel when timestamp newer than binary, otherwise use binary

- only broadcast kernel if needed (each node keeps a binary)


## Project TODOs

- build 64 bit

- port to MPICH (Linux)

- support more compilers