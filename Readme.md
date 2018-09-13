# clustermcml

Implementation of MCML in MPI and OpenCL.



## Build

Currently only Windows.
1. Get MSVC, MSMPI and OpenCL SDK.
2. Set paths in Makefile.
3. Run nmake.

Alternative make targets:
- clustermcml-windows-debug.exe: output windows debug symbols
- clustermcml-gl-windows.exe: launch a GL shader for output buffer visualization after the kernel has run



## Run

mpiexec [hosts] clustermcml-windows.exe kernel.cl "-Werror -D DEBUG" sample.mci

[hosts] will be listing some network host addresses that MPI should use as computing nodes.

"-Werror" is given to OpenCL compiler.
Without this option you won't see warnings because the program prints only errors (currently).
For multiple OpenCL compiler options separate them by spaces and wrap the whole string in "".
[A list of all options is found in the spec](https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/clBuildProgram.html#notes).

Pass the DEBUG define to the kernel to trigger the creation of a debug buffer.
There is for example an assert macro in kernel.cl that will print to the debug buffer.
The host code will print the debug buffer contents to the console if "error" string is found.
The debug buffer can also be visualized with a GL shader.



## Changing code

MPI and OpenCL boilerplate setup is done in "main.cpp".
Real work is done in the "run\*.cpp" files and "kernel.cl".

//TODO document interfaces

Currently there are two other kernel functions in kernel.cl besides "mcml" for testing purposes:
- "simpson" integrates the function "simpson_f" with the simpson method (integrate from 0 to 1 to approximate PI)
- "mcpi" approximates PI with a monte carlo method

Replace "runMCML" in the Makefile against "runSimpson" or "runMonteCarloPi" to run these kernels.
This works by simply linking different "getCLKernelName" and "runCLKernel" functions into the main program.
As you can see in the "run\*.cpp" files, the latter function sets the kernel arguments,
places it in a command queue, waits for it to finish and accumulates the results from all threads and processes.

//TODO The interfaces of runSimpson and runMonteCarloPi need to be updated, because it was changed for runMCML!

//TODO split kernel file



## Code TODOs

- use all available devices: multiple command queues for multiple GPUs, multiple contexts when adding CPUs

- consider using -cl-mad-enable, native_log() and other kernel optimizations

- only compile kernel when timestamp newer than binary, otherwise use binary

- only broadcast kernel if needed (each node keeps a binary)

- Hash by Dave Hoskins: https://www.shadertoy.com/view/4djSRW


## Project TODOs

- Turn main.cpp into module ("cluster.cpp")
  - nicer interface functions (fewer arguments, fewer functions, wrappers?)
  - header files for each module

- add rule to Makefile to bake kernel into cpp source to make the exe self-contained
  - so you can run "clustermcml.exe sample.mci"

- build 64 bit

- port to MPICH (Linux)

- support more compilers



## References

[Original MCML](https://omlc.org/software/mc/)

[CUDA MCML](http://www.atomic.physics.lu.se/biophotonics/research/monte-carlo-simulations/gpu-monte-carlo/)