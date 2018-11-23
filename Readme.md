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
- clustermcpi-windows.exe: very simple example program that approximates PI using a Monte Carlo method

## Run
<details>
<summary>mpiexec [hosts] clustermcml-windows.exe kernel.cl "-Werror -D DEBUG" sample.mci</summary>
<br>
[hosts] will be listing some network host addresses that MPI should use as computing nodes.
mpiexec will communicate with hosts by connecting to a service process (called smpd when using MSMPI).
Also make sure when using MSMPI that all Windows PCs have the same username and password and you are logged in,
otherwise the authentication with smpd fails.
<br><br>
"-Werror" is given to OpenCL compiler.
Without this option you won't see warnings because the program prints only errors (currently).
For multiple OpenCL compiler options separate them by spaces and wrap the whole string in "".
[A list of all options is found in the spec](https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/clBuildProgram.html#notes).
<br><br>
Pass the DEBUG define to the kernel to trigger the creation of a debug buffer.
There is for example an assert macro in kernel.cl that will print to the debug buffer.
The host code will print the debug buffer contents to the console if "error" string is found.
The debug buffer can also be visualized with a GL shader.
</p>
</details>


## Changing code

MPI and OpenCL boilerplate setup is done in "main.cpp" and "clusterlib.cpp".
Real work is done in the "run\*.cpp" and "\*Kernel.c" files.

Currently there are two other kernels for testing purposes:
- "simpson" integrates the function "simpson_f" with the simpson method (integrate from 0 to 1 to approximate PI)
- "mcpi" approximates PI with a monte carlo method

Here are the interfaces of the "run\*.cpp" files, which are called by main:
```c
/**
* Returns the name of the kernel function.
*/
void getCLKernelName();
```
```c
/**
* Allocates host memory and reports sizes for device buffers.
*/
void allocCLKernelResources(size_t totalThreadCount, char* kernelOptions, char* otherOptions,
int* inputBufferCount, size_t* inputBufferSizes,
int* outputBufferCount, size_t* outputBufferSizes, int maxBufferCount)
```
```c
/**
* Sets the kernel arguments,
* places it in a command queue, 
* waits for it to finish 
* and accumulates the results from all threads and processes.
*/
void runCLKernel(cl_context context, cl_command_queue cmdQueue, cl_kernel kernel, cl_mem* inputBuffers, cl_mem* outputBuffers,
size_t totalThreadCount, size_t simdThreadCount, int processCount, int rank);
```



## Code TODOs

- Use all available devices: multiple command queues for multiple GPUs, multiple contexts when adding CPUs

- Consider using -cl-mad-enable, native_log() and other kernel optimizations

- Only compile kernel when a) timestamp newer than binary or b) compiler options changed, otherwise use binary
    - Only broadcast kernel if needed (each node keeps a binary)

- Hash by Dave Hoskins: https://www.shadertoy.com/view/4djSRW

- Try quasirandom Monte Carlo (blue noise instead of true random)
    - http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/

- OpenCL 2.0: [work_group_reduce](https://www.khronos.org/registry/OpenCL/sdk/2.0/docs/man/xhtml/work_group_reduce.html)
instead of atomics may help to speed up binning

## Project Modules

### clusterlib

Provides convenience functions for managing a cluster setup using MPI and OpenCL.

### randomlib

Provides GPU-friendly random number generators and hashes.

### clmem

Provides macros for struct alignment, memory allocation and data layout.

Useful for [data](https://www.youtube.com/watch?v=rX0ItVEVjHc)-[orie](https://twitter.com/aras_p/status/1044656885100675072)[nted](https://twitter.com/BrookeHodgman/status/1049278775215570944) [des](https://www.youtube.com/watch?v=yy8jQgmhbAU)[ign](http://www.asawicki.info/news_1422_data-oriented_design_-_links_and_thoughts.html).

### cl2cpp

OpenCL C to C++ transpiler for enabling better debuggable "No GPU" builds of the program.

### gl-windows

Provides convenience methods for running an OpenGL shader with some uniform buffers on Windows.
The shader runs on a viewport-filling quad and can be edited "live" in the console.

Useful for visualizing simulation outputs.

## Project TODOs

- Header files for each module

- The "*run.cpp" modules should have the main method because only one can be used for execution

- Add rule to Makefile to bake kernel into cpp source to make the exe self-contained
  - so you can run "clustermcml.exe sample.mci"

- Build 64 bit

- Port to MPICH (Linux)

- Support more compilers



## References

[Original MCML](https://omlc.org/software/mc/)

[CUDA MCML](http://www.atomic.physics.lu.se/biophotonics/research/monte-carlo-simulations/gpu-monte-carlo/)