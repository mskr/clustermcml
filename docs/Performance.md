Some things to do, when one wants to optimize performance:

1.  Print time spent on GPU resp. on CPU for each round (std::chrono::high_resolution_clock::time_point, std::chrono::duration<double, std::milli>)
2.  A first simple improvement: OpenMP (#pragma omp parallel for) around the finished photon checking loop
3.  Verify totalThreadCount: increment by one should nearly double GPU time
4. Turn off OpenCL command queue profiling mode
5. Look at results of sampling profiler for more insight [1][2]
6. Get total number of instructions, number of used registers and number of local/global/constant memory load/stores (there may be LLVM Tools for that and CL compilers sometimes based on LLVM)
7. Global memory reads/writes are likely to be bottlenecks. Alternatives are image memory and constant memory.
8. CUDAMCML accumulates weights in local variable before performing atomic add on global output array
9. CUDAMCML restarts finished photons immediately from the kernel, to avoid inactive threads, but at the cost of needing another atomic counter to keep track of total finished photons

[1] http://jamie-wong.com/post/speedscope/

[2] https://en.wikipedia.org/wiki/List_of_performance_analysis_tools