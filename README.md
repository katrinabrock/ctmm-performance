# CTMM Performance Improvement Project

These are some scripts I used in an attempt to profile, and improve runtime of Chris Flemming's [ctmm package](https://github.com/ctmm-initiative/ctmm). These scripts are not intended for re-use, but might be a helpful reference for similar projects. (Improving performance of an existing R statistical package.)

This work was conducted at and therefore funded by the Max Planck Institute of Animal Behavior.

# Lab Notebook

## TODO
- Look for vectorization opportunities in Kalman
- Go back to try to map ctmm kalman filter function to other kalman filter libraries (read the paper carefully). Validate the idea that this would cause issues with Inf and NaNs.
- Try nvblass on version with bigger matrices.

## 29062024

Sometime before this, I worked throught the vignettes in the package.

### nvblas attempt

I have a GPU machine, so I tried running CTMM with [nvblas](https://docs.nvidia.com/cuda/nvblas/index.html). I ran a small test script [see bench_gpu.R](R/bench_gpu.R) and benchmarked using [hyperfine](https://github.com/sharkdp/hyperfine). I was able to get R to use nvblas. I could see the process on `nvtop`. However, this made the script slower, not faster. Talking to Chris later, his suspicion was that this package is not acting on big enough matrices to make the GPU trade off worth it. So far, I haven't found a case (even outside of ctmm) where R+nvblas does better than openblas with my setup. Might try later to contrive some example to make this the case.

### profiling 

Profiled some of code from one of the ctmm vignettes using visprof.

- [Code here](R/profile_buf.R)
- [Results here](profiles/profile_all3.html). 

High level, Kalman filter portion is where all the CPU time is going. There are some small functions (notably sinch() and pd solvers) that can likely be optimized to save time.

I though it might be worth while to try popping in an existing kalman filter function. On my own I wasn't able to figure out how to map the code (and didn't find the right paper supplement). [Attempt here](R/kalman_mapping.R)

### dexp optimization attempt

Targetted the dexp1 and dexp2 functions for improvements.

- [Benchmarks Here](R/bench_dexp.R)
- [Resulting PR Here](https://github.com/ctmm-initiative/ctmm/pull/58#issue-2380386672) (with the outputs from above script printed)
- [Early attempts to user rcpp](R/cpp_dexp.R)