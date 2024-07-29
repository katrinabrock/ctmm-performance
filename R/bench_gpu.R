library(ctmm)
data("buffalo")
Cilla <- buffalo$Cilla
m.ouf <- ctmm(sigma=23 %#% "km^2",tau=c(6 %#% "day",1 %#% "hour"))
M.OUF <- ctmm.fit(Cilla,m.ouf)

# This script was  run with
# hyperfine \
# 'LD_PRELOAD=/usr/local/cuda-12/lib64/libnvblas.so NVBLAS_CONFIG_FILE=/etc/nvblas.conf Rscript bench_gpu.R' \
# 'Rscript bench_gpu.R' 


# Results:
 
# Summary
#   'Rscript bench_gpu.R' ran
#     1.18 ± 0.18 times faster than 'LD_PRELOAD=/usr/local/cuda-12/lib64/libnvblas.so NVBLAS_CONFIG_FILE=/etc/nvblas.conf Rscript bench_gpu.R'
    
#     Benchmark 1: LD_PRELOAD=/usr/local/cuda-12/lib64/libnvblas.so NVBLAS_CONFIG_FILE=/etc/nvblas.conf Rscript bench_gpu.R
#   Time (mean ± σ):     362.302 s ± 63.907 s    [User: 358.800 s, System: 11.398 s]
#   Range (min … max):   284.585 s … 462.310 s    10 runs
 
# Benchmark 2: Rscript bench_gpu.R
#   Time (mean ± σ):     383.261 s ± 64.066 s    [User: 381.764 s, System: 8.471 s]
#   Range (min … max):   283.340 s … 467.294 s    10 runs
 
# Summary
#   'LD_PRELOAD=/usr/local/cuda-12/lib64/libnvblas.so NVBLAS_CONFIG_FILE=/etc/nvblas.conf Rscript bench_gpu.R' ran
#     1.06 ± 0.26 times faster than 'Rscript bench_gpu.R'


# Conclusion: using nvblas did not help

