# Cuda profile driver
if [ -f ./saxpy ]
then
  export CUDA_PROFILE=1
  ./saxpy 50000000
  export CUDA_PROFILE=0

  nvprof ./saxpy 50000000
  nvprof --print-gpu-trace ./saxpy 50000000
  nvprof -‐metrics ipc,flops_sp,inst_executed,l2_read_throughput,dram_read_throughput ./saxpy 50000000
fi

# Textual nvidia profiler
# nvprof EXEC
