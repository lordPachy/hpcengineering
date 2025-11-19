# CUDA and Thrust
We can specify, in thrust, which things goes on the cpu and which on device. This is done through CUDA decorators.\

When using GPUs, we need to specify a "compute capability". This is specified through --arch_sm89 (for example) flag.\

We can compile once and obtain .cpu and .gpu files.