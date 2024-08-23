# NEXUS-CUDA
CUDA implementation of NEXUS, the first non-interactive protocol for secure transformer inference.

<br/>

> [!IMPORTANT]  
> This is a research project and is not intended for production use.

<br/>

## Prerequisites
- CUDA 11
- Tested on A100 GPUs

<br/>

## Compiling NEXUS-CUDA
NEXUS-CUDA uses a modified version of [PhantomFHE](https://github.com/encryptorion-lab/phantom-fhe/tree/5988c9c0a82ef86934c34a044e54032b94fd5a16), whose code is included in full under the directory `thirdparty`. There is no need to compile and install PhantomFHE as a dedicated library prior to compiling NEXUS-CUDA.

Use the following commands to build NEXUS-CUDA:

```bash
mkdir build
cd build
cmake ..
make
```

This should produce two binary executables inside the `build` directory:
  - `bin/main`
  - `bin/bootstrapping`
