# NEXUS-CUDA
CUDA implementation of NEXUS

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

This should produce a binary file `bin/main` under `build` directory.

<br>

## TO-DO
- Bootstrapping
  - `Bootstrapper::modraise_inplace`
- Matrix Multiplication
  - `transform_from_ntt_inplace` (see `moddown_from_NTT`)
  - `transform_to_ntt_inplace`
