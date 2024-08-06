# NEXUS
*NEXUS* is the first non-interactive protocol for secure transformer inference.

<br/>

## Installing SEAL
NEXUS uses a modified (bootstrapping-friendly) version of [Microsoft SEAL version 4.1](https://github.com/microsoft/SEAL/tree/4.1.2).

Please install the SEAL library `thirdparty/SEAL-4.1-bs` included as part of this repository (follow the instructions in the above link) before compiling NEXUS.

<br/>

## Compiling NEXUS
Once Microsoft SEAL 4.1 is installed, to build NEXUS simply run:

```bash
mkdir build
cd build
cmake ..
make
```

This should produce two binary executables inside the `build` directory:
- `bin/main`
- `bin/bootstrapping`

<br/>

## NEXUS-CUDA
We also provide a GPU accelerated implementation of NEXUS under the directory `cuda`.

For detailed instructions on how to compile and run NEXUS-CUDA, see [README.md](https://github.com/Kevin-Zh-CS/NEXUS/tree/main/cuda/README.md)
