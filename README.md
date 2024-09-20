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

<br/>

## Preprint
```tex
@misc{cryptoeprint:2024/136,
      author = {Jiawen Zhang and Xinpeng Yang and Lipeng He and Kejia Chen and Wen-jie Lu and Yinghao Wang and Xiaoyang Hou and Jian Liu and Kui Ren and Xiaohu Yang},
      title = {Secure Transformer Inference Made Non-interactive},
      howpublished = {Cryptology {ePrint} Archive, Paper 2024/136},
      year = {2024},
      url = {https://eprint.iacr.org/2024/136}
}
```
