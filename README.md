# NEXUS
NEXUS is the first non-interactive protocol for secure transformer inference.

<br/>

## Compiling SEAL
NEXUS depends on [Microsoft SEAL version 4.1](https://github.com/microsoft/SEAL/tree/4.1).
Download and install SEAL (follow the instructions in the above link) before compiling NEXUS.

<br/>

## Compiling NEXUS
Once Microsoft SEAL 4.1 is installed, to build NEXUS simply run:

```bash
cmake .
make
```

This should produce a binary file `bin/main`.
