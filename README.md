# SecureGPT
SecureGPT is an efficient secure transformer framework for chatGPT inference


# Compiling SEAL

SecureGPT depends on [Microsoft SEAL version 4.1](https://github.com/microsoft/SEAL/tree/4.1).

Download and install SEAL (follow the instructions in the above link) before before compiling SecureGPT.

# Compiling SecureGPT

Once Microsoft SEAL 4.1 is installed, to build SecureGPT simply run:

```
cmake .
make
```

This should produce a binary file ``bin/main``.