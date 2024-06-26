# NEXUS-CUDA
CUDA implementation of NEXUS

<br/>

## TO-DO
- [ ] PhantomFHE bug: multiplications result in incorrect scale, reported: [Github Issue](https://github.com/encryptorion-lab/phantom-fhe/issues/5).
- [ ] Reimplement `mod_switch_to_next_inplace` (`mod_switch_drop_to_next_inplace`) and `rescale_to_next_inplace` (`mod_switch_scale_to_next_inplace`)
  - ref: [evaluate.cuh](https://github.com/encryptorion-lab/phantom-fhe/blob/1b4d892de7155285ac655bb7280256e3ae3c5e3b/include/evaluate.cuh), [evaluate.cu](https://github.com/encryptorion-lab/phantom-fhe/blob/1b4d892de7155285ac655bb7280256e3ae3c5e3b/src/evaluate.cu)
