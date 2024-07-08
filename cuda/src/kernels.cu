#include "kernels.cuh"

__global__ void negacyclic_shift_poly_coeffmod_kernel(
    const uint64_t *d_poly, size_t poly_degree, size_t shift, DModulus *modulus, size_t coeff_mod_size, uint64_t *d_result) {
  for (size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
       tid < poly_degree * coeff_mod_size;
       tid += blockDim.x * gridDim.x) {
    size_t twr = tid / poly_degree;
    DModulus mod = modulus[twr];

    if (tid < poly_degree) {
      uint64_t index_raw = shift + tid;
      uint64_t coeff_count_mod_mask = static_cast<uint64_t>(poly_degree) - 1;
      uint64_t index = index_raw & coeff_count_mod_mask;

      if (!(index_raw & static_cast<uint64_t>(poly_degree)) || !*(d_poly + tid)) {
        *(d_result + index) = *(d_poly + tid);
      } else {
        *(d_result + index) = mod.value() - *(d_poly + tid);
      }
    }
  }
}

__global__ void expand_encode_kernel(const double *d_val, size_t poly_modulus_degree, DModulus *modulus, uint64_t *d_p) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < poly_modulus_degree) {
    auto coeffd = round(d_val[i] * 10000000000);
    bool is_negative = signbit(coeffd);
    auto coeffu = static_cast<uint64_t>(fabs(coeffd));

    if (is_negative) {
      for (size_t j = 0; j < 2; j++) {
        d_p[i + (j * poly_modulus_degree)] = negate_uint64_mod(
            barrett_reduce_uint64_uint64(coeffu, modulus[j].value(), modulus[j].const_ratio()[1]), modulus[j].value());
      }
    } else {
      for (size_t j = 0; j < 2; j++) {
        d_p[i + (j * poly_modulus_degree)] = barrett_reduce_uint64_uint64(coeffu, modulus[j].value(), modulus[j].const_ratio()[1]);
      }
    }
  }
}
