#include "kernels.cuh"

__global__ void negacyclic_shift_poly_coeffmod_kernel(
    const uint64_t *d_poly, size_t coeff_count, size_t shift, uint64_t modulus_value, uint64_t *d_result) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < coeff_count) {
    uint64_t index_raw = shift + i;
    uint64_t coeff_count_mod_mask = static_cast<uint64_t>(coeff_count) - 1;
    uint64_t index = index_raw & coeff_count_mod_mask;

    if (!(index_raw & static_cast<uint64_t>(coeff_count)) || !d_poly[i]) {
      *(d_result + index) = *(d_poly + i);
    } else {
      *(d_result + index) = modulus_value - *(d_poly + i);
    }
  }
}

// __global__ void expand_encode_kernel(const double *d_val, size_t poly_modulus_degree, uint64_t modulus_value, uint64_t modulus_const_ratio, uint64_t *d_p) {
//   size_t i = blockIdx.x * blockDim.x + threadIdx.x;

//   if (i < poly_modulus_degree) {
//     auto coeffd = round(d_val[i] * 10000000000);
//     bool is_negative = signbit(coeffd);
//     auto coeffu = static_cast<uint64_t>(fabs(coeffd));

//     if (is_negative) {
//       for (size_t j = 0; j < 2; j++) {
//         uint64_t result;
//         barrett_reduce_64_kernel<<<1, 1>>>(coeffu, modulus_value, modulus_const_ratio, &result);
//         negate_uint_mod_kernel<<<1, 1>>>(result, modulus_value, d_p + (i + (j * poly_modulus_degree)));
//       }
//     } else {
//       for (size_t j = 0; j < 2; j++) {
//         barrett_reduce_64_kernel<<<1, 1>>>(coeffu, modulus_value, modulus_const_ratio, d_p + (i + (j * poly_modulus_degree)));
//       }
//     }
//   }
// }

// __global__ void barrett_reduce_64_kernel(uint64_t input, uint64_t modulus_value, uint64_t const_ratio, uint64_t *result) {
//   // Reduces input using base 2^64 Barrett reduction
//   // floor(2^64 / mod) == floor( floor(2^128 / mod) )
//   uint64_t tmp[2];
//   const std::uint64_t *const_ratio = modulus->const_ratio().data();
//   multiply_uint64_hw64_kernel(input, const_ratio[1], tmp + 1);

//   // Barrett subtraction
//   tmp[0] = input - tmp[1] * modulus->value();

//   // One more subtraction is enough
//   *result = tmp[0] >= modulus->value() ? tmp[0] - modulus->value() : tmp[0];
// }

// __global__ void negate_uint_mod_kernel(uint64_t operand, uint64_t modulus_value, uint64_t *result) {
//   int64_t non_zero = static_cast<int64_t>(operand != 0);
//   *result = (modulus_value - operand) & static_cast<uint64_t>(-non_zero);
// }

// __global__ void multiply_uint64_hw64_kernel(uint64_t operand1, uint64_t operand2, uint64_t *hw64) {
//   auto operand1_coeff_right = operand1 & 0x00000000FFFFFFFF;
//   auto operand2_coeff_right = operand2 & 0x00000000FFFFFFFF;
//   operand1 >>= 32;
//   operand2 >>= 32;

//   auto middle1 = operand1 * operand2_coeff_right;
//   uint64_t middle;
//   add_uint64_kernel(middle1, operand2 * operand1_coeff_right, &middle);
//   auto left = operand1 * operand2 +
//               (static_cast<uint64_t>(static_cast<unsigned char>(middle < middle1)) << 32);
//   auto right = operand1_coeff_right * operand2_coeff_right;
//   auto temp_sum = (right >> 32) + (middle & 0x00000000FFFFFFFF);

//   *hw64 = static_cast<uint64_t>(left + (middle >> 32) + (temp_sum >> 32));
// }
