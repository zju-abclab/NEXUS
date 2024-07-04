#pragma once

#include "phantom.h"

using namespace std;
using namespace phantom;
using namespace phantom::arith;

__global__ void negacyclic_shift_poly_coeffmod_kernel(
    const uint64_t *d_poly, size_t coeff_count, size_t shift, uint64_t modulus_value, uint64_t *d_result);

// __global__ void expand_encode_kernel(const double *d_val, size_t poly_modulus_degree, Modulus *d_coeff_modulus, uint64_t *d_p);

// __global__ void negate_uint_mod_kernel(uint64_t operand, uint64_t modulus_value, uint64_t *result);

// __global__ void barrett_reduce_64_kernel(uint64_t input, Modulus *modulus, uint64_t *result);

// __global__ void multiply_uint64_hw64_kernel(uint64_t operand1, uint64_t operand2, uint64_t *hw64);

// __global__ void add_uint64_kernel(uint64_t operand1, uint64_t operand2, uint64_t *result) {
//   *result = operand1 + operand2;
// }
