#pragma once
#include <complex>

#include "kernels.cuh"
#include "phantom.h"

namespace nexus {
using namespace std;
using namespace phantom;
using namespace phantom::util;

class Encoder {
 private:
  PhantomContext *context;
  PhantomCKKSEncoder *encoder;

 public:
  Encoder() = default;

  Encoder(PhantomContext &context, PhantomCKKSEncoder &encoder) {
    this->context = &context;
    this->encoder = &encoder;
  }

  inline size_t message_length() { return encoder->message_length(); }

  // Vector (of doubles or complexes) inputs
  inline void encode(vector<double> values, size_t chain_index, double scale, PhantomPlaintext &plain) {
    if (values.size() == 1) {
      encode(values[0], chain_index, scale, plain);
      return;
    }
    values.resize(encoder->message_length(), 0.0);
    encoder->encode(*context, values, scale, plain, chain_index);
  }

  inline void encode(vector<double> values, double scale, PhantomPlaintext &plain) {
    if (values.size() == 1) {
      encode(values[0], scale, plain);
      return;
    }
    values.resize(encoder->message_length(), 0.0);
    encoder->encode(*context, values, scale, plain);
  }

  inline void encode(vector<complex<double>> complex_values, double scale, PhantomPlaintext &plain) {
    if (complex_values.size() == 1) {
      encode(complex_values[0], scale, plain);
      return;
    }
    complex_values.resize(encoder->message_length(), 0.0 + 0.0i);
    encoder->encode(*context, complex_values, scale, plain);
  }

  // Value inputs (fill all slots with that value)
  inline void encode(double value, size_t chain_index, double scale, PhantomPlaintext &plain) {
    vector<double> values(encoder->message_length(), value);
    encoder->encode(*context, values, scale, plain, chain_index);
  }

  inline void encode(double value, double scale, PhantomPlaintext &plain) {
    vector<double> values(encoder->message_length(), value);
    encoder->encode(*context, values, scale, plain);
  }

  inline void encode(complex<double> complex_value, double scale, PhantomPlaintext &plain) {
    vector<complex<double>> complex_values(encoder->message_length(), complex_value);
    encoder->encode(*context, complex_values, scale, plain);
  }

  inline void decode(PhantomPlaintext &plain, vector<double> &values) {
    encoder->decode(*context, plain, values);
  }
};

class Encryptor {
 private:
  PhantomContext *context;
  PhantomPublicKey *encryptor;

 public:
  Encryptor() = default;

  Encryptor(PhantomContext &context, PhantomPublicKey &encryptor) {
    this->context = &context;
    this->encryptor = &encryptor;
  }

  inline void encrypt(PhantomPlaintext &plain, PhantomCiphertext &ct) {
    encryptor->encrypt_asymmetric(*context, plain, ct);
  }
};

class Evaluator {
 private:
  PhantomContext *context;
  PhantomCKKSEncoder *encoder;

 public:
  Evaluator() = default;
  Evaluator(PhantomContext &context, PhantomCKKSEncoder &encoder) {
    this->context = &context;
    this->encoder = &encoder;
  }

  // Mod switch
  inline void mod_switch_to_next_inplace(PhantomCiphertext &ct) {
    ::mod_switch_to_next_inplace(*context, ct);
  }

  inline void mod_switch_to_inplace(PhantomCiphertext &ct, size_t chain_index) {
    ::mod_switch_to_inplace(*context, ct, chain_index);
  }

  inline void mod_switch_to_inplace(PhantomPlaintext &pt, size_t chain_index) {
    ::mod_switch_to_inplace(*context, pt, chain_index);
  }

  inline void rescale_to_next_inplace(PhantomCiphertext &ct) {
    ::rescale_to_next_inplace(*context, ct);
  }

  // Relinearization
  inline void relinearize_inplace(PhantomCiphertext &ct, const PhantomRelinKey &relin_keys) {
    ::relinearize_inplace(*context, ct, relin_keys);
  }

  // Multiplication
  inline void square(PhantomCiphertext &ct, PhantomCiphertext &dest) {
    multiply(ct, ct, dest);
  }

  inline void square_inplace(PhantomCiphertext &ct) {
    multiply_inplace(ct, ct);
  }

  inline void multiply(PhantomCiphertext &ct1, const PhantomCiphertext &ct2, PhantomCiphertext &dest) {
    if (&ct2 == &dest) {
      multiply_inplace(dest, ct1);
    } else {
      dest = ct1;
      multiply_inplace(dest, ct2);
    }
  }

  inline void multiply_inplace(PhantomCiphertext &ct1, const PhantomCiphertext &ct2) {
    ::multiply_inplace(*context, ct1, ct2);
  }

  inline void multiply_plain(PhantomCiphertext &ct, PhantomPlaintext &plain, PhantomCiphertext &dest) {
    dest = ct;
    multiply_plain_inplace(dest, plain);
  }

  inline void multiply_plain_inplace(PhantomCiphertext &ct, PhantomPlaintext &plain) {
    ::multiply_plain_inplace(*context, ct, plain);
  }

  // Addition
  inline void add_plain(PhantomCiphertext &ct, PhantomPlaintext &plain, PhantomCiphertext &dest) {
    dest = ct;
    add_plain_inplace(dest, plain);
  }

  inline void add_plain_inplace(PhantomCiphertext &ct, PhantomPlaintext &plain) {
    ::add_plain_inplace(*context, ct, plain);
  }

  inline void add(PhantomCiphertext &ct1, const PhantomCiphertext &ct2, PhantomCiphertext &dest) {
    if (&ct2 == &dest) {
      add_inplace(dest, ct1);
    } else {
      dest = ct1;
      add_inplace(dest, ct2);
    }
  }

  inline void add_inplace(PhantomCiphertext &ct1, const PhantomCiphertext &ct2) {
    ::add_inplace(*context, ct1, ct2);
  }

  inline void add_many(vector<PhantomCiphertext> &cts, PhantomCiphertext &dest) {
    ::add_many(*context, cts, dest);
  }

  // Subtraction
  inline void sub_plain(PhantomCiphertext &ct, PhantomPlaintext &plain, PhantomCiphertext &dest) {
    dest = ct;
    sub_plain_inplace(dest, plain);
  }

  inline void sub_plain_inplace(PhantomCiphertext &ct, PhantomPlaintext &plain) {
    ::sub_plain_inplace(*context, ct, plain);
  }

  inline void sub(PhantomCiphertext &ct1, const PhantomCiphertext &ct2, PhantomCiphertext &dest) {
    if (&ct2 == &dest) {
      sub_inplace(dest, ct1);
      negate_inplace(dest);
    } else {
      dest = ct1;
      sub_inplace(dest, ct2);
    }
  }

  inline void sub_inplace(PhantomCiphertext &ct1, const PhantomCiphertext &ct2) {
    ::sub_inplace(*context, ct1, ct2);
  }

  // Rotation
  inline void rotate_vector(PhantomCiphertext &ct, int steps, PhantomGaloisKey &galois_keys, PhantomCiphertext &dest) {
    dest = ct;
    rotate_vector_inplace(dest, steps, galois_keys);
  }

  inline void rotate_vector_inplace(PhantomCiphertext &ct, int steps, PhantomGaloisKey &galois_keys) {
    ::rotate_vector_inplace(*context, ct, steps, galois_keys);
  }

  // Negation
  inline void negate(PhantomCiphertext &ct, PhantomCiphertext &dest) {
    dest = ct;
    negate_inplace(dest);
  }

  inline void negate_inplace(PhantomCiphertext &ct) {
    ::negate_inplace(*context, ct);
  }

  // Galois
  inline void apply_galois(PhantomCiphertext &ct, int step, PhantomGaloisKey &galois_keys, PhantomCiphertext &dest) {
    dest = ct;
    apply_galois_inplace(dest, step, galois_keys);
  }

  inline void apply_galois_inplace(PhantomCiphertext &ct, int step, PhantomGaloisKey &galois_keys) {
    ::apply_galois_inplace(*context, ct, step, galois_keys);
  }

  // Matrix Multiplication
  inline void transform_from_ntt(const PhantomCiphertext &ct, PhantomCiphertext &dest) {
    dest = ct;
    transform_from_ntt_inplace(dest);
    cudaStreamSynchronize(phantom::util::global_variables::default_stream->get_stream());
  }

  inline void transform_from_ntt_inplace(PhantomCiphertext &ct) {
    auto &context_data = context->get_context_data(ct.chain_index());
    auto &parms = context_data.parms();
    auto &coeff_modulus = parms.coeff_modulus();
    const size_t coeff_modulus_size = coeff_modulus.size();
    auto &rns_tool = context_data.gpu_rns_tool();

    const auto &stream = phantom::util::global_variables::default_stream->get_stream();

    for (size_t i = 0; i < ct.size(); i++) {
      nwt_2d_radix8_backward_inplace(
          ct.data() + i * ct.poly_modulus_degree() * ct.coeff_modulus_size(),
          context->gpu_rns_tables(), coeff_modulus_size, 0, stream);
      rns_tool.base_Ql().compose_array(ct.data() + i * ct.poly_modulus_degree() * ct.coeff_modulus_size(),
                                       ct.data() + i * ct.poly_modulus_degree() * ct.coeff_modulus_size(),
                                       ct.poly_modulus_degree(), 1, stream);
    }

    ct.set_ntt_form(false);
    cudaStreamSynchronize(stream);
  }

  inline void transform_to_ntt(const PhantomCiphertext &ct, PhantomCiphertext &dest) {
    dest = ct;
    transform_to_ntt_inplace(dest);
    cudaStreamSynchronize(phantom::util::global_variables::default_stream->get_stream());
  }

  inline void transform_to_ntt_inplace(PhantomCiphertext &ct) {
    auto &context_data = context->get_context_data(ct.chain_index());
    auto &parms = context_data.parms();
    auto &coeff_modulus = parms.coeff_modulus();
    const size_t coeff_modulus_size = coeff_modulus.size();

    const auto &stream = phantom::util::global_variables::default_stream->get_stream();

    for (size_t i = 0; i < ct.size(); i++) {
      nwt_2d_radix8_forward_inplace(
          ct.data() + i * ct.poly_modulus_degree() * ct.coeff_modulus_size(),
          context->gpu_rns_tables(), coeff_modulus_size, 0, stream);
    }

    ct.set_ntt_form(true);
    cudaStreamSynchronize(stream);
  }

  inline void negacyclic_shift_poly_coeffmod(const std::uint64_t *poly, size_t coeff_count, size_t shift, const Modulus &modulus, std::uint64_t *result) {
    int block_size = blockDimGlb.x;
    int num_blocks = coeff_count / block_size;

    // cout << block_size << endl;

    // if (shift == 0) {
    //   set_uint(poly, coeff_count, result);
    //   return;
    // }

    // uint64_t index_raw = shift;
    // uint64_t coeff_count_mod_mask = static_cast<uint64_t>(coeff_count) - 1;
    // for (size_t i = 0; i < coeff_count; i++, poly++, index_raw++) {
    //   uint64_t index = index_raw & coeff_count_mod_mask;
    //   if (!(index_raw & static_cast<uint64_t>(coeff_count)) || !*poly) {
    //     *(result + index) = *poly;
    //   } else {
    //     *(result + index) = modulus.value() - *poly;
    //   }
    // }

    negacyclic_shift_poly_coeffmod_kernel<<<num_blocks, block_size>>>(poly, coeff_count, shift, modulus.value(), result);
    cudaStreamSynchronize(phantom::util::global_variables::default_stream->get_stream());
  }

  // {
  //   if (shift == 0) {
  //     set_uint(poly, coeff_count, result);
  //     return;
  //   }

  //   uint64_t index_raw = shift;
  //   uint64_t coeff_count_mod_mask = static_cast<uint64_t>(coeff_count) - 1;
  //   for (size_t i = 0; i < coeff_count; i++, poly++, index_raw++) {
  //     uint64_t index = index_raw & coeff_count_mod_mask;
  //     if (!(index_raw & static_cast<uint64_t>(coeff_count)) || !*poly) {
  //       result[index] = *poly;
  //     } else {
  //       result[index] = modulus.value() - *poly;
  //     }
  //   }
  // }

  // Bootstrapping
  inline void multiply_const(const PhantomCiphertext &ct, double value, PhantomCiphertext &dest) {
    dest = ct;
    multiply_const_inplace(dest, value);
  }

  inline void multiply_const_inplace(PhantomCiphertext &ct, double value) {
    PhantomPlaintext const_plain;

    vector<double> values(encoder->message_length(), value);
    encoder->encode(*context, values, ct.scale(), const_plain);
    mod_switch_to_inplace(const_plain, ct.params_id());
    multiply_plain_inplace(ct, const_plain);
  }

  inline void add_const(PhantomCiphertext &ct, double value, PhantomCiphertext &dest) {
    dest = ct;
    add_const_inplace(dest, value);
  }

  inline void add_const_inplace(PhantomCiphertext &ct, double value) {
    PhantomPlaintext const_plain;

    vector<double> values(encoder->message_length(), value);
    encoder->encode(*context, values, ct.scale(), const_plain);
    mod_switch_to_inplace(const_plain, ct.params_id());
    add_plain_inplace(ct, const_plain);
  }

  inline void add_reduced_error(const PhantomCiphertext &ct1, const PhantomCiphertext &ct2, PhantomCiphertext &dest) {
    if (&ct1 == &dest) {
      add_inplace_reduced_error(dest, ct1);
    } else {
      dest = ct1;
      add_inplace_reduced_error(dest, ct2);
    }
  }

  void add_inplace_reduced_error(PhantomCiphertext &ct1, const PhantomCiphertext &ct2);

  inline void sub_reduced_error(const PhantomCiphertext &ct1, const PhantomCiphertext &ct2, PhantomCiphertext &dest) {
    if (&ct1 == &dest) {
      sub_inplace_reduced_error(dest, ct1);
    } else {
      dest = ct1;
      sub_inplace_reduced_error(dest, ct2);
    }
  }

  void sub_inplace_reduced_error(PhantomCiphertext &ct1, const PhantomCiphertext &ct2);

  inline void multiply_reduced_error(const PhantomCiphertext &ct1, const PhantomCiphertext &ct2, const PhantomRelinKey &relin_keys, PhantomCiphertext &dest) {
    if (&ct2 == &dest) {
      multiply_inplace_reduced_error(dest, ct1, relin_keys);
    } else {
      dest = ct1;
      multiply_inplace_reduced_error(dest, ct2, relin_keys);
    }
  }

  void multiply_inplace_reduced_error(PhantomCiphertext &ct1, const PhantomCiphertext &ct2, const PhantomRelinKey &relin_keys);

  inline void double_inplace(PhantomCiphertext &ct) const {
    ::add_inplace(*context, ct, ct);
  }

  template <typename T, typename = std::enable_if_t<std::is_same<std::remove_cv_t<T>, double>::value || std::is_same<std::remove_cv_t<T>, std::complex<double>>::value>>
  inline void multiply_vector_reduced_error(PhantomCiphertext &ct, const std::vector<T> &value, PhantomCiphertext &dest) {
    dest = ct;
    multiply_vector_inplace_reduced_error(dest, value);
  }

  template <typename T, typename = std::enable_if_t<std::is_same<std::remove_cv_t<T>, double>::value || std::is_same<std::remove_cv_t<T>, std::complex<double>>::value>>
  inline void multiply_vector_inplace_reduced_error(PhantomCiphertext &ct, const std::vector<T> &value) {
    PhantomPlaintext plain;

    encoder->encode(*context, value, ct.scale(), plain);
    mod_switch_to_inplace(plain, ct.params_id());
    multiply_plain_inplace(ct, plain);
  }

  inline void complex_conjugate(PhantomCiphertext &ct, const PhantomGaloisKey &galois_keys, PhantomCiphertext &dest) {
    dest = ct;
    complex_conjugate_inplace(dest, galois_keys);
  }

  inline void complex_conjugate_inplace(PhantomCiphertext &ct, const PhantomGaloisKey &galois_keys) {
    ::complex_conjugate_inplace(*context, ct, galois_keys);
  }
};

class Decryptor {
 private:
  PhantomContext *context;
  PhantomSecretKey *decryptor;

 public:
  Decryptor() = default;
  Decryptor(PhantomContext &context, PhantomSecretKey &decryptor) {
    this->context = &context;
    this->decryptor = &decryptor;
  }

  inline void decrypt(PhantomCiphertext &ct, PhantomPlaintext &plain) {
    decryptor->decrypt(*context, ct, plain);
  }

  inline void create_galois_keys_from_steps(vector<int> &steps, PhantomGaloisKey &galois_keys) {
    galois_keys = decryptor->create_galois_keys_from_steps(*context, steps);
  }
};

class CKKSEvaluator {
 private:
  // Sign function coefficients
  double sgn_factor = 0.5;

  int g4_scale = (1 << 10);
  vector<double> g4_coeffs = {0, 5850, 0, -34974, 0, 97015, 0, -113492, 0, 46623};
  vector<double> g4_coeffs_last;

  int f4_scale = (1 << 7);
  vector<double> f4_coeffs = {0, 315, 0, -420, 0, 378, 0, -180, 0, 35};
  vector<double> f4_coeffs_last;

  // Helper functions
  uint64_t get_modulus(PhantomCiphertext &x, int k);

  PhantomCiphertext init_guess(PhantomCiphertext x);
  PhantomCiphertext eval_line(PhantomCiphertext x, PhantomPlaintext m, PhantomPlaintext c);

  // Evaluation functions
  PhantomCiphertext newton_iter(PhantomCiphertext x, PhantomCiphertext res, int iter);
  pair<PhantomCiphertext, PhantomCiphertext> goldschmidt_iter(PhantomCiphertext v, PhantomCiphertext y, int d = 1);
  void eval_odd_deg9_poly(vector<double> &a, PhantomCiphertext &x, PhantomCiphertext &dest);

 public:
  PhantomContext *context;
  Encryptor encryptor;
  Decryptor decryptor;
  Encoder encoder;
  Evaluator evaluator;
  PhantomRelinKey *relin_keys;
  PhantomGaloisKey *galois_keys;
  std::vector<std::uint32_t> rots;

  size_t degree;
  double scale;
  size_t slot_count;

  CKKSEvaluator(PhantomContext &context, PhantomPublicKey &encryptor, PhantomSecretKey &decryptor,
                PhantomCKKSEncoder &encoder, PhantomRelinKey &relin_keys, PhantomGaloisKey &galois_keys,
                double scale) {
    this->context = &context;
    this->relin_keys = &relin_keys;
    this->galois_keys = &galois_keys;

    this->scale = scale;
    this->slot_count = encoder.slot_count();
    this->degree = this->slot_count * 2;

    // Instantiate the component classes
    Encoder ckks_encoder(context, encoder);
    this->encoder = ckks_encoder;

    Encryptor ckks_encryptor(context, encryptor);
    this->encryptor = ckks_encryptor;

    Evaluator ckks_evaluator(context, encoder);
    this->evaluator = ckks_evaluator;

    Decryptor ckks_decryptor(context, decryptor);
    this->decryptor = ckks_decryptor;

    // Compute rotation steps
    for (int i = 0; i < uint(std::ceil(log2(degree))); i++) {
      rots.push_back((degree + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
    }

    // Compute sign function coefficients
    f4_coeffs_last.resize(10, 0);
    g4_coeffs_last.resize(10, 0);
    for (int i = 0; i <= 9; i++) {
      f4_coeffs[i] /= f4_scale;
      f4_coeffs_last[i] = f4_coeffs[i] * sgn_factor;

      g4_coeffs[i] /= g4_scale;
      g4_coeffs_last[i] = g4_coeffs[i] * sgn_factor;
    }
  }

  // Helper functions
  vector<double> init_vec_with_value(double value);
  PhantomPlaintext init_plain_power_of_x(size_t exponent);

  void re_encrypt(PhantomCiphertext &ct);
  void print_decrypted_ct(PhantomCiphertext &ct, int num);

  // Evaluation functions
  PhantomCiphertext sgn_eval(PhantomCiphertext x, int d_g, int d_f);
  PhantomCiphertext invert_sqrt(PhantomCiphertext x, int d_newt = 20, int d_gold = 1);
  PhantomCiphertext exp(PhantomCiphertext x);
  PhantomCiphertext inverse(PhantomCiphertext x, int iter = 4);

  // Metrics calcuation functions
  double calculate_MAE(vector<double> &y_true, PhantomCiphertext &ct, int N);
};

}  // namespace nexus
