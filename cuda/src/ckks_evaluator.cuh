#pragma once
#include <complex>

#include "phantom.h"
#include "troy.cuh"

namespace nexus {
using namespace std;
using namespace phantom;
using namespace troy;

class Encoder {
 private:
  PhantomContext *context;
  PhantomCKKSEncoder *encoder;

 public:
  Encoder() = default;

  Encoder(PhantomContext *context, PhantomCKKSEncoder *encoder) {
    this->context = context;
    this->encoder = encoder;
  }

  inline size_t message_length() { return encoder->message_length(); }

  inline void reset_sparse_slots() { encoder->reset_sparse_slots(); }

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

  template <typename T, typename = std::enable_if_t<std::is_same<std::remove_cv_t<T>, double>::value || std::is_same<std::remove_cv_t<T>, std::complex<double>>::value>>
  inline void decode(PhantomPlaintext &plain, vector<T> &values) {
    encoder->decode(*context, plain, values);
  }
};

class Encryptor {
 private:
  PhantomContext *context;
  PhantomPublicKey *encryptor;

 public:
  Encryptor() = default;

  Encryptor(PhantomContext *context, PhantomPublicKey *encryptor) {
    this->context = context;
    this->encryptor = encryptor;
  }

  inline void encrypt(PhantomPlaintext &plain, PhantomCiphertext &ct) {
    encryptor->encrypt_asymmetric(*context, plain, ct);
  }

  inline void encrypt_zero(PhantomCiphertext &ct, size_t chain_index) {
    const phantom::util::cuda_stream_wrapper &stream_wrapper = *phantom::util::global_variables::default_stream;
    const auto &stream = stream_wrapper.get_stream();

    ct.set_correction_factor(1);
    ct.set_scale(1.0);

    encryptor->encrypt_zero_asymmetric_internal(*context, ct, chain_index, stream);
  }
};

class Evaluator {
 private:
  PhantomContext *context;
  PhantomCKKSEncoder *encoder;

 public:
  Evaluator() = default;
  Evaluator(PhantomContext *context, PhantomCKKSEncoder *encoder) {
    this->context = context;
    this->encoder = encoder;
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
  inline void add_plain(const PhantomCiphertext &ct, PhantomPlaintext &plain, PhantomCiphertext &dest) {
    dest = ct;
    add_plain_inplace(dest, plain);
  }

  inline void add_plain_inplace(PhantomCiphertext &ct, PhantomPlaintext &plain) {
    ::add_plain_inplace(*context, ct, plain);
  }

  inline void add(PhantomCiphertext &ct1, const PhantomCiphertext &ct2, PhantomCiphertext &dest) {
    dest = ::add(*context, ct1, ct2);
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
  inline void rotate_vector(const PhantomCiphertext &ct, int steps, PhantomGaloisKey &galois_keys, PhantomCiphertext &dest) {
    dest = ::rotate_vector(*context, ct, steps, galois_keys);
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
  inline void apply_galois(PhantomCiphertext &ct, uint32_t elt, PhantomGaloisKey &galois_keys, PhantomCiphertext &dest) {
    auto &galois_elts = context->key_galois_tool_->galois_elts();
    auto iter = find(galois_elts.begin(), galois_elts.end(), elt);
    auto galois_elt_index = iter - galois_elts.begin();

    dest = ::apply_galois(*context, ct, galois_elt_index, galois_keys);
  }

  inline void apply_galois_inplace(PhantomCiphertext &ct, int step, PhantomGaloisKey &galois_keys) {
    ::apply_galois_inplace(*context, ct, step, galois_keys);
  }

  // Complex Conjugate
  inline void complex_conjugate(PhantomCiphertext &ct, const PhantomGaloisKey &galois_keys, PhantomCiphertext &dest) {
    dest = ct;
    complex_conjugate_inplace(dest, galois_keys);
  }

  inline void complex_conjugate_inplace(PhantomCiphertext &ct, const PhantomGaloisKey &galois_keys) {
    ::complex_conjugate_inplace(*context, ct, galois_keys);
  }

  // Matrix Multiplication
  inline void transform_from_ntt(const PhantomCiphertext &ct, PhantomCiphertext &dest) {
    dest = ct;
    transform_from_ntt_inplace(dest);
  }

  inline void transform_from_ntt_inplace(PhantomCiphertext &ct) {
    auto rns_coeff_count = ct.poly_modulus_degree() * ct.coeff_modulus_size();

    const auto &stream = phantom::util::global_variables::default_stream->get_stream();

    for (size_t i = 0; i < ct.size(); i++) {
      uint64_t *ci = ct.data() + i * rns_coeff_count;
      nwt_2d_radix8_backward_inplace(ci, context->gpu_rns_tables(), ct.coeff_modulus_size(), 0, stream);
    }

    ct.set_ntt_form(false);
    cudaStreamSynchronize(stream);
  }

  inline void transform_to_ntt(const PhantomCiphertext &ct, PhantomCiphertext &dest) {
    dest = ct;
    transform_to_ntt_inplace(dest);
  }

  inline void transform_to_ntt_inplace(PhantomCiphertext &ct) {
    auto rns_coeff_count = ct.poly_modulus_degree() * ct.coeff_modulus_size();
    const auto &stream = phantom::util::global_variables::default_stream->get_stream();

    for (size_t i = 0; i < ct.size(); i++) {
      uint64_t *ci = ct.data() + i * rns_coeff_count;
      nwt_2d_radix8_forward_inplace(ci, context->gpu_rns_tables(), ct.coeff_modulus_size(), 0, stream);
    }

    ct.set_ntt_form(true);
    cudaStreamSynchronize(stream);
  }

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
  inline void multiply_vector_reduced_error(PhantomCiphertext &ct, std::vector<T> &values, PhantomCiphertext &dest) {
    dest = ct;
    multiply_vector_inplace_reduced_error(dest, values);
  }

  inline void multiply_vector_inplace_reduced_error(PhantomCiphertext &ct, vector<double> &values) {
    PhantomPlaintext plain;

    values.resize(encoder->message_length(), 0.0);
    encoder->encode(*context, values, ct.scale(), plain);
    mod_switch_to_inplace(plain, ct.params_id());
    multiply_plain_inplace(ct, plain);
  }

  inline void multiply_vector_inplace_reduced_error(PhantomCiphertext &ct, vector<complex<double>> &values) {
    PhantomPlaintext plain;

    values.resize(encoder->message_length(), 0.0 + 0.0i);
    encoder->encode(*context, values, ct.scale(), plain);
    mod_switch_to_inplace(plain, ct.params_id());
    multiply_plain_inplace(ct, plain);
  }
};

class Decryptor {
 private:
  PhantomContext *context;
  PhantomSecretKey *decryptor;

 public:
  Decryptor() = default;
  Decryptor(PhantomContext *context, PhantomSecretKey *decryptor) {
    this->context = context;
    this->decryptor = decryptor;
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
  // PhantomFHE
  PhantomContext *context;
  PhantomRelinKey *relin_keys;
  PhantomGaloisKey *galois_keys;
  std::vector<std::uint32_t> galois_elts;

  Encoder encoder;
  Encryptor encryptor;
  Evaluator evaluator;
  Decryptor decryptor;

  // TroyNova
  HeContextPointer troy_context;
  GaloisKeys *troy_galois_keys;
  std::vector<std::uint64_t> troy_galois_elts;

  troy::CKKSEncoder *troy_encoder;
  troy::Encryptor *troy_encryptor;
  troy::Evaluator *troy_evaluator;
  troy::Decryptor *troy_decryptor;

  size_t degree;
  double scale;
  size_t slot_count;

  CKKSEvaluator(PhantomContext *context, PhantomPublicKey *encryptor, PhantomSecretKey *decryptor,
                PhantomCKKSEncoder *encoder, PhantomRelinKey *relin_keys, PhantomGaloisKey *galois_keys,
                double scale, vector<uint32_t> galois_elts = {}) {
    this->context = context;
    this->relin_keys = relin_keys;
    this->galois_keys = galois_keys;
    this->galois_elts = galois_elts;

    this->scale = scale;
    this->slot_count = encoder->slot_count();
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

  CKKSEvaluator(HeContextPointer context, troy::Encryptor *encryptor, troy::Decryptor *decryptor,
                troy::Evaluator *evaluator, troy::CKKSEncoder *encoder, troy::GaloisKeys *galois_keys,
                double scale, vector<uint64_t> galois_elts = {}) {
    this->troy_context = context;
    this->troy_galois_keys = galois_keys;
    this->troy_galois_elts = galois_elts;

    this->scale = scale;
    this->slot_count = encoder->slot_count();
    this->degree = this->slot_count * 2;

    this->troy_encoder = encoder;
    this->troy_encryptor = encryptor;
    this->troy_evaluator = evaluator;
    this->troy_decryptor = decryptor;
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
