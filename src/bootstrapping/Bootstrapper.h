#pragma once

#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>

#include "ModularReducer.h"
// #include "ScaleInvEvaluator.h"

using namespace std;
using namespace seal;
using namespace seal::util;

class Bootstrapper {
 public:
  long loge;
  long logn;
  long n;
  long logNh;
  long Nh;
  long L;

  double initial_scale;
  double final_scale;

  long boundary_K;
  long sin_cos_deg;
  long scale_factor;
  long inverse_deg;

  SEALContext &context;
  KeyGenerator &keygen;
  CKKSEncoder &encoder;
  Encryptor &encryptor;
  Decryptor &decryptor;
  Evaluator &evaluator;
  RelinKeys &relin_keys;
  GaloisKeys &gal_keys;

  vector<long> slot_vec;
  long slot_index = 0;
  vector<vector<vector<vector<complex<double>>>>> orig_coeffvec, orig_invcoeffvec;
  vector<vector<vector<complex<double>>>> fftcoeff1, fftcoeff2, fftcoeff3;
  vector<vector<vector<complex<double>>>> invfftcoeff1, invfftcoeff2, invfftcoeff3;

  vector<Plaintext> fftcoeff_plain1, fftcoeff_plain2, invfftcoeff_plain1, invfftcoeff_plain2;
  vector<RNSIter> fftcoeff_iter1, fftcoeff_iter2, invfftcoeff_iter1, invfftcoeff_iter2;

  ModularReducer *mod_reducer;

  Bootstrapper(
      long _loge,
      long _logn,
      long _logNh,
      long _L,
      double _final_scale,
      long _boundary_K,
      long _sin_cos_deg,
      long _scale_factor,
      long _inverse_deg,
      SEALContext &_context,
      KeyGenerator &_keygen,
      CKKSEncoder &_encoder,
      Encryptor &_encryptor,
      Decryptor &_decryptor,
      Evaluator &_evaluator,
      RelinKeys &_relin_keys,
      GaloisKeys &_gal_keys);

  inline void set_final_scale(double _final_scale) {
    final_scale = _final_scale;
  }

  void print_decrypted_ct(Ciphertext &ct, int num);

  // Add rotation keys needed in bootstrapping (private function)
  void addLeftRotKeys_Linear_to_vector(vector<int> &gal_steps_vector);
  void addLeftRotKeys_Linear_to_vector_3(vector<int> &gal_steps_vector);

  void addLeftRotKeys_Linear_to_vector_other_slots(vector<int> &gal_steps_vector, long other_logn);
  void addLeftRotKeys_Linear_to_vector_3_other_slots(vector<int> &gal_steps_vector, long other_logn);

  void addLeftRotKeys_Linear_to_vector_one_depth(vector<int> &gal_steps_vector);
  void addLeftRotKeys_Linear_to_vector_one_depth_more_depth(vector<int> &gal_steps_vector);

  // Add rotation keys needed in bootstrapping (public function)
  void addBootKeys(GaloisKeys &gal_keys);
  void addBootKeys_3(GaloisKeys &gal_keys);
  void addBootKeys_other_keys(GaloisKeys &gal_keys, vector<int> &other_keys);
  void addBootKeys_3_other_keys(GaloisKeys &gal_keys, vector<int> &other_keys);

  void addBootKeys_other_slots(GaloisKeys &gal_keys, vector<long> &other_logn_vec);
  void addBootKeys_3_other_slots(GaloisKeys &gal_keys, vector<long> &other_logn_vec);
  void addBootKeys_3_other_slots_keys(GaloisKeys &gal_keys, vector<long> &other_logn_vec, vector<int> &other_keys);

  void addBootKeys_hoisting(GaloisKeys &gal_keys);
  void addBootKeys_one_depth(GaloisKeys &gal_keys);
  void addBootKeys_one_depth_more_depth(GaloisKeys &gal_keys);

  void change_logn(long new_logn);

  // Prepare the FFT coefficients
  void genorigcoeff();
  void merge_coeff(vector<vector<complex<double>>> merged_coeff, vector<vector<vector<complex<double>>>> orig_coeff);
  void rotated_merge_coeff(complex<double> **merged_coeff, complex<double> ***orig_coeff);
  void genfftcoeff_one_depth();
  void genfftcoeff_full_one_depth();
  void geninvfftcoeff_one_depth();
  void generate_LT_coefficient_one_depth();

  void genfftcoeff();
  void genfftcoeff_full();
  void geninvfftcoeff();
  void geninvfftcoeff_full();

  void genfftcoeff_3();
  void genfftcoeff_full_3();
  void geninvfftcoeff_3();
  void geninvfftcoeff_full_3();
  void generate_LT_coefficient();
  void generate_LT_coefficient_3();

  // Prepare the approximate polynomial
  void prepare_mod_polynomial();

  void subsum(double scale, Ciphertext &cipher);

  void bsgs_linear_transform(
      Ciphertext &rtncipher, Ciphertext &cipher, int totlen, int basicstep, int coeff_logn, const vector<vector<complex<double>>> &fftcoeff);
  void rotated_bsgs_linear_transform(
      Ciphertext &rtncipher, Ciphertext &cipher, int totlen, int basicstep, int coeff_logn, const vector<vector<complex<double>>> &fftcoeff);
  void rotated_nobsgs_linear_transform(Ciphertext &rtncipher, Ciphertext &cipher, int totlen, int coeff_logn, vector<vector<complex<double>>> fftcoeff);

  void bsgs_linear_transform_hoisting(
      Ciphertext &rtncipher, Ciphertext &cipher, int totlen, int basicstep, int coeff_logn, vector<vector<complex<double>>> fftcoeff);
  void rotated_bsgs_linear_transform_hoisting(
      Ciphertext &rtncipher, Ciphertext &cipher, int totlen, int basicstep, int coeff_logn, vector<vector<complex<double>>> fftcoeff);

  void sfl_one_depth(Ciphertext &rtncipher, Ciphertext &cipher);
  void sfl_full_one_depth(Ciphertext &rtncipher, Ciphertext &cipher);
  void sflinv_one_depth(Ciphertext &rtncipher, Ciphertext &cipher);
  void sflinv_one_depth_more_depth(Ciphertext &rtncipher, Ciphertext &cipher);

  void sfl(Ciphertext &rtncipher, Ciphertext &cipher);
  void sfl_full(Ciphertext &rtncipher, Ciphertext &cipher);
  void sflinv(Ciphertext &rtncipher, Ciphertext &cipher);
  void sflinv_full(Ciphertext &rtncipher, Ciphertext &cipher);

  void sfl_3(Ciphertext &rtncipher, Ciphertext &cipher);
  void sfl_full_3(Ciphertext &rtncipher, Ciphertext &cipher);
  void sfl_half_3(Ciphertext &rtncipher, Ciphertext &cipher);
  void sfl_full_half_3(Ciphertext &rtncipher, Ciphertext &cipher);
  void sflinv_3(Ciphertext &rtncipher, Ciphertext &cipher);
  void sflinv_full_3(Ciphertext &rtncipher, Ciphertext &cipher);

  void sfl_hoisting(Ciphertext &rtncipher, Ciphertext &cipher);
  void sfl_full_hoisting(Ciphertext &rtncipher, Ciphertext &cipher);
  void sflinv_hoisting(Ciphertext &rtncipher, Ciphertext &cipher);
  void sflinv_full_hoisting(Ciphertext &rtncipher, Ciphertext &cipher);

  // original bootstrapping
  void coefftoslot(Ciphertext &rtncipher, Ciphertext &cipher);
  void slottocoeff(Ciphertext &rtncipher, Ciphertext &cipher);

  void coefftoslot_full(Ciphertext &rtncipher1, Ciphertext &rtncipher2, Ciphertext &cipher);
  void slottocoeff_full(Ciphertext &rtncipher, Ciphertext &cipher1, Ciphertext &cipher2);

  // level-3 LT
  void coefftoslot_3(Ciphertext &rtncipher, Ciphertext &cipher);
  void slottocoeff_3(Ciphertext &rtncipher, Ciphertext &cipher);
  void slottocoeff_half_3(Ciphertext &rtncipher, Ciphertext &cipher);

  void coefftoslot_full_3(Ciphertext &rtncipher1, Ciphertext &rtncipher2, Ciphertext &cipher);
  void slottocoeff_full_3(Ciphertext &rtncipher, Ciphertext &cipher1, Ciphertext &cipher2);
  void slottocoeff_full_half_3(Ciphertext &rtncipher, Ciphertext &cipher1, Ciphertext &cipher2);
  // original bootstrapping hoisting version
  void coefftoslot_hoisting(Ciphertext &rtncipher, Ciphertext &cipher);
  void slottocoeff_hoisting(Ciphertext &rtncipher, Ciphertext &cipher);

  void coefftoslot_full_hoisting(Ciphertext &rtncipher1, Ciphertext &rtncipher2, Ciphertext &cipher);
  void slottocoeff_full_hoisting(Ciphertext &rtncipher, Ciphertext &cipher1, Ciphertext &cipher2);

  // one depth bootstrapping
  void coefftoslot_one_depth(Ciphertext &rtncipher, Ciphertext &cipher);
  void slottocoeff_one_depth(Ciphertext &rtncipher, Ciphertext &cipher);

  void coefftoslot_full_one_depth(Ciphertext &rtncipher1, Ciphertext &rtncipher2, Ciphertext &cipher);
  void slottocoeff_full_one_depth(Ciphertext &rtncipher, Ciphertext &cipher1, Ciphertext &cipher2);

  // mul_first bootstrapping
  void coefftoslot_full_mul_first(Ciphertext &rtncipher1, Ciphertext &rtncipher2, Ciphertext &cipher);
  void slottocoeff_full_mul_first(Ciphertext &rtncipher, Ciphertext &cipher1, Ciphertext &cipher2);
  void modraise_inplace(Ciphertext &cipher);

  // API bootstrapping
  void bootstrap_sparse(Ciphertext &rtncipher, Ciphertext &cipher);
  void bootstrap_full(Ciphertext &rtncipher, Ciphertext &cipher);

  void bootstrap_sparse_3(Ciphertext &rtncipher, Ciphertext &cipher);
  void bootstrap_full_3(Ciphertext &rtncipher, Ciphertext &cipher);

  void bootstrap_sparse_real_3(Ciphertext &rtncipher, Ciphertext &cipher);
  void bootstrap_full_real_3(Ciphertext &rtncipher, Ciphertext &cipher);

  void bootstrap_sparse_hoisting(Ciphertext &rtncipher, Ciphertext &cipher);
  void bootstrap_full_hoisting(Ciphertext &rtncipher, Ciphertext &cipher);
  void bootstrap_one_depth(Ciphertext &rtncipher, Ciphertext &cipher);
  void bootstrap_more_depth(Ciphertext &rtncipher, Ciphertext &cipher);
  void bootstrap(Ciphertext &rtncipher, Ciphertext &cipher);
  void bootstrap_inplace(Ciphertext &cipher);

  void bootstrap_3(Ciphertext &rtncipher, Ciphertext &cipher);
  void bootstrap_inplace_3(Ciphertext &cipher);

  void bootstrap_real_3(Ciphertext &rtncipher, Ciphertext &cipher);
  void bootstrap_inplace_real_3(Ciphertext &cipher);

  void bootstrap_hoisting(Ciphertext &rtncipher, Ciphertext &cipher);
  void bootstrap_inplace_hoisting(Ciphertext &cipher);
};
