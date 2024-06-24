#include <chrono>
#include <iostream>

#include "nexus_util.cuh"
#include "phantom.h"
#include "util.cuh"

using namespace std;
using namespace phantom;
using namespace phantom::arith;
using namespace phantom::util;
using namespace std::chrono;

size_t N = 1ULL << 16;
double SCALE = pow(2.0, 40);
double TOLERANCE = 1e-2;

int main() {
  EncryptionParameters params(scheme_type::ckks);

  params.set_poly_modulus_degree(N);
  params.set_coeff_modulus(CoeffModulus::Create(
      N,
      {60, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
       50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 53}  // 1763-bit
      ));

  PhantomContext context(params);
  print_parameters(context);

  PhantomSecretKey secret_key(context);
  PhantomPublicKey public_key = secret_key.gen_publickey(context);
  PhantomRelinKey relin_keys = secret_key.gen_relinkey(context);

  PhantomCKKSEncoder encoder(context);

  size_t slot_count = encoder.slot_count();
  // cout << "Number of slots: " << slot_count << endl;

  // Generate messages
  vector<cuDoubleComplex> x_msg, y_msg, result;
  double rand_real, rand_imag;

  size_t x_size = slot_count;
  size_t y_size = slot_count;
  x_msg.reserve(x_size);
  for (size_t i = 0; i < x_size; i++) {
    rand_real = (double)rand() / RAND_MAX;
    rand_imag = 0.0;
    x_msg.push_back(make_cuDoubleComplex(rand_real, rand_imag));
  }
  // cout << "Message vector 1: " << endl;
  // print_vector(x_msg, 3, 7);

  y_msg.reserve(y_size);
  for (size_t i = 0; i < y_size; i++) {
    rand_real = (double)rand() / RAND_MAX;
    rand_imag = 0.0;
    y_msg.push_back(make_cuDoubleComplex(rand_real, rand_imag));
  }
  // cout << "Message vector 2: " << endl;
  // print_vector(y_msg, 3, 7);

  // Encode & encrypt
  PhantomPlaintext x_plain;
  PhantomPlaintext y_plain;
  PhantomPlaintext xy_plain;

  auto start = high_resolution_clock::now();
  encoder.encode(context, x_msg, SCALE, x_plain);
  auto end = high_resolution_clock::now();
  cout << "Encoding took: " << duration_cast<microseconds>(end - start).count() << " μs" << endl;

  encoder.encode(context, y_msg, SCALE, y_plain);

  PhantomCiphertext x_cipher;
  PhantomCiphertext y_cipher;

  start = high_resolution_clock::now();
  public_key.encrypt_asymmetric(context, x_plain, x_cipher);
  end = high_resolution_clock::now();
  cout << "Encryption took: " << duration_cast<microseconds>(end - start).count() << " μs" << endl;

  public_key.encrypt_asymmetric(context, y_plain, y_cipher);

  // Compute ciphertext-ciphertext multiplication
  // cout << "Compute, relinearize, and rescale x*y." << endl;

  start = high_resolution_clock::now();
  multiply_inplace(context, x_cipher, y_cipher);
  end = high_resolution_clock::now();
  cout << "Multiplication took: " << duration_cast<microseconds>(end - start).count() << " μs" << endl;

  start = high_resolution_clock::now();
  relinearize_inplace(context, x_cipher, relin_keys);
  end = high_resolution_clock::now();
  cout << "Relinerization took: " << duration_cast<microseconds>(end - start).count() << " μs" << endl;

  start = high_resolution_clock::now();
  mod_switch_to_next_inplace(context, x_cipher);
  end = high_resolution_clock::now();
  cout << "Mod switch took: " << duration_cast<microseconds>(end - start).count() << " μs" << endl;

  // cout << "Scale of x*y after rescale: " << log2(x_cipher.scale()) << " bits" << endl;

  // Decrypt & decode
  start = high_resolution_clock::now();
  secret_key.decrypt(context, x_cipher, xy_plain);
  end = high_resolution_clock::now();
  cout << "Decryption took: " << duration_cast<microseconds>(end - start).count() << " μs" << endl;

  start = high_resolution_clock::now();
  encoder.decode(context, xy_plain, result);
  end = high_resolution_clock::now();
  cout << "Decoding took: " << duration_cast<microseconds>(end - start).count() << " μs" << endl;

  // cout << "Result vector: " << endl;
  // print_vector(result, 3, 7);

  // Verify correctness
  bool correctness = true;
  for (size_t i = 0; i < x_size; i++) {
    cuDoubleComplex res = result[i];
    cuDoubleComplex msg = cuCmul(x_msg[i], y_msg[i]);
    correctness &= near_equal(res, msg, TOLERANCE);
  }
  if (!correctness)
    throw std::logic_error("Homomorphic multiplication error");

  cout << "Homomorphic multiplication successful." << endl;

  result.clear();
  x_msg.clear();
  y_msg.clear();
}
