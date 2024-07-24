#include <fstream>
#include <iostream>

#include "argmax.cuh"
#include "ckks_evaluator.cuh"
#include "gelu.cuh"
#include "layer_norm.cuh"
#include "matrix_mul.cuh"
#include "phantom.h"
#include "softmax.cuh"
#include "utils.cuh"

using namespace std;
using namespace phantom;
using namespace phantom::arith;
using namespace phantom::util;
using namespace troy;
using namespace troy::utils;
using namespace nexus;

// Choose test target here:
int TEST_TARGET_IDX = 2;

size_t N = 1ULL << 16;
size_t MM_LOG_N = 13;
size_t MM_N = 1ULL << MM_LOG_N;

double SCALE = pow(2.0, 40);

vector<string> TEST_TARGETS = {"MatMul", "MatMul_Phantom", "Argmax", "SoftMax", "LayerNorm", "GELU"};
vector<vector<int>> COEFF_MODULI =
    {
        {60, 40, 60},                                                                      // MatMul (0)
        {60, 40, 60},                                                                      // MatMul_Phantom (1)
        {17},                                                                              // Argmax (2) - Number of Moduli
        {58, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 58},          // SoftMax (3)
        {58, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 58},  // LayerNorm (4)
        {58, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 58}   // GELU (5)
};

string TEST_TARGET = TEST_TARGETS[TEST_TARGET_IDX];
vector<int> TEST_COEFF_MODULI = COEFF_MODULI[TEST_TARGET_IDX];

void MM_test() {
  // Convert vector<int> to vector<size_t>
  vector<size_t> MM_TEST_COEFF_MODULI(TEST_COEFF_MODULI.size());
  for (size_t i = 0; i < MM_TEST_COEFF_MODULI.size(); ++i) {
    MM_TEST_COEFF_MODULI[i] = static_cast<size_t>(TEST_COEFF_MODULI[i]);
  }

  troy::EncryptionParameters parms(SchemeType::CKKS);

  parms.set_poly_modulus_degree(MM_N);
  parms.set_coeff_modulus(troy::CoeffModulus::create(MM_N, MM_TEST_COEFF_MODULI));

  auto context = HeContext::create(parms, true, SecurityLevel::Nil);
  troy::CKKSEncoder encoder(context);

  context->to_device_inplace();
  encoder.to_device_inplace();

  KeyGenerator keygen(context);
  PublicKey public_key = keygen.create_public_key(false);

  troy::Encryptor encryptor(context);
  encryptor.set_public_key(public_key);
  encryptor.to_device_inplace();

  troy::Evaluator evaluator(context);
  troy::Decryptor decryptor(context, keygen.secret_key());

  std::vector<std::uint64_t> galois_elts;
  for (int i = 0; i < MM_LOG_N; i++) {
    galois_elts.push_back((MM_N + pow(2, i)) / pow(2, i));
  }
  GaloisKeys galois_keys = keygen.create_galois_keys_from_elements(galois_elts, false);
  galois_keys.to_device_inplace();

  CKKSEvaluator ckks_evaluator(context, &encryptor, &decryptor, &evaluator, &encoder, &galois_keys, SCALE, galois_elts);
  MMEvaluator mme(ckks_evaluator);

  std::vector<std::vector<double>> matrix_4096x768 = mme.read_matrix("../../data/input/matrixmul_input_m_128_n_768_k_64_batch_128.txt", 4096, 768);
  std::vector<std::vector<double>> matrix_768x64 = mme.read_matrix("../../data/input/matrix_input_n_768_k_64.txt", 768, 64);

  vector<Ciphertext> res;

  auto matrix_4096x768_T = mme.transpose_matrix(matrix_4096x768);
  auto matrix_768x64_T = mme.transpose_matrix(matrix_768x64);

  std::vector<std::vector<double>> row_pack;

  std::vector<double> row_ct(4096, 0.0);
  for (auto i = 0; i < 64 * 768; i++) {
    int row = i / 768;
    int col = i % 768;
    row_ct[i % 4096] = matrix_768x64_T[row][col];
    if (i % 4096 == 4095) {
      row_pack.push_back(row_ct);
    }
  }

  auto timer = Timer();

  mme.matrix_mul(matrix_4096x768_T, row_pack, res);

  timer.stop();
  cout << "[MatMul] 4096x768 x 768x64 takes: " << timer.duration<milliseconds>() << " milliseconds" << endl;

  std::vector<std::vector<double>> matrix_4096x64 = mme.read_matrix("../../data/calibration/matrix_output_m_128_k_64_batch_128.txt", 4096, 64);
  auto matrix_4096x64_T = mme.transpose_matrix(matrix_4096x64);

  // Calculate the error of the first col
  Plaintext res_pt;
  vector<complex<double>> mm_res;
  decryptor.decrypt(res[0], res_pt);
  encoder.decode_complex64_simd(res_pt, mm_res);

  double average_err = 0.0;
  for (auto i = 0; i < 4096; i++) {
    average_err += fabs(mm_res[i].real() / 2.0 - matrix_4096x64_T[0][i]);
    if (i < 10) printf("%+.10lf %+.10lf\n", mm_res[i].real() / 2.0, matrix_4096x64_T[0][i]);
  }
  std::cout << "Average error: " << average_err / 4096.0 << std::endl;

  MemoryPool::Destroy();
}

void MM_test_p() {
  phantom::EncryptionParameters parms(scheme_type::ckks);

  parms.set_poly_modulus_degree(MM_N);
  parms.set_coeff_modulus(phantom::arith::CoeffModulus::Create(MM_N, TEST_COEFF_MODULI));

  PhantomContext context(parms);
  PhantomCKKSEncoder encoder(context);

  PhantomSecretKey secret_key(context);
  std::ifstream sk_in("../../sk.txt");
  secret_key.load_secret_key(context, sk_in);

  PhantomPublicKey public_key = secret_key.gen_publickey(context);
  PhantomRelinKey relin_keys = secret_key.gen_relinkey(context);
  PhantomGaloisKey galois_keys;

  std::vector<std::uint32_t> galois_elts;
  for (int i = 0; i < MM_LOG_N; i++) {
    galois_elts.push_back((MM_N + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
  }

  CKKSEvaluator ckks_evaluator(&context, &public_key, &secret_key, &encoder, &relin_keys, &galois_keys, SCALE, galois_elts);
  MMEvaluator mme(ckks_evaluator);

  ckks_evaluator.decryptor.create_galois_keys_from_elts(galois_elts, *(ckks_evaluator.galois_keys));

  std::vector<std::vector<double>> matrix_4096x768 = mme.read_matrix("../../data/input/matrixmul_input_m_128_n_768_k_64_batch_128.txt", 4096, 768);
  std::vector<std::vector<double>> matrix_768x64 = mme.read_matrix("../../data/input/matrix_input_n_768_k_64.txt", 768, 64);

  vector<PhantomCiphertext> res;

  auto matrix_4096x768_T = mme.transpose_matrix(matrix_4096x768);
  auto matrix_768x64_T = mme.transpose_matrix(matrix_768x64);

  std::vector<std::vector<double>> row_pack;

  std::vector<double> row_ct(MM_N, 0.0);
  for (auto i = 0; i < 64 * 768; i++) {
    int row = i / 768;
    int col = i % 768;
    row_ct[i % MM_N] = matrix_768x64_T[row][col];
    if ((i % MM_N) == (MM_N - 1)) {
      row_pack.push_back(row_ct);
    }
  }

  auto timer = Timer();

  mme.matrix_mul(matrix_4096x768_T, row_pack, res);

  timer.stop();
  cout << "[MatMul] 4096x768 x 768x64 takes: " << timer.duration<milliseconds>() << " milliseconds" << endl;

  std::vector<std::vector<double>> matrix_4096x64 = mme.read_matrix("../../data/calibration/matrix_output_m_128_k_64_batch_128.txt", 4096, 64);
  auto matrix_4096x64_T = mme.transpose_matrix(matrix_4096x64);

  // err of the first col
  PhantomPlaintext res_pt;
  vector<double> mm_res;
  ckks_evaluator.decryptor.decrypt(res[0], res_pt);
  ckks_evaluator.encoder.decode(res_pt, mm_res);

  double average_err = 0.0;
  for (auto i = 0; i < 4096; i++) {
    average_err += fabs(mm_res[i] / 2.0 - matrix_4096x64_T[0][i]);
    if (i < 10) printf("%+.10lf <- %+.10lf\n", mm_res[i] / 2.0, matrix_4096x64_T[0][i]);
  }
  std::cout << "Average Error: " << average_err / 4096.0 << std::endl;
}

void argmax_test() {
  long logN = 15;
  long logn = logN - 2;
  long sparse_slots = (1 << logn);

  int logp = 46;
  int logq = 51;
  int log_special_prime = 51;

  // QuickMax: 17
  int main_mod_count = TEST_COEFF_MODULI[0];

  // Must be greater than 14: subsum 1 + coefftoslot 2 + ModReduction 9 + slottocoeff 2
  int bs_mod_count = 14;

  int secret_key_hamming_weight = 192;

  vector<int> coeff_bit_vec;

  coeff_bit_vec.push_back(logq);

  for (int i = 0; i < main_mod_count; i++) {
    coeff_bit_vec.push_back(logp);
  }

  for (int i = 0; i < bs_mod_count; i++) {
    coeff_bit_vec.push_back(logq);
  }

  coeff_bit_vec.push_back(log_special_prime);

  phantom::EncryptionParameters parms(scheme_type::ckks);
  size_t poly_modulus_degree = (size_t)(1 << logN);
  double scale = pow(2.0, logp);

  parms.set_poly_modulus_degree(poly_modulus_degree);
  parms.set_coeff_modulus(phantom::arith::CoeffModulus::Create(poly_modulus_degree, coeff_bit_vec));
  parms.set_secret_key_hamming_weight(secret_key_hamming_weight);
  parms.set_sparse_slots(sparse_slots);

  PhantomContext context(parms);

  PhantomSecretKey secret_key(context);
  PhantomPublicKey public_key = secret_key.gen_publickey(context);
  PhantomRelinKey relin_keys = secret_key.gen_relinkey(context);
  PhantomGaloisKey galois_keys = secret_key.create_galois_keys(context);

  PhantomCKKSEncoder encoder(context);

  CKKSEvaluator ckks_evaluator(&context, &public_key, &secret_key, &encoder, &relin_keys, &galois_keys, scale);
  ArgmaxEvaluator argmax_evaluator(ckks_evaluator, main_mod_count);

  size_t slot_count = encoder.slot_count();

  PhantomPlaintext plain_input;
  PhantomCiphertext cipher_input;
  PhantomCiphertext cipher_output;
  vector<double> input(slot_count, 0.0);

  double num;
  int argmax_input_size = 0;
  vector<double> argmax_input(sparse_slots, 0.0), argmax_calibration;

  ifstream input_file("../../data/input/argmax_input_8.txt");
  while (input_file >> num) {
    argmax_input[argmax_input_size] = num;
    argmax_input_size++;
  }
  input_file.close();

  ifstream calibration_file("../../data/calibration/argmax_calibration_8.txt");
  while (calibration_file >> num) {
    argmax_calibration.push_back(num);
  }
  calibration_file.close();

  // Spare input
  for (size_t i = 0; i < slot_count; i++) {
    input[i] = argmax_input[i % sparse_slots];
  }

  ckks_evaluator.encoder.encode(input, scale, plain_input);
  ckks_evaluator.encryptor.encrypt(plain_input, cipher_input);

  // Mod switch to remaining level
  for (int i = 0; i < bs_mod_count; i++) {
    ckks_evaluator.evaluator.mod_switch_to_next_inplace(cipher_input);
  }

  auto timer = Timer();
  argmax_evaluator.argmax(cipher_input, cipher_output, argmax_input_size);
  timer.stop();

  cout << "[Argmax] " << argmax_input_size << " takes: "
       << timer.duration<milliseconds>() << " milliseconds" << endl;
  ckks_evaluator.calculate_MAE(argmax_calibration, cipher_output, argmax_input_size);
}

int main() {
  /*
    MatMul
  */
  if (TEST_TARGET == "MatMul") {
    MM_test();
    return 0;
  }

  if (TEST_TARGET == "MatMul_Phantom") {
    MM_test_p();
    return 0;
  }

  if (TEST_TARGET == "Argmax") {
    argmax_test();
    return 0;
  }

  phantom::EncryptionParameters params(scheme_type::ckks);

  params.set_poly_modulus_degree(N);
  params.set_coeff_modulus(phantom::arith::CoeffModulus::Create(N, TEST_COEFF_MODULI));

  PhantomContext context(params);

  PhantomSecretKey secret_key(context);
  PhantomPublicKey public_key = secret_key.gen_publickey(context);
  PhantomRelinKey relin_keys = secret_key.gen_relinkey(context);
  PhantomGaloisKey galois_keys = secret_key.create_galois_keys(context);

  PhantomCKKSEncoder encoder(context);

  CKKSEvaluator ckks_evaluator(&context, &public_key, &secret_key, &encoder, &relin_keys, &galois_keys, SCALE);

  vector<double> input;
  PhantomPlaintext plain_input;
  PhantomCiphertext cipher_input;
  PhantomCiphertext cipher_output;
  PhantomPlaintext plain_output;
  vector<double> output;

  /*
    GELU
  */
  if (TEST_TARGET == "GELU") {
    GELUEvaluator gelu_evaluator(ckks_evaluator);

    double num;
    vector<double> gelu_calibration;
    ifstream input_file("../../data/input/gelu_input_32768.txt");
    while (input_file >> num) {
      input.push_back(num);
    }
    input_file.close();

    ifstream calibration_file("../../data/calibration/gelu_calibration_32768.txt");
    while (calibration_file >> num) {
      gelu_calibration.push_back(num);
    }
    calibration_file.close();

    ckks_evaluator.encoder.encode(input, SCALE, plain_input);
    ckks_evaluator.encryptor.encrypt(plain_input, cipher_input);

    auto timer = Timer();
    gelu_evaluator.gelu(cipher_input, cipher_output);
    timer.stop();

    cout << "[GELU] 32768 takes: " << timer.duration() << " milliseconds" << endl;
    cout << "Mean Absolute Error: " << ckks_evaluator.calculate_MAE(gelu_calibration, cipher_output, N / 2) << endl;
  }

  /*
    LayerNorm
  */
  if (TEST_TARGET == "LayerNorm") {
    LNEvaluator ln_evaluator(ckks_evaluator);

    double num;
    vector<double> input, layernorm_calibration;
    ifstream input_file("../../data/input/layernorm_input_16_768.txt");
    while (input_file >> num) {
      input.push_back(num);
    }
    input_file.close();

    ifstream calibration_file("../../data/calibration/layernorm_calibration_16_768.txt");
    while (calibration_file >> num) {
      layernorm_calibration.push_back(num);
    }
    calibration_file.close();

    ckks_evaluator.encoder.encode(input, SCALE, plain_input);
    ckks_evaluator.encryptor.encrypt(plain_input, cipher_input);

    auto timer = Timer();
    ln_evaluator.layer_norm(cipher_input, cipher_output, 1024);
    timer.stop();

    cout << "[LayerNorm] 16 x 768 takes: " << timer.duration() << " milliseconds" << endl;
    cout << "Mean Absolute Error: " << ckks_evaluator.calculate_MAE(layernorm_calibration, cipher_output, 768) << endl;
  }

  /*
    Softmax
  */
  if (TEST_TARGET == "SoftMax") {
    SoftmaxEvaluator softmax_evaluator(ckks_evaluator);

    double num;
    vector<double> input, softmax_calibration;
    ifstream input_file("../../data/input/softmax_input_128_128.txt");
    while (input_file >> num) {
      input.push_back(num);
    }
    input_file.close();

    ifstream calibration_file("../../data/calibration/softmax_calibration_128_128.txt");
    while (calibration_file >> num) {
      softmax_calibration.push_back(num);
    }
    calibration_file.close();

    ckks_evaluator.encoder.encode(input, SCALE, plain_input);
    ckks_evaluator.encryptor.encrypt(plain_input, cipher_input);

    auto timer = Timer();
    softmax_evaluator.softmax(cipher_input, cipher_output, 128);
    timer.stop();

    cout << "[Softmax] 128 x 128 takes: " << timer.duration() << " milliseconds" << endl;
    cout << "Mean Absolute Error: " << ckks_evaluator.calculate_MAE(softmax_calibration, cipher_output, 128) << endl;
  }
}
