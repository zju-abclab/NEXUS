#include <fstream>
#include <iostream>

#include "ckks_evaluator.cuh"
#include "gelu.cuh"
#include "phantom.h"
#include "utils.cuh"

using namespace std;
using namespace phantom;
using namespace phantom::arith;
using namespace phantom::util;
using namespace nexus;

size_t N = 1ULL << 16;
double SCALE = pow(2.0, 40);

int main() {
  EncryptionParameters params(scheme_type::ckks);

  params.set_poly_modulus_degree(N);
  params.set_coeff_modulus(CoeffModulus::Create(
      N,
      {58, 40, 40, 40, 40, 40, 40, 40, 40, 58}
      // {60, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
      //  50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 53}  // 1763-bit
      ));
  params.set_special_modulus_size(1);

  PhantomContext context(params);

  PhantomSecretKey secret_key(context);
  PhantomPublicKey public_key = secret_key.gen_publickey(context);
  PhantomRelinKey relin_keys = secret_key.gen_relinkey(context);
  PhantomGaloisKey galois_keys = secret_key.create_galois_keys(context);

  PhantomCKKSEncoder encoder(context);

  CKKSEvaluator ckks_evaluator(context, public_key, secret_key, encoder, relin_keys, galois_keys, SCALE);
  GELUEvaluator gelu_evaluator(ckks_evaluator);

  vector<double> input;
  PhantomPlaintext plain_input;
  PhantomCiphertext cipher_input;
  PhantomCiphertext cipher_output;
  vector<double> output;

  /*
    GELU
  */
  double num;

  vector<double> gelu_calibration;
  ifstream input_file("../data/input/gelu_input_32768.txt");
  while (input_file >> num) {
    input.push_back(num);
  }
  input_file.close();

  ifstream calibration_file("../data/calibration/gelu_calibration_32768.txt");
  while (calibration_file >> num) {
    gelu_calibration.push_back(num);
  }
  calibration_file.close();

  ckks_evaluator.encoder.encode(input, SCALE, plain_input);
  ckks_evaluator.encryptor.encrypt(plain_input, cipher_input);

  auto timer = Timer();
  cout << "cipher_input scale: " << cipher_input.scale() << endl;
  cout << "depth = " << context.get_context_data_from_params_id(cipher_input.params_id()).chain_index() << endl;
  gelu_evaluator.gelu(cipher_input, cipher_output);
  timer.stop();
  cout << N / 2 << " times gelu() takes: " << timer.duration() << " milliseconds" << endl;

  cout << "Mean Absolute Error: " << ckks_evaluator.calculate_MAE(gelu_calibration, cipher_output);
}
