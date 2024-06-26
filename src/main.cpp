#include <chrono>
#include <fstream>
#include <iostream>
#include <math.h>
#include <random>
#include <seal/ciphertext.h>
#include <seal/plaintext.h>
#include <seal/seal.h>
#include <string>
#include <vector>
#include "gelu.h"
#include "layer_norm.h"
#include "matrix_mul.h"
#include "softmax.h"

using namespace std;
using namespace seal;
using namespace seal::util;
using namespace std::chrono;

void MM_test();

int main()
{
    EncryptionParameters parms(scheme_type::ckks);
    long logN = 16;
    size_t poly_modulus_degree = 1 << logN;
    double scale = pow(2.0, 40);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 58, 40, 40, 40, 40, 40, 40, 40, 40, 58 }));
    SEALContext context(parms, true, sec_level_type::none);

    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);

    Encryptor encryptor(context, public_key);
    CKKSEncoder encoder(context);
    Evaluator evaluator(context, encoder);
    Decryptor decryptor(context, secret_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    GaloisKeys galois_keys;

    // std::vector<std::uint32_t> rots;
    // for (int i = 0; i < 12; i++) {
    //     rots.push_back((poly_modulus_degree + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
    // }
    // keygen.create_galois_keys(rots, galois_keys);
    keygen.create_galois_keys(galois_keys);

    CKKSEvaluator ckks_evaluator(context, encryptor, decryptor, encoder, evaluator, scale, relin_keys, galois_keys);
    GeLUEvaluator gelu_evaluator(ckks_evaluator);
    LNEvaluator ln_evaluator(ckks_evaluator);
    SoftmaxEvaluator softmax_evaluator(ckks_evaluator);

    // vector<double> input = {-7, -6, -5, -4, -3, -2, -1, 0};
    vector<double> input;
    Plaintext plain_input;
    Ciphertext cipher_input;
    Ciphertext cipher_output;
    vector<double> output;

    /*
        GELU
    */
    double num;
    vector<double> gelu_calibration;
    ifstream input_file("data/input/gelu_input_32768.txt");
    while (input_file >> num) {
        input.push_back(num);
    }
    input_file.close();
    ifstream calibration_file("data/calibration/gelu_calibration_32768.txt");
    while (calibration_file >> num) {
        gelu_calibration.push_back(num);
    }
    calibration_file.close();
    ckks_evaluator.encoder->encode(input, scale, plain_input);
    ckks_evaluator.encryptor->encrypt(plain_input, cipher_input);
    auto start = high_resolution_clock::now();
    gelu_evaluator.gelu(cipher_input, cipher_output);
    auto end = high_resolution_clock::now();
    cout << poly_modulus_degree / 2 << " times gelu() takes: " << duration_cast<milliseconds>(end - start).count() / 1.0 << " milliseconds" << endl;

    // auto start = high_resolution_clock::now(); int size = input.size();
    // ln_evaluator.layer_norm(cipher_input, cipher_output, size);
    // auto end = high_resolution_clock::now();
    // cout << poly_modulus_degree/4 << " times LN() takes: " << duration_cast<milliseconds>(end - start).count() / 1.0
    // << " milliseconds" << endl;

    /*
        Softmax
    */
    // auto start = high_resolution_clock::now();
    // int size = input.size();
    // softmax_evaluator.softmax2(cipher_input, cipher_output, size);
    // auto end = high_resolution_clock::now();
    // cout << poly_modulus_degree / 4
    //      << " times softmax() takes: " << duration_cast<milliseconds>(end - start).count() / 1.0 << " milliseconds"
    //      << endl;
    // ckks_evaluator.print_decrypted_ct(cipher_output, 32768);

    cout << "Mean Absolute Error: " << ckks_evaluator.calculateMAE(gelu_calibration, cipher_output);
    // cout << "communication cost: " << ckks_evaluator.comm << " bytes" << endl;
    // cout << "communication round: " << ckks_evaluator.round << endl;
    // MM_test();
}

void MM_test()
{
    EncryptionParameters parms(scheme_type::ckks);
    long logN = 13;
    size_t poly_modulus_degree = 1 << logN;
    double scale = pow(2.0, 40);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 60, 40, 60 }));
    SEALContext context(parms, true, sec_level_type::none);

    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);

    Encryptor encryptor(context, public_key, secret_key);
    CKKSEncoder encoder(context);
    Evaluator evaluator(context, encoder);
    Decryptor decryptor(context, secret_key);

    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    GaloisKeys galois_keys;

    std::vector<std::uint32_t> rots;
    for (int i = 0; i < logN; i++) {
        rots.push_back((poly_modulus_degree + exponentiate_uint(2, i)) / exponentiate_uint(2, i));
    }
    keygen.create_galois_keys(rots, galois_keys);

    CKKSEvaluator ckks_evaluator(context, encryptor, decryptor, encoder, evaluator, scale, relin_keys, galois_keys);

    MMEvaluator mme(ckks_evaluator);

    vector<vector<double>> X(768);
    vector<vector<double>> Y(72, vector<double>(poly_modulus_degree, 0.0));
    for (auto i = 0; i < 768; i++) {
        vector<double> val(poly_modulus_degree / 2);
        for (auto j = 0; j < poly_modulus_degree / 2; j++) {
            val[j] = 10.0 * 2.0 * (1.0 * rand() / RAND_MAX - 0.5);
        }
        X[i] = val;
    }
    vector<Ciphertext> res;

    mme.matrix_mul(X, Y, res);

    // for (auto i = 0; i < 10; i++) {
    //     printf("%+.10lf\n", -3.5153774 * X[0][i]);
    // }

    Plaintext res_pt;
    vector<double> mm_res;
    ckks_evaluator.decryptor->decrypt(res[0], res_pt);
    ckks_evaluator.encoder->decode(res_pt, mm_res);
    for (auto i = 0; i < 10; i++) {
        printf("%+.10lf\n", mm_res[i]);
    }
}