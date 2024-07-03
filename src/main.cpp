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

    MM_test();
    exit(0);


    EncryptionParameters parms(scheme_type::ckks);
    long logN = 16;
    size_t poly_modulus_degree = 1 << logN;
    double scale = pow(2.0, 40);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(
        poly_modulus_degree, { 58, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 58 }));

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

    keygen.create_galois_keys(galois_keys);

    CKKSEvaluator ckks_evaluator(context, encryptor, decryptor, encoder, evaluator, scale, relin_keys, galois_keys);
    GeLUEvaluator gelu_evaluator(ckks_evaluator);
    LNEvaluator ln_evaluator(ckks_evaluator);
    SoftmaxEvaluator softmax_evaluator(ckks_evaluator);

    Plaintext plain_input;
    Ciphertext cipher_input;
    Ciphertext cipher_output;
    vector<double> output;

    /*
        GELU
    */
    // double num;
    // vector<double> input, gelu_calibration;
    // ifstream input_file("data/input/gelu_input_32768.txt");
    // while (input_file >> num) {
    //     input.push_back(num);
    // }
    // input_file.close();
    // ifstream calibration_file("data/calibration/gelu_calibration_32768.txt");
    // while (calibration_file >> num) {
    //     gelu_calibration.push_back(num);
    // }
    // calibration_file.close();
    // ckks_evaluator.encoder->encode(input, scale, plain_input);
    // ckks_evaluator.encryptor->encrypt(plain_input, cipher_input);
    // auto start = high_resolution_clock::now();
    // gelu_evaluator.gelu(cipher_input, cipher_output);
    // auto end = high_resolution_clock::now();
    // cout << "[GELU] 32768 takes:" << duration_cast<milliseconds>(end -
    // start).count() << " milliseconds" << endl; cout << "Mean Absolute Error: "
    // << ckks_evaluator.calculateMAE(gelu_calibration, cipher_output,
    // poly_modulus_degree/2) << endl;

    /*
        LayerNorm
    */
    // double num;
    // vector<double> input, layernorm_calibration;
    // ifstream input_file("data/input/layernorm_input_16_768.txt");
    // while (input_file >> num) {
    //     input.push_back(num);
    // }
    // input_file.close();
    // ifstream
    // calibration_file("data/calibration/layernorm_calibration_16_768.txt");
    // while (calibration_file >> num) {
    //     layernorm_calibration.push_back(num);
    // }
    // calibration_file.close();
    // ckks_evaluator.encoder->encode(input, scale, plain_input);
    // ckks_evaluator.encryptor->encrypt(plain_input, cipher_input);
    // auto start = high_resolution_clock::now();
    // ln_evaluator.layer_norm(cipher_input, cipher_output, 1024);
    // auto end = high_resolution_clock::now();
    // cout << "[LayerNorm] 16 x 768 takes: " << duration_cast<milliseconds>(end -
    // start).count() << " milliseconds" << endl; cout << "Mean Absolute Error: "
    // << ckks_evaluator.calculateMAE(layernorm_calibration, cipher_output, 768)
    // << endl;

    /*
        Softmax
    */
    double num;
    vector<double> input, softmax_calibration;
    ifstream input_file("data/input/softmax_input_128_128.txt");
    while (input_file >> num) {
        input.push_back(num);
    }
    input_file.close();
    ifstream calibration_file("data/calibration/softmax_calibration_128_128.txt");
    while (calibration_file >> num) {
        softmax_calibration.push_back(num);
    }
    calibration_file.close();
    ckks_evaluator.encoder->encode(input, scale, plain_input);
    ckks_evaluator.encryptor->encrypt(plain_input, cipher_input);
    auto start = high_resolution_clock::now();
    softmax_evaluator.softmax(cipher_input, cipher_output, 128);
    auto end = high_resolution_clock::now();
    cout << "[Softmax] 128 x 128 takes: " << duration_cast<milliseconds>(end - start).count() << " milliseconds"
         << endl;
    cout << "Mean Absolute Error: " << ckks_evaluator.calculateMAE(softmax_calibration, cipher_output, 128) << endl;
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

    std::vector<std::vector<double>> matrix_4096x768 = mme.readMatrix("data/input/matrixmul_input_m_128_n_768_k_64_batch_128.txt", 4096, 768);

    std::vector<std::vector<double>> matrix_768x64 = mme.readMatrix("data/input/matrix_input_n_768_k_64.txt", 768, 64);

    vector<Ciphertext> res;

    auto matrix_4096x768_T = mme.transposeMatrix(matrix_4096x768);
    auto matrix_768x64_T = mme.transposeMatrix(matrix_768x64);

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
    mme.matrix_mul(matrix_4096x768_T, row_pack, res);

    std::vector<std::vector<double>> matrix_4096x64 = mme.readMatrix("data/calibration/matrix_output_m_128_k_64_batch_128.txt", 4096, 64);
    auto matrix_4096x64_T = mme.transposeMatrix(matrix_4096x64);

    double average_err = 0.0;

    // err of the first col
    for (auto col = 0; col < 64; col++) {
        Plaintext res_pt;
        vector<double> mm_res;
        ckks_evaluator.decryptor->decrypt(res[col], res_pt);
        ckks_evaluator.encoder->decode(res_pt, mm_res);
        for (auto i = 0; i < 4096; i++) {
            average_err += fabs(mm_res[i] / 2.0 - matrix_4096x64_T[col][i]);
        }
    }

    std::cout << "average_err: " << average_err / 4096 / 64 << std::endl;
}
