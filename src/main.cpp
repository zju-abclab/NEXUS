#include <seal/ciphertext.h>
#include<seal/seal.h>
#include<iostream>
#include<vector>
#include<string>
#include<random>
#include<chrono>
#include<math.h>
//#include "gelu.h"
//#include "layer_norm.h"
#include "softmax.h"

using namespace std;
using namespace seal;
using namespace std::chrono;

int main()
{
    
    EncryptionParameters parms(scheme_type::ckks);
	long logN = 14;
	size_t poly_modulus_degree = 1 << logN;
	double scale = pow(2.0, 40);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, { 58, 40, 40, 40, 40, 40, 40, 40, 40, 58 }));
    SEALContext context(parms);

	KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context);
    Decryptor decryptor(context, secret_key);
    CKKSEncoder encoder(context);		
    RelinKeys relin_keys;
	keygen.create_relin_keys(relin_keys);
	GaloisKeys galois_keys;
    keygen.create_galois_keys(galois_keys);

    CKKSEvaluator ckks_evaluator(context, encryptor, decryptor, encoder, evaluator, scale, relin_keys, galois_keys);
    //GeLUEvaluator gelu_evaluator(ckks_evaluator);
    //LNEvaluator ln_evaluator(ckks_evaluator);
    SoftmaxEvaluator softmax_evaluator(ckks_evaluator);
    double bound = 1.0 / (1 << 16);
    vector<double> input = {-0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4};
    Plaintext plain_input;
    Ciphertext cipher_input;
    Ciphertext cipher_output;
    vector<double> output;
    ckks_evaluator.encoder->encode(input, scale, plain_input);
    ckks_evaluator.encryptor->encrypt(plain_input, cipher_input);

    // auto start = high_resolution_clock::now();
    // gelu_evaluator.gelu(cipher_input, cipher_output);
    // auto end = high_resolution_clock::now();
    // cout << poly_modulus_degree/2 << " times gelu() takes: " << duration_cast<milliseconds>(end - start).count() / 2.5 << " milliseconds" << endl;
    // auto start = high_resolution_clock::now();
    // int size = input.size();
    // ln_evaluator.layer_norm(cipher_input, cipher_output, size);
    // //ckks_evaluator.sgn_eval(cipher_input, 7, 3, 0.5);
    // auto end = high_resolution_clock::now();
    // cout << poly_modulus_degree/4 << " times LN() takes: " << duration_cast<milliseconds>(end - start).count() / 2.5 << " milliseconds" << endl;
    // ckks_evaluator.print_decrypted_ct(cipher_output, 8);
    auto start = high_resolution_clock::now();
    int size = input.size();
    softmax_evaluator.softmax2(cipher_input, cipher_output, size);
    auto end = high_resolution_clock::now();
    cout << poly_modulus_degree/4 << " times softmax() takes: " << duration_cast<milliseconds>(end - start).count() / 2.5 << " milliseconds" << endl;

    
    ckks_evaluator.print_decrypted_ct(cipher_output, 8);
    cout << "communication cost: " << ckks_evaluator.comm  << " bytes" << endl;
    cout << "communication round: " << ckks_evaluator.round  <<  endl;
}