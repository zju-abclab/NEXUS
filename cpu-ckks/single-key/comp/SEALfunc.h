#pragma	once

#include "seal/seal.h"
#include "func.h"
#include "MinicompFunc.h"
#include "PolyUpdate.h"
#include <memory>

using namespace seal;

namespace seal{
    void MultipleAdd_SEAL(Evaluator &evaluator, Ciphertext &cipher, Ciphertext& result, long long n);
	void test_evaluation(Evaluator &evaluator, CKKSEncoder &encoder, PublicKey &public_key, SecretKey &secret_key, RelinKeys &relin_keys, const Ciphertext& cipher_in, Ciphertext& cipher_out);
    void geneT0T1(Encryptor &encryptor, Evaluator &evaluator, CKKSEncoder &encoder, PublicKey &public_key, SecretKey &secret_key, RelinKeys &relin_keys, Ciphertext& T0, Ciphertext& T1, Ciphertext& cipher);
    void evalT(Evaluator &evaluator, PublicKey &public_key, SecretKey &secret_key, RelinKeys &relin_keys, Ciphertext& Tmplusn, const Ciphertext& Tm, const Ciphertext& Tn, const Ciphertext& Tmminusn);
	void eval_polynomial_integrate(Encryptor &encryptor, Evaluator &evaluator, Decryptor &decryptor, CKKSEncoder &encoder, PublicKey &public_key, SecretKey &secret_key, RelinKeys &relin_keys, Ciphertext& res, Ciphertext& cipher, long deg, const vector<double> &decomp_coeff, Tree &tree);
    long coeff_number(long deg, Tree& tree);
    void coeff_change(long comp_no, long deg[], double* coeff[], long type[], vector<Tree> &tree);
	long ShowFailure_ReLU(Decryptor &decryptor, CKKSEncoder &encoder, Ciphertext& cipher, vector<double>& x, long precision, long n);
}   