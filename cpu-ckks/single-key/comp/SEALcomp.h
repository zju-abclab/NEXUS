#pragma	once

#include "seal/seal.h"
#include "SEALfunc.h"
//#include "func.h"
#include "MinicompFunc.h"

using namespace seal;
using namespace seal::util;

void minimax_ReLU_seal(long comp_no, vector<int> deg, long alpha, vector<Tree> &tree, double scaled_val, long scalingfactor, Encryptor &encryptor, Evaluator &evaluator, Decryptor &decryptor, CKKSEncoder &encoder, PublicKey &public_key, SecretKey &secret_key, RelinKeys &relin_keys, Ciphertext& cipher_x, Ciphertext& cipher_res);