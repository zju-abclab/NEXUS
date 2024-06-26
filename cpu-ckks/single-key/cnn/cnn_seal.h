#pragma	once
#include "seal/seal.h"
#include "SEALcomp.h"
#include "MinicompFunc.h"
#include "func.h"
#include "PolyUpdate.h"
#include "program.h"
#include "Bootstrapper.h"
#include <omp.h>
#include <NTL/RR.h>
#include <fstream>
#include <vector>
#include <chrono>


using namespace std;
using namespace seal;
using namespace minicomp;

class TensorCipher
{
private:
	int k_;		// k: gap
	int h_;		// w: height
	int w_;		// w: width
	int c_;		// c: number of channels
	int t_;		// t: \lfloor c/k^2 \rfloor
	int p_;		// p: 2^log2(nt/k^2hwt)
	int logn_;
	Ciphertext cipher_;

public:
	TensorCipher();
	TensorCipher(int logn, int k, int h, int w, int c, int t, int p, vector<double> data, Encryptor &encryptor, CKKSEncoder &encoder, int logp); 	// data vector contains hxwxc real numbers. 
	TensorCipher(int logn, int k, int h, int w, int c, int t, int p, Ciphertext cipher);
	int k() const;
    int h() const;
    int w() const;
	int c() const;
	int t() const;
	int p() const;
    int logn() const;
	Ciphertext cipher() const;
	void set_ciphertext(Ciphertext cipher);
	void print_parms();
};

void multiplexed_parallel_convolution_print(const TensorCipher &cnn_in, TensorCipher &cnn_out, int co, int st, int fh, int fw, const vector<double> &data, vector<double> running_var, vector<double> constant_weight, double epsilon, CKKSEncoder &encoder, Encryptor &encryptor, Evaluator &evaluator, GaloisKeys &gal_keys, vector<Ciphertext> &cipher_pool, ofstream &output, Decryptor &decryptor, SEALContext &context, size_t stage, bool end = false);
void multiplexed_parallel_batch_norm_seal_print(const TensorCipher &cnn_in, TensorCipher &cnn_out, vector<double> bias, vector<double> running_mean, vector<double> running_var, vector<double> weight, double epsilon, CKKSEncoder &encoder, Encryptor &encryptor, Evaluator &evaluator, double B, ofstream &output, Decryptor &decryptor, SEALContext &context, size_t stage, bool end = false);
void approx_ReLU_seal_print(const TensorCipher &cnn_in, TensorCipher &cnn_out, long comp_no, vector<int> deg, long alpha, vector<Tree> &tree, double scaled_val, long scalingfactor, Encryptor &encryptor, Evaluator &evaluator, Decryptor &decryptor, CKKSEncoder &encoder, PublicKey &public_key, SecretKey &secret_key, RelinKeys &relin_keys, double B, ofstream &output, SEALContext &context, GaloisKeys &gal_keys, size_t stage);
void bootstrap_print(const TensorCipher &cnn_in, TensorCipher &cnn_out, Bootstrapper &bootstrapper, ofstream &output, Decryptor &decryptor, CKKSEncoder &encoder, SEALContext &context, size_t stage);
void cipher_add_seal_print(const TensorCipher &cnn1, const TensorCipher &cnn2, TensorCipher &destination, Evaluator &evaluator, ofstream &output, Decryptor &decryptor, CKKSEncoder &encoder, SEALContext &context);
void multiplexed_parallel_downsampling_seal_print(const TensorCipher &cnn_in, TensorCipher &cnn_out, Evaluator &evaluator, Decryptor &decryptor, CKKSEncoder &encoder, SEALContext &context, GaloisKeys &gal_keys, ofstream &output);
void averagepooling_seal_scale_print(const TensorCipher &cnn_in, TensorCipher &cnn_out, Evaluator &evaluator, GaloisKeys &gal_keys, double B, ofstream &output, Decryptor &decryptor, CKKSEncoder &encoder, SEALContext &context);
void fully_connected_seal_print(const TensorCipher &cnn_in, TensorCipher &cnn_out, vector<double> matrix, vector<double> bias, int q, int r, Evaluator &evaluator, GaloisKeys &gal_keys, ofstream &output, Decryptor &decryptor, CKKSEncoder &encoder, SEALContext &context);

void multiplexed_parallel_convolution_seal(const TensorCipher &cnn_in, TensorCipher &cnn_out, int co, int st, int fh, int fw, const vector<double> &data, vector<double> running_var, vector<double> constant_weight, double epsilon, CKKSEncoder &encoder, Encryptor &encryptor, Evaluator &evaluator, GaloisKeys &gal_keys, vector<Ciphertext> &cipher_pool, bool end = false);
void multiplexed_parallel_batch_norm_seal(const TensorCipher &cnn_in, TensorCipher &cnn_out, vector<double> bias, vector<double> running_mean, vector<double> running_var, vector<double> weight, double epsilon, CKKSEncoder &encoder, Encryptor &encryptor, Evaluator &evaluator, double B, bool end = false);
void ReLU_seal(const TensorCipher &cnn_in, TensorCipher &cnn_out, long comp_no, vector<int> deg, long alpha, vector<Tree> &tree, double scaled_val, long scalingfactor, Encryptor &encryptor, Evaluator &evaluator, Decryptor &decryptor, CKKSEncoder &encoder, PublicKey &public_key, SecretKey &secret_key, RelinKeys &relin_keys, double scale = 1.0);
void cnn_add_seal(const TensorCipher &cnn1, const TensorCipher &cnn2, TensorCipher &destination, Evaluator &evaluator);
void multiplexed_parallel_downsampling_seal(const TensorCipher &cnn_in, TensorCipher &cnn_out, Evaluator &evaluator, GaloisKeys &gal_keys);
void averagepooling_seal_scale(const TensorCipher &cnn_in, TensorCipher &cnn_out, Evaluator &evaluator, GaloisKeys &gal_keys, double B, CKKSEncoder &encoder, Decryptor &decryptor, ofstream &output);
void matrix_multiplication_seal(const TensorCipher &cnn_in, TensorCipher &cnn_out, vector<double> matrix, vector<double> bias, int q, int r, Evaluator &evaluator, GaloisKeys &gal_keys);
void memory_save_rotate(const Ciphertext &cipher_in, Ciphertext &cipher_out, int steps, Evaluator &evaluator, GaloisKeys &gal_keys);

