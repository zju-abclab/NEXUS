#include "softmax.h"

using namespace std;
using namespace seal;


void SoftmaxEvaluator::softmax(Ciphertext &x, Ciphertext &res, int len) {
    Ciphertext tmp, exp_x;
    int log_step = log2(len);
    ckks->evaluator->rotate_vector(x, -len, *ckks->galois_keys, tmp);
    ckks->evaluator->add_inplace(x, tmp);
    exp_x = ckks->exp(x);
    tmp = exp_x;
    for (int i = 0; i < log_step; ++i) {
        ckks->evaluator->rotate_vector(tmp, pow(2, i), *ckks->galois_keys, res);
        ckks->evaluator->add_inplace(res, tmp);
        tmp = res;
    }
    ckks->evaluator->square_inplace(res);
    ckks->evaluator->relinearize_inplace(res, *ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(res);
    ckks->re_encrypt(res);
    res = ckks->invert_sqrt(res, 15, 5);
    ckks->evaluator->mod_switch_to_inplace(exp_x, res.parms_id());
    ckks->evaluator->multiply(res, exp_x, res);
    ckks->evaluator->relinearize_inplace(res, *ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(res);
}