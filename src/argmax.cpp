#include "argmax.h"

using namespace std;
using namespace seal;

void ArgmaxEvaluator::argmax(Ciphertext &x, Ciphertext &res, int len) {
    Ciphertext tmp, a, b, sign, a_plus_b, a_minus_b;
    Plaintext one, zero_point_five;
    int log_step = log2(len);
    ckks->evaluator->rotate_vector(x, -len, *ckks->galois_keys, tmp);
    ckks->evaluator->add_inplace(x, tmp);
    a = x;
    for (int i = 0; i < log_step; ++i) {
        ckks->evaluator->rotate_vector(a, pow(2, i), *ckks->galois_keys, b);

        //a = max(a,b)
        ckks->evaluator->add(a, b, a_plus_b);
        ckks->evaluator->sub(a, b, a_minus_b);
        sign = ckks->sgn_eval2(a_minus_b, 2, 2);
        ckks->encoder->encode(0.5, a.parms_id(), a.scale(), zero_point_five);
        ckks->evaluator->multiply_plain_inplace(a_plus_b, zero_point_five);
        ckks->evaluator->rescale_to_next_inplace(a_plus_b);
        ckks->evaluator->mod_switch_to_inplace(a_minus_b, sign.parms_id());
        ckks->evaluator->multiply_inplace(a_minus_b, sign);
        ckks->evaluator->relinearize_inplace(a_minus_b, *ckks->relin_keys);
        ckks->evaluator->rescale_to_next_inplace(a_minus_b);
        a_plus_b.scale() = ckks->scale;
        a_minus_b.scale() = ckks->scale;
        ckks->evaluator->mod_switch_to_inplace(a_plus_b, a_minus_b.parms_id());
        ckks->evaluator->add(a_plus_b, a_minus_b, a);

        
        ckks->re_encrypt(a);
    }
    a.scale() = ckks->scale;
    ckks->evaluator->mod_switch_to_inplace(x, a.parms_id());
    ckks->evaluator->sub(x, a, res);
    res = ckks->sgn(res, 2, 2);
    ckks->encoder->encode(1.0, res.parms_id(),res.scale(), one);
    ckks->evaluator->add_plain_inplace(res, one);
    
}