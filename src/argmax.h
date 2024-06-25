#include <iostream>
#include <vector>
#include <seal/seal.h>
#include "ckks_evaluator.h"

class ArgmaxEvaluator
{
private:
    /* data */

public:
    CKKSEvaluator *ckks = nullptr;
    ArgmaxEvaluator(CKKSEvaluator &ckks) {
        this->ckks = &ckks;
    }
    void argmax(Ciphertext &x, Ciphertext &res, int len);
};