#include <iostream>
#include <vector>
#include <seal/seal.h>
#include "ckks_evaluator.h"

class SoftmaxEvaluator
{
private:
    /* data */

public:
    CKKSEvaluator *ckks = nullptr;
    SoftmaxEvaluator(CKKSEvaluator &ckks) {
        this->ckks = &ckks;
    }
    void softmax(Ciphertext &x, Ciphertext &res, int len);
};