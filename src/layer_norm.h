#include <iostream>
#include <vector>
#include <seal/seal.h>
#include "ckks_evaluator.h"

class LNEvaluator
{
private:
    /* data */

public:
    CKKSEvaluator *ckks = nullptr;
    LNEvaluator(CKKSEvaluator &ckks) {
        this->ckks = &ckks;
    }
    void layer_norm(Ciphertext &x, Ciphertext &res, int len);
};