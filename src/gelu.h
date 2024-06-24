#include <iostream>
#include <vector>
#include <seal/seal.h>
#include "ckks_evaluator.h"

class GeLUEvaluator
{
private:
    /* data */

public:
    CKKSEvaluator *ckks = nullptr;
    GeLUEvaluator(CKKSEvaluator &ckks) {
        this->ckks = &ckks;
    }
    void gelu(Ciphertext &x, Ciphertext &res);
    vector<double> gelu_plain(vector<double>& input);
};