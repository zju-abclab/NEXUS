#include <seal/ciphertext.h>
#include <vector>
#include "ckks_evaluator.h"
class MMEvaluator {
private:
    /* data */
    vector<Ciphertext> expand_ciphertext(
        const Ciphertext &encrypted, uint32_t m, GaloisKeys &galkey, vector<uint32_t> &galois_elts);

    void multiply_power_of_X(Ciphertext &encrypted, Ciphertext &destination, int index);

    void expandEncode(vector<double> &vec, Ciphertext &ct);

public:
    CKKSEvaluator *ckks = nullptr;
    MMEvaluator(CKKSEvaluator &ckks)
    {
        this->ckks = &ckks;
    }

    void matrix_mul(vector<vector<double>> &x, vector<vector<double>> &y, vector<Ciphertext> &res);
};