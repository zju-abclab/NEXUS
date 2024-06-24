#include "gelu.h"

using namespace std;
using namespace seal;


void GeLUEvaluator::gelu(Ciphertext &x, Ciphertext &res)
{
    Ciphertext b0, b1, b2;
    Plaintext p0, p1, p2, delta;
    vector<double> dest;

    ckks->encoder->encode(ckks->init_vec_with_value(ckks->slot_count, -4.0), x.parms_id(), x.scale(), p0);
    ckks->encoder->encode(ckks->init_vec_with_value(ckks->slot_count, -1.95), x.parms_id(), x.scale(), p1);
    ckks->encoder->encode(ckks->init_vec_with_value(ckks->slot_count, 3.0), x.parms_id(), x.scale(), p2);
    ckks->encoder->encode(
        ckks->init_vec_with_value(ckks->slot_count, 1.0 / 8), x.parms_id(), x.scale(), delta);

    ckks->evaluator->sub_plain(x, p0, b0);
    ckks->evaluator->multiply_plain_inplace(b0, delta);
    ckks->evaluator->rescale_to_next_inplace(b0);
    ckks->evaluator->sub_plain(x, p1, b1);
    ckks->evaluator->multiply_plain_inplace(b1, delta);
    ckks->evaluator->rescale_to_next_inplace(b1);
    ckks->evaluator->sub_plain(x, p2, b2);
    ckks->evaluator->multiply_plain_inplace(b2, delta);
    ckks->evaluator->rescale_to_next_inplace(b2);

    // time_start = high_resolution_clock::now();
    b0 = ckks->sgn_eval2(b0, 2, 2);
    b1 = ckks->sgn_eval2(b1, 2, 2);
    b2 = ckks->sgn_eval2(b2, 2, 2);

    
    Plaintext zero_point_five;
    ckks->encoder->encode(
        ckks->init_vec_with_value(ckks->slot_count, 0.5), b2.parms_id(), b2.scale(), zero_point_five);
    Ciphertext a0, a1, a2, a3;

    ckks->evaluator->sub(b0, b1, a1);                    // a1 = b0 - b1
    ckks->evaluator->sub(b1, b2, a2);                    // a2 = b1 - b2
    ckks->evaluator->add_plain(b2, zero_point_five, a3); // a3 = b2 + 0.5

    double A[] = {-0.5054031199708174, -0.42226581151983866, -0.11807612951181953, -0.011034134030615728};
    double B[] = {0.008526321541038084, 0.5, 0.3603292692789629, 0, -0.037688200365904236, 0, 0.0018067462606141187};

    Ciphertext x_2;
    ckks->evaluator->square(x, x_2);
    ckks->evaluator->relinearize_inplace(x_2, *ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(x_2);

    Ciphertext x_3;
    ckks->evaluator->mod_switch_to_inplace(x, x_2.parms_id());
    ckks->evaluator->multiply(x_2, x, x_3);
    ckks->evaluator->relinearize_inplace(x_3, *ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(x_3);

    Ciphertext x_4;
    ckks->evaluator->square(x_2, x_4);
    ckks->evaluator->relinearize_inplace(x_4, *ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(x_4);

    Ciphertext x_6;
    ckks->evaluator->square(x_3, x_6);
    ckks->evaluator->relinearize_inplace(x_6, *ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(x_6);

    vector<Plaintext> coeff_A(4), coeff_B(7);
    for (size_t i = 0; i < coeff_A.size(); i++)
    {
        ckks->encoder->encode(A[i], ckks->scale, coeff_A[i]);
    }
    for (size_t i = 0; i < coeff_B.size(); i++)
    {
        ckks->encoder->encode(B[i], ckks->scale, coeff_B[i]);
    }
    
    //Ax = A[0]+A[1]x+A[2]x^2+A[3]x^3
    Ciphertext Ax, a1x, a2x2, a3x3;
    ckks->evaluator->mod_switch_to_inplace(coeff_A[1], x.parms_id());
    ckks->evaluator->multiply_plain(x, coeff_A[1], a1x);
    ckks->evaluator->rescale_to_next_inplace(a1x);

    ckks->evaluator->mod_switch_to_inplace(coeff_A[2], x_2.parms_id());
    ckks->evaluator->multiply_plain(x_2, coeff_A[2], a2x2);
    ckks->evaluator->rescale_to_next_inplace(a2x2);

    ckks->evaluator->mod_switch_to_inplace(coeff_A[3], x_3.parms_id());
    ckks->evaluator->multiply_plain(x_3, coeff_A[3], a3x3);
    ckks->evaluator->rescale_to_next_inplace(a3x3);

    x.scale() = ckks->scale;    
    a1x.scale() = ckks->scale;
    a2x2.scale() = ckks->scale;
    a3x3.scale() = ckks->scale;
    
    ckks->evaluator->mod_switch_to_inplace(coeff_A[0], a3x3.parms_id());
    ckks->evaluator->mod_switch_to_inplace(a1x, a3x3.parms_id());
    ckks->evaluator->mod_switch_to_inplace(a2x2, a3x3.parms_id());
    ckks->evaluator->add_inplace(a3x3, a2x2);
    ckks->evaluator->add_inplace(a3x3, a1x);
    ckks->evaluator->add_plain(a3x3, coeff_A[0], Ax);

    //Bx = B[0]+B[1]x+B[2]x^2+B[4]x^4+B[6]x^6
    Ciphertext Bx, b1x, b2x2, b4x4, b6x6;
    ckks->evaluator->mod_switch_to_inplace(coeff_B[1], x.parms_id());
    ckks->evaluator->multiply_plain(x, coeff_B[1], b1x);
    ckks->evaluator->rescale_to_next_inplace(b1x);

    ckks->evaluator->mod_switch_to_inplace(coeff_B[2], x_2.parms_id());
    ckks->evaluator->multiply_plain(x_2, coeff_B[2], b2x2);
    ckks->evaluator->rescale_to_next_inplace(b2x2);

    ckks->evaluator->mod_switch_to_inplace(coeff_B[4], x_4.parms_id());
    ckks->evaluator->multiply_plain(x_4, coeff_B[4], b4x4);
    ckks->evaluator->rescale_to_next_inplace(b4x4);

    ckks->evaluator->mod_switch_to_inplace(coeff_B[6], x_6.parms_id());
    ckks->evaluator->multiply_plain(x_6, coeff_B[6], b6x6);
    ckks->evaluator->rescale_to_next_inplace(b6x6);
    
    x.scale() = ckks->scale;    
    b1x.scale() = ckks->scale;
    b2x2.scale() = ckks->scale;
    b4x4.scale() = ckks->scale;
    b6x6.scale() = ckks->scale;
    
    ckks->evaluator->mod_switch_to_inplace(coeff_B[0], b6x6.parms_id());
    ckks->evaluator->mod_switch_to_inplace(b1x, b6x6.parms_id());
    ckks->evaluator->mod_switch_to_inplace(b2x2, b6x6.parms_id());
    ckks->evaluator->mod_switch_to_inplace(b4x4, b6x6.parms_id());
    ckks->evaluator->mod_switch_to_inplace(b6x6, b6x6.parms_id());
    ckks->evaluator->add_inplace(b6x6, b4x4);
    ckks->evaluator->add_inplace(b6x6, b2x2);
    ckks->evaluator->add_inplace(b6x6, b1x);
    ckks->evaluator->add_plain(b6x6, coeff_B[0], Bx);

    Ciphertext s1, s2, s3;
    // cout << Ax.scale() << " " << Bx.scale() << " " << a1.scale() << endl;
    ckks->evaluator->mod_switch_to_inplace(Ax, a1.parms_id());
    ckks->evaluator->multiply(Ax, a1, s1);
    ckks->evaluator->relinearize_inplace(s1, *ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(s1);

    ckks->evaluator->mod_switch_to_inplace(Bx, a2.parms_id());
    ckks->evaluator->multiply(Bx, a2, s2);
    ckks->evaluator->relinearize_inplace(s2, *ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(s2);

    ckks->evaluator->mod_switch_to_inplace(x, a3.parms_id());
    ckks->evaluator->multiply(x, a3, s3);
    ckks->evaluator->relinearize_inplace(s3, *ckks->relin_keys);
    ckks->evaluator->rescale_to_next_inplace(s3);

    s1.scale() = ckks->scale;
    s2.scale() = ckks->scale;
    s3.scale() = ckks->scale;
    ckks->evaluator->mod_switch_to_inplace(s3, s2.parms_id());
    ckks->evaluator->mod_switch_to_inplace(s1, s2.parms_id());

    ckks->evaluator->add_inplace(s2, s1);
    ckks->evaluator->add(s3, s2, res);
    // ckks->print_decrypted_ct(s3, 2);
    // ckks->evaluator->add_inplace(s3, s2);

    // time_end = high_resolution_clock::now();
    // cout << N << " times GELU" << endl;
    // cout << "Computation cost:  " << duration_cast<milliseconds>(time_end - time_start).count() << " ms"
    //      << endl; // TODO: HEXL Acc
}


vector<double> GeLUEvaluator::gelu_plain(vector<double>& input) {
    vector<double> output;
    output.reserve(input.size());
    
    for (double x : input) {
        double gelu_x = 0.5 * x * (1.0 + std::tanh(std::sqrt(2.0 / M_PI) * (x + 0.044715 * x * x * x)));
        output.push_back(gelu_x);
    }

    return output;
}