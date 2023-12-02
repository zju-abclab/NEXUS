#include "ckks_evaluator.h"

void CKKSEvaluator::re_encrypt(Ciphertext &ct)
{
    Plaintext temp;
    vector<double> v;
    decryptor->decrypt(ct, temp);
    encoder->decode(temp, v);
    encoder->encode(v, scale, temp);
    encryptor->encrypt(temp, ct);
    // cout << "depth = " << context->get_context_data(ct.parms_id())->chain_index() << "\n";
}

void CKKSEvaluator::print_decrypted_ct(Ciphertext &ct, int nums)
{
    Plaintext temp;
    vector<double> v;
    decryptor->decrypt(ct, temp);
    encoder->decode(temp, v);
    for (int i = 0; i < nums; i++)
    {
        cout << v[i] << " ";
    }
    cout << "\n";
}

vector<double> CKKSEvaluator::init_vec_with_value(int N, double init_value)
{
    std::vector<double> v(N);

    for (int i = 0; i < N; ++i)
    {
        v[i] = init_value;
    }

    return v;
}

Ciphertext CKKSEvaluator::poly_eval(Ciphertext x, vector<Plaintext> coeff)
{
    // cout << "Initial depth " << context->get_context_data(x.parms_id())->chain_index() << "\n";
    // x^2
    Ciphertext x_2;
    evaluator->square(x, x_2);
    evaluator->relinearize_inplace(x_2, *relin_keys);
    evaluator->rescale_to_next_inplace(x_2);
    // cout << "x: " << x.scale() << " x_2: " << x_2.scale() << "\n";

    // x^3
    Ciphertext x_3;
    evaluator->mod_switch_to_inplace(x, x_2.parms_id());
    evaluator->multiply(x_2, x, x_3);
    evaluator->relinearize_inplace(x_3, *relin_keys);
    evaluator->rescale_to_next_inplace(x_3);
    // cout << "x_2: " << x_2.scale() << " x_3: " << x_3.scale() << "\n";

    // x^4
    Ciphertext x_4;
    evaluator->square(x_2, x_4);
    evaluator->relinearize_inplace(x_4, *relin_keys);
    evaluator->rescale_to_next_inplace(x_4);
    // cout << "x_4: " << x_4.scale() << " x_2: " << x_2.scale() << "\n";

    // x^5
    Ciphertext x_5;
    evaluator->mod_switch_to_inplace(x_2, x_3.parms_id());
    evaluator->multiply(x_2, x_3, x_5);
    evaluator->relinearize_inplace(x_5, *relin_keys);
    evaluator->rescale_to_next_inplace(x_5);
    // cout << "x_5: " << x_5.scale() << " x_3: " << x_3.scale() << "\n";

    // x^7
    Ciphertext x_7;
    evaluator->multiply(x_3, x_4, x_7);
    evaluator->relinearize_inplace(x_7, *relin_keys);
    evaluator->rescale_to_next_inplace(x_7);
    // cout << "x_7: " << x_7.scale() << "\n";

    // cout << depth(x_7) <<"\n";
    // Multiply constants
    evaluator->mod_switch_to_inplace(coeff[1], x.parms_id());
    evaluator->multiply_plain_inplace(x, coeff[1]);
    evaluator->rescale_to_next_inplace(x);

    evaluator->mod_switch_to_inplace(coeff[3], x_3.parms_id());
    evaluator->multiply_plain_inplace(x_3, coeff[3]);
    evaluator->rescale_to_next_inplace(x_3);

    evaluator->mod_switch_to_inplace(coeff[5], x_5.parms_id());
    evaluator->multiply_plain_inplace(x_5, coeff[5]);
    evaluator->rescale_to_next_inplace(x_5);

    evaluator->mod_switch_to_inplace(coeff[7], x_7.parms_id());
    evaluator->multiply_plain_inplace(x_7, coeff[7]);
    evaluator->rescale_to_next_inplace(x_7);

    // c7*x^7 + c5*x^5 + c3*x^3 + c1*x
    x.scale() = scale;
    x_3.scale() = scale;
    x_5.scale() = scale;
    x_7.scale() = scale;
    evaluator->mod_switch_to_inplace(x, x_7.parms_id());
    evaluator->mod_switch_to_inplace(x_3, x_7.parms_id());
    evaluator->mod_switch_to_inplace(x_5, x_7.parms_id());
    evaluator->add_inplace(x, x_3);
    evaluator->add_inplace(x, x_5);
    evaluator->add_inplace(x, x_7);

    // cout << "Final depth " << context->get_context_data(x.parms_id())->chain_index() << "\n";

    return x;
}

Ciphertext CKKSEvaluator::sgn_eval(Ciphertext x, int d_g, int d_f, double factor)
{
    vector<double> df_coeff = {0, 35.0 / 16, 0, -35.0 / 16, 0, 21.0 / 16, 0, -5.0 / 16};
    vector<double> dg_coeff = {0, 4589.0 / 1024, 0, -16577.0 / 1024, 0, 25614.0 / 1024, 0, -12860.0 / 1024};

    vector<double> df_coeff_last;
    vector<double> dg_coeff_last;
    for (size_t i = 0; i < df_coeff.size(); i++)
    {
        df_coeff_last.push_back(df_coeff[i] * factor);
    }
    for (size_t i = 0; i < dg_coeff.size(); i++)
    {
        dg_coeff_last.push_back(dg_coeff[i] * factor);
    }

    vector<Plaintext> encode_df(8), encode_dg(8), encode_df_last(8), encode_dg_last(8);
    for (int i = 0; i < 8; i++)
    {
        encoder->encode(df_coeff[i], scale, encode_df[i]);
        encoder->encode(dg_coeff[i], scale, encode_dg[i]);
        encoder->encode(df_coeff_last[i], scale, encode_df_last[i]);
        encoder->encode(dg_coeff_last[i], scale, encode_dg_last[i]);
    }
    for (int i = 0; i < d_g; i++)
    {
        if (context->get_context_data(x.parms_id())->chain_index() < 4)
        {
            re_encrypt(x);
        }
        if (i == d_g - 1)
        {
            x = poly_eval(x, encode_dg_last);
        }
        else
        {
            x = poly_eval(x, encode_dg);
        }
    }
    for (int i = 0; i < d_f; i++)
    {
        if (context->get_context_data(x.parms_id())->chain_index() < 4)
        {
            re_encrypt(x);
        }
        if (i == d_f - 1)
        {
            x = poly_eval(x, encode_df_last);
        }
        else
        {
            x = poly_eval(x, encode_df);
        }
    }
    re_encrypt(x);
    return x;
}

Ciphertext CKKSEvaluator::newtonIter(Ciphertext x, Ciphertext res, int iter)
{

    for (int i = 0; i < iter; i++)
    {
        if (context->get_context_data(res.parms_id())->chain_index() < 4)
            re_encrypt(res);
        // cout << i << " " << depth(res) << "\n";
        Plaintext three_half, neg_half;
        encoder->encode(1.5, scale, three_half);
        encoder->encode(-0.5, scale, neg_half);

        // x^2
        Ciphertext res_sq;
        evaluator->square(res, res_sq);
        evaluator->relinearize_inplace(res_sq, *relin_keys);
        evaluator->rescale_to_next_inplace(res_sq);
        // printVector(res_sq);
        // evaluator->negate_inplace(res_sq);
        // printVector(res_sq, 3);
        // cout << "square\n";

        //-0.5*x*b
        Ciphertext res_x;
        evaluator->mod_switch_to_inplace(neg_half, x.parms_id());
        evaluator->multiply_plain(x, neg_half, res_x);
        evaluator->rescale_to_next_inplace(res_x);
        if (context->get_context_data(res.parms_id())->chain_index() < context->get_context_data(res_x.parms_id())->chain_index())
            evaluator->mod_switch_to_inplace(res_x, res.parms_id());
        else
            evaluator->mod_switch_to_inplace(res, res_x.parms_id());

        evaluator->multiply_inplace(res_x, res);
        evaluator->relinearize_inplace(res_x, *relin_keys);
        evaluator->rescale_to_next_inplace(res_x);
        // cout << "negate\n";

        //-0.5*b*x^3
        evaluator->mod_switch_to_inplace(res_sq, res_x.parms_id());
        evaluator->multiply_inplace(res_x, res_sq);
        evaluator->relinearize_inplace(res_x, *relin_keys);
        evaluator->rescale_to_next_inplace(res_x);
        // cout << "res_x\n";
        // printVector(res_x, 3);
        // 1.5*x
        evaluator->mod_switch_to_inplace(three_half, res.parms_id());
        evaluator->multiply_plain_inplace(res, three_half);
        evaluator->rescale_to_next_inplace(res);
        // cout << "constant\n";

        //-0.5*b*x^3 + 1.5*x
        evaluator->mod_switch_to_inplace(res, res_x.parms_id());
        res_x.scale() = scale;
        res.scale() = scale;
        evaluator->add_inplace(res, res_x);
        // cout << "final\n";
    }
    re_encrypt(res);
    return res;
}

pair<Ciphertext, Ciphertext> CKKSEvaluator::goldSchmidtIter(Ciphertext v, Ciphertext y, int d)
{
    Ciphertext x, h, r, temp;
    Plaintext constant;
    encoder->encode(0.5, scale, constant);

    // GoldSchmidt's algorithm
    evaluator->mod_switch_to_inplace(y, v.parms_id());
    evaluator->multiply(v, y, x);
    evaluator->relinearize_inplace(x, *relin_keys);
    evaluator->rescale_to_next_inplace(x);
    evaluator->mod_switch_to_inplace(constant, y.parms_id());
    evaluator->multiply_plain(y, constant, h);
    evaluator->rescale_to_next_inplace(h);
    // cout << "gold\n";

    for (int i = 0; i < d; i++)
    {
        encoder->encode(0.5, scale, constant);
        // r = 0.5 - xh
        if (context->get_context_data(x.parms_id())->chain_index() < 3)
        {
            re_encrypt(x);
            re_encrypt(h);
        }
        evaluator->multiply(x, h, r);
        evaluator->relinearize_inplace(r, *relin_keys);
        evaluator->rescale_to_next_inplace(r);
        r.scale() = scale;
        evaluator->negate(r, temp);
        evaluator->mod_switch_to_inplace(constant, temp.parms_id());
        evaluator->add_plain(temp, constant, r);
        // cout << "r\n";

        // x = x + x*r
        evaluator->mod_switch_to_inplace(x, r.parms_id());
        evaluator->multiply(x, r, temp);
        evaluator->relinearize_inplace(temp, *relin_keys);
        evaluator->rescale_to_next_inplace(temp);
        x.scale() = scale;
        temp.scale() = scale;
        evaluator->mod_switch_to_inplace(x, temp.parms_id());
        evaluator->add_inplace(x, temp);
        // cout << "x\n";

        // h = h + h*r
        evaluator->mod_switch_to_inplace(h, r.parms_id());
        evaluator->multiply(h, r, temp);
        evaluator->relinearize_inplace(temp, *relin_keys);
        evaluator->rescale_to_next_inplace(temp);
        h.scale() = scale;
        temp.scale() = scale;
        evaluator->mod_switch_to_inplace(h, temp.parms_id());
        evaluator->add_inplace(h, temp);
        // cout << "h\n";
    }
    encoder->encode(2.0, scale, constant);
    evaluator->mod_switch_to_inplace(constant, h.parms_id());
    evaluator->multiply_plain_inplace(h, constant);
    evaluator->rescale_to_next_inplace(h);

    return make_pair(x, h);
}

Ciphertext CKKSEvaluator::invert_sqrt(Ciphertext x, int d_newt, int d_gold)
{
    Ciphertext res = initGuess(x);
    Ciphertext y = newtonIter(x, res, d_newt);
    pair<Ciphertext, Ciphertext> sqrt_inv_sqrt = goldSchmidtIter(x, y, d_gold);
    // printVector(sqrt_inv_sqrt.first, 1);
    // printVector(sqrt_inv_sqrt.second, 1);
    return sqrt_inv_sqrt.second;
}

uint64_t CKKSEvaluator::get_modulus(Ciphertext &x, int k)
{
    const vector<Modulus> &modulus = context->get_context_data(x.parms_id())->parms().coeff_modulus();
    int sz = modulus.size();
    return modulus[sz - k].value();
}

Ciphertext CKKSEvaluator::initGuess(Ciphertext x)
{
    Plaintext A, B;
    // a = 1e-3; b = 750
    // encoder->encode(-0.00019703, scale, A);
    // encoder->encode(0.14777278, scale, B);

    // a = 1e-4; b = 1000
    encoder->encode(-1.29054537e-04, scale, A);
    encoder->encode(1.29054537e-01, scale, B);

    return evalLine(x, A, B);
}

Ciphertext CKKSEvaluator::evalLine(Ciphertext x, Plaintext m, Plaintext c)
{
    // cout << "line\n";
    evaluator->mod_switch_to_inplace(m, x.parms_id());
    evaluator->multiply_plain_inplace(x, m);
    evaluator->rescale_to_next_inplace(x);
    evaluator->mod_switch_to_inplace(c, x.parms_id());
    x.scale() = scale;
    evaluator->add_plain_inplace(x, c);
    return x;
}

Ciphertext CKKSEvaluator::exp(Ciphertext x)
{
    Ciphertext b0, b1;
    Plaintext p0, p1, delta;
    vector<double> dest;

    encoder->encode(init_vec_with_value(N, -8.0), x.parms_id(), x.scale(), p0);
    encoder->encode(init_vec_with_value(N, 1.5), x.parms_id(), x.scale(), p1);
    encoder->encode(
        init_vec_with_value(N, 1.0 / 32), x.parms_id(), x.scale(), delta);

    evaluator->sub_plain(x, p0, b0);
    evaluator->multiply_plain_inplace(b0, delta);
    evaluator->rescale_to_next_inplace(b0);
    evaluator->sub_plain(x, p1, b1);
    evaluator->multiply_plain_inplace(b1, delta);
    evaluator->rescale_to_next_inplace(b1);

    b0 = sgn_eval(b0, 7, 3, 0.5);
    b1 = sgn_eval(b1, 7, 3, 0.5);

    Plaintext zero_point_five;
    encoder->encode(init_vec_with_value(N, 0.5), b1.parms_id(), b1.scale(), zero_point_five);
    Ciphertext a1, a2;

    evaluator->sub(b0, b1, a1);                    // a1 = b0 - b1
    evaluator->add_plain(b1, zero_point_five, a2); // a2 = b1 + 0.5

    double A[] = {0.999469891622, 0.998104199650, 0.501415542413, 0.169660297661, 0.042133244334, 0.007501312598, 0.000879175634, 0.000059258169, 0.000001716078};
    double B[] = {6.943979090878, -16.061172554433, 21.461821805218, -14.267410003218, 6.156317208726, -1.632082712072, 0.275766518989, -0.026342660111, 0.001204185268};

    vector<Plaintext> coeff_A(9), coeff_B(9);
    for (size_t i = 0; i < 9; i++)
    {
        encoder->encode(A[i], scale, coeff_A[i]);
        encoder->encode(B[i], scale, coeff_B[i]);
    }

    vector<Ciphertext> x_pow(9);
    x_pow[1] = x;
    evaluator->square(x, x_pow[2]);
    evaluator->relinearize_inplace(x_pow[2], *relin_keys);
    evaluator->rescale_to_next_inplace(x_pow[2]);

    evaluator->mod_switch_to_inplace(x, x_pow[2].parms_id());
    evaluator->multiply(x_pow[2], x, x_pow[3]);
    evaluator->relinearize_inplace(x_pow[3], *relin_keys);
    evaluator->rescale_to_next_inplace(x_pow[3]);

    evaluator->square(x_pow[2], x_pow[4]);
    evaluator->relinearize_inplace(x_pow[4], *relin_keys);
    evaluator->rescale_to_next_inplace(x_pow[4]);

    evaluator->mod_switch_to_inplace(x_pow[2], x_pow[3].parms_id());
    evaluator->multiply(x_pow[2], x_pow[3], x_pow[5]);
    evaluator->relinearize_inplace(x_pow[5], *relin_keys);
    evaluator->rescale_to_next_inplace(x_pow[5]);

    evaluator->square(x_pow[3], x_pow[6]);
    evaluator->relinearize_inplace(x_pow[6], *relin_keys);
    evaluator->rescale_to_next_inplace(x_pow[6]);

    evaluator->mod_switch_to_inplace(x_pow[3], x_pow[4].parms_id());
    evaluator->multiply(x_pow[3], x_pow[4], x_pow[7]);
    evaluator->relinearize_inplace(x_pow[7], *relin_keys);
    evaluator->rescale_to_next_inplace(x_pow[7]);

    evaluator->square(x_pow[4], x_pow[8]);
    evaluator->relinearize_inplace(x_pow[8], *relin_keys);
    evaluator->rescale_to_next_inplace(x_pow[8]);

    Ciphertext s1 = exp_poly_eval(x_pow, coeff_A);
    Ciphertext s2 = exp_poly_eval(x_pow, coeff_B);
    evaluator->mod_switch_to_inplace(a1, s1.parms_id());
    evaluator->multiply(a1, s1, s1);
    evaluator->relinearize_inplace(s1, *relin_keys);
    evaluator->rescale_to_next_inplace(s1);
    evaluator->mod_switch_to_inplace(a2, s2.parms_id());
    evaluator->multiply(a2, s2, s2);
    evaluator->relinearize_inplace(s2, *relin_keys);
    evaluator->rescale_to_next_inplace(s2);
    s1.scale() = scale;
    s2.scale() = scale;
    evaluator->mod_switch_to_inplace(s1, s2.parms_id());
    Ciphertext res;
    evaluator->add(s1, s2, res);
    return res;
}

Ciphertext CKKSEvaluator::exp_poly_eval(vector<Ciphertext> x_pow, vector<Plaintext> coeff) {
    vector<Ciphertext> coeff_x(coeff.size());
    for (size_t i = 1; i < coeff.size(); i++) {
        evaluator->mod_switch_to_inplace(coeff[i], x_pow[i].parms_id());
        evaluator->multiply_plain(x_pow[i], coeff[i], coeff_x[i]);
        evaluator->rescale_to_next_inplace(coeff_x[i]);
    }
    coeff_x[coeff.size() - 1].scale() = scale;
    for (size_t i = 1; i < coeff.size() - 1; i++) {
        coeff_x[i].scale() = scale;
        evaluator->mod_switch_to_inplace(coeff_x[i], coeff_x[coeff.size() - 1].parms_id());
        evaluator->add_inplace(coeff_x[coeff.size() - 1], coeff_x[i]);
    }
    evaluator->mod_switch_to_inplace(coeff[0], coeff_x[coeff.size() - 1].parms_id());
    evaluator->add_plain_inplace(coeff_x[coeff.size() - 1], coeff[0]);
    return coeff_x[coeff.size() - 1];
}