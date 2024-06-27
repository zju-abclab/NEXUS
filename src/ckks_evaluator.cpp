#include "ckks_evaluator.h"
#include <seal/util/defines.h>
using namespace std::chrono;

void CKKSEvaluator::re_encrypt(Ciphertext &ct)
{
    auto start = high_resolution_clock::now();
    while (ct.coeff_modulus_size() > 1) {
        evaluator->mod_switch_to_next_inplace(ct);
    }
    vector<seal_byte> data;
    data.resize(ct.save_size(compr_mode_type::zstd));
    comm += ct.save(data.data(), data.size(), compr_mode_type::zstd);
    round++;
    // cout << "Communication cost:  " <<
    // ct.save(data.data(),data.size(),compr_mode_type::zstd) << " bytes" << endl;
    Plaintext temp;
    vector<double> v;
    decryptor->decrypt(ct, temp);
    encoder->decode(temp, v);
    encoder->encode(v, scale, temp);
    encryptor->encrypt(temp, ct);
    data.resize(ct.save_size(compr_mode_type::zstd));
    comm += ct.save(data.data(), data.size(), compr_mode_type::zstd);
    auto end = high_resolution_clock::now();
    cout << duration_cast<milliseconds>(end - start).count() / 2 << " milliseconds" << endl;
    // cout << "depth = " <<
    // context->get_context_data(ct.parms_id())->chain_index() << "\n";
}

void CKKSEvaluator::print_decrypted_ct(Ciphertext &ct, int nums)
{
    Plaintext temp;
    vector<double> v;
    decryptor->decrypt(ct, temp);
    encoder->decode(temp, v);
    for (int i = 0; i < nums; i++) {
        cout << v[i] << " ";
    }
    cout << "\n";
}

double CKKSEvaluator::calculateMAE(vector<double> &y_true, Ciphertext &ct)
{
    Plaintext temp;
    vector<double> y_pred;
    decryptor->decrypt(ct, temp);
    encoder->decode(temp, y_pred);

    double sum_absolute_errors = 0.0;
    for (size_t i = 0; i < y_true.size(); ++i) {
        sum_absolute_errors += abs(y_true[i] - y_pred[i]);
    }

    return sum_absolute_errors / y_true.size();
}

vector<double> CKKSEvaluator::init_vec_with_value(int N, double init_value)
{
    std::vector<double> v(N);

    for (int i = 0; i < N; ++i) {
        v[i] = init_value;
    }

    return v;
}

vector<double> CKKSEvaluator::init_mask(int N, int m)
{
    std::vector<double> v(N);

    for (int i = 0; i < N / m; ++i) {
        if (i % 2 == 0) {
            for (int j = 0; j < m; ++j) {
                v[i * m + j] = 1;
            }
        } else {
            for (int j = 0; j < m; ++j) {
                v[i * m + j] = 0;
            }
        }
    }

    return v;
}

Ciphertext CKKSEvaluator::poly_eval(Ciphertext x, vector<Plaintext> coeff)
{
    // cout << "Initial depth " <<
    // context->get_context_data(x.parms_id())->chain_index() << "\n"; x^2
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

    // x^5
    Ciphertext x_5;

    // x^7
    Ciphertext x_7;

    Ciphertext new_x, new_x3, sum_5_7;

    evaluator->mod_switch_to_inplace(coeff[1], x.parms_id());
    evaluator->multiply_plain(x, coeff[1], new_x);
    evaluator->rescale_to_next_inplace(new_x);

    evaluator->mod_switch_to_inplace(coeff[3], x_3.parms_id());
    evaluator->multiply_plain(x_3, coeff[3], new_x3);
    evaluator->rescale_to_next_inplace(new_x3);

    // cout << x.scale() << " " << x_3.scale() << " " << scale << endl;

    evaluator->mod_switch_to_inplace(coeff[5], x_3.parms_id());
    evaluator->mod_switch_to_inplace(x, x_3.parms_id());
    evaluator->multiply_plain(x, coeff[5], x_5);
    evaluator->rescale_to_next_inplace(x_5);

    evaluator->mod_switch_to_inplace(coeff[7], x_3.parms_id());
    evaluator->multiply_plain(x_3, coeff[7], x_7);
    evaluator->rescale_to_next_inplace(x_7);

    // cout << x_5.scale() << " " << x_7.scale() << " " << scale << endl;

    x_5.scale() = scale;
    x_7.scale() = scale;

    evaluator->add(x_5, x_7, sum_5_7);

    evaluator->mod_switch_to_next_inplace(x_4);
    evaluator->multiply_inplace(sum_5_7, x_4);
    evaluator->relinearize_inplace(sum_5_7, *relin_keys);
    evaluator->rescale_to_next_inplace(sum_5_7);

    // c7*x^7 + c5*x^5 + c3*x^3 + c1*x
    // x.scale() = scale;
    // x_3.scale() = scale;
    // x_5.scale() = scale;
    // x_7.scale() = scale;

    new_x.scale() = scale;
    new_x3.scale() = scale;
    sum_5_7.scale() = scale;

    evaluator->mod_switch_to_inplace(new_x, sum_5_7.parms_id());
    evaluator->mod_switch_to_inplace(new_x3, sum_5_7.parms_id());

    evaluator->add_inplace(sum_5_7, new_x);
    evaluator->add_inplace(sum_5_7, new_x3);

    // evaluator->mod_switch_to_inplace(x, x_7.parms_id());
    // evaluator->mod_switch_to_inplace(x_3, x_7.parms_id());
    // evaluator->mod_switch_to_inplace(x_5, x_7.parms_id());
    // evaluator->add_inplace(x, x_3);
    // evaluator->add_inplace(x, x_5);
    // evaluator->add_inplace(x, x_7);

    // cout << "Final depth " <<
    // context->get_context_data(x.parms_id())->chain_index() << "\n";

    return sum_5_7;
}

Ciphertext CKKSEvaluator::sgn_eval2(Ciphertext x, int d_g, int d_f)
{
    Ciphertext dest = x;
    for (int i = 0; i < d_g; i++) {
        // cout << "depth: " << context->get_context_data(dest.parms_id())->chain_index() << endl;
        // if (context->get_context_data(dest.parms_id())->chain_index() < 4) {
        //     re_encrypt(dest);
        // }
        if (i == d_g - 1) {
            eval_odd_deg9_poly(g4_coeffs_last, dest, dest);
        } else {
            eval_odd_deg9_poly(g4_coeffs, dest, dest);
        }
    }
    for (int i = 0; i < d_f; i++) {
        // cout << "depth: " << context->get_context_data(dest.parms_id())->chain_index() << endl;

        // if (context->get_context_data(dest.parms_id())->chain_index() < 4) {
        //     re_encrypt(dest);
        // }
        if (i == d_f - 1) {
            eval_odd_deg9_poly(f4_coeffs_last, dest, dest);
        } else {
            eval_odd_deg9_poly(f4_coeffs, dest, dest);
        }
    }
    // re_encrypt(x);
    return dest;
}

void CKKSEvaluator::eval_odd_deg9_poly(vector<double> &a, Ciphertext &x, Ciphertext &dest)
{
    /*
        (polyeval/odd9.h)
        P(x) = a9 x^9 + a7 x^7 + a5 x^5 + a3 x^3 + a1 x

        T1 = (a3 + a5 x^2) x^3
        T2 = (a7 x + a9 x^3) x^6
        T3 = a1 x
        P(x) = T1 + T2 + T3

        Depth=4, #Muls=5

        Exactly what babystep_giantstep would do, but written explicitly to optimize

        ###

        -> Errorless Polynomial Evaluation (3.2. of https://eprint.iacr.org/2020/1203)
        GOAL: evaluate a polynomial exactly so no need to stabilize and lose precision
        (x at level L and scale D --> P(x) at level L-4 and scale D)
        it's possible to do this exactly for polyeval as (x,x2,x3,x6) determine the scale D_L for each involved level L:
        (assume the primes at levels L to L-4 are p, q, r, s)

        level       ctx       scale (D_l)
        ==================================
          L          x          D
          L-1        x2         D^2 / p
          L-2        x3         D^3 / pq
          L-3        x6         D^6 / p^2 q^2 r

        Just by encoding constants at different scales we can make every ctx at level l be at scale D_l
        (not possible in general, e.g. rescale(x2*x2) produces L-2 ciphertext with scale D^4/ppq)
        (to fix this we would use the Adjust op. that multiplies ctx by constants and Algo 3 for primes from https://eprint.iacr.org/2020/1118)

        Now we know that sc(P(x)) should be D, so we recursively go back to compute the scales for each coefficient
        sc(T1)=sc(T2)=sc(T3)=sc(P(x))=D

        T3:
            sc(a1) = q (should be p but it gets multiplied with modswitched x)

        T2:
            sc(x^6) = D^6 / p^2 q^2 r, so sc(a7*x) = sc(a9*x^3) = p^2 q^2 r s / D^5
            next, sc(a7) = p^2 q^3 r s / D^6
            similarly, sc(a9) = p^3 q^3 r^2 s / D^8

        T1:
            sc(x^3) = D^3 / pq
            implying sc(a3) = pqr / D^2 and also sc(a5*x^2) = pqr / D^2
            as sc(x^2) = D^2 / p this implies sc(a5) = p^2 q^2 r / D^4
    */
    // chrono::high_resolution_clock::time_point time_start, time_end;
    // time_start = high_resolution_clock::now();
    double D = x.scale(); // maybe not init_scale but preserved

    uint64_t p = get_modulus(x, 1);
    uint64_t q = get_modulus(x, 2);
    uint64_t r = get_modulus(x, 3);
    uint64_t s = get_modulus(x, 4);
    uint64_t t = get_modulus(x, 5);

    p = q;
    q = r;
    r = s;
    s = t;

    vector<double> a_scales(10);
    a_scales[1] = q;
    a_scales[3] = (double)p / D * q / D * r;
    a_scales[5] = (double)p / D * p / D * q / D * q / D * r;
    a_scales[7] = (double)p / D * p / D * q / D * q / D * q / D * r / D * s;
    a_scales[9] = (double)p / D * p / D * p / D * q / D * q / D * q / D * r / D * r / D * s;

    ///////////////////////////////////////////////
    Ciphertext x2, x3, x6;

    evaluator->square(x, x2);
    evaluator->relinearize_inplace(x2, *relin_keys);
    evaluator->rescale_to_next_inplace(x2); // L-1

    evaluator->mod_switch_to_next_inplace(x); // L-1
    evaluator->multiply(x2, x, x3);
    evaluator->relinearize_inplace(x3, *relin_keys);
    evaluator->rescale_to_next_inplace(x3); // L-2

    evaluator->square(x3, x6);
    evaluator->relinearize_inplace(x6, *relin_keys);
    evaluator->rescale_to_next_inplace(x6); // L-3

    Plaintext a1, a3, a5, a7, a9;

    // Build T1
    Ciphertext T1;
    double a5_scale = D / x2.scale() * p / x3.scale() * q;
    encoder->encode(a[5], x2.parms_id(), a5_scale, a5); // L-1
    evaluator->multiply_plain(x2, a5, T1);
    evaluator->rescale_to_next_inplace(T1); // L-2

    // Update: using a_scales[3] is only approx. correct, so we directly use T1.scale()
    encoder->encode(a[3], T1.parms_id(), T1.scale(), a3); // L-2

    evaluator->add_plain_inplace(T1, a3); // L-2

    evaluator->multiply_inplace(T1, x3);
    evaluator->relinearize_inplace(T1, *relin_keys);
    evaluator->rescale_to_next_inplace(T1); // L-3

    // Build T2
    Ciphertext T2;
    Plaintext a9_switched;
    double a9_scale = D / x3.scale() * r / x6.scale() * q;
    encoder->encode(a[9], x3.parms_id(), a9_scale, a9); // L-2
    evaluator->multiply_plain(x3, a9, T2);
    evaluator->rescale_to_next_inplace(T2); // L-3

    Ciphertext a7x;
    double a7_scale = T2.scale() / x.scale() * p;
    encoder->encode(a[7], x.parms_id(), a7_scale, a7); // L-1 (x was modswitched)
    evaluator->multiply_plain(x, a7, a7x);
    evaluator->rescale_to_next_inplace(a7x);              // L-2
    evaluator->mod_switch_to_inplace(a7x, T2.parms_id()); // L-3

    double mid_scale = (T2.scale() + a7x.scale()) / 2;
    T2.scale() = a7x.scale() = mid_scale; // this is the correct scale now, need to set it still to avoid SEAL assert
    evaluator->add_inplace(T2, a7x);      // L-3
    evaluator->multiply_inplace(T2, x6);
    evaluator->relinearize_inplace(T2, *relin_keys);
    evaluator->rescale_to_next_inplace(T2); // L-4

    // Build T3
    Ciphertext T3;
    encoder->encode(a[1], x.parms_id(), p, a1); // L-1 (x was modswitched)
    evaluator->multiply_plain(x, a1, T3);
    evaluator->rescale_to_next_inplace(T3); // L-2

    // T1, T2 and T3 should be on the same scale up to floating point
    // but we still need to set them manually to avoid SEAL assert
    double mid3_scale = (T1.scale() + T2.scale() + T3.scale()) / 3;
    T1.scale() = T2.scale() = T3.scale() = mid3_scale;

    dest = T2;
    evaluator->mod_switch_to_inplace(T1, dest.parms_id()); // L-4
    evaluator->add_inplace(dest, T1);
    evaluator->mod_switch_to_inplace(T3, dest.parms_id()); // L-4
    evaluator->add_inplace(dest, T3);

    /////////////////////////////////////////
    // it should be ==D but we don't stabilize if it's not, D' != D is ok
    // the goal was to make T1+T2+T3 work with minimal loss in precision
    // time_end = high_resolution_clock::now();
    // cout << "Poly eval took " << duration_cast<milliseconds>(time_end - time_start).count() << " ms" << endl;
}

Ciphertext CKKSEvaluator::sgn_eval(Ciphertext x, int d_g, int d_f, double factor)
{
    vector<double> df_coeff = { 0, 35.0 / 16, 0, -35.0 / 16, 0, 21.0 / 16, 0, -5.0 / 16 };
    vector<double> dg_coeff = { 0, 4589.0 / 1024, 0, -16577.0 / 1024, 0, 25614.0 / 1024, 0, -12860.0 / 1024 };

    vector<double> df_coeff_last;
    vector<double> dg_coeff_last;
    for (size_t i = 0; i < df_coeff.size(); i++) {
        df_coeff_last.push_back(df_coeff[i] * factor);
    }
    for (size_t i = 0; i < dg_coeff.size(); i++) {
        dg_coeff_last.push_back(dg_coeff[i] * factor);
    }

    vector<Plaintext> encode_df(8), encode_dg(8), encode_df_last(8), encode_dg_last(8);
    for (int i = 0; i < 8; i++) {
        encoder->encode(df_coeff[i], scale, encode_df[i]);
        encoder->encode(dg_coeff[i], scale, encode_dg[i]);
        encoder->encode(df_coeff_last[i], scale, encode_df_last[i]);
        encoder->encode(dg_coeff_last[i], scale, encode_dg_last[i]);
    }
    for (int i = 0; i < d_g; i++) {
        if (context->get_context_data(x.parms_id())->chain_index() < 4) {
            re_encrypt(x);
        }
        if (i == d_g - 1) {
            x = poly_eval(x, encode_dg_last);
        } else {
            x = poly_eval(x, encode_dg);
        }
    }
    for (int i = 0; i < d_f; i++) {
        if (context->get_context_data(x.parms_id())->chain_index() < 4) {
            re_encrypt(x);
        }
        if (i == d_f - 1) {
            x = poly_eval(x, encode_df_last);
        } else {
            x = poly_eval(x, encode_df);
        }
    }
    re_encrypt(x);
    return x;
}

Ciphertext CKKSEvaluator::newtonIter(Ciphertext x, Ciphertext res, int iter)
{
    for (int i = 0; i < iter; i++) {
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
    return res;
}

pair<Ciphertext, Ciphertext> CKKSEvaluator::goldSchmidtIter(Ciphertext v, Ciphertext y, int d)
{
    Ciphertext x, h, r, temp;
    Plaintext constant;
    encoder->encode(0.5, scale, constant);

    // GoldSchmidt's algorithm
    evaluator->mod_switch_to_inplace(v, y.parms_id());
    evaluator->multiply(v, y, x);
    evaluator->relinearize_inplace(x, *relin_keys);
    evaluator->rescale_to_next_inplace(x);
    evaluator->mod_switch_to_inplace(constant, y.parms_id());
    evaluator->multiply_plain(y, constant, h);
    evaluator->rescale_to_next_inplace(h);

    for (int i = 0; i < d; i++) {
        encoder->encode(0.5, scale, constant);
        // r = 0.5 - xh
        if (context->get_context_data(x.parms_id())->chain_index() < 3) {
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
    cout << "depth = " << context->get_context_data(y.parms_id())->chain_index() << "\n";
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

    encoder->encode(init_vec_with_value(slot_count, -8.0), x.parms_id(), x.scale(), p0);
    encoder->encode(init_vec_with_value(slot_count, 1.5), x.parms_id(), x.scale(), p1);
    encoder->encode(init_vec_with_value(slot_count, 1.0 / 32), x.parms_id(), x.scale(), delta);

    evaluator->sub_plain(x, p0, b0);
    evaluator->multiply_plain_inplace(b0, delta);
    evaluator->rescale_to_next_inplace(b0);
    evaluator->sub_plain(x, p1, b1);
    evaluator->multiply_plain_inplace(b1, delta);
    evaluator->rescale_to_next_inplace(b1);

    b0 = sgn_eval2(b0, 2, 2);
    b1 = sgn_eval2(b1, 2, 2);

    Plaintext zero_point_five;
    encoder->encode(init_vec_with_value(slot_count, 0.5), b1.parms_id(), b1.scale(), zero_point_five);
    Ciphertext a1, a2;

    evaluator->sub(b0, b1, a1);                    // a1 = b0 - b1
    evaluator->add_plain(b1, zero_point_five, a2); // a2 = b1 + 0.5

    double A[] = { 0.999469891622, 0.998104199650, 0.501415542413, 0.169660297661, 0.042133244334,
                   0.007501312598, 0.000879175634, 0.000059258169, 0.000001716078 };
    double B[] = { 6.943979090878,  -16.061172554433, 21.461821805218, -14.267410003218, 6.156317208726,
                   -1.632082712072, 0.275766518989,   -0.026342660111, 0.001204185268 };

    vector<Plaintext> coeff_A(9), coeff_B(9);
    for (size_t i = 0; i < 9; i++) {
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

Ciphertext CKKSEvaluator::exp_poly_eval(vector<Ciphertext> x_pow, vector<Plaintext> coeff)
{
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