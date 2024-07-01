#include "matrix_mul.cuh"
#include "utils.cuh"

using namespace nexus;

vector<PhantomCiphertext> MMEvaluator::expand_ciphertext(
    const PhantomCiphertext &encrypted, uint32_t m, PhantomGaloisKey &galkey, vector<uint32_t> &galois_elts) {
  uint32_t logm = ceil(log2(m));
  auto n = ckks->degree;

  vector<PhantomCiphertext> temp;
  temp.push_back(encrypted);

  PhantomCiphertext tempctxt;
  PhantomCiphertext tempctxt_rotated;
  PhantomCiphertext tempctxt_shifted;
  PhantomCiphertext tempctxt_rotatedshifted;

  for (uint32_t i = 0; i < logm; i++) {
    vector<PhantomCiphertext> newtemp(temp.size() << 1);
    int index_raw = (n << 1) - (1 << i);
    int index = (index_raw * galois_elts[i]) % (n << 1);
    for (uint32_t a = 0; a < temp.size(); a++) {
      ckks->evaluator.apply_galois(temp[a], ckks->rots[i], *(ckks->galois_keys), tempctxt_rotated);  // sub
      ckks->evaluator.add(temp[a], tempctxt_rotated, newtemp[a]);
      ckks->multiply_power_of_x(temp[a], tempctxt_shifted, index_raw);  // x**-1
      ckks->multiply_power_of_x(tempctxt_rotated, tempctxt_rotatedshifted, index);
      ckks->evaluator.add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
    }
    temp = newtemp;
  }
  return temp;
}

void MMEvaluator::expand_encode(vector<double> &val, PhantomCiphertext &ct) {
  PhantomPlaintext zero_pt;
  ckks->encoder.encode(std::vector<double>(ckks->degree / 2, 0.0), ckks->scale, zero_pt);
  PhantomCiphertext zero;
  ckks->encryptor.encrypt(zero_pt, zero);

  auto &context_data = ckks->context->first_context_data();
  auto param = context_data.parms();

  auto poly_modulus_degree = ckks->degree;
  const phantom::util::cuda_stream_wrapper &stream_wrapper = *phantom::util::global_variables::default_stream;
  const auto &stream = stream_wrapper.get_stream();

  PhantomPlaintext p;
  p.resize(param.coeff_modulus().size(), poly_modulus_degree * 2, stream_wrapper.get_stream());

  for (auto i = 0; i < poly_modulus_degree; i++) {
    val[i] = 10.0 * 2.0 * (1.0 * rand() / RAND_MAX - 0.5);
  }
  for (auto i = 0; i < poly_modulus_degree; i++) {
    auto coeffd = std::round(val[i] * 10000000000);
    bool is_negative = std::signbit(coeffd);
    auto coeffu = static_cast<std::uint64_t>(std::fabs(coeffd));
    if (is_negative) {
      for (std::size_t j = 0; j < 2; j++) {
        // Does not work, need to implement to_host, to_device for plaintext data
        p.data()[i + (j * poly_modulus_degree)] = negate_uint_mod(
            barrett_reduce_64(coeffu, param.coeff_modulus()[j]), param.coeff_modulus()[j]);
      }
    } else {
      for (std::size_t j = 0; j < 2; j++) {
        p.data()[i + (j * poly_modulus_degree)] = barrett_reduce_64(coeffu, param.coeff_modulus()[j]);
      }
    }
  }
  for (std::size_t i = 0; i < 2; i++) {
    nwt_2d_radix8_forward_inplace(&p.data()[i * poly_modulus_degree], ckks->context->gpu_rns_tables(), param.coeff_modulus().size(), 0, stream);
  }
  p.set_chain_index(context_data.chain_index());
  p.scale() = 10000000000;

  zero.scale() = p.scale();

  ckks->evaluator.add_plain(zero, p, ct);
}

void MMEvaluator::matrix_mul(vector<vector<double>> &x, vector<vector<double>> &y, vector<PhantomCiphertext> &res) {
  chrono::high_resolution_clock::time_point time_start, time_end;

  vector<PhantomPlaintext> a_pts;
  a_pts.reserve(768);
  for (int i = 0; i < 768; i++) {
    PhantomPlaintext pt;
    ckks->encoder.encode(x[i], ckks->scale, pt);
    a_pts.emplace_back(pt);
  }

  vector<PhantomCiphertext> b_compressed_cts;
  for (int i = 0; i < 768 * 768 / ckks->degree; i++) {
    PhantomPlaintext pt;
    PhantomCiphertext ct;
    expand_encode(y[i], ct);
    b_compressed_cts.push_back(ct);
  }

  vector<PhantomCiphertext> b_expanded_cts;

  auto timer = Timer();
  for (auto i = 0; i < b_compressed_cts.size(); i++) {
    vector<PhantomCiphertext> temp_cts =
        expand_ciphertext(b_compressed_cts[i], ckks->degree, *ckks->galois_keys, ckks->rots);
    cout << "Expanding..." << endl;
    b_expanded_cts.insert(
        b_expanded_cts.end(), make_move_iterator(temp_cts.begin()), make_move_iterator(temp_cts.end()));
  }
  timer.stop();
  cout << "Expanding time: " << timer.duration<seconds>() << " seconds" << endl;

  PhantomPlaintext pt;
  PhantomCiphertext zero;
  ckks->encoder.encode(std::vector<double>(ckks->degree / 2, 0.0), ckks->scale, pt);
  ckks->encryptor.encrypt(pt, zero);

  time_start = high_resolution_clock::now();
  PhantomCiphertext temp;

  for (int i = 0; i < 768; i++) {
    PhantomCiphertext res_col_ct = zero;
    vector<PhantomCiphertext> temp_cts(768);
    for (int j = 0; j < 768; j++) {
      temp_cts[j] = b_expanded_cts[i * 768 + j];
      ckks->evaluator.multiply_plain_inplace(temp_cts[j], a_pts[j]);
    }
    res_col_ct.scale() = temp_cts[0].scale();
    ckks->evaluator.add_many(temp_cts, res_col_ct);
    res_col_ct.scale() *= 4096;
    res.push_back(res_col_ct);
  }

  for (auto &ct : res) {
    while (ct.coeff_modulus_size() > 1) {
      ckks->evaluator.rescale_to_next_inplace(ct);
    }
  }

  time_end = high_resolution_clock::now();
  cout << "calculating res time: " << duration_cast<seconds>(time_end - time_start).count() << " seconds" << endl;
}
