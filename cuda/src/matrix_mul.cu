#include "matrix_mul.cuh"
#include "utils.cuh"

using namespace troy::utils;
using namespace nexus;

__global__ void kernel_compress_ciphertext(Slice<uint64_t> plain_poly, size_t degree, ConstSlice<troy::Modulus> moduli, ConstSlice<double> values) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < degree) {
    auto coeffd = std::round(values[i] * 10000000000);
    bool is_negative = std::signbit(coeffd);
    auto coeffu = static_cast<std::uint64_t>(std::fabs(coeffd));

    if (is_negative) {
      for (std::size_t j = 0; j < 2; j++) {
        plain_poly[i + (j * degree)] = negate_uint64_mod(
            barrett_reduce_uint64(coeffu, moduli[j]), moduli[j]);
      }
    } else {
      for (std::size_t j = 0; j < 2; j++) {
        plain_poly[i + (j * degree)] = barrett_reduce_uint64(coeffu, moduli[j]);
      }
    }
  }
}

void MMEvaluator::multiply_power_of_x(Ciphertext &encrypted, Ciphertext &destination, int index) {
  destination = encrypted;
  ckks->troy_evaluator->transform_from_ntt_inplace(destination);
  ckks->troy_evaluator->negacyclic_shift_inplace(destination, index);
  ckks->troy_evaluator->transform_to_ntt_inplace(destination);
}

void MMEvaluator::enc_compress_ciphertext(vector<double> &values, Ciphertext &ct) {
  Plaintext zero_pt;
  Ciphertext zero;
  std::vector<complex<double>> zero_input(ckks->slot_count, 0.0 + 0.0i);

  ckks->troy_encoder->encode_complex64_simd(zero_input, std::nullopt, ckks->scale, zero_pt);
  ckks->troy_encryptor->encrypt_asymmetric(zero_pt, zero);

  auto context_data = ckks->troy_context->get_context_data(ckks->troy_context->first_parms_id()).value();
  auto param = context_data->parms();
  auto ntt_tables = context_data->small_ntt_tables();
  auto poly_modulus_degree = ckks->degree;

  Plaintext p;
  p.to_device_inplace();
  p.resize(poly_modulus_degree * 2);

  ConstSlice<double> values_slice(values.data(), values.size(), false, nullptr);
  Array<double> values_device = Array<double>::create_and_copy_from_slice(values_slice, true, MemoryPool::GlobalPool());
  kernel_compress_ciphertext<<<ceil(poly_modulus_degree / KERNEL_THREAD_COUNT), KERNEL_THREAD_COUNT>>>(
      p.poly(), poly_modulus_degree, param.coeff_modulus(), values_device.const_reference());

  ntt_negacyclic_harvey_p(p.poly(), poly_modulus_degree, ntt_tables);

  p.parms_id() = ckks->troy_context->first_parms_id();
  p.coeff_modulus_size() = param.coeff_modulus().size();
  p.poly_modulus_degree() = poly_modulus_degree;
  p.is_ntt_form() = true;
  p.scale() = 10000000000;

  zero.scale() = p.scale();

  ckks->troy_evaluator->add_plain(zero, p, ct);
}

vector<Ciphertext> MMEvaluator::decompress_ciphertext(const Ciphertext &encrypted) {
  auto timer = Timer();

  auto N = ckks->degree;
  uint32_t logN = ceil(log2(N));

  cout << logN << endl;

  vector<Ciphertext> temp;
  temp.push_back(encrypted);

  Ciphertext tempctxt_rotated;
  Ciphertext tempctxt_shifted;
  Ciphertext tempctxt_rotatedshifted;

  for (uint32_t i = 0; i < logN; i++) {
    vector<Ciphertext> newtemp(temp.size() << 1);

    uint64_t galois_elt = ckks->troy_galois_elts[i];
    int index_raw = (N << 1) - (1 << i);
    int index = (index_raw * galois_elt) % (N << 1);

    cout << i << " => " << temp.size() << endl;

    for (uint32_t a = 0; a < temp.size(); a++) {
      ckks->troy_evaluator->apply_galois(temp[a], galois_elt, *(ckks->troy_galois_keys), tempctxt_rotated);  // sub
      ckks->troy_evaluator->add(temp[a], tempctxt_rotated, newtemp[a]);
      multiply_power_of_x(temp[a], tempctxt_shifted, index_raw);  // x**-1
      multiply_power_of_x(tempctxt_rotated, tempctxt_rotatedshifted, index);
      ckks->troy_evaluator->add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
    }

    temp = newtemp;
  }

  return temp;
}

void MMEvaluator::matrix_mul(vector<vector<double>> &x, vector<vector<double>> &y, vector<Ciphertext> &res) {
  auto timer = Timer();

  // Encode plaintext
  vector<Plaintext> a_pts;
  a_pts.reserve(768);

  for (int i = 0; i < 768; i++) {
    Plaintext pt;
    vector<complex<double>> complex_x_i;
    for (auto &v : x[i]) {
      complex_x_i.emplace_back(v, 0.0);
    }
    ckks->troy_encoder->encode_complex64_simd(complex_x_i, std::nullopt, ckks->scale, pt);
    a_pts.emplace_back(pt);
  }

  // Ciphertext encoding & compression
  timer.start();

  int b_cts_count = 768 * 64 / ckks->degree;
  vector<Ciphertext> b_compressed_cts;
  b_compressed_cts.reserve(b_cts_count);

  for (int i = 0; i < b_cts_count; i++) {
    Ciphertext ct;
    enc_compress_ciphertext(y[i], ct);
    b_compressed_cts.push_back(ct);
  }

  timer.stop();
  cout << "Compression took: " << timer.duration<milliseconds>() << " milliseconds" << endl;

  // Ciphertext decompression
  timer.start();

  vector<Ciphertext> b_expanded_cts;

  for (auto i = 0; i < b_compressed_cts.size(); i++) {
    vector<Ciphertext> temp_cts = decompress_ciphertext(b_compressed_cts[i]);
    cout << "Expanded ciphertext #" << i + 1 << endl;
    b_expanded_cts.insert(b_expanded_cts.end(), make_move_iterator(temp_cts.begin()), make_move_iterator(temp_cts.end()));
  }

  timer.stop();
  cout << "Decompression took: " << timer.duration<seconds>() << " seconds" << endl;

  // Perform plain-cipher matrix multiplication
  timer.start();

  Plaintext zero_pt;
  Ciphertext zero;
  vector<complex<double>> zero_input(ckks->slot_count, 0.0 + 0.0i);

  ckks->troy_encoder->encode_complex64_simd(zero_input, std::nullopt, ckks->scale, zero_pt);
  ckks->troy_encryptor->encrypt_asymmetric(zero_pt, zero);

  Ciphertext temp;

  for (int i = 0; i < 64; i++) {
    Ciphertext res_col_ct = zero;
    vector<Ciphertext> temp_cts(768);

    for (int j = 0; j < 768; j++) {
      ckks->troy_evaluator->multiply_plain(b_expanded_cts[i * 768 + j], a_pts[j], temp_cts[j]);
    }

    res_col_ct.scale() = temp_cts[0].scale();
    for (size_t i = 0; i < temp_cts.size(); i++) {
      ckks->troy_evaluator->add_inplace(res_col_ct, temp_cts[i]);
    }

    res_col_ct.scale() *= 4096;
    res.push_back(res_col_ct);
  }

  for (auto &ct : res) {
    while (ct.coeff_modulus_size() > 1) {
      ckks->troy_evaluator->rescale_to_next_inplace(ct);
    }
  }

  timer.stop();
  cout << "Result calculation time: " << timer.duration<milliseconds>() << " milliseconds" << endl;
}
