#include <fstream>

#include "matrix_mul.cuh"
#include "utils.cuh"

using namespace std;
using namespace phantom::util;
using namespace phantom::arith;
using namespace nexus;

__global__ void kernel_compress_ciphertext(uint64_t *plain_data, size_t plain_scale, size_t degree,
                                           const DModulus *moduli, const double *values) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < degree) {
    auto coeffd = std::round(values[i] * plain_scale);
    bool is_negative = std::signbit(coeffd);
    auto coeffu = static_cast<std::uint64_t>(std::fabs(coeffd));

    if (is_negative) {
      for (std::size_t j = 0; j < 2; j++) {
        plain_data[i + (j * degree)] = negate_uint64_mod(
            barrett_reduce_uint64_uint64(coeffu, moduli[j].value(), moduli[j].const_ratio()[1]), moduli[j].value());
      }
    } else {
      for (std::size_t j = 0; j < 2; j++) {
        plain_data[i + (j * degree)] = barrett_reduce_uint64_uint64(coeffu, moduli[j].value(), moduli[j].const_ratio()[1]);
      }
    }
  }
}

// FIXME:
__global__ void kernel_negacyclic_shift(const uint64_t *cipher_data, size_t cipher_count, uint64_t coeff_count, size_t mod_count,
                                        int shift, const DModulus *moduli, uint64_t *result) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx < cipher_count * mod_count * coeff_count) {
    int i = idx / (mod_count * coeff_count);
    int j = (idx / coeff_count) % mod_count;
    size_t k = idx % coeff_count;
    uint64_t mask = coeff_count - 1;

    uint64_t index = (shift + k) & mask;
    uint64_t result_index = i * mod_count * coeff_count + j * coeff_count + index;
    if (!cipher_data[idx] || !(static_cast<uint64_t>(shift + k) & coeff_count)) {
      result[result_index] = cipher_data[idx];
    } else {
      result[result_index] = moduli[j].value() - cipher_data[idx];
    }
  }
}

void MMEvaluator::multiply_power_of_x(PhantomCiphertext &encrypted, PhantomCiphertext &destination, int index) {
  auto context = ckks->context;
  auto coeff_count = ckks->degree;
  auto param = context->get_context_data(encrypted.params_id()).parms();
  auto moduli = param.coeff_modulus();
  // auto moduli = context->gpu_rns_tables().modulus();
  auto coeff_mod_count = param.coeff_modulus().size();
  auto encrypted_count = encrypted.size();
  auto rns_coeff_count = coeff_count * coeff_mod_count;

  const auto &stream = phantom::util::global_variables::default_stream->get_stream();

  ckks->evaluator.transform_from_ntt_inplace(encrypted);

  destination = encrypted;

  // uint64_t gridDimGlb = ceil_div<uint64_t>(encrypted_count * rns_coeff_count, blockDimGlb.x);
  // kernel_negacyclic_shift<<<gridDimGlb, blockDimGlb, 0, stream>>>(
  //     encrypted.data(), encrypted_count, coeff_count, coeff_mod_count, index, moduli, destination.data());

  auto dest_data = new uint64_t[rns_coeff_count * encrypted_count];
  auto ct_data = new uint64_t[rns_coeff_count * encrypted_count];
  cudaMemcpyAsync(ct_data, encrypted.data(), encrypted_count * rns_coeff_count * sizeof(uint64_t), cudaMemcpyDeviceToHost, stream);

  for (int i = 0; i < encrypted_count; i++) {
    for (int j = 0; j < coeff_mod_count; j++) {
      uint64_t *poly = ct_data + i * rns_coeff_count + j * coeff_count;
      uint64_t *result = dest_data + i * rns_coeff_count + j * coeff_count;

      uint64_t index_raw = index;
      uint64_t coeff_count_mod_mask = static_cast<uint64_t>(coeff_count) - 1;
      for (size_t k = 0; k < coeff_count; k++, poly++, index_raw++) {
        uint64_t index = index_raw & coeff_count_mod_mask;
        if (!(index_raw & static_cast<uint64_t>(coeff_count)) || !*poly) {
          result[index] = *poly;
        } else {
          result[index] = moduli[j].value() - *poly;
        }
      }
    }
  }

  cudaMemcpyAsync(destination.data(), dest_data, encrypted_count * rns_coeff_count * sizeof(uint64_t), cudaMemcpyHostToDevice, stream);

  ckks->evaluator.transform_to_ntt_inplace(encrypted);
  ckks->evaluator.transform_to_ntt_inplace(destination);
}

void MMEvaluator::enc_compress_ciphertext(vector<double> &values, PhantomCiphertext &ct) {
  size_t plain_scale = 10000000000;

  auto &context_data = ckks->context->first_context_data();
  auto param = context_data.parms();
  auto moduli = ckks->context->gpu_rns_tables().modulus();
  auto coeff_modulus_size = param.coeff_modulus().size();
  auto poly_modulus_degree = param.poly_modulus_degree();

  const phantom::util::cuda_stream_wrapper &stream_wrapper = *phantom::util::global_variables::default_stream;
  const auto &stream = stream_wrapper.get_stream();

  // mod_count = 2, poly_degree = 8192 => p coeff_count = 2 * 8192 = 16384
  PhantomPlaintext p;
  p.resize(coeff_modulus_size, poly_modulus_degree, stream);

  auto gpu_values = make_cuda_auto_ptr<double>(values.size(), stream);
  cudaMemcpyAsync(gpu_values.get(), values.data(), values.size() * sizeof(double), cudaMemcpyHostToDevice, stream);

  kernel_compress_ciphertext<<<poly_modulus_degree / blockDimGlb.x, blockDimGlb, 0, stream>>>(
      p.data(), plain_scale, poly_modulus_degree, moduli, gpu_values.get());

  // Transform polynomials to the NTT domain
  nwt_2d_radix8_forward_inplace(p.data(), ckks->context->gpu_rns_tables(), coeff_modulus_size, 0, stream);

  // Update plaintext parameters
  p.parms_id() = context_data.parms().parms_id();
  p.set_chain_index(context_data.chain_index());
  p.scale() = plain_scale;

  // Create a ciphertext encrypting zero
  PhantomPlaintext zero_pt;
  PhantomCiphertext zero;
  ckks->encoder.encode(0.0, ckks->scale, zero_pt);
  ckks->encryptor.encrypt(zero_pt, zero);

  zero.scale() = p.scale();
  ckks->evaluator.add_plain(zero, p, ct);
}

vector<PhantomCiphertext> MMEvaluator::decompress_ciphertext(const PhantomCiphertext &encrypted) {
  auto N = ckks->degree;
  uint32_t logN = ceil(log2(N));

  vector<PhantomCiphertext> temp;
  temp.push_back(encrypted);

  PhantomCiphertext tempctxt_rotated;
  PhantomCiphertext tempctxt_shifted;
  PhantomCiphertext tempctxt_rotatedshifted;

  for (uint32_t i = 0; i < logN; i++) {
    vector<PhantomCiphertext> newtemp(temp.size() << 1);

    uint32_t galois_elt = ckks->galois_elts[i];
    int index_raw = (N << 1) - (1 << i);
    int index = (index_raw * galois_elt) % (N << 1);

    for (uint32_t a = 0; a < temp.size(); a++) {
      ckks->evaluator.apply_galois(temp[a], galois_elt, *(ckks->galois_keys), tempctxt_rotated);  // sub
      ckks->evaluator.add(temp[a], tempctxt_rotated, newtemp[a]);
      multiply_power_of_x(temp[a], tempctxt_shifted, index_raw);  // x**-1
      multiply_power_of_x(tempctxt_rotated, tempctxt_rotatedshifted, index);
      ckks->evaluator.add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
    }

    temp = newtemp;
  }

  return temp;
}

void MMEvaluator::matrix_mul(vector<vector<double>> &x, vector<vector<double>> &y, vector<PhantomCiphertext> &res) {
  const phantom::util::cuda_stream_wrapper &stream_wrapper = *phantom::util::global_variables::default_stream;
  const auto &stream = stream_wrapper.get_stream();

  auto timer = Timer();

  // Encode plaintext
  vector<PhantomPlaintext> a_pts;
  a_pts.reserve(768);

  for (int i = 0; i < 768; i++) {
    PhantomPlaintext pt;
    ckks->encoder.encode(x[i], ckks->scale, pt);
    a_pts.emplace_back(pt);
  }

  // Ciphertext encoding & compression
  timer.start();

  int b_cts_count = 768 * 64 / ckks->degree;
  vector<PhantomCiphertext> b_compressed_cts;
  b_compressed_cts.reserve(b_cts_count);

  for (int i = 0; i < b_cts_count; i++) {
    PhantomCiphertext ct;
    enc_compress_ciphertext(y[i], ct);
    b_compressed_cts.emplace_back(ct);
  }

  timer.stop();
  cout << "Compression took: " << timer.duration<milliseconds>() << " milliseconds" << endl;

  // Ciphertext decompression
  timer.start();

  vector<PhantomCiphertext> b_expanded_cts;

  for (auto i = 0; i < b_compressed_cts.size(); i++) {
    vector<PhantomCiphertext> temp_cts = decompress_ciphertext(b_compressed_cts[i]);
    cout << "Expanded ciphertext #" << i + 1 << endl;
    b_expanded_cts.insert(b_expanded_cts.end(), make_move_iterator(temp_cts.begin()), make_move_iterator(temp_cts.end()));
  }

  timer.stop();
  cout << "Decompression took: " << timer.duration<seconds>() << " seconds" << endl;

  // Perform plain-cipher matrix multiplication
  timer.start();

  for (int i = 0; i < 64; i++) {
    PhantomCiphertext res_col_ct;
    vector<PhantomCiphertext> temp_cts(768);

    for (int j = 0; j < 768; j++) {
      ckks->evaluator.multiply_plain(b_expanded_cts[i * 768 + j], a_pts[j], temp_cts[j]);
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

  timer.stop();
  cout << "Result calculation time: " << timer.duration<milliseconds>() << " milliseconds" << endl;
}
