#include <fstream>

#include "matrix_mul.cuh"
#include "utils.cuh"

using namespace phantom::util;
using namespace phantom::arith;
using namespace nexus;

__global__ void negacyclic_shift_poly_coeffmod_kernel(
    const uint64_t *d_poly, size_t poly_degree, size_t shift, DModulus *modulus, size_t coeff_mod_size, uint64_t *d_result) {
  for (size_t tid = blockIdx.x * blockDim.x + threadIdx.x;
       tid < poly_degree * coeff_mod_size;
       tid += blockDim.x * gridDim.x) {
    size_t twr = tid / poly_degree;
    DModulus mod = modulus[twr];

    uint64_t index_raw = shift + tid;
    uint64_t coeff_count_mod_mask = static_cast<uint64_t>(poly_degree) - 1;
    uint64_t index = index_raw & coeff_count_mod_mask;

    if (!(index_raw & static_cast<uint64_t>(poly_degree)) || !*(d_poly + tid)) {
      *(d_result + index) = *(d_poly + tid);
    } else {
      *(d_result + index) = mod.value() - *(d_poly + tid);
    }
  }
}

__global__ void expand_encode_kernel(const double *d_val, size_t poly_modulus_degree, DModulus *modulus, uint64_t *d_p) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < poly_modulus_degree) {
    auto coeffd = round(d_val[i] * 10000000000);
    bool is_negative = signbit(coeffd);
    auto coeffu = static_cast<uint64_t>(fabs(coeffd));

    if (is_negative) {
      for (size_t j = 0; j < 2; j++) {
        d_p[i + (j * poly_modulus_degree)] = negate_uint64_mod(
            barrett_reduce_uint64_uint64(coeffu, modulus[j].value(), modulus[j].const_ratio()[1]), modulus[j].value());
      }
    } else {
      for (size_t j = 0; j < 2; j++) {
        d_p[i + (j * poly_modulus_degree)] = barrett_reduce_uint64_uint64(coeffu, modulus[j].value(), modulus[j].const_ratio()[1]);
      }
    }
  }
}

void MMEvaluator::multiply_power_of_x(PhantomCiphertext &encrypted, PhantomCiphertext &destination, int index) {
  auto context = ckks->context;
  auto coeff_count = ckks->degree;
  // auto moduli_gpu = context->gpu_rns_tables().modulus();
  auto param = context->get_context_data(encrypted.params_id()).parms();
  auto moduli = param.coeff_modulus();
  auto coeff_mod_count = param.coeff_modulus().size();
  auto encrypted_count = encrypted.size();
  auto rns_coeff_count = coeff_count * coeff_mod_count;

  const auto &stream = phantom::util::global_variables::default_stream->get_stream();
  // uint64_t gridDimGlb = coeff_count * coeff_mod_count / blockDimGlb.x;

  ckks->evaluator.transform_from_ntt_inplace(encrypted);

  destination = encrypted;

  auto dest_data = new uint64_t[rns_coeff_count * encrypted_count];
  auto ct_data = new uint64_t[rns_coeff_count * encrypted_count];
  cudaMemcpy(ct_data, encrypted.data(), rns_coeff_count * encrypted_count * sizeof(uint64_t), cudaMemcpyDeviceToHost);

  if (print_count == 0) {
    for (int i = 0; i < 10; i++) {
      cout << ct_data[i] << " ";
    }
    cout << endl;
  }

  for (int i = 0; i < encrypted_count; i++) {
    // negacyclic_shift_poly_coeffmod_kernel<<<gridDimGlb, blockDimGlb, 0, stream>>>(
    //     encrypted.data(i),
    //     coeff_count,
    //     index,
    //     moduli,
    //     coeff_mod_count,
    //     destination.data(i));
    // -----------------------------------------------------------------------------------
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

  if (print_count == 0) {
    for (int i = 0; i < 10; i++) {
      cout << dest_data[i] << " ";
    }
    cout << endl;
    print_count += 1;
  }

  cudaMemcpy(destination.data(), dest_data, rns_coeff_count * encrypted_count * sizeof(uint64_t), cudaMemcpyHostToDevice);

  ckks->evaluator.transform_to_ntt_inplace(encrypted);
  ckks->evaluator.transform_to_ntt_inplace(destination);

  cudaStreamSynchronize(stream);
}

void MMEvaluator::enc_compress_ciphertext(vector<double> &values, PhantomCiphertext &ct) {
  // PhantomPlaintext zero_pt;
  // ckks->encoder.encode(0.0, ckks->scale, zero_pt);
  // PhantomCiphertext zero;
  // ckks->encryptor.encrypt(zero_pt, zero);

  auto &context_data = ckks->context->first_context_data();
  auto param = context_data.parms();

  auto coeff_modulus_size = param.coeff_modulus().size();
  auto poly_modulus_degree = ckks->degree;
  auto rns_coeff_count = poly_modulus_degree * coeff_modulus_size;

  const phantom::util::cuda_stream_wrapper &stream_wrapper = *phantom::util::global_variables::default_stream;
  const auto &stream = stream_wrapper.get_stream();

  PhantomPlaintext p;                                         // empty plaintext
  p.resize(coeff_modulus_size, poly_modulus_degree, stream);  // mod_count = 2, poly_degree = 8192 => p coeff_count = 2 * 8192 = 16384

  // // Copy vector values to GPU
  // auto val_gpu = make_cuda_auto_ptr<double>(val.size(), stream);
  // cudaMemcpy(val_gpu.get(), val.data(), val.size(), cudaMemcpyHostToDevice);

  // // Execute expand kernel
  auto gpu_moduli = ckks->context->gpu_rns_tables().modulus();
  // uint64_t block_size = blockDimGlb.x;
  // int num_blocks = (poly_modulus_degree + block_size - 1) / block_size;
  // expand_encode_kernel<<<num_blocks, block_size, 0, stream>>>(val_gpu.get(), poly_modulus_degree, gpu_moduli, p.data());

  cudaStreamSynchronize(stream);

  auto p_data = new uint64_t[rns_coeff_count];

  // Coefficients of the two RNS polynomails should be the same except with different mod
  for (auto i = 0; i < poly_modulus_degree; i++) {
    auto coeffd = std::round(values[i] * 10000000000);
    bool is_negative = std::signbit(coeffd);
    auto coeffu = static_cast<std::uint64_t>(std::fabs(coeffd));
    if (is_negative) {
      for (std::size_t j = 0; j < 2; j++) {
        p_data[i + (j * poly_modulus_degree)] = negate_uint_mod(
            barrett_reduce_64(coeffu, param.coeff_modulus()[j]), param.coeff_modulus()[j]);
      }
    } else {
      for (std::size_t j = 0; j < 2; j++) {
        p_data[i + (j * poly_modulus_degree)] = barrett_reduce_64(coeffu, param.coeff_modulus()[j]);
      }
    }
  }

  // Debug
  cout << "Plaintext RNS polynomial coefficients:" << endl;
  for (auto i = 0; i < 10; i++) {
    cout << p_data[i] << " ";
  }
  cout << endl;

  // Copy plaintext data to GPU (synchronous)
  cudaMemcpy(p.data(), p_data, rns_coeff_count * sizeof(uint64_t), cudaMemcpyHostToDevice);

  // Transform all 2 RNS polynomials to the NTT domain
  nwt_2d_radix8_forward_inplace(p.data(), ckks->context->gpu_rns_tables(), coeff_modulus_size, 0, stream);

  // Debug
  cudaStreamSynchronize(stream);
  cudaMemcpy(p_data, p.data(), rns_coeff_count * sizeof(uint64_t), cudaMemcpyDeviceToHost);
  cout << "Plaintext RNS polynomial coefficients in NTT:" << endl;
  for (auto i = 0; i < 10; i++) {
    cout << p_data[i] << " ";
  }
  cout << endl;

  // Update plaintext parameters
  p.set_chain_index(context_data.chain_index());
  p.scale() = 10000000000;

  // Encrypt zero asymmetrically
  ckks->encryptor.encrypt_zero(ct, p.chain_index());

  // (c0, c1) => (c0 + p, c1)
  uint64_t gridDimGlb = rns_coeff_count / blockDimGlb.x;
  add_rns_poly<<<gridDimGlb, blockDimGlb, 0, stream>>>(
      ct.data(), p.data(), gpu_moduli, ct.data(), poly_modulus_degree, coeff_modulus_size);

  ct.scale() = p.scale();
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
    vector<PhantomCiphertext> newtemp(temp.size() * 2);
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
  auto timer = Timer();

  vector<PhantomPlaintext> a_pts;
  a_pts.reserve(768);

  for (int i = 0; i < 768; i++) {
    PhantomPlaintext pt;
    ckks->encoder.encode(x[i], ckks->scale, pt);
    a_pts.emplace_back(pt);
  }

  timer.start();
  int b_cts_count = 768 * 64 / ckks->degree;
  vector<PhantomCiphertext> b_compressed_cts;
  b_compressed_cts.reserve(b_cts_count);

  for (int i = 0; i < b_cts_count; i++) {
    PhantomCiphertext ct;
    enc_compress_ciphertext(y[i], ct);
    b_compressed_cts.push_back(ct);
  }
  timer.stop();
  cout << "Compression took: " << timer.duration<milliseconds>() << " milliseconds" << endl;

  timer.start();
  vector<PhantomCiphertext> b_expanded_cts;
  for (auto i = 0; i < b_compressed_cts.size(); i++) {
    vector<PhantomCiphertext> temp_cts = decompress_ciphertext(b_compressed_cts[i]);
    cout << "Expanded ciphertext #" << i << endl;

    // cudaStreamSynchronize(phantom::util::global_variables::default_stream->get_stream());
    // auto temp_ct0_data = new uint64_t[ckks->degree * 2];
    // cudaMemcpy(temp_ct0_data, temp_cts[0].data(), ckks->degree * 2 * sizeof(uint64_t), cudaMemcpyDeviceToHost);

    // for (int j = 0; j < 10; j++) {
    //   cout << temp_ct0_data[j] << " ";
    // }
    // cout << endl;

    ckks->print_decrypted_ct(temp_cts[0], 10);

    b_expanded_cts.insert(
        b_expanded_cts.end(), make_move_iterator(temp_cts.begin()), make_move_iterator(temp_cts.end()));
  }
  timer.stop();
  cout << "Decompression took: " << timer.duration<seconds>() << " seconds" << endl;

  PhantomPlaintext zero_pt;
  PhantomCiphertext zero;
  ckks->encoder.encode(0.0, ckks->scale, zero_pt);
  ckks->encryptor.encrypt(zero_pt, zero);

  timer.start();
  for (int i = 0; i < 64; i++) {
    PhantomCiphertext res_col_ct = zero;
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
