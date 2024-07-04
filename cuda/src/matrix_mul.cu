#include <fstream>

#include "matrix_mul.cuh"
#include "utils.cuh"

using namespace nexus;

std::vector<std::vector<double>> MMEvaluator::transpose_matrix(const std::vector<std::vector<double>> &matrix) {
  if (matrix.empty()) {
    return {};
  }
  int rows = matrix.size();
  int cols = matrix[0].size();
  std::vector<std::vector<double>> transposedMatrix(cols, std::vector<double>(rows));

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      transposedMatrix[j][i] = matrix[i][j];
    }
  }

  return transposedMatrix;
}

std::vector<std::vector<double>> MMEvaluator::read_matrix(const std::string &filename, int rows, int cols) {
  std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Can not open file: " << filename << std::endl;
    return matrix;
  }

  std::string line;
  for (int i = 0; i < rows; ++i) {
    if (std::getline(file, line)) {
      std::istringstream iss(line);
      for (int j = 0; j < cols; ++j) {
        if (!(iss >> matrix[i][j])) {
          std::cerr << "read error: " << filename << " (row: " << i << ", column: " << j << ")" << std::endl;
        }
      }
    }
  }

  file.close();
  return matrix;
}

void MMEvaluator::multiply_power_of_x(PhantomCiphertext &encrypted, PhantomCiphertext &destination, int index) {
  auto context = ckks->context;
  auto param = context->first_context_data().parms();

  auto coeff_mod_count = param.coeff_modulus().size();
  auto coeff_count = ckks->degree;
  auto encrypted_count = encrypted.size();

  // // After NTT
  // auto encrypted_data_1 = new uint64_t[encrypted_count * coeff_count * encrypted.coeff_modulus_size()];
  // cudaMemcpy(encrypted_data_1, encrypted.data(), encrypted_count * coeff_count * coeff_mod_count * sizeof(uint64_t), cudaMemcpyDeviceToHost);

  // for (auto i = 0; i < 7; i++) {
  //   cout << fixed << setprecision(5) << encrypted_data_1[i] << ", ";
  // }
  // cout << endl;

  ckks->evaluator.transform_from_ntt_inplace(encrypted);

  destination = encrypted;

  // auto encrypted_data = new uint64_t[encrypted_count * coeff_count * encrypted.coeff_modulus_size()];
  // auto destination_data = new uint64_t[encrypted_count * coeff_count * destination.coeff_modulus_size()];

  // // cudaMemcpy(encrypted_data, encrypted.data(), encrypted_count * coeff_count * coeff_mod_count * sizeof(uint64_t), cudaMemcpyDeviceToHost);
  // cudaMemcpy(destination_data, destination.data(), encrypted_count * coeff_count * coeff_mod_count * sizeof(uint64_t), cudaMemcpyDeviceToHost);

  // // Before NTT
  // for (auto i = 0; i < 7; i++) {
  //   cout << fixed << setprecision(5) << destination_data[i] << ", ";
  // }
  // cout << endl;

  for (int i = 0; i < encrypted_count; i++) {
    for (int j = 0; j < coeff_mod_count; j++) {
      ckks->evaluator.negacyclic_shift_poly_coeffmod(
          encrypted.data(i) + (j * coeff_count),
          coeff_count,
          index,
          param.coeff_modulus()[j],
          destination.data(i) + (j * coeff_count));
    }
  }

  // // After NTT
  // auto destination_data_2 = new uint64_t[encrypted_count * coeff_count * encrypted.coeff_modulus_size()];
  // auto encrypted_data_2 = new uint64_t[encrypted_count * coeff_count * encrypted.coeff_modulus_size()];
  // cudaMemcpy(destination_data_2, destination.data(), encrypted_count * coeff_count * coeff_mod_count * sizeof(uint64_t), cudaMemcpyDeviceToHost);
  // cudaMemcpy(encrypted_data_2, encrypted.data(), encrypted_count * coeff_count * coeff_mod_count * sizeof(uint64_t), cudaMemcpyDeviceToHost);

  // for (auto i = 0; i < coeff_mod_count; i++) {
  //   cout << fixed << setprecision(5) << destination_data_2[i * coeff_count] << ", ";
  // }
  // cout << endl;
  // for (auto i = 0; i < coeff_mod_count; i++) {
  //   cout << fixed << setprecision(5) << encrypted_data_2[i * coeff_count + coeff_count - 1] << ", ";
  // }
  // cout << endl;

  ckks->evaluator.transform_to_ntt_inplace(encrypted);
  ckks->evaluator.transform_to_ntt_inplace(destination);
}

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
      ckks->evaluator.apply_galois(temp[a], i, *(ckks->galois_keys), tempctxt_rotated);  // sub
      ckks->evaluator.add(temp[a], tempctxt_rotated, newtemp[a]);
      multiply_power_of_x(temp[a], tempctxt_shifted, index_raw);  // x**-1
      multiply_power_of_x(tempctxt_rotated, tempctxt_rotatedshifted, index);
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
  auto coeff_modulus_size = param.coeff_modulus().size();
  auto poly_modulus_degree = ckks->degree;

  const phantom::util::cuda_stream_wrapper &stream_wrapper = *phantom::util::global_variables::default_stream;
  const auto &stream = stream_wrapper.get_stream();

  PhantomPlaintext p;
  p.resize(coeff_modulus_size, poly_modulus_degree * 2, stream);

  // int block_size = blockDimGlb.x;
  // int num_blocks = (poly_modulus_degree + block_size - 1) / block_size;

  // auto gpu_val = make_cuda_auto_ptr<double>(val.size(), stream);

  // cudaMemcpyAsync(gpu_val.get(), val.data(), val.size() * sizeof(double), cudaMemcpyHostToDevice, stream);

  // expand_encode_kernel<<<num_blocks, block_size>>>(
  //     gpu_val.get(), poly_modulus_degree, , p.data());

  cudaStreamSynchronize(stream);

  uint64_t *p_data = new uint64_t[poly_modulus_degree * 2 * coeff_modulus_size];
  cudaMemcpy(p_data, p.data(), poly_modulus_degree * 2 * coeff_modulus_size * sizeof(uint64_t), cudaMemcpyDeviceToHost);

  for (auto i = 0; i < poly_modulus_degree; i++) {
    auto coeffd = std::round(val[i] * 10000000000);
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

  cudaMemcpy(p.data(), p_data, poly_modulus_degree * 2 * coeff_modulus_size * sizeof(uint64_t), cudaMemcpyHostToDevice);

  for (std::size_t i = 0; i < 2; i++) {
    nwt_2d_radix8_forward_inplace(p.data(i * poly_modulus_degree), ckks->context->gpu_rns_tables(), coeff_modulus_size, i, stream);
  }

  p.set_chain_index(ckks->context->get_first_index());
  p.scale() = 10000000000;

  zero.scale() = p.scale();

  ckks->evaluator.add_plain(zero, p, ct);
}

void MMEvaluator::matrix_mul(vector<vector<double>> &x, vector<vector<double>> &y, vector<PhantomCiphertext> &res) {
  vector<PhantomPlaintext> a_pts;
  a_pts.reserve(768);
  for (int i = 0; i < 768; i++) {
    PhantomPlaintext pt;
    ckks->encoder.encode(x[i], ckks->scale, pt);
    a_pts.emplace_back(pt);
  }

  vector<PhantomCiphertext> b_compressed_cts;
  for (int i = 0; i < 768 * 64 / ckks->degree; i++) {
    PhantomPlaintext pt;
    PhantomCiphertext ct;
    expand_encode(y[i], ct);
    b_compressed_cts.push_back(ct);
  }

  cout << "1" << endl;

  auto timer = Timer();
  vector<PhantomCiphertext> b_expanded_cts;
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

  timer = Timer();
  PhantomCiphertext temp;

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
  cout << "calculating res time: " << timer.duration<seconds>() << " seconds" << endl;
}
