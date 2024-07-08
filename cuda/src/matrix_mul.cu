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
  auto base_rns = context->gpu_rns_tables().modulus();
  auto encrypted_count = encrypted.size();

  ckks->evaluator.transform_from_ntt_inplace(encrypted);

  destination = encrypted;

  // auto destination_data = new uint64_t[encrypted_count * coeff_count * destination.coeff_modulus_size()];
  // cudaMemcpy(destination_data, destination.data(), encrypted_count * coeff_count * coeff_mod_count * sizeof(uint64_t), cudaMemcpyDeviceToHost);
  // for (auto i = 0; i < 4; i++) {
  //   cout << fixed << setprecision(5) << destination_data[i] << ", ";
  // }
  // cout << "..., ";
  // for (auto i = coeff_count - 4; i < coeff_count; i++) {
  //   cout << fixed << setprecision(5) << destination_data[i] << ", ";
  // }
  // cout << endl;

  for (int i = 0; i < encrypted_count; i++) {
    // TODO: fix me
    ckks->evaluator.negacyclic_shift_poly_coeffmod(
        encrypted.data(i),
        coeff_count,
        index,
        base_rns,
        coeff_mod_count,
        destination.data(i));
  }

  // auto destination_data_2 = new uint64_t[encrypted_count * coeff_count * encrypted.coeff_modulus_size()];
  // cudaMemcpy(destination_data_2, destination.data(), encrypted_count * coeff_count * coeff_mod_count * sizeof(uint64_t), cudaMemcpyDeviceToHost);
  // for (auto i = 0; i < 4; i++) {
  //   cout << fixed << setprecision(5) << destination_data_2[i] << ", ";
  // }
  // cout << "..., ";
  // for (auto i = coeff_count - 4; i < coeff_count; i++) {
  //   cout << fixed << setprecision(5) << destination_data_2[i] << ", ";
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
  auto moduli = ckks->context->gpu_rns_tables().modulus();

  const phantom::util::cuda_stream_wrapper &stream_wrapper = *phantom::util::global_variables::default_stream;
  const auto &stream = stream_wrapper.get_stream();

  PhantomPlaintext p;
  p.resize(coeff_modulus_size, poly_modulus_degree * 2, stream);

  // Copy vector values to the GPU
  auto val_gpu = make_cuda_auto_ptr<double>(val.size(), stream);
  cudaMemcpy(val_gpu.get(), val.data(), val.size(), cudaMemcpyHostToDevice);

  // Execute expand kernel
  uint64_t block_size = blockDimGlb.x;
  int num_blocks = (poly_modulus_degree + block_size - 1) / block_size;
  expand_encode_kernel<<<num_blocks, block_size, 0, stream>>>(val_gpu.get(), poly_modulus_degree, moduli, p.data());

  // To NTT
  for (std::size_t i = 0; i < 2; i++) {
    nwt_2d_radix8_forward_inplace(p.data(i * poly_modulus_degree), ckks->context->gpu_rns_tables(), coeff_modulus_size, 0, stream);
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

  auto timer = Timer();
  vector<PhantomCiphertext> b_expanded_cts;
  for (auto i = 0; i < b_compressed_cts.size(); i++) {
    vector<PhantomCiphertext> temp_cts =
        expand_ciphertext(b_compressed_cts[i], ckks->degree, *ckks->galois_keys, ckks->rots);
    // cout << "Expanding..." << endl;
    b_expanded_cts.insert(
        b_expanded_cts.end(), make_move_iterator(temp_cts.begin()), make_move_iterator(temp_cts.end()));
  }
  timer.stop();
  cout << "Expanding time: " << timer.duration<seconds>() << " seconds" << endl;

  vector<double> b0_expanded_output;
  PhantomPlaintext b0_expanded;

  ckks->decryptor.decrypt(b_expanded_cts[0], b0_expanded);
  ckks->encoder.decode(b0_expanded, b0_expanded_output);

  for (auto i = 0; i < 10; i++) {
    cout << b0_expanded_output[i] << ", ";
  }
  cout << endl;

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
