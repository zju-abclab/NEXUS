#include "matrix_mul.h"

#include <seal/ciphertext.h>
#include <seal/plaintext.h>
#include <seal/util/defines.h>
#include <seal/valcheck.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "seal/util/polyarithsmallmod.h"

using namespace std;
using namespace seal;
using namespace std::chrono;
using namespace seal::util;

std::vector<std::vector<double>> MMEvaluator::transposeMatrix(const std::vector<std::vector<double>> &matrix) {
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

std::vector<std::vector<double>> MMEvaluator::readMatrix(const std::string &filename, int rows, int cols) {
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

void MMEvaluator::enc_compress_ciphertext(vector<double> &val, Ciphertext &ct) {
  Plaintext zero_pt;
  ckks->encoder->encode(std::vector<double>(ckks->N / 2, 0.0), ckks->scale, zero_pt);
  Ciphertext zero;
  ckks->encryptor->encrypt(zero_pt, zero);

  auto context = *ckks->context;
  auto context_data = context.get_context_data(context.first_parms_id());
  auto param = context_data->parms();
  auto ntt_tables = context_data->small_ntt_tables();

  auto poly_modulus_degree = ckks->degree;

  Plaintext p(poly_modulus_degree * 2);

  for (auto i = 0; i < poly_modulus_degree; i++) {
    auto coeffd = std::round(val[i] * 10000000000);
    bool is_negative = std::signbit(coeffd);
    auto coeffu = static_cast<std::uint64_t>(std::fabs(coeffd));
    if (is_negative) {
      for (std::size_t j = 0; j < 2; j++) {
        p[i + (j * poly_modulus_degree)] = util::negate_uint_mod(
            util::barrett_reduce_64(coeffu, param.coeff_modulus()[j]), param.coeff_modulus()[j]);
      }
    } else {
      for (std::size_t j = 0; j < 2; j++) {
        p[i + (j * poly_modulus_degree)] = util::barrett_reduce_64(coeffu, param.coeff_modulus()[j]);
      }
    }
  }

  for (std::size_t i = 0; i < 2; i++) {
    util::ntt_negacyclic_harvey(p.data(i * poly_modulus_degree), ntt_tables[i]);
  }

  p.parms_id() = context.first_parms_id();
  p.scale() = 10000000000;

  zero.scale() = p.scale();

  ckks->evaluator->add_plain(zero, p, ct);
}

vector<Ciphertext> MMEvaluator::expand_ciphertext(const Ciphertext &encrypted, uint32_t m, GaloisKeys &galkey, vector<uint32_t> &galois_elts) {
  uint32_t logm = ceil(log2(m));
  auto n = ckks->N;

  vector<Ciphertext> temp;
  temp.push_back(encrypted);

  Ciphertext tempctxt_rotated;
  Ciphertext tempctxt_shifted;
  Ciphertext tempctxt_rotatedshifted;

  for (uint32_t i = 0; i < logm; i++) {
    vector<Ciphertext> newtemp(temp.size() << 1);
    int index_raw = (n << 1) - (1 << i);
    int index = (index_raw * galois_elts[i]) % (n << 1);

    for (uint32_t a = 0; a < temp.size(); a++) {
      ckks->evaluator->apply_galois(temp[a], ckks->rots[i], *(ckks->galois_keys), tempctxt_rotated);  // sub
      ckks->evaluator->add(temp[a], tempctxt_rotated, newtemp[a]);
      multiply_power_of_x(temp[a], tempctxt_shifted, index_raw);  // x**-1
      // if (temp.size() == 1) ckks->print_decrypted_ct(tempctxt_shifted, 10);
      multiply_power_of_x(tempctxt_rotated, tempctxt_rotatedshifted, index);
      // if (temp.size() == 1) ckks->print_decrypted_ct(tempctxt_rotatedshifted, 10);
      ckks->evaluator->add(tempctxt_shifted, tempctxt_rotatedshifted, newtemp[a + temp.size()]);
    }
    temp = newtemp;
  }
  return temp;
}

void MMEvaluator::multiply_power_of_x(Ciphertext &encrypted, Ciphertext &destination, int index) {
  // Method 1:
  // string s = "";
  // destination = encrypted;
  // while (index >= ckks->N - 1) {
  //     s = "1x^" + to_string(ckks->N - 1);
  //     Plaintext p(s);
  //     ckks->evaluator->multiply_plain(destination, p, destination);
  //     index -= ckks->N - 1;
  // }

  // s = "1x^" + to_string(index);

  // Plaintext p(s);
  // ckks->evaluator->multiply_plain(destination, p, destination);

  // Method 2:
  auto context = *ckks->context;
  auto context_data = context.get_context_data(context.first_parms_id());
  auto param = context_data->parms();

  ckks->evaluator->transform_from_ntt_inplace(encrypted);
  auto coeff_mod_count = param.coeff_modulus().size();
  auto coeff_count = ckks->degree;
  auto encrypted_count = encrypted.size();

  destination = encrypted;

  for (int i = 0; i < encrypted_count; i++) {
    for (int j = 0; j < coeff_mod_count; j++) {
      negacyclic_shift_poly_coeffmod(
          encrypted.data(i) + (j * coeff_count),
          coeff_count,
          index,
          param.coeff_modulus()[j],
          destination.data(i) + (j * coeff_count));
    }
  }

  ckks->evaluator->transform_to_ntt_inplace(encrypted);
  ckks->evaluator->transform_to_ntt_inplace(destination);
}

void MMEvaluator::matrix_mul(vector<vector<double>> &x, vector<vector<double>> &y, vector<Ciphertext> &res) {
  chrono::high_resolution_clock::time_point time_start, time_end;

  vector<Plaintext> a_pts;
  a_pts.reserve(768);
  for (int i = 0; i < 768; i++) {
    Plaintext pt;
    ckks->encoder->encode(x[i], ckks->scale, pt);
    a_pts.emplace_back(pt);
  }

  vector<Ciphertext> b_compressed_cts;
  for (int i = 0; i < 768 * 64 / ckks->degree; i++) {
    Ciphertext ct;
    enc_compress_ciphertext(y[i], ct);
    b_compressed_cts.push_back(ct);
  }

  vector<seal::seal_byte> ct_bytes(b_compressed_cts[0].save_size());
  auto send_size = 0;
  for (auto &ct : b_compressed_cts) {
    auto ctt = ckks->encryptor->encrypt(a_pts[0]);
    send_size += ctt.save(ct_bytes.data(), ct_bytes.size());
  }
  cout << send_size / 1024.0 / 1024.0 << " MB" << endl;

  time_start = high_resolution_clock::now();
  vector<Ciphertext> b_expanded_cts;

  for (auto i = 0; i < b_compressed_cts.size(); i++) {
    vector<Ciphertext> temp_cts = expand_ciphertext(b_compressed_cts[i], ckks->degree, *ckks->galois_keys, ckks->rots);
    cout << "Expanded ciphertext #" << i + 1 << endl;
    b_expanded_cts.insert(b_expanded_cts.end(), make_move_iterator(temp_cts.begin()), make_move_iterator(temp_cts.end()));
  }

  time_end = high_resolution_clock::now();
  cout << "Expanding time: " << duration_cast<std::chrono::seconds>(time_end - time_start).count() << " seconds"
       << endl;

  time_start = high_resolution_clock::now();
  Ciphertext temp;

  for (int i = 0; i < 64; i++) {
    Ciphertext res_col_ct;
    vector<Ciphertext> temp_cts(768);

    for (int j = 0; j < 768; j++) {
      ckks->evaluator->multiply_plain(b_expanded_cts[i * 768 + j], a_pts[j], temp_cts[j]);
    }

    res_col_ct.scale() = temp_cts[0].scale();
    ckks->evaluator->add_many(temp_cts, res_col_ct);

    res_col_ct.scale() *= 4096;
    res.push_back(res_col_ct);
  }

  for (auto &ct : res) {
    while (ct.coeff_modulus_size() > 1) {
      ckks->evaluator->rescale_to_next_inplace(ct);
    }
  }

  vector<seal::seal_byte> rece_bytes(res[0].save_size());
  auto rece_size = 0;
  for (auto &ct : res) {
    rece_size += ct.save(rece_bytes.data(), rece_bytes.size());
  }
  cout << rece_size / 1024.0 / 1024.0 << " MB" << endl;

  time_end = high_resolution_clock::now();
  cout << "Result calculation time: " << duration_cast<seconds>(time_end - time_start).count() << " seconds" << endl;
}
