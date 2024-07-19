#include "Polynomial.cuh"

namespace boot {
Polynomial::Polynomial(long _deg) {
  deg = _deg;
  coeff = new RR[deg + 1];
  chebcoeff = new RR[deg + 1];
  for (int i = 0; i < deg + 1; i++) {
    coeff[i] = RR(0);
    chebcoeff[i] = RR(0);
  }
}

Polynomial::Polynomial(long _deg, RR *_coeff, string tag) {
  deg = _deg;
  coeff = new RR[deg + 1];
  chebcoeff = new RR[deg + 1];
  if (tag == "power") {
    for (int i = 0; i < deg + 1; i++) {
      coeff[i] = _coeff[i];
    }
    power_to_cheb();
  }

  else if (tag == "cheb") {
    for (int i = 0; i < deg + 1; i++) {
      chebcoeff[i] = _coeff[i];
    }
    cheb_to_power();
  }
}

Polynomial::~Polynomial() {
  if (coeff) delete[] coeff;
  if (chebcoeff) delete[] chebcoeff;
  if (poly_heap) {
    for (int i = 0; i < heaplen; i++) {
      if (poly_heap[i]) delete poly_heap[i];
    }
  }
}

void Polynomial::set_polynomial(long _deg, RR *_coeff, string tag) {
  deg = _deg;
  coeff = new RR[deg + 1];
  chebcoeff = new RR[deg + 1];
  if (tag == "power") {
    for (int i = 0; i < deg + 1; i++) {
      coeff[i] = _coeff[i];
    }
    power_to_cheb();
  }

  else if (tag == "cheb") {
    for (int i = 0; i < deg + 1; i++) {
      chebcoeff[i] = _coeff[i];
    }
    cheb_to_power();
  }
}

void Polynomial::set_zero_polynomial(long _deg) {
  deg = _deg;
  coeff = new RR[deg + 1]{};
  chebcoeff = new RR[deg + 1]{};
}

void Polynomial::showcoeff() {
  for (int i = 0; i < deg + 1; i++) {
    cout << "term " << i << " : " << coeff[i] << endl;
  }
  cout << endl;
}

void Polynomial::showchebcoeff() {
  for (int i = 0; i < deg + 1; i++) {
    cout << "chebterm " << i << " : " << chebcoeff[i] << endl;
  }
  cout << endl;
}
void Polynomial::copy(Polynomial &poly) {
  deg = poly.deg;
  coeff = new RR[deg + 1];
  chebcoeff = new RR[deg + 1];
  for (int i = 0; i < deg + 1; i++) {
    coeff[i] = poly.coeff[i];
    chebcoeff[i] = poly.chebcoeff[i];
  }
}

void Polynomial::power_to_cheb() {
  Polynomial chebbasis, tmp(deg);
  for (int i = 0; i <= deg; i++) {
    tmp.coeff[i] = coeff[i];
  }

  for (int i = deg; i >= 0; i--) {
    chebyshev(chebbasis, i);
    chebcoeff[i] = tmp.coeff[i] / chebbasis.coeff[i];
    for (int j = 0; j <= i; j++) {
      chebbasis.coeff[j] *= chebcoeff[i];
    }
    tmp.subtinplace(chebbasis);
  }
}

void Polynomial::cheb_to_power() {
  Polynomial chebbasis, tmp(deg);
  for (int i = 0; i <= deg; i++) {
    chebyshev(chebbasis, i);
    for (int j = 0; j <= i; j++) {
      chebbasis.coeff[j] *= chebcoeff[i];
    }
    tmp.addinplace(chebbasis);
  }

  for (int i = 0; i <= deg; i++) {
    coeff[i] = tmp.coeff[i];
  }
}

RR Polynomial::evaluate(RR value) {
  RR rtn = RR(0), term = RR(1);
  for (int i = 0; i <= deg; i++) {
    rtn += coeff[i] * term;
    term *= value;
  }

  return rtn;
}

void Polynomial::constmul(RR constant) {
  for (int i = 0; i <= deg; i++) {
    coeff[i] *= constant;
    chebcoeff[i] *= constant;
  }
}

void Polynomial::mulinplace(Polynomial &poly) {
  Polynomial rtn;
  mul(rtn, (*this), poly);
  copy(rtn);
}

void Polynomial::addinplace(Polynomial &poly) {
  Polynomial rtn;
  add(rtn, (*this), poly);
  copy(rtn);
}

void Polynomial::subtinplace(Polynomial &poly) {
  Polynomial rtn;
  subt(rtn, (*this), poly);
  copy(rtn);
}

void Polynomial::change_variable_scale(RR scale) {
  RR div = RR(1);
  for (int i = 0; i <= deg; i++) {
    coeff[i] /= div;
    div *= scale;
  }
  power_to_cheb();
}

void Polynomial::generate_poly_heap_manual(long k, long m) {
  heap_k = k;
  heap_m = m;
  heaplen = (1 << (heap_m + 1)) - 1;
  poly_heap = new Polynomial *[heaplen];
  for (long i = 0; i < heaplen; i++) {
    poly_heap[i] = 0;
  }

  poly_heap[0] = new Polynomial();
  poly_heap[0]->copy((*this));

  long chebdeg = heap_k << heap_m;
  long first, last;
  Polynomial chebgiant;
  for (int i = 0; i < heap_m; i++) {
    chebdeg = chebdeg >> 1;
    first = (1 << i) - 1;
    last = (1 << (i + 1)) - 1;
    chebyshev(chebgiant, chebdeg);
    for (int j = first; j < last; j++) {
      if (poly_heap[j]) {
        if (poly_heap[j]->deg < chebdeg) {
          poly_heap[2 * (j + 1)] = new Polynomial();
          poly_heap[2 * (j + 1)]->copy(*poly_heap[j]);
        }

        else {
          poly_heap[2 * (j + 1) - 1] = new Polynomial();
          poly_heap[2 * (j + 1)] = new Polynomial();
          divide_poly(*poly_heap[2 * (j + 1) - 1], *poly_heap[2 * (j + 1)], *poly_heap[j], chebgiant);
          poly_heap[2 * (j + 1) - 1]->power_to_cheb();
          poly_heap[2 * (j + 1)]->power_to_cheb();
        }
      }
    }
  }
}

void Polynomial::generate_poly_heap_odd() {
  oddbabycount(heap_k, heap_m, deg);
  generate_poly_heap_manual(heap_k, heap_m);
}

void Polynomial::generate_poly_heap() {
  babycount(heap_k, heap_m, deg);
  generate_poly_heap_manual(heap_k, heap_m);
}

void Polynomial::write_heap_to_file(ofstream &out) {
  RR::SetOutputPrecision(17);
  out << heaplen << endl;

  for (int index = 0; index < heaplen; index++) {
    if (poly_heap[index]) {
      out << index << " " << poly_heap[index]->deg << endl;
      for (int i = 0; i <= poly_heap[index]->deg; i++) {
        out << poly_heap[index]->chebcoeff[i] << endl;
      }
      out << endl;
    }
  }
}

void Polynomial::read_heap_from_file(ifstream &in) {
  long index = 0, in_deg;

  in >> heaplen;

  poly_heap = new Polynomial *[heaplen];

  for (int i = 0; i < heaplen; i++) {
    poly_heap[i] = 0;
  }

  while (index < heaplen - 1) {
    in >> index >> in_deg;
    poly_heap[index] = new Polynomial(in_deg);
    for (int i = 0; i <= in_deg; i++) {
      in >> poly_heap[index]->chebcoeff[i];
    }
    poly_heap[index]->cheb_to_power();
  }

  copy(*poly_heap[0]);
}

void Polynomial::homomorphic_poly_evaluation(CKKSEvaluator *ckks, PhantomCiphertext &rtn, PhantomCiphertext &cipher) {
  double zero = 1.0 / cipher.scale();

  if (deg == 1) {
    ckks->evaluator.multiply_const(cipher, to_double(coeff[1]), rtn);
    ckks->evaluator.rescale_to_next_inplace(rtn);
    ckks->evaluator.add_const(rtn, to_double(coeff[0]), rtn);
    return;
  }

  else if (deg == 2) {
    PhantomCiphertext squared;
    ckks->evaluator.square(cipher, squared);
    ckks->evaluator.relinearize_inplace(squared, *(ckks->relin_keys));
    ckks->evaluator.rescale_to_next_inplace(squared);

    ckks->evaluator.multiply_const_inplace(squared, to_double(coeff[2]));
    ckks->evaluator.rescale_to_next_inplace(squared);

    if (abs(to_double(coeff[1])) >= zero) {
      ckks->evaluator.multiply_const(cipher, to_double(coeff[1]), rtn);
      ckks->evaluator.rescale_to_next_inplace(rtn);
      ckks->evaluator.add_reduced_error(rtn, squared, rtn);
    } else {
      rtn = squared;
    }

    ckks->evaluator.add_const_inplace(rtn, to_double(coeff[0]));
    return;
  }

  else if (deg == 3) {
    PhantomCiphertext squared, cubic;
    ckks->evaluator.square(cipher, squared);
    ckks->evaluator.relinearize_inplace(squared, *(ckks->relin_keys));
    ckks->evaluator.rescale_to_next_inplace(squared);

    ckks->evaluator.multiply_const(cipher, to_double(coeff[3]), cubic);
    ckks->evaluator.rescale_to_next_inplace(cubic);

    ckks->evaluator.multiply_inplace_reduced_error(cubic, squared, *(ckks->relin_keys));
    ckks->evaluator.rescale_to_next_inplace(cubic);

    if (abs(to_double(coeff[1])) >= zero) {
      ckks->evaluator.multiply_const(cipher, to_double(coeff[1]), rtn);
      ckks->evaluator.rescale_to_next_inplace(rtn);
      ckks->evaluator.add_reduced_error(rtn, cubic, rtn);
    } else {
      rtn = cubic;
    }

    if (abs(to_double(coeff[2])) >= zero) {
      ckks->evaluator.multiply_const_inplace(squared, to_double(coeff[2]));
      ckks->evaluator.rescale_to_next_inplace(squared);
      ckks->evaluator.add_reduced_error(rtn, squared, rtn);
    }

    ckks->evaluator.add_const_inplace(rtn, to_double(coeff[0]));
    return;
  }

  else {
    vector<PhantomCiphertext> baby(heap_k, PhantomCiphertext());
    vector<bool> babybool(heap_k, false);

    PhantomCiphertext addtmp1, addtmp2;

    PhantomPlaintext tmpplain;

    baby[1] = cipher;
    babybool[1] = true;

    // cout << endl;
    // cout << "heap_k: " << heap_k << endl;  // 8

    // i = 2, 4, 8
    for (int i = 2; i < heap_k; i *= 2) {
      ckks->evaluator.square(baby[i / 2], baby[i]);
      ckks->evaluator.relinearize_inplace(baby[i], *(ckks->relin_keys));
      ckks->evaluator.rescale_to_next_inplace(baby[i]);
      ckks->evaluator.double_inplace(baby[i]);

      ckks->evaluator.add_const(baby[i], -1.0, baby[i]);
      babybool[i] = true;
    }

    // cout << "Baby step 1: " << endl;
    // for (int i = 0; i < baby.size(); i++) {
    //   cout << i << " ";
    //   ckks->print_decrypted_ct(baby[i], 10);
    // }
    // cout << endl;

    long lpow2, res, diff;
    PhantomCiphertext tmp;

    for (int i = 1; i < heap_k; i++) {
      if (!babybool[i]) {
        // i = 3, 5, 6, 7
        lpow2 = (1 << (int)floor(log(i) / log(2)));
        res = i - lpow2;
        diff = abs(lpow2 - res);

        // cout << "lpow2: " << lpow2 << " res: " << res << " diff: " << diff << endl;

        ckks->evaluator.multiply_reduced_error(baby[lpow2], baby[res], *(ckks->relin_keys), baby[i]);
        ckks->evaluator.rescale_to_next_inplace(baby[i]);
        ckks->evaluator.double_inplace(baby[i]);

        ckks->evaluator.sub_reduced_error(baby[i], baby[diff], baby[i]);
        babybool[i] = true;
      }
    }

    // cout << "Baby step 2: " << endl;
    // for (int i = 0; i < baby.size(); i++) {
    //   cout << i << " ";
    //   ckks->print_decrypted_ct(baby[i], 10);
    // }
    // cout << endl;

    vector<PhantomCiphertext> giant(heap_m, PhantomCiphertext());
    vector<bool> giantbool(heap_m, false);

    giantbool[0] = true;
    lpow2 = (1 << ((int)ceil(log(heap_k) / log(2)) - 1));
    res = heap_k - lpow2;
    diff = abs(lpow2 - res);

    if (res == 0) {
      giant[0] = baby[lpow2];
    }

    else if (diff == 0) {
      ckks->evaluator.square(baby[lpow2], giant[0]);
      ckks->evaluator.relinearize_inplace(giant[0], *(ckks->relin_keys));
      ckks->evaluator.rescale_to_next_inplace(giant[0]);
      ckks->evaluator.double_inplace(giant[0]);
      ckks->evaluator.add_const(giant[0], -1.0, giant[0]);
    }

    else {
      ckks->evaluator.multiply_reduced_error(baby[lpow2], baby[res], *(ckks->relin_keys), giant[0]);
      ckks->evaluator.rescale_to_next_inplace(giant[0]);
      ckks->evaluator.double_inplace(giant[0]);
      ckks->evaluator.sub_reduced_error(giant[0], baby[diff], giant[0]);
    }

    for (int i = 1; i < heap_m; i++) {
      giantbool[i] = true;

      ckks->evaluator.square(giant[i - 1], giant[i]);
      ckks->evaluator.relinearize_inplace(giant[i], *(ckks->relin_keys));
      ckks->evaluator.rescale_to_next_inplace(giant[i]);
      ckks->evaluator.double_inplace(giant[i]);

      ckks->evaluator.add_const_inplace(giant[i], -1.0);
    }

    // cout << "giant: " << endl;
    // for (auto &giant_ct : giant) {
    //   ckks->print_decrypted_ct(giant_ct, 10);
    // }
    // cout << endl;

    vector<PhantomCiphertext> cipherheap((1 << (heap_m + 1)) - 1, PhantomCiphertext());
    vector<bool> cipherheapbool((1 << (heap_m + 1)) - 1, false);

    long heapfirst = (1 << heap_m) - 1;
    long heaplast = (1 << (heap_m + 1)) - 1;
    double zero = 1.0 / cipher.scale();

    for (int i = heapfirst; i < heaplast; i++) {
      if (poly_heap[i]) {
        cipherheapbool[i] = true;

        ckks->evaluator.multiply_const(baby[1], to_double(poly_heap[i]->chebcoeff[1]), cipherheap[i]);
        ckks->evaluator.rescale_to_next_inplace(cipherheap[i]);

        if (!(abs(to_double(poly_heap[i]->chebcoeff[1])) <= zero))
          ckks->evaluator.add_const_inplace(cipherheap[i], to_double(poly_heap[i]->chebcoeff[0]));

        for (int j = 2; j <= poly_heap[i]->deg; j++) {
          if (abs(to_double(poly_heap[i]->chebcoeff[j])) <= zero)
            continue;

          if (j < heap_k) {
            ckks->evaluator.multiply_const(baby[j], to_double(poly_heap[i]->chebcoeff[j]), tmp);
            ckks->evaluator.rescale_to_next_inplace(tmp);
          } else {
            ckks->evaluator.multiply_const(giant[0], to_double(poly_heap[i]->chebcoeff[j]), tmp);
            ckks->evaluator.rescale_to_next_inplace(tmp);
          }

          ckks->evaluator.add_reduced_error(cipherheap[i], tmp, cipherheap[i]);
        }
      }
    }

    long depth = heap_m;
    long gindex = 0;

    while (depth != 0) {
      depth--;
      heapfirst = (1 << depth) - 1;
      heaplast = (1 << (depth + 1)) - 1;
      for (int i = heapfirst; i < heaplast; i++) {
        if (poly_heap[i]) {
          cipherheapbool[i] = true;
          if (!cipherheapbool[2 * (i + 1) - 1]) {
            cipherheap[i] = cipherheap[2 * (i + 1)];
          } else {
            ckks->evaluator.multiply_reduced_error(cipherheap[2 * (i + 1) - 1], giant[gindex], *(ckks->relin_keys), cipherheap[i]);
            ckks->evaluator.rescale_to_next_inplace(cipherheap[i]);
            ckks->evaluator.add_reduced_error(cipherheap[i], cipherheap[2 * (i + 1)], cipherheap[i]);
          }
        }
      }
      gindex++;
    }

    rtn = cipherheap[0];
  }
}

void mul(Polynomial &rtn, Polynomial &a, Polynomial &b) {
  rtn.set_zero_polynomial(a.deg + b.deg);
  for (long i = 0; i <= rtn.deg; i++) {
    for (long j = 0; j <= i; j++) {
      if (j <= a.deg && i - j <= b.deg) {
        rtn.coeff[i] += a.coeff[j] * b.coeff[i - j];
      }
    }
  }
}

void add(Polynomial &rtn, Polynomial &a, Polynomial &b) {
  if (a.deg >= b.deg) {
    rtn.copy(a);
    for (int i = 0; i <= b.deg; i++) {
      rtn.coeff[i] += b.coeff[i];
    }
  }

  else {
    rtn.copy(b);
    for (int i = 0; i < a.deg; i++) {
      rtn.coeff[i] += a.coeff[i];
    }
  }
}

void subt(Polynomial &rtn, Polynomial &a, Polynomial &b) {
  if (a.deg >= b.deg) {
    rtn.copy(a);
    for (int i = 0; i <= b.deg; i++) {
      rtn.coeff[i] -= b.coeff[i];
    }
  }

  else {
    rtn.copy(b);
    for (int i = 0; i < b.deg; i++) {
      rtn.coeff[i] *= RR(-1);
    }

    for (int i = 0; i < a.deg; i++) {
      rtn.coeff[i] += a.coeff[i];
    }
  }
}

void divide_poly(Polynomial &quotient, Polynomial &remainder, Polynomial &target, Polynomial &divider) {
  if (target.deg < divider.deg) {
    quotient.set_zero_polynomial(0);
    quotient.coeff[0] = 0;

    remainder.copy(target);
  }

  else {
    quotient.set_zero_polynomial(target.deg - divider.deg);

    Polynomial tmp(target.deg, target.coeff, "power");
    RR currcoeff;
    long currdeg = target.deg;
    for (long i = quotient.deg; i >= 0; i--) {
      currcoeff = tmp.coeff[currdeg] / divider.coeff[divider.deg];
      quotient.coeff[i] = currcoeff;
      tmp.coeff[currdeg] = 0;
      for (int j = 0; j < divider.deg; j++) {
        tmp.coeff[j + i] -= divider.coeff[j] * currcoeff;
      }
      currdeg--;
    }
    remainder.set_polynomial(divider.deg - 1, tmp.coeff, "power");
  }
}

void chebyshev(Polynomial &rtn, long deg) {
  if (deg == 0) {
    rtn.set_zero_polynomial(0);
    rtn.coeff[0] = RR(1);
  }

  else if (deg == 1) {
    rtn.set_zero_polynomial(1);
    rtn.coeff[1] = RR(1);
  }

  else {
    Polynomial iden2(1);
    iden2.coeff[1] = RR(2);

    Polynomial tmp1(0), tmp2(1), tmp3;

    tmp1.coeff[0] = RR(1);
    tmp2.coeff[1] = RR(1);

    for (int i = 2; i <= deg; i++) {
      mul(tmp3, iden2, tmp2);
      subt(rtn, tmp3, tmp1);
      tmp1.copy(tmp2);
      tmp2.copy(rtn);
    }
  }
}

void second_chebyshev_times_x_for_sine(Polynomial &rtn, long deg) {
  if (deg == 1) {
    rtn.set_zero_polynomial(1);
    rtn.coeff[1] = RR(1);
  }

  else if (deg == 3) {
    rtn.set_zero_polynomial(3);
    rtn.coeff[1] = RR(3);
    rtn.coeff[3] = RR(-4);
  }

  else if (deg % 2 == 1 && deg > 0) {
    Polynomial mul_fac(2);
    mul_fac.coeff[0] = RR(2);
    mul_fac.coeff[2] = RR(-4);

    Polynomial iden(1);
    iden.coeff[1] = RR(1);

    Polynomial tmp1(0), tmp2(2), tmp3;

    tmp1.coeff[0] = RR(1);
    tmp2.coeff[0] = RR(3);
    tmp2.coeff[2] = RR(-4);

    for (int i = 4; i < deg; i += 2) {
      mul(tmp3, mul_fac, tmp2);
      subt(rtn, tmp3, tmp1);
      tmp1.copy(tmp2);
      tmp2.copy(rtn);
    }
    rtn.mulinplace(iden);
  }
}
}  // namespace boot
