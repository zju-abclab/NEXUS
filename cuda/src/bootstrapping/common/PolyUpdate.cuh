#pragma once

#include <NTL/RR.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "MinicompFunc.cuh"

using namespace std;
using namespace NTL;
using namespace minicomp;

enum class evaltype : int {
  // no evaluation method set; cannot be used for evaluating
  none = 0,

  // upgrade odd-baby + optimal level consumption
  oddbaby = 1,

  // upgrade baby + optimal level consumption
  baby = 2,
};

namespace minicomp {
class Tree {
 public:
  int depth;
  evaltype type;
  vector<int> tree;
  int m;
  int l;
  int b;

 public:
  Tree();
  Tree(evaltype ty);
  Tree(Tree a, Tree b, int g);
  void clear();
  void print();
  void merge(Tree a, Tree b, int g);
};

class Polynomial {
 public:
  vector<RR> coeff;
  long deg;
  vector<RR> chebcoeff;

  Polynomial();
  Polynomial(long _deg);
  Polynomial(long _deg, vector<RR> _coeff, string tag);
  Polynomial(long _deg, RR *_coeff, string tag);

  void get_coeff(vector<RR> &_coeff);
  void showcoeff();
  void showchebcoeff();
  void chebround(Polynomial &poly, long n);
  void copy(Polynomial poly);
  void power_to_cheb();
  void cheb_to_power();
  void power_to_cheb_scale(RR scale);
  void cheb_to_power_scale(RR scale);
  RR evaluate(RR input);
  RR evaluate_cheb(RR input);
};

void mul(Polynomial &rtn, Polynomial &a, Polynomial &b);
void mul(Polynomial &rtn, Polynomial &a, RR b);
void add(Polynomial &rtn, Polynomial &a, Polynomial &b);
void subt(Polynomial &rtn, Polynomial &a, Polynomial &b);

void mulinplace(Polynomial &a, Polynomial &b);
void addinplace(Polynomial &a, Polynomial &b);
void subtinplace(Polynomial &a, Polynomial &b);

void divide_poly(Polynomial &quotient, Polynomial &remainder, Polynomial &target, Polynomial &divider);
void chebyshev(Polynomial &rtn, long deg);
void chebyshev_scale(Polynomial &rtn, long deg, RR scale);

void power_to_cheb_int(Polynomial &q, int type, RR scale);
void cheb_to_power_int(Polynomial &qround, int type, RR scale);
void eval_divide(Polynomial &pround, Polynomial &qround, Polynomial &rround, Polynomial &Ti, int type, RR scale);
void geneTi(Polynomial &Ti, int deg, int type, RR scale);

// void poly_decompose_integrate(int deg, int type, RR scale, vector<RR> &coeff, Tree& tree, ofstream &out, evaltype eval);
/**
@param[in] deg The degree of root polynomial to evaluate
@param[in] input_type The input type. 0: power basis, 1: Cheb basis, 2: scaled Cheb basis
@param[in] input_scale The input scale K. If original Chebshev polynomial is used, K should be 1.
@param[in] coeff The coefficients of root polynomial based on basis polynomials (power, Cheb, scaled Cheb)
@param[in] tree The T-tree
@param[in] eval_type evaluation type. evaltype::none : not set, evaltype::oddbaby : level consumption optimization odd baby-giant, evaltype::baby : level consumption optimization baby-giant
@param[in] output_type The output type.
@param[in] output_scale The output scale K.
@param[in] out_decomp_coeff The tree-decomposed coefficients are stored in this ofstream.
*/
void poly_decompose_integrate(int deg, int input_type, RR input_scale, vector<RR> &coeff, Tree &tree, evaltype eval_type, int output_type, RR output_scale, ofstream &out_decomp_coeff);
void print_text_chebcoeff(ofstream &out, Polynomial &q);

}  // namespace minicomp
