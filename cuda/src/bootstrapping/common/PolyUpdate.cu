#include <NTL/LLL.h>
#include <NTL/RR.h>
#include <NTL/mat_RR.h>

#include "PolyUpdate.cuh"

namespace minicomp {
// n : odd
Tree::Tree() {
  depth = 0;
  type = evaltype::none;
  m = 0;
  l = 0;
  b = 0;
  tree.resize(2);
  tree[0] = -1;
  tree[1] = 0;
}

Tree::Tree(evaltype ty) {
  depth = 0;
  type = ty;
  m = 0;
  l = 0;
  b = 0;
  tree.resize(2);
  tree[0] = -1;
  tree[1] = 0;
}

Tree::Tree(Tree a, Tree b, int g) {
  if (a.type != b.type) throw std::invalid_argument("the types of two trees are not the same");

  if (a.depth > b.depth)
    depth = a.depth + 1;
  else
    depth = b.depth + 1;
  type = a.type;
  tree.resize(pow2(depth + 1), -1);
  tree[1] = g;

  for (int i = 1; i <= pow2(a.depth + 1) - 1; i++) {
    int temp = pow2(static_cast<int>(log(static_cast<double>(i)) / log(2.0)));
    tree[i + temp] = a.tree[i];
  }
  for (int i = 1; i <= pow2(b.depth + 1) - 1; i++) {
    int temp = pow2(static_cast<int>(log(static_cast<double>(i)) / log(2.0)));
    tree[i + 2 * temp] = b.tree[i];
  }
}

void Tree::clear() {
  depth = 0;
  type = evaltype::none;
  tree.resize(2);
  tree[0] = -1;
  tree[1] = 0;
}

void Tree::print() {
  cout << "depth of tree: " << depth << endl;
  for (int i = 0; i <= depth; i++) {
    for (int j = pow2(i); j < pow2(i + 1); j++) {
      cout << tree[j] << " ";
    }
    cout << endl;
  }

  if (type == evaltype::oddbaby) {
    cout << "m: " << m << endl;
    cout << "l: " << l << endl;

    // number of nonscalar multiplications
    int nonscalar = m - 1 + pow2(l - 1) - 1;
    for (size_t i = 0; i < pow2(depth + 1); i++)
      if (tree[i] > 0) nonscalar++;
    cout << "nonscalar: " << nonscalar << endl;
    cout << endl;
  } else if (type == evaltype::baby) {
    cout << "m: " << m << endl;
    cout << "b: " << b << endl;

    // number of nonscalar multiplications
    int nonscalar = m + b - 2;

    for (size_t i = 0; i < pow2(depth + 1); i++)
      if (tree[i] > 0) nonscalar++;
    cout << "nonscalar: " << nonscalar << endl;
    cout << endl;
  }
}

// void Tree::babyprint()
// {
// 	cout << "depth of tree: " << depth << endl;
// 	for(int i=0; i<=depth; i++)
// 	{
// 		for(int j=pow2(i); j<pow2(i+1); j++)
// 		{
// 			cout << tree[j] << " ";
// 		}
// 		cout << endl;
// 	}
// 	cout << "m: " << m << endl;
// 	cout << "b: " << b << endl;

// 	// number of nonscalar multiplications
// 	int nonscalar = m+b-2;

// 	for(size_t i=0; i<pow2(depth+1); i++) if(tree[i]>0) nonscalar++;
// 	cout << "nonscalar: " << nonscalar << endl;
// 	cout << endl;
// }
void Tree::merge(Tree a, Tree b, int g) {
  clear();
  if (a.type != b.type) throw std::invalid_argument("the types of two trees are not the same");
  type = a.type;

  if (a.depth > b.depth)
    depth = a.depth + 1;
  else
    depth = b.depth + 1;
  tree.resize(pow2(depth + 1), -1);
  tree[1] = g;

  for (int i = 1; i <= pow2(a.depth + 1) - 1; i++) {
    int temp = pow2(static_cast<int>(log(static_cast<double>(i)) / log(2.0)));
    tree[i + temp] = a.tree[i];
  }
  for (int i = 1; i <= pow2(b.depth + 1) - 1; i++) {
    int temp = pow2(static_cast<int>(log(static_cast<double>(i)) / log(2.0)));
    tree[i + 2 * temp] = b.tree[i];
  }
}

Polynomial::Polynomial() {
  deg = -1;
}

Polynomial::Polynomial(long _deg) {
  deg = _deg;
  for (int i = 0; i < deg + 1; i++) {
    coeff.emplace_back(RR(0));
    chebcoeff.emplace_back(RR(0));
  }
}

Polynomial::Polynomial(long _deg, vector<RR> _coeff, string tag) {
  deg = _deg;
  if (tag == "power") {
    coeff = _coeff;
    for (int i = 0; i < deg + 1; i++) chebcoeff.emplace_back(RR(0));
  } else if (tag == "cheb") {
    chebcoeff = _coeff;
    for (int i = 0; i < deg + 1; i++) coeff.emplace_back(RR(0));
  }
}

Polynomial::Polynomial(long _deg, RR *_coeff, string tag) {
  deg = _deg;
  if (tag == "power") {
    for (int i = 0; i < deg + 1; i++) coeff.emplace_back(_coeff[i]);
    for (int i = 0; i < deg + 1; i++) chebcoeff.emplace_back(RR(0));
  } else if (tag == "cheb") {
    for (int i = 0; i < deg + 1; i++) chebcoeff.emplace_back(_coeff[i]);
    for (int i = 0; i < deg + 1; i++) coeff.emplace_back(RR(0));
  }
}

void Polynomial::get_coeff(vector<RR> &_coeff) {
  _coeff = coeff;
}

void Polynomial::showcoeff() {
  for (int i = 0; i < deg + 1; i++) {
    //	cout << "x^" << i << " : " << coeff[i] << endl;
    cout << coeff[i] << endl;
  }
  cout << endl;
}

void Polynomial::showchebcoeff() {
  for (int i = 0; i < deg + 1; i++) {
    cout << "term " << i << " : " << chebcoeff[i] << endl;
  }
  cout << endl;
}

void Polynomial::chebround(Polynomial &poly, long bitprec) {
  RR expprec;
  pow(expprec, static_cast<RR>(2.0), static_cast<RR>(bitprec));

  deg = poly.deg;
  coeff.clear();
  chebcoeff.clear();
  for (int i = 0; i < deg + 1; i++) {
    coeff.emplace_back(0);
    chebcoeff.emplace_back(floor(poly.chebcoeff[i] * expprec + static_cast<RR>(0.5)) / expprec);
  }
}

void Polynomial::copy(Polynomial poly) {
  deg = poly.deg;
  coeff.clear();
  chebcoeff.clear();
  for (int i = 0; i < deg + 1; i++) {
    coeff.emplace_back(poly.coeff[i]);
    chebcoeff.emplace_back(poly.chebcoeff[i]);
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
    subtinplace(tmp, chebbasis);
  }
}

void Polynomial::cheb_to_power() {
  Polynomial chebbasis, tmp(deg);
  for (int i = 0; i <= deg; i++) {
    chebyshev(chebbasis, i);
    for (int j = 0; j <= i; j++) {
      chebbasis.coeff[j] *= chebcoeff[i];
    }
    addinplace(tmp, chebbasis);
  }

  for (int i = 0; i <= deg; i++) {
    coeff[i] = tmp.coeff[i];
  }
}

void Polynomial::power_to_cheb_scale(RR scale) {
  Polynomial chebbasis, tmp(deg);
  for (int i = 0; i <= deg; i++) {
    tmp.coeff[i] = coeff[i];
  }

  for (int i = deg; i >= 0; i--) {
    chebyshev_scale(chebbasis, i, scale);
    chebcoeff[i] = tmp.coeff[i] / chebbasis.coeff[i];
    for (int j = 0; j <= i; j++) {
      chebbasis.coeff[j] *= chebcoeff[i];
    }
    subtinplace(tmp, chebbasis);
  }
}

void Polynomial::cheb_to_power_scale(RR scale) {
  Polynomial chebbasis, tmp(deg);
  for (int i = 0; i <= deg; i++) {
    chebyshev_scale(chebbasis, i, scale);
    for (int j = 0; j <= i; j++) {
      chebbasis.coeff[j] *= chebcoeff[i];
    }
    addinplace(tmp, chebbasis);
  }

  for (int i = 0; i <= deg; i++) {
    coeff[i] = tmp.coeff[i];
  }
}

// coeff should be set
RR Polynomial::evaluate(RR input) {
  RR rtn = RR(0);
  RR pow = RR(1);
  for (int i = 0; i <= deg; i++) {
    rtn += coeff[i] * pow;
    pow *= input;
  }
  return rtn;
}

RR Polynomial::evaluate_cheb(RR input) {
  RR rtn = RR(0), val = input;
  if (deg == 0)
    return RR(1);
  else {
    RR tmp1, tmp2, tmp3;
    tmp1 = 1;
    tmp2 = val;
    rtn = chebcoeff[0] * tmp1 + chebcoeff[1] * tmp2;
    for (int i = 2; i <= deg; i++) {
      tmp3 = RR(2.0) * val * tmp2 - tmp1;
      tmp1 = tmp2;
      tmp2 = tmp3;
      rtn += chebcoeff[i] * tmp3;
    }
  }
  return rtn;
}

void mul(Polynomial &rtn, Polynomial &a, Polynomial &b) {
  rtn.deg = a.deg + b.deg;
  rtn.coeff.resize(a.deg + b.deg + 1);
  rtn.chebcoeff.resize(a.deg + b.deg + 1, RR(0));
  for (long i = 0; i <= rtn.deg; i++) rtn.coeff[i] = RR(0);
  for (long i = 0; i <= rtn.deg; i++) {
    for (long j = 0; j <= i; j++) {
      if (j <= a.deg && i - j <= b.deg) rtn.coeff[i] += a.coeff[j] * b.coeff[i - j];
    }
  }
}

void mul(Polynomial &rtn, Polynomial &a, RR b) {
  rtn.deg = a.deg;
  rtn.coeff.resize(a.deg + 1);
  rtn.chebcoeff.resize(a.deg + 1, RR(0));
  for (long i = 0; i <= rtn.deg; i++) rtn.coeff[i] = a.coeff[i] * b;
}

void add(Polynomial &rtn, Polynomial &a, Polynomial &b) {
  if (a.deg >= b.deg) {
    rtn.copy(Polynomial(a.deg, a.coeff, "power"));
    for (int i = 0; i <= b.deg; i++) rtn.coeff[i] += b.coeff[i];
  } else {
    rtn.copy(Polynomial(b.deg, b.coeff, "power"));
    for (int i = 0; i <= a.deg; i++) rtn.coeff[i] += a.coeff[i];
  }
}

void subt(Polynomial &rtn, Polynomial &a, Polynomial &b) {
  if (a.deg >= b.deg) {
    rtn.copy(Polynomial(a.deg, a.coeff, "power"));
    for (int i = 0; i <= b.deg; i++) rtn.coeff[i] -= b.coeff[i];
  } else {
    rtn.copy(Polynomial(b.deg, b.coeff, "power"));
    for (int i = 0; i <= b.deg; i++) rtn.coeff[i] *= RR(-1);
    for (int i = 0; i <= a.deg; i++) rtn.coeff[i] += a.coeff[i];
  }
}

void mulinplace(Polynomial &a, Polynomial &b) {
  Polynomial rtn;
  mul(rtn, a, b);
  a.copy(rtn);
}

void addinplace(Polynomial &a, Polynomial &b) {
  Polynomial rtn;
  add(rtn, a, b);
  a.copy(rtn);
}

void subtinplace(Polynomial &a, Polynomial &b) {
  Polynomial rtn;
  subt(rtn, a, b);
  a.copy(rtn);
}

void divide_poly(Polynomial &quotient, Polynomial &remainder, Polynomial &target, Polynomial &divider) {
  if (target.deg < divider.deg) {
    quotient.copy(Polynomial(0));
    quotient.coeff[0] = 0;

    remainder.copy(Polynomial(target.deg, target.coeff, "power"));
  } else {
    quotient.copy(Polynomial(target.deg - divider.deg));
    Polynomial tmp(target.deg, target.coeff, "power");
    RR currcoeff;
    long currdeg = target.deg;
    for (long i = quotient.deg; i >= 0; i--) {
      currcoeff = tmp.coeff[currdeg] / divider.coeff[divider.deg];
      quotient.coeff[i] = currcoeff;
      tmp.coeff[currdeg] = 0;
      for (int j = 0; j < divider.deg; j++) tmp.coeff[j + i] -= divider.coeff[j] * currcoeff;
      currdeg--;
    }
    remainder.copy(Polynomial(divider.deg - 1, tmp.coeff, "power"));
  }
}

// coeff를 업데이트한다. 즉, power basis의 계수들을 구해준다. chebcoeff는 따로 구해주어야한다.
void chebyshev(Polynomial &rtn, long deg) {
  if (deg == 0) {
    rtn.copy(Polynomial(0));
    rtn.coeff[0] = RR(1);
  } else if (deg == 1) {
    rtn.copy(Polynomial(1));
    rtn.coeff[0] = RR(0);
    rtn.coeff[1] = RR(1);
  } else {
    Polynomial iden2(1);
    iden2.coeff[0] = RR(0);
    iden2.coeff[1] = RR(2);

    Polynomial tmp1(0), tmp2(1), tmp3;
    tmp1.coeff[0] = RR(1);
    tmp2.coeff[0] = RR(0);
    tmp2.coeff[1] = RR(1);

    for (int i = 2; i <= deg; i++) {
      mul(tmp3, iden2, tmp2);
      subt(rtn, tmp3, tmp1);
      tmp1.copy(tmp2);
      tmp2.copy(rtn);
    }
  }
}

// scale = K
void chebyshev_scale(Polynomial &rtn, long deg, RR scale) {
  if (deg == 0) {
    rtn.copy(Polynomial(0));
    rtn.coeff[0] = RR(1);
  } else if (deg == 1) {
    rtn.copy(Polynomial(1));
    rtn.coeff[0] = RR(0);
    rtn.coeff[1] = RR(1.0) / RR(scale);
  } else {
    Polynomial iden2(1);
    iden2.coeff[0] = RR(0);
    iden2.coeff[1] = RR(2.0) / RR(scale);

    Polynomial tmp1(0), tmp2(1), tmp3;
    tmp1.coeff[0] = RR(1);
    tmp2.coeff[0] = RR(0);
    tmp2.coeff[1] = RR(1) / RR(scale);

    for (int i = 2; i <= deg; i++) {
      mul(tmp3, iden2, tmp2);
      subt(rtn, tmp3, tmp1);
      tmp1.copy(tmp2);
      tmp2.copy(rtn);
    }
  }
}

void power_to_cheb_int(Polynomial &q, int type, RR scale) {
  if (type == 1)
    q.power_to_cheb();
  else if (type == 2)
    q.power_to_cheb_scale(scale);
}

void cheb_to_power_int(Polynomial &qround, int type, RR scale) {
  if (type == 1)
    qround.cheb_to_power();
  else if (type == 2)
    qround.cheb_to_power_scale(scale);
}

void eval_divide(Polynomial &pround, Polynomial &qround, Polynomial &rround, Polynomial &Ti, int type, RR scale) {
  Polynomial temp1;
  mul(temp1, qround, Ti);
  add(pround, temp1, rround);
  if (type == 1)
    pround.power_to_cheb();
  else if (type == 2)
    pround.power_to_cheb_scale(scale);
}

void geneTi(Polynomial &Ti, int deg, int type, RR scale) {
  if (type == 1) {
    chebyshev(Ti, deg);
    Ti.power_to_cheb();
  } else if (type == 2) {
    chebyshev_scale(Ti, deg, scale);
    Ti.power_to_cheb_scale(scale);
  }
}

void print_text_chebcoeff(ofstream &out, Polynomial &q) {
  for (int i = 0; i <= q.deg; i++) out << q.chebcoeff[i] << endl;
}
// void poly_decompose_integrate(int deg, int type, RR scale, vector<RR> &coeff, Tree& tree, ofstream &out, evaltype eval_type)
// {
// 	vector<Polynomial> T;
// 	for(int i=0; i<100; i++) T.emplace_back(Polynomial(0));

// 	long m = tree.m, l = tree.l, b = tree.b;
// 	if(eval_type == evaltype::baby)
// 	{
// 		for(int i=1; i<=b; i++)
// 		{
// 			T[i].copy(Polynomial(i));
// 			geneTi(T[i], i, type, scale);
// 		}
// 		for(int i=1; i<=m-1; i++)
// 		{
// 			T[pow2(i)*b].copy(Polynomial(pow2(i)*b));
// 			geneTi(T[pow2(i)*b], pow2(i)*b, type, scale);
// 		}
// 	}
// 	else if(eval_type == evaltype::oddbaby)	// comptype = 3
// 	{
// 		for(int i=0; i<=m-1; i++)
// 		{
// 			T[pow2(i)].copy(Polynomial(pow2(i)));
// 			geneTi(T[pow2(i)], pow2(i), type, scale);
// 		}
// 		for(int i=0; i<=l; i++)
// 		{
// 			for(int j=pow2(i-1)+1; j<=pow2(i)-1; j+=2)
// 			{
// 				T[j].copy(Polynomial(j));
// 				geneTi(T[j], j, type, scale);
// 			}
// 		}
// 	}
// 	else throw std::invalid_argument("evaluation type is not correct");

// //	long m = tree.m, l = tree.l;
// 	// generate Tis.  compute all Tis including needless Tis.

// 	// decompose setting
// 	vector<Polynomial> pt;
// 	for(int i=0; i<pow2(tree.depth+1); i++) pt.emplace_back(Polynomial(0));
// 	pt[1].copy(Polynomial(deg, coeff, "cheb"));
// 	if(type == 1) pt[1].cheb_to_power();
// 	else if(type == 2) pt[1].cheb_to_power_scale(scale);

// 	// decompose
// 	for(int i=0; i<= tree.depth; i++)
// 	{
// 		for(int j=pow2(i); j<pow2(i+1); j++)
// 		{
// 			if(tree.tree[j]>0) divide_poly(pt[2*j+1], pt[2*j], pt[j], T[tree.tree[j]]);		// pt[j]'s chebcoeff not set
// 		}
// 	}

// 	for(int i=0; i<pow2(tree.depth+1); i++) if(tree.tree[i] == 0) power_to_cheb_int(pt[i], type, scale);
// 	for(int i=0; i<pow2(tree.depth+1); i++) if(tree.tree[i] == 0) print_text_chebcoeff(out, pt[i]);
// }

void poly_decompose_integrate(int deg, int input_type, RR input_scale, vector<RR> &coeff, Tree &tree, evaltype eval_type, int output_type, RR output_scale, ofstream &out_decomp_coeff) {
  vector<Polynomial> T;
  for (int i = 0; i < 100; i++) T.emplace_back(Polynomial(0));

  long m = tree.m, l = tree.l, b = tree.b;
  if (eval_type == evaltype::baby) {
    for (int i = 1; i <= b; i++) {
      T[i].copy(Polynomial(i));
      geneTi(T[i], i, output_type, output_scale);
    }
    for (int i = 1; i <= m - 1; i++) {
      T[pow2(i) * b].copy(Polynomial(pow2(i) * b));
      geneTi(T[pow2(i) * b], pow2(i) * b, output_type, output_scale);
    }
  } else if (eval_type == evaltype::oddbaby)  // comptype = 3
  {
    for (int i = 0; i <= m - 1; i++) {
      T[pow2(i)].copy(Polynomial(pow2(i)));
      geneTi(T[pow2(i)], pow2(i), output_type, output_scale);
    }
    for (int i = 0; i <= l; i++) {
      for (int j = pow2(i - 1) + 1; j <= pow2(i) - 1; j += 2) {
        T[j].copy(Polynomial(j));
        geneTi(T[j], j, output_type, output_scale);
      }
    }
  } else
    throw std::invalid_argument("evaluation type is not correct");

  // decompose setting
  vector<Polynomial> pt;
  for (int i = 0; i < pow2(tree.depth + 1); i++) pt.emplace_back(Polynomial(0));

  // input polynomial setting
  if (input_type == 0)
    pt[1].copy(Polynomial(deg, coeff, "power"));
  else
    pt[1].copy(Polynomial(deg, coeff, "cheb"));
  if (input_type == 1)
    pt[1].cheb_to_power();  // input_type = 0 means power basis, no need to perform cheb_to_power
  else if (input_type == 2)
    pt[1].cheb_to_power_scale(input_scale);

  cout << "coeff print" << endl;
  for (auto num : pt[1].coeff) cout << num << " ";
  cout << endl;

  // decompose (power basis operation)
  for (int i = 0; i <= tree.depth; i++) {
    for (int j = pow2(i); j < pow2(i + 1); j++) {
      if (tree.tree[j] > 0) divide_poly(pt[2 * j + 1], pt[2 * j], pt[j], T[tree.tree[j]]);  // pt[j]'s chebcoeff not set
    }
  }

  for (int i = 0; i < pow2(tree.depth + 1); i++)
    if (tree.tree[i] == 0) power_to_cheb_int(pt[i], output_type, output_scale);
  for (int i = 0; i < pow2(tree.depth + 1); i++)
    if (tree.tree[i] == 0) print_text_chebcoeff(out_decomp_coeff, pt[i]);
}
}  // namespace minicomp
