#include <NTL/RR.h>

#include "Remez.cuh"

namespace boot {
Remez::Remez(RemezParam _params, long _boundary_K, double _log_width, long _deg)
    : params(_params), boundary_K(_boundary_K), log_width(_log_width), deg(_deg) {
  width = to_RR(pow(2.0, -_log_width));
  sc = width / pow(2.0, params.log_scan_step_diff);
  approx = power2_RR(-params.log_approx_degree);

  max_err = 1000;
  min_err = 1;

  sample_point = new Point[deg + 2];
  extreme_point = new Point[2 * deg + 2 * boundary_K];

  coeff = new RR[deg + 1];
}

void Remez::better_initialize() {
  int* nodecount = new int[boundary_K];
  int deg_bdd = deg + 2;

  for (int i = 0; i < boundary_K; i++)
    nodecount[i] = 1;
  int tot_deg = 2 * boundary_K - 1;
  double err = pow(2.0, -log_width);

  RR::SetPrecision(params.RR_prec);

  double* bdd = new double[boundary_K];

  double temp = 0;
  for (int i = 1; i <= (2 * boundary_K - 1); i++)
    temp -= log2((double)i);
  temp += (2 * boundary_K - 1) * log2(2 * M_PI);
  temp += log2(err);

  for (int i = 0; i < boundary_K; i++) {
    bdd[i] = temp;
    for (int j = 1; j <= boundary_K - 1 - i; j++)
      bdd[i] += log2((double)j + err);
    for (int j = 1; j <= boundary_K - 1 + i; j++)
      bdd[i] += log2((double)j + err);
  }

  int max_iter = 200;
  int iter;

  for (iter = 0; iter < max_iter; iter++) {
    if (tot_deg >= deg_bdd)
      break;
    int maxi = max_index(bdd, boundary_K);

    if (maxi != 0) {
      if ((tot_deg + 2) > deg_bdd)
        break;

      for (int i = 0; i < boundary_K; i++) {
        bdd[i] -= log2(tot_deg + 1);
        bdd[i] -= log2(tot_deg + 2);
        bdd[i] += 2.0 * log2(2.0 * M_PI);

        if (i != maxi) {
          bdd[i] += log2(abs((double)(i - maxi)) + err);
          bdd[i] += log2((double)(i + maxi) + err);
        } else {
          bdd[i] += (log2(err) - 1.0);
          bdd[i] += log2(2.0 * (double)i + err);
        }
      }

      tot_deg += 2;
    } else {
      bdd[0] -= log2(tot_deg + 1);
      bdd[0] += (log2(err) - 1.0);
      bdd[0] += log2(2.0 * M_PI);
      for (int i = 1; i < boundary_K; i++) {
        bdd[i] -= log2(tot_deg + 1);
        bdd[i] += log2(2.0 * M_PI);
        bdd[i] += log2((double)i + err);
      }

      tot_deg += 1;
    }

    nodecount[maxi] += 1;
  }

  delete[] bdd;

  if (tot_deg == deg_bdd - 1) {
    nodecount[0]++;
    tot_deg++;
  }

  RR inter_size = RR(pow(2.0, -log_width));

  int cnt = 0;
  if ((nodecount[0] % 2) != 0)
    sample_point[cnt++].x = RR(0);

  for (int i = boundary_K - 1; i > 0; i--) {
    for (int j = 1; j <= nodecount[i]; j++) {
      RR temp = ((RR(2 * j - 1)) * ComputePi_RR()) / (RR(2 * nodecount[i]));
      sample_point[cnt++].x = RR(i) + inter_size * cos(temp);
      sample_point[cnt++].x = RR(-i) - inter_size * cos(temp);
    }
  }

  for (int j = 1; j <= (nodecount[0] / 2); j++) {
    RR temp = ((RR(2 * j - 1)) * ComputePi_RR()) / (RR(2 * nodecount[0]));
    sample_point[cnt++].x = inter_size * cos(temp);
    sample_point[cnt++].x = -inter_size * cos(temp);
  }

  sort(sample_point, sample_point + (deg + 2), xcompare);

  for (int i = 0; i < deg + 2; i++) {
    sample_point[i].y = function_value(sample_point[i].x);
    // std::cout << sample_point[i].x  << std::endl;
  }

  delete[] nodecount;
}

void Remez::initialize() {
  RR::SetPrecision(params.RR_prec);
  long* nodecount = new long[boundary_K];
  long ind = boundary_K - 1;
  bool alter = true;
  for (long j = 0; j < boundary_K; j++) {
    nodecount[j] = 0;
  }
  for (long j = 0; j < (deg + 2) / 2; j++) {
    if (ind > 0) {
      nodecount[ind]++;
      ind--;
    } else {
      if (alter) {
        nodecount[ind]++;
        ind = boundary_K - 1;
        alter = false;
      } else {
        nodecount[boundary_K - 1]++;
        ind = boundary_K - 2;
        alter = true;
      }
    }
  }
  //	for(long j = 0; j < K; j++) {
  //		cout << nodecount[j] << " ";
  //	}
  //	cout << endl;
  ind = 0;
  for (long j = boundary_K - 1; j >= 0; j--) {
    for (long k = 0; k < nodecount[j]; k++) {
      if (j == 0)
        sample_point[ind].x = (-1) * width + width / (nodecount[j] + 1) * (k + 1);
      else
        sample_point[ind].x = (-1) * j - width + 2 * width / (nodecount[j] + 1) * (k + 1);
      sample_point[deg + 1 - ind].x = (-1) * sample_point[ind].x;
      sample_point[ind].y = function_value(sample_point[ind].x);
      sample_point[deg + 1 - ind].y = function_value(sample_point[deg + 1 - ind].x);
      // cout << sam[ind].x << endl;
      ind++;
    }
    //		cout << endl;
  }
  // for(long j = 0; j < deg + 2; j++) {
  //         cout << sample_point[j].x << endl;
  // }
}

void Remez::getcoeffwitherr() {
  RR::SetPrecision(params.RR_prec);

  vec_RR v, w, v_0;
  mat_RR m;

  v.SetLength(deg + 2);
  w.SetLength(deg + 2);
  v_0.SetLength(deg + 2);
  m.SetDims(deg + 2, deg + 2);

  for (long i = 0; i < deg + 2; i++) {
    v[i] = sample_point[i].x;
    w[i] = sample_point[i].y;
  }

  RR var;
  for (long i = 0; i < deg + 2; i++) {
    var = v[i];
    m[i][0] = 1;
    m[i][1] = var / boundary_K;
    for (long j = 2; j < deg + 1; j++) {
      m[i][j] = 2 * (var / boundary_K) * m[i][j - 1] - m[i][j - 2];
    }
    m[i][deg + 1] = 2 * (i % 2) - 1;
  }
  RR deter;
  mat_RR mtrns;
  mat_RR minv;
  transpose(mtrns, m);
  inv(deter, minv, mtrns);
  v_0 = w * minv;
  for (int i = 0; i <= deg; i++) {
    coeff[i] = v_0[i];
  }
  current_err = floor(power2_RR(params.log_round_prec) * abs(v_0[deg + 1])) / power2_RR(params.log_round_prec);
  // cout << current_err << endl;
}

void Remez::getextreme_local(Point* local_extreme_point, long& local_extreme_count, long k) {
  RR::SetPrecision(params.RR_prec);
  RR::SetOutputPrecision(20);
  // cout << "getextreme start" << endl;
  long inc_1 = 0;
  long inc_2 = 0;
  long tmpinc;
  RR scan_1;
  RR scan_2;

  scan_2 = k - width;
  // std::cout << "k = " << k << " width = " << width << std::endl;
  RR scan_y1;
  RR scan_y2;
  scan_y2 = chebeval(deg, coeff, scan_2 / boundary_K) - function_value(scan_2);
  local_extreme_count = 0;

  RR detail[3];
  RR prec_sc, prec_ext, prec_x, tmp;
  long prec_iter, prec_ind;
  bool prec_end;
  long tmp_inc;

  while (scan_2 < k + width + sc) {
    scan_1 = scan_2;
    scan_2 = scan_1 + sc;
    // std::cout << scan_2 << std::endl;
    if (fracpart(scan_2) > width + sc / 2) {
      // cout << "boundary detected" << endl;
      // cout << "inc : " << inc_2 << endl;
      scan_1 = round(scan_1) + width;
      scan_2 = round(scan_2) + 1 - width;
      scan_y1 = chebeval(deg, coeff, scan_1 / boundary_K) - function_value(scan_1);
      scan_y2 = chebeval(deg, coeff, scan_2 / boundary_K) - function_value(scan_2);

      prec_end = false;
      prec_x = scan_1 - sc;
      while (!prec_end) {
        prec_sc = (scan_1 - prec_x) / 2;
        prec_end = true;
        for (long k = 0; k < 3; k++) {
          detail[k] = scan_1 - 2 * prec_sc + prec_sc * k;
        }
        prec_iter = 0;
        while (prec_iter < params.binary_prec) {
          prec_ext = chebeval(deg, coeff, detail[0] / boundary_K) - function_value(detail[0]);
          prec_ind = 0;
          for (long k = 1; k < 3; k++) {
            tmp = chebeval(deg, coeff, detail[k] / boundary_K) - function_value(detail[k]);
            if ((inc_2 == 1 && prec_ext < tmp) || (inc_2 == -1 && prec_ext > tmp)) {
              prec_ext = tmp;
              prec_ind = k;
            }
          }
          if (prec_ind != 2) prec_end = false;
          prec_x = detail[prec_ind];
          prec_sc = prec_sc / 2;
          for (long k = 0; k < 3; k++) {
            detail[k] = std::min(prec_x + prec_sc, scan_1) - 2 * prec_sc + prec_sc * k;
          }
          prec_iter++;
          // cout << "prec_ind : " << prec_ind << " prec_x : " << prec_x << endl;
        }
        if (inc_2 == 1)
          tmpinc = 1;
        else
          tmpinc = -1;

        // cout << "###" << inc_2 << " " << prec_ext << " " << tmpinc * prec_ext << endl;
        if (tmpinc * prec_ext >= current_err) {
          local_extreme_point[local_extreme_count].x = prec_x;
          local_extreme_point[local_extreme_count].y = prec_ext;
          if (inc_2 == 1)
            local_extreme_point[local_extreme_count].locmm = 1;
          else
            local_extreme_point[local_extreme_count].locmm = -1;
          // cout << extreme_point[extreme_count].x << " " << extreme_point[extreme_count].y << " " << extreme_point[extreme_count].locmm << endl;
          local_extreme_count++;
          // cout << "prec_end : " << prec_end << endl;
        }
        if (!prec_end) inc_2 = (-1) * inc_2;
      }
      inc_2 = 0;
    } else {
      inc_1 = inc_2;
      scan_y1 = scan_y2;
      scan_y2 = chebeval(deg, coeff, scan_2 / boundary_K) - function_value(scan_2);
      if (scan_y1 < scan_y2)
        inc_2 = 1;
      else if (scan_y1 > scan_y2)
        inc_2 = -1;
      else
        inc_2 = 0;
      // cout << "inc : " << inc_2 << endl;
      if ((inc_1 == 1 && inc_2 != 1) || (inc_1 == -1 && inc_2 != -1) || inc_1 == 0) {
        prec_end = false;
        tmp_inc = inc_2;
        prec_x = scan_2;
        while (!prec_end) {
          prec_sc = (prec_x - scan_1) / 2;
          prec_end = true;
          if (inc_1 != 0) {
            // cout << "extreme detected" << endl;
            for (long k = 0; k < 3; k++) {
              detail[k] = scan_1 - prec_sc + prec_sc * k;
            }
          } else {
            // cout << "boundary start" << endl;
            for (long k = 0; k < 3; k++) {
              detail[k] = scan_1 + prec_sc * k;
            }
          }
          prec_iter = 0;
          while (prec_iter < params.binary_prec) {
            prec_ext = chebeval(deg, coeff, detail[0] / boundary_K) - function_value(detail[0]);
            prec_ind = 0;
            for (long k = 1; k < 3; k++) {
              tmp = chebeval(deg, coeff, detail[k] / boundary_K) - function_value(detail[k]);
              if ((inc_1 == 1 && prec_ext < tmp) || (inc_1 == -1 && prec_ext > tmp)) {
                prec_ext = tmp;
                prec_ind = k;
              } else if (inc_1 == 0) {
                if ((inc_2 == 1 && prec_ext > tmp) || (inc_2 == -1 && prec_ext < tmp)) {
                  prec_ext = tmp;
                  prec_ind = k;
                }
              }
            }

            if (inc_1 == 0 && prec_ind != 0) prec_end = false;
            prec_x = detail[prec_ind];
            prec_sc = prec_sc / 2;
            // cout << "prec_ind : " << prec_ind << " prec_x : " << prec_x << endl;
            for (long k = 0; k < 3; k++) {
              if (inc_1 == 0)
                detail[k] = std::max(prec_x - prec_sc, scan_1) + prec_sc * k;
              else
                detail[k] = prec_x - prec_sc + prec_sc * k;
            }
            prec_iter++;
          }
          if (inc_2 == 1)
            tmpinc = -1;
          else
            tmpinc = 1;
          // cout << "###" << inc_2 << " " << prec_ext << " " << tmpinc * prec_ext << endl;
          if (tmpinc * prec_ext >= current_err) {
            local_extreme_point[local_extreme_count].x = prec_x;
            local_extreme_point[local_extreme_count].y = prec_ext;
            if (inc_2 == 1)
              local_extreme_point[local_extreme_count].locmm = -1;
            else
              local_extreme_point[local_extreme_count].locmm = +1;

            // cout << extreme_point[local_extreme_count].x << " " << extreme_point[local_extreme_count].y << " " << extreme_point[local_extreme_count].locmm << endl;
            local_extreme_count++;
            // cout << "prec_end : " << prec_end << endl;
          }
          if (!prec_end) {
            inc_2 = (-1) * inc_2;
            // cout << "inc : " << inc_2 << endl;
          }
        }
        inc_2 = tmp_inc;
      }
    }
  }
}

void Remez::getextreme() {
  Point** local_extreme_point_array = new Point*[2 * boundary_K - 1];
  for (int i = 0; i < 2 * boundary_K - 1; i++) {
    local_extreme_point_array[i] = new Point[4 * deg];
  }
  long* local_extreme_count_array = new long[2 * boundary_K - 1];

  vector<thread> vec_thr;

  for (int i = 0; i < 2 * boundary_K - 1; i++) {
    vec_thr.emplace_back(&Remez::getextreme_local, this, local_extreme_point_array[i], ref(local_extreme_count_array[i]), i - boundary_K + 1);
    // getextreme_local(local_extreme_point_array[i], local_extreme_count_array[i], i - boundary_K + 1);
  }

  for (auto& t : vec_thr)
    t.join();

  extreme_count = 0;
  for (int i = 0; i < 2 * boundary_K - 1; i++) {
    for (int j = 0; j < local_extreme_count_array[i]; j++) {
      extreme_point[extreme_count].x = local_extreme_point_array[i][j].x;
      extreme_point[extreme_count].y = local_extreme_point_array[i][j].y;
      extreme_point[extreme_count].locmm = local_extreme_point_array[i][j].locmm;

      extreme_count++;
    }
  }

  // cout << "ok" << endl;
  max_err = 0;
  for (long i = 0; i < extreme_count; i++) {
    // cout << extreme_point[i].x << " " << extreme_point[i].y << " " << extreme_point[i].locmm << endl;
    if (max_err < abs(extreme_point[i].y)) max_err = abs(extreme_point[i].y);
  }
  sort(extreme_point, extreme_point + extreme_count, xcompare);
  // cout << "**" << extreme_count << endl;

  for (int i = 0; i < 2 * boundary_K - 1; i++) {
    delete[] local_extreme_point_array[i];
  }

  delete[] local_extreme_point_array;
  delete[] local_extreme_count_array;
}

void Remez::choosemaxs() {
  RR::SetPrecision(params.RR_prec);
  Point* extract = new Point[extreme_count];
  long count = 0, ind = 0;
  max_err = RR(0);
  min_err = RR(1000);

  long* temparray = new long[extreme_count];
  long tempcount = 0;
  RR maxtempvalue;
  long maxtemp;
  while (ind < extreme_count) {
    if (tempcount == 0) {
      temparray[tempcount] = ind;
      tempcount++;
      ind++;
    }

    else {
      if (extreme_point[ind].locmm * extreme_point[ind - 1].locmm == 1) {
        temparray[tempcount] = ind;
        tempcount++;
        ind++;
      } else {
        maxtempvalue = RR(0);
        for (long i = 0; i < tempcount; i++) {
          if (maxtempvalue < abs(extreme_point[temparray[i]].y)) {
            maxtempvalue = abs(extreme_point[temparray[i]].y);
            maxtemp = temparray[i];
          }
        }
        extract[count].x = extreme_point[maxtemp].x;
        extract[count].y = extreme_point[maxtemp].y;
        extract[count].locmm = extreme_point[maxtemp].locmm;
        count++;
        tempcount = 0;
      }
    }
  }
  maxtempvalue = RR(0);
  for (long i = 0; i < tempcount; i++) {
    if (maxtempvalue < abs(extreme_point[temparray[i]].y)) {
      maxtempvalue = abs(extreme_point[temparray[i]].y);
      maxtemp = temparray[i];
    }
  }
  extract[count].x = extreme_point[maxtemp].x;
  extract[count].y = extreme_point[maxtemp].y;
  extract[count].locmm = extreme_point[maxtemp].locmm;
  count++;
  tempcount = 0;

  // cout << count << " " << deg + 2 << endl;
  RR minsum;
  long minindex = 0;
  while (count > deg + 2) {
    minsum = 100000;
    if (count == deg + 3) {
      if (abs(extract[0].y) > abs(extract[count - 1].y)) {
        count--;
      } else {
        for (long i = 0; i < count - 1; i++) {
          extract[i].x = extract[i + 1].x;
          extract[i].y = extract[i + 1].y;
          extract[i].locmm = extract[i + 1].locmm;
        }
        count--;
      }
    }

    else if (count == deg + 4) {
      for (long i = 0; i < count; i++) {
        if (minsum > abs(extract[i].y) + abs(extract[(i + 1) % count].y)) {
          minsum = abs(extract[i].y) + abs(extract[(i + 1) % count].y);
          minindex = i;
        }
      }
      if (minindex == count - 1) {
        for (long i = 0; i < count - 2; i++) {
          extract[i].x = extract[i + 1].x;
          extract[i].y = extract[i + 1].y;
          extract[i].locmm = extract[i + 1].locmm;
        }
        count -= 2;
      } else {
        for (long i = minindex; i < count - 2; i++) {
          extract[i].x = extract[i + 2].x;
          extract[i].y = extract[i + 2].y;
          extract[i].locmm = extract[i + 2].locmm;
        }
        count -= 2;
      }
    }

    else {
      for (long i = 0; i < count - 1; i++) {
        if (minsum > abs(extract[i].y) + abs(extract[i + 1].y)) {
          minsum = abs(extract[i].y) + abs(extract[i + 1].y);
          minindex = i;
        }
      }
      if (minindex == 0) {
        for (long i = 0; i < count - 1; i++) {
          extract[i].x = extract[i + 1].x;
          extract[i].y = extract[i + 1].y;
          extract[i].locmm = extract[i + 1].locmm;
        }
        count--;
      } else if (minindex == count - 2) {
        count--;
      } else {
        for (long i = minindex; i < count - 2; i++) {
          extract[i].x = extract[i + 2].x;
          extract[i].y = extract[i + 2].y;
          extract[i].locmm = extract[i + 2].locmm;
        }
        count -= 2;
      }
    }
  }
  for (long i = 0; i < (deg + 2); i++) {
    sample_point[i].x = extract[i].x;
    sample_point[i].y = function_value(sample_point[i].x);
    // cout << extract[i].y << endl;
    if (max_err < abs(extract[i].y)) {
      //			cout << maxerr << endl;
      max_err = abs(extract[i].y);
    }
    if (min_err > abs(extract[i].y)) min_err = abs(extract[i].y);
  }

  delete[] extract;
  delete[] temparray;
}

void Remez::generate_optimal_poly(Polynomial& poly) {
  RR::SetOutputPrecision(17);

  better_initialize();
  long it = 0;

  while ((max_err - min_err) / min_err > approx) {
    getcoeffwitherr();
    getextreme();
    choosemaxs();
    it++;
    // for(long i = 0; i < deg + 2; i++) {
    //         cout << sample_point[i].x << endl;
    // }
    // cout << it << "th end" << endl;
    // cout << max_err << " " << min_err << endl;
  }

  // for(long i = 0; i < extreme_count; i++) {
  //         cout << extreme_point[i].x << " " << extreme_point[i].y << " " << extreme_point[i].locmm << endl;
  // }
  // cout << endl;

  // showcoeff();
  poly.set_polynomial(deg, coeff, "cheb");
  //	showgraph(out, coeff, deg, K, sc);
}

void Remez::showcoeff() {
  for (int i = 0; i < deg + 1; i++) {
    cout << i << " : " << coeff[i] << endl;
  }
}
}  // namespace boot
