#include "Bootstrapper.cuh"

Bootstrapper::Bootstrapper(
    long _loge,
    long _logn,
    long _logNh,
    long _L,
    double _final_scale,
    long _boundary_K,
    long _sin_cos_deg,
    long _scale_factor,
    long _inverse_deg,
    CKKSEvaluator *_ckks)
    : loge(_loge), logn(_logn), logNh(_logNh), L(_L), final_scale(_final_scale), boundary_K(_boundary_K), sin_cos_deg(_sin_cos_deg), scale_factor(_scale_factor), inverse_deg(_inverse_deg), ckks(_ckks) {
  n = 1 << logn;
  Nh = 1 << logNh;
  mod_reducer =
      new ModularReducer(boundary_K, (double)loge, sin_cos_deg, scale_factor, inverse_deg, ckks);
}

void Bootstrapper::addLeftRotKeys_Linear_to_vector(vector<int> &gal_steps_vector) {
  int split_point = floor(logn / 2.0);
  int totlen1 = (1 << split_point) - 1;
  int totlen2 = (1 << (logn - split_point)) - 1;
  int gs1 = giantstep(2 * totlen1 + 1);
  int gs2 = giantstep(totlen2 + 1);
  int gs2_e = giantstep(2 * totlen2 + 1);
  int basicstart1 = -totlen1 + gs1 * floor((totlen1 + 0.0) / (gs1 + 0.0));
  int giantfirst1 = -floor((totlen1 + 0.0) / (gs1 + 0.0));
  int giantlast1 = floor((2 * totlen1 + 0.0) / (gs1 + 0.0)) + giantfirst1;
  int basicstart2_e = -totlen2 + gs2_e * floor((totlen2 + 0.0) / (gs2_e + 0.0));
  int giantfirst2_e = -floor((totlen2 + 0.0) / (gs2_e + 0.0));
  int giantlast2_e = floor((2 * totlen2 + 0.0) / (gs2_e + 0.0)) + giantfirst2_e;
  int giantlast2 = floor((totlen2 + 0.0) / (gs2 + 0.0));
  int basicstep = (1 << split_point);

  for (int i = basicstart1; i < basicstart1 + gs1; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i) % Nh);
      }
    }
  }

  for (int i = basicstart2_e; i < basicstart2_e + gs2_e; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * basicstep) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * basicstep) % Nh);
      }
    }
  }

  for (int i = 1; i < gs2; i++) {
    if (find(gal_steps_vector.begin(), gal_steps_vector.end(), i * basicstep) == gal_steps_vector.end()) {
      gal_steps_vector.push_back(i * basicstep);
    }
  }

  for (int i = giantfirst1; i <= giantlast1; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * gs1) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * gs1) % Nh);
      }
    }
  }

  for (int i = giantfirst2_e; i <= giantlast2_e; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * gs2_e * basicstep) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * gs2_e * basicstep) % Nh);
      }
    }
  }

  for (int i = 1; i <= giantlast2; i++) {
    if (find(gal_steps_vector.begin(), gal_steps_vector.end(), i * gs2 * basicstep) == gal_steps_vector.end()) {
      gal_steps_vector.push_back(i * gs2 * basicstep);
    }
  }
}

void Bootstrapper::addLeftRotKeys_Linear_to_vector_3(vector<int> &gal_steps_vector) {
  int div_part1 = floor(logn / 3.0);
  int div_part2 = floor((logn - div_part1) / 2.0);
  int div_part3 = logn - div_part1 - div_part2;

  int totlen1 = (1 << div_part1) - 1;
  int totlen2 = (1 << div_part2) - 1;
  int totlen3 = (1 << div_part3) - 1;

  int basicstep1 = (1 << (logn - div_part1));
  int basicstep2 = (1 << (logn - div_part1 - div_part2));
  int basicstep3 = 1;

  int gs1, gs1_e = 0;
  gs1 = giantstep(totlen1 + 1);
  if (logn != logNh)
    gs1_e = giantstep(2 * totlen1 + 1);

  int gs2 = giantstep(2 * totlen2 + 1);
  int gs3 = giantstep(2 * totlen3 + 1);

  int basicstart1 = -totlen1 + gs1 * floor((totlen1 + 0.0) / (gs1 + 0.0));
  int giantfirst1 = -floor((totlen1 + 0.0) / (gs1 + 0.0));
  int giantlast1 = floor((2 * totlen1 + 0.0) / (gs1 + 0.0)) + giantfirst1;

  int giantlast1_e = 0;
  if (logn != logNh)
    giantlast1_e = floor((totlen1 + 0.0) / (gs1 + 0.0));

  int basicstart2 = -totlen2 + gs2 * floor((totlen2 + 0.0) / (gs2 + 0.0));
  int giantfirst2 = -floor((totlen2 + 0.0) / (gs2 + 0.0));
  int giantlast2 = floor((2 * totlen2 + 0.0) / (gs2 + 0.0)) + giantfirst2;

  int basicstart3 = -totlen3 + gs3 * floor((totlen3 + 0.0) / (gs3 + 0.0));
  int giantfirst3 = -floor((totlen3 + 0.0) / (gs3 + 0.0));
  int giantlast3 = floor((2 * totlen3 + 0.0) / (gs3 + 0.0)) + giantfirst3;

  for (int i = basicstart1; i < basicstart1 + gs1; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * basicstep1) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * basicstep1) % Nh);
      }
    }
  }

  for (int i = 1; i < gs1_e; i++) {
    if (find(gal_steps_vector.begin(), gal_steps_vector.end(), i * basicstep1) == gal_steps_vector.end()) {
      gal_steps_vector.push_back(i * basicstep1);
    }
  }

  for (int i = basicstart2; i < basicstart2 + gs2; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * basicstep2) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * basicstep2) % Nh);
      }
    }
  }

  for (int i = basicstart3; i < basicstart3 + gs3; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * basicstep3) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * basicstep3) % Nh);
      }
    }
  }

  for (int i = giantfirst1; i <= giantlast1; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * gs1 * basicstep1) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * gs1 * basicstep1) % Nh);
      }
    }
  }

  for (int i = 1; i <= giantlast1_e; i++) {
    if (find(gal_steps_vector.begin(), gal_steps_vector.end(), i * gs1_e * basicstep1) == gal_steps_vector.end()) {
      gal_steps_vector.push_back(i * gs1_e * basicstep1);
    }
  }

  for (int i = giantfirst2; i <= giantlast2; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * gs2 * basicstep2) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * gs2 * basicstep2) % Nh);
      }
    }
  }

  for (int i = giantfirst3; i <= giantlast3; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * gs3 * basicstep3) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * gs3 * basicstep3) % Nh);
      }
    }
  }
}
void Bootstrapper::addLeftRotKeys_Linear_to_vector_3_other_slots(vector<int> &gal_steps_vector, long other_logn) {
  int div_part1 = floor(logn / 3.0);
  int div_part2 = floor((logn - div_part1) / 2.0);
  int div_part3 = logn - div_part1 - div_part2;

  int totlen1 = (1 << div_part1) - 1;
  int totlen2 = (1 << div_part2) - 1;
  int totlen3 = (1 << div_part3) - 1;

  int basicstep1 = (1 << (logn - div_part1));
  int basicstep2 = (1 << (logn - div_part1 - div_part2));
  int basicstep3 = 1;

  int gs1, gs1_e = 0;
  gs1 = giantstep(totlen1 + 1);
  if (logn != logNh)
    gs1_e = giantstep(2 * totlen1 + 1);

  int gs2 = giantstep(2 * totlen2 + 1);
  int gs3 = giantstep(2 * totlen3 + 1);

  int basicstart1 = -totlen1 + gs1 * floor((totlen1 + 0.0) / (gs1 + 0.0));
  int giantfirst1 = -floor((totlen1 + 0.0) / (gs1 + 0.0));
  int giantlast1 = floor((2 * totlen1 + 0.0) / (gs1 + 0.0)) + giantfirst1;

  int giantlast1_e = 0;
  if (logn != logNh)
    giantlast1_e = floor((totlen1 + 0.0) / (gs1 + 0.0));

  int basicstart2 = -totlen2 + gs2 * floor((totlen2 + 0.0) / (gs2 + 0.0));
  int giantfirst2 = -floor((totlen2 + 0.0) / (gs2 + 0.0));
  int giantlast2 = floor((2 * totlen2 + 0.0) / (gs2 + 0.0)) + giantfirst2;

  int basicstart3 = -totlen3 + gs3 * floor((totlen3 + 0.0) / (gs3 + 0.0));
  int giantfirst3 = -floor((totlen3 + 0.0) / (gs3 + 0.0));
  int giantlast3 = floor((2 * totlen3 + 0.0) / (gs3 + 0.0)) + giantfirst3;

  for (int i = basicstart1; i < basicstart1 + gs1; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * basicstep1) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * basicstep1) % Nh);
      }
    }
  }

  for (int i = 1; i < gs1_e; i++) {
    if (find(gal_steps_vector.begin(), gal_steps_vector.end(), i * basicstep1) == gal_steps_vector.end()) {
      gal_steps_vector.push_back(i * basicstep1);
    }
  }

  for (int i = basicstart2; i < basicstart2 + gs2; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * basicstep2) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * basicstep2) % Nh);
      }
    }
  }

  for (int i = basicstart3; i < basicstart3 + gs3; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * basicstep3) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * basicstep3) % Nh);
      }
    }
  }

  for (int i = giantfirst1; i <= giantlast1; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * gs1 * basicstep1) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * gs1 * basicstep1) % Nh);
      }
    }
  }

  for (int i = 1; i <= giantlast1_e; i++) {
    if (find(gal_steps_vector.begin(), gal_steps_vector.end(), i * gs1_e * basicstep1) == gal_steps_vector.end()) {
      gal_steps_vector.push_back(i * gs1_e * basicstep1);
    }
  }

  for (int i = giantfirst2; i <= giantlast2; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * gs2 * basicstep2) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * gs2 * basicstep2) % Nh);
      }
    }
  }

  for (int i = giantfirst3; i <= giantlast3; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * gs3 * basicstep3) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * gs3 * basicstep3) % Nh);
      }
    }
  }
}
void Bootstrapper::addLeftRotKeys_Linear_to_vector_one_depth(vector<int> &gal_steps_vector) {
  int totlen1 = (1 << logn) - 1;
  int totlen2 = (1 << logn) - 1;
  int gs1 = giantstep(2 * totlen1 + 1);
  int gs2 = giantstep(totlen2 + 1);
  int gs2_e = giantstep(2 * totlen2 + 1);
  int basicstart1 = -totlen1 + gs1 * floor((totlen1 + 0.0) / (gs1 + 0.0));
  int giantfirst1 = -floor((totlen1 + 0.0) / (gs1 + 0.0));
  int giantlast1 = floor((2 * totlen1 + 0.0) / (gs1 + 0.0)) + giantfirst1;
  int basicstart2_e = -totlen2 + gs2_e * floor((totlen2 + 0.0) / (gs2_e + 0.0));
  int giantfirst2_e = -floor((totlen2 + 0.0) / (gs2_e + 0.0));
  int giantlast2_e = floor((2 * totlen2 + 0.0) / (gs2_e + 0.0)) + giantfirst2_e;
  int giantlast2 = floor((totlen2 + 0.0) / (gs2 + 0.0));

  for (int i = basicstart1; i < basicstart1 + gs1; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i) % Nh);
      }
    }
  }

  for (int i = basicstart2_e; i < basicstart2_e + gs2_e; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i) % Nh);
      }
    }
  }

  for (int i = 1; i < gs2; i++) {
    if (find(gal_steps_vector.begin(), gal_steps_vector.end(), i) == gal_steps_vector.end()) {
      gal_steps_vector.push_back(i);
    }
  }

  for (int i = giantfirst1; i <= giantlast1; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * gs1) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * gs1) % Nh);
      }
    }
  }

  for (int i = giantfirst2_e; i <= giantlast2_e; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), (Nh + i * gs2_e) % Nh) == gal_steps_vector.end()) {
        gal_steps_vector.push_back((Nh + i * gs2_e) % Nh);
      }
    }
  }

  for (int i = 1; i <= giantlast2; i++) {
    if (find(gal_steps_vector.begin(), gal_steps_vector.end(), i * gs2) == gal_steps_vector.end()) {
      gal_steps_vector.push_back(i * gs2);
    }
  }
}

void Bootstrapper::addLeftRotKeys_Linear_to_vector_one_depth_more_depth(vector<int> &gal_steps_vector) {
  for (int i = 0; i < Nh; i++) {
    if (i != 0) {
      if (find(gal_steps_vector.begin(), gal_steps_vector.end(), i) == gal_steps_vector.end()) {
        gal_steps_vector.push_back(i);
      }
    }
  }
}

void Bootstrapper::addBootKeys(PhantomGaloisKey &gal_keys) {
  vector<int> gal_steps_vector;
  gal_steps_vector.push_back(0);
  for (int i = 0; i < logNh; i++) {
    gal_steps_vector.push_back((1 << i));
  }
  addLeftRotKeys_Linear_to_vector(gal_steps_vector);
  ckks->decryptor.create_galois_keys_from_steps(gal_steps_vector, *(ckks->galois_keys));

  slot_vec.push_back(logn);

  slot_index = -1;
  for (int i = 0; i < slot_vec.size(); i++) {
    if (slot_vec[i] == logn) {
      slot_index = i;
      break;
    }
  }
  if (slot_index == -1)
    throw("LT coefficients were not generated for this logn");
}

void Bootstrapper::addBootKeys_3(PhantomGaloisKey &gal_keys) {
  vector<int> gal_steps_vector;
  gal_steps_vector.push_back(0);
  for (int i = 0; i < logNh; i++) {
    gal_steps_vector.push_back((1 << i));
  }
  addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
  ckks->decryptor.create_galois_keys_from_steps(gal_steps_vector, *(ckks->galois_keys));

  slot_vec.push_back(logn);
  slot_index = -1;
  for (int i = 0; i < slot_vec.size(); i++) {
    if (slot_vec[i] == logn) {
      slot_index = i;
      break;
    }
  }
  if (slot_index == -1)
    throw("LT coefficients were not generated for this logn");
}

void Bootstrapper::addBootKeys_3_other_slots(PhantomGaloisKey &gal_keys, vector<long> &other_logn_vec) {
  vector<int> gal_steps_vector;
  gal_steps_vector.push_back(0);
  for (int i = 0; i < logNh; i++) {
    gal_steps_vector.push_back((1 << i));
  }
  for (auto i : other_logn_vec) {
    addLeftRotKeys_Linear_to_vector_3_other_slots(gal_steps_vector, i);
    slot_vec.push_back(i);
  }
  ckks->decryptor.create_galois_keys_from_steps(gal_steps_vector, *(ckks->galois_keys));

  slot_index = -1;
  for (int i = 0; i < slot_vec.size(); i++) {
    if (slot_vec[i] == logn) {
      slot_index = i;
      break;
    }
  }
  if (slot_index == -1)
    throw("LT coefficients were not generated for this logn");
}

void Bootstrapper::addBootKeys_3_other_slots_keys(PhantomGaloisKey &gal_keys, vector<long> &other_logn_vec, vector<int> &other_keys) {
  vector<int> gal_steps_vector;
  gal_steps_vector.push_back(0);
  for (int i = 0; i < logNh; i++) {
    gal_steps_vector.push_back((1 << i));
  }
  for (auto rot : other_keys) {
    if (find(gal_steps_vector.begin(), gal_steps_vector.end(), rot) == gal_steps_vector.end())
      gal_steps_vector.push_back(rot);
  }
  for (auto i : other_logn_vec) {
    addLeftRotKeys_Linear_to_vector_3_other_slots(gal_steps_vector, i);
    slot_vec.push_back(i);
  }
  ckks->decryptor.create_galois_keys_from_steps(gal_steps_vector, *(ckks->galois_keys));

  slot_index = -1;
  for (int i = 0; i < slot_vec.size(); i++) {
    if (slot_vec[i] == logn) {
      slot_index = i;
      break;
    }
  }
  if (slot_index == -1)
    throw("LT coefficients were not generated for this logn");
}

void Bootstrapper::addBootKeys_other_keys(PhantomGaloisKey &gal_keys, vector<int> &other_keys) {
  vector<int> gal_steps_vector;
  gal_steps_vector.push_back(0);
  for (int i = 0; i < logNh; i++) {
    gal_steps_vector.push_back((1 << i));
  }

  for (auto rot : other_keys) {
    //	if(find(gal_steps_vector.begin(), gal_steps_vector.end(), rot) == gal_steps_vector.end()) gal_steps_vector.emplace_back(rot);
    if (find(gal_steps_vector.begin(), gal_steps_vector.end(), rot) == gal_steps_vector.end())
      gal_steps_vector.push_back(rot);
  }

  addLeftRotKeys_Linear_to_vector(gal_steps_vector);
  //	for(auto num : gal_steps_vector) std::cout << num << " ";
  //	cout << endl;
  ckks->decryptor.create_galois_keys_from_steps(gal_steps_vector, *(ckks->galois_keys));
}
void Bootstrapper::addBootKeys_3_other_keys(PhantomGaloisKey &gal_keys, vector<int> &other_keys) {
  vector<int> gal_steps_vector;
  gal_steps_vector.push_back(0);
  for (int i = 0; i < logNh; i++) {
    gal_steps_vector.push_back((1 << i));
  }
  for (auto rot : other_keys) {
    if (find(gal_steps_vector.begin(), gal_steps_vector.end(), rot) == gal_steps_vector.end())
      gal_steps_vector.push_back(rot);
  }

  addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
  ckks->decryptor.create_galois_keys_from_steps(gal_steps_vector, *(ckks->galois_keys));
}

void Bootstrapper::addBootKeys_hoisting(PhantomGaloisKey &gal_keys) {
  vector<int> gal_steps_vector;
  gal_steps_vector.push_back(0);
  for (int i = 0; i < logNh; i++) {
    gal_steps_vector.push_back((1 << i));
  }
  addLeftRotKeys_Linear_to_vector(gal_steps_vector);
  ckks->decryptor.create_galois_keys_from_steps(gal_steps_vector, *(ckks->galois_keys));
}

void Bootstrapper::addBootKeys_one_depth(PhantomGaloisKey &gal_keys) {
  vector<int> gal_steps_vector;
  gal_steps_vector.push_back(0);
  for (int i = 0; i < logNh; i++) {
    gal_steps_vector.push_back((1 << i));
  }
  addLeftRotKeys_Linear_to_vector_one_depth(gal_steps_vector);
  ckks->decryptor.create_galois_keys_from_steps(gal_steps_vector, *(ckks->galois_keys));
}

void Bootstrapper::addBootKeys_one_depth_more_depth(PhantomGaloisKey &gal_keys) {
  vector<int> gal_steps_vector;
  gal_steps_vector.push_back(0);
  for (int i = 0; i < logNh; i++) {
    gal_steps_vector.push_back((1 << i));
  }
  addLeftRotKeys_Linear_to_vector_one_depth_more_depth(gal_steps_vector);
  ckks->decryptor.create_galois_keys_from_steps(gal_steps_vector, *(ckks->galois_keys));
}

void Bootstrapper::change_logn(long new_logn) {
  logn = new_logn;
  n = (1 << logn);
  slot_index = -1;
  for (int i = 0; i < slot_vec.size(); i++) {
    if (slot_vec[i] == logn) {
      slot_index = i;
      break;
    }
  }
  if (slot_index == -1)
    throw("LT coefficients were not generated for this logn");
}

void Bootstrapper::genorigcoeff() {
  orig_coeffvec.resize(slot_vec.size());
  orig_invcoeffvec.resize(slot_vec.size());
  for (int u = 0; u < slot_vec.size(); u++) {
    long new_logn = slot_vec[u];
    long new_n = (1 << new_logn);
    double theta_0 = M_PI / (2 * new_n);
    double current_theta;
    complex<double> current_zeta;
    int current_power;
    int blocklen = 1;
    int blockcount = new_n;

    orig_coeffvec[u].resize(new_logn);

    for (int i = 0; i < new_logn; i++) {
      orig_coeffvec[u][i].resize(3);
      for (int j = 0; j < 3; j++) {
        orig_coeffvec[u][i][j].resize(new_n);
      }
    }

    for (int i = 0; i < new_logn; i++) {
      blocklen = blocklen << 1;
      blockcount = blockcount >> 1;
      current_theta = theta_0 * (1 << (new_logn - 1 - i));
      current_power = 1;
      current_zeta = std::polar(1.0, current_theta * current_power);
      for (int j = 0; j < blocklen / 2; j++) {
        for (int k = 0; k < blockcount; k++) {
          orig_coeffvec[u][i][1][k * blocklen + j] = 1;
          orig_coeffvec[u][i][1][k * blocklen + j + blocklen / 2] = -current_zeta;

          orig_coeffvec[u][i][0][k * blocklen + j] = 0;                 // index 2 -> 0
          orig_coeffvec[u][i][0][k * blocklen + j + blocklen / 2] = 1;  // index 2 -> 0

          orig_coeffvec[u][i][2][k * blocklen + j] = current_zeta;      // index 0 -> 2
          orig_coeffvec[u][i][2][k * blocklen + j + blocklen / 2] = 0;  // index 0 -> 2
        }
        current_power = (5 * current_power) % (1 << (i + 3));
        current_zeta = std::polar(1.0, current_theta * current_power);
      }
    }

    theta_0 = -M_PI / (2 * new_n);
    blocklen = new_n;
    blockcount = 1;

    orig_invcoeffvec.clear();
    orig_invcoeffvec[u].resize(new_logn);

    for (int i = 0; i < new_logn; i++) {
      orig_invcoeffvec[u][i].resize(3);
      for (int j = 0; j < 3; j++) {
        orig_invcoeffvec[u][i][j].resize(new_n);
      }
    }

    for (int i = 0; i < new_logn; i++) {
      current_theta = theta_0 * (1 << i);
      current_power = 1;
      current_zeta = std::polar(1.0, current_theta * current_power);
      for (int j = 0; j < blocklen / 2; j++) {
        for (int k = 0; k < blockcount; k++) {
          orig_invcoeffvec[u][i][1][k * blocklen + j] = 0.5;
          orig_invcoeffvec[u][i][1][k * blocklen + j + blocklen / 2] = -0.5 * current_zeta;

          orig_invcoeffvec[u][i][0][k * blocklen + j] = 0;
          orig_invcoeffvec[u][i][0][k * blocklen + j + blocklen / 2] = 0.5 * current_zeta;

          orig_invcoeffvec[u][i][2][k * blocklen + j] = 0.5;
          orig_invcoeffvec[u][i][2][k * blocklen + j + blocklen / 2] = 0;
        }
        current_power = (5 * current_power) % (1 << ((new_logn - 1 - i) + 3));
        current_zeta = std::polar(1.0, current_theta * current_power);
      }
      blocklen = blocklen >> 1;
      blockcount = blockcount << 1;
    }
  }
}

void Bootstrapper::genfftcoeff_one_depth() {
  fftcoeff1.resize(slot_vec.size());
  vector<complex<double>> tmpvec;
  vector<int> tmpcount;
  int totlen1, all_case_count, pos, current_pos, ind, ind_res, curr_logn, curr_n;

  for (int u = 0; u < slot_vec.size(); u++) {
    curr_logn = slot_vec[u];
    curr_n = (1 << curr_logn);

    totlen1 = (1 << curr_logn) - 1;
    fftcoeff1[u].resize(2 * totlen1 + 1);

    for (int i = 0; i < 2 * totlen1 + 1; i++)
      fftcoeff1[u][i].resize(curr_n, 0);

    tmpvec.clear();
    tmpvec.resize(curr_n);

    tmpcount.clear();
    tmpcount.resize(curr_logn);

    all_case_count = pow(3, curr_logn);

    for (int j = 0; j < all_case_count; j++) {
      ind = j;
      pos = 0;
      for (int p = 0; p < curr_logn; p++) {
        ind_res = ind % 3;
        pos += (ind_res - 1) * (1 << p);
        tmpcount[p] = ind_res;
        ind = (ind - ind_res) / 3;
      }

      current_pos = pos;
      for (int k = 0; k < curr_n; k++) {
        tmpvec[k] = 1.0;
      }

      for (int p = 0; p < curr_logn; p++) {
        current_pos = current_pos - (tmpcount[p] - 1) * (1 << p);
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = tmpvec[k] * orig_coeffvec[u][p][tmpcount[p]][((k + current_pos + curr_n) % curr_n)];
        }
      }

      for (int k = 0; k < curr_n; k++) {
        fftcoeff1[u][pos + totlen1][k] += tmpvec[k];
      }
    }
  }
}

void Bootstrapper::genfftcoeff_full_one_depth() {
  fftcoeff1.resize(slot_vec.size());
  vector<complex<double>> tmpvec;
  vector<int> tmpcount;
  int totlen1, all_case_count, pos, current_pos, ind, ind_res, curr_logn, curr_n;
  for (int u = 0; u < slot_vec.size(); u++) {
    curr_logn = slot_vec[u];
    curr_n = (1 << curr_logn);
    totlen1 = (1 << curr_logn) - 1;

    fftcoeff1[u].resize(totlen1 + 1);

    for (int i = 0; i < totlen1 + 1; i++)
      fftcoeff1[u][i].resize(curr_n, 0);

    tmpvec.clear();
    tmpvec.resize(curr_n);

    tmpcount.clear();
    tmpcount.resize(curr_logn);
    all_case_count = pow(3, curr_logn);

    for (int j = 0; j < all_case_count; j++) {
      ind = j;
      pos = 0;
      for (int p = 0; p < curr_logn; p++) {
        ind_res = ind % 3;
        pos += (ind_res - 1) * (1 << p);
        tmpcount[p] = ind_res;
        ind = (ind - ind_res) / 3;
      }

      current_pos = pos;
      for (int k = 0; k < curr_n; k++) {
        tmpvec[k] = 1.0;
      }

      for (int p = 0; p < curr_logn; p++) {
        current_pos = current_pos - (tmpcount[p] - 1) * (1 << p);
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = tmpvec[k] * orig_coeffvec[u][p][tmpcount[p]][((k + (curr_n + current_pos)) % curr_n)];
        }
      }

      for (int k = 0; k < curr_n; k++) {
        fftcoeff1[u][(pos + totlen1 + 1) % (totlen1 + 1)][k] += tmpvec[k];
      }
    }
  }
}
void Bootstrapper::geninvfftcoeff_one_depth() {
  invfftcoeff1.resize(slot_vec.size());
  vector<complex<double>> tmpvec;
  vector<int> tmpcount;
  int totlen1, all_case_count, pos, current_pos, ind, ind_res, curr_logn, curr_n;
  for (int u = 0; u < slot_vec.size(); u++) {
    curr_logn = slot_vec[u];
    curr_n = (1 << curr_logn);
    totlen1 = (1 << curr_logn) - 1;

    invfftcoeff1[u].resize(totlen1 + 1);

    for (int i = 0; i < totlen1 + 1; i++)
      invfftcoeff1[u][i].resize(curr_n, 0);

    tmpvec.clear();
    tmpvec.resize(curr_n);

    tmpcount.clear();
    tmpcount.resize(curr_logn);
    all_case_count = pow(3, curr_logn);

    for (int j = 0; j < all_case_count; j++) {
      ind = j;
      pos = 0;
      for (int p = 0; p < curr_logn; p++) {
        ind_res = ind % 3;
        pos += (ind_res - 1) * (1 << (curr_logn - 1 - p));
        tmpcount[p] = ind_res;
        ind = (ind - ind_res) / 3;
      }
      current_pos = pos;
      for (int k = 0; k < curr_n; k++) {
        tmpvec[k] = 1.0;
      }

      for (int p = 0; p < curr_logn; p++) {
        current_pos = current_pos - (tmpcount[p] - 1) * (1 << (curr_logn - 1 - p));
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = tmpvec[k] * orig_invcoeffvec[u][p][tmpcount[p]][((k + (curr_n + current_pos)) % curr_n)];
        }
      }

      for (int k = 0; k < curr_n; k++) {
        invfftcoeff1[u][(pos + totlen1 + 1) % (totlen1 + 1)][k] += tmpvec[k];
      }
    }
  }
}

void Bootstrapper::generate_LT_coefficient_one_depth() {
  genorigcoeff();
  if (logn == logNh)
    genfftcoeff_full_one_depth();
  else
    genfftcoeff_one_depth();
  geninvfftcoeff_one_depth();
}

void Bootstrapper::genfftcoeff() {
  fftcoeff1.resize(slot_vec.size());
  fftcoeff2.resize(slot_vec.size());
  vector<complex<double>> tmpvec;
  vector<int> tmpcount;
  int totlen1, totlen2, split_point, basicstep, all_case_count, pos, current_pos, ind, ind_res, curr_logn, curr_n;
  for (int u = 0; u < slot_vec.size(); u++) {
    curr_logn = slot_vec[u];
    curr_n = (1 << curr_logn);
    split_point = floor(curr_logn / 2.0);
    totlen1 = (1 << split_point) - 1;
    totlen2 = (1 << (curr_logn - split_point)) - 1;

    fftcoeff1[u].resize(2 * totlen1 + 1);
    fftcoeff2[u].resize(2 * totlen2 + 1);

    for (int i = 0; i < 2 * totlen1 + 1; i++)
      fftcoeff1[u][i].resize(2 * curr_n, 0);
    for (int i = 0; i < 2 * totlen2 + 1; i++)
      fftcoeff2[u][i].resize(2 * curr_n, 0);

    tmpvec.clear();
    tmpvec.resize(curr_n);

    tmpcount.clear();
    tmpcount.resize(split_point);
    all_case_count = pow(3, split_point);

    for (int j = 0; j < all_case_count; j++) {
      ind = j;
      pos = 0;
      for (int p = 0; p < split_point; p++) {
        ind_res = ind % 3;
        pos += (ind_res - 1) * (1 << p);
        tmpcount[p] = ind_res;
        ind = (ind - ind_res) / 3;
      }

      current_pos = pos;
      for (int k = 0; k < curr_n; k++) {
        tmpvec[k] = 1.0;
      }

      for (int p = 0; p < split_point; p++) {
        current_pos = current_pos - (tmpcount[p] - 1) * (1 << p);
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = tmpvec[k] * orig_coeffvec[u][p][tmpcount[p]][((k + current_pos + curr_n) % curr_n)];
        }
      }

      for (int k = 0; k < curr_n; k++) {
        fftcoeff1[u][pos + totlen1][k] += tmpvec[k];
      }
    }

    for (int i = 0; i < 2 * totlen1 + 1; i++) {
      for (int j = 0; j < curr_n; j++)
        fftcoeff1[u][i][j + curr_n] = fftcoeff1[u][i][j];
    }

    basicstep = (1 << split_point);
    tmpcount.clear();
    tmpcount.resize(curr_logn - split_point);
    all_case_count = pow(3, curr_logn - split_point);

    for (int j = 0; j < all_case_count; j++) {
      ind = j;
      pos = 0;
      for (int p = 0; p < curr_logn - split_point; p++) {
        ind_res = ind % 3;
        pos += (ind_res - 1) * (1 << p);
        tmpcount[p] = ind_res;
        ind = (ind - ind_res) / 3;
      }

      current_pos = pos;
      for (int k = 0; k < curr_n; k++) {
        tmpvec[k] = 1.0;
      }

      for (int p = 0; p < curr_logn - split_point; p++) {
        current_pos = current_pos - (tmpcount[p] - 1) * (1 << p);
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = tmpvec[k] * orig_coeffvec[u][p + split_point][tmpcount[p]][((k + basicstep * (curr_n + current_pos)) % curr_n)];
        }
      }

      for (int k = 0; k < curr_n; k++) {
        fftcoeff2[u][pos + totlen2][k] += tmpvec[k];
      }
    }

    for (int i = 0; i < 2 * totlen2 + 1; i++) {
      for (int j = 0; j < curr_n; j++)
        fftcoeff2[u][i][j + curr_n] = complex<double>(0, 1) * fftcoeff2[u][i][j];
    }
  }
}

void Bootstrapper::genfftcoeff_full() {
  fftcoeff1.resize(slot_vec.size());
  fftcoeff2.resize(slot_vec.size());
  vector<complex<double>> tmpvec;
  vector<int> tmpcount;
  int totlen1, totlen2, split_point, basicstep, all_case_count, pos, current_pos, ind, ind_res, curr_logn, curr_n;
  for (int u = 0; u < slot_vec.size(); u++) {
    curr_logn = slot_vec[u];
    curr_n = (1 << curr_logn);
    split_point = floor(curr_logn / 2.0);
    totlen1 = (1 << split_point) - 1;
    totlen2 = (1 << (curr_logn - split_point)) - 1;

    fftcoeff1[u].resize(2 * totlen1 + 1);
    fftcoeff2[u].resize(totlen2 + 1);

    for (int i = 0; i < 2 * totlen1 + 1; i++)
      fftcoeff1[u][i].resize(curr_n, 0);
    for (int i = 0; i < totlen2 + 1; i++)
      fftcoeff2[u][i].resize(curr_n, 0);

    tmpvec.clear();
    tmpvec.resize(curr_n);

    tmpcount.clear();
    tmpcount.resize(split_point);
    all_case_count = pow(3, split_point);

    for (int j = 0; j < all_case_count; j++) {
      ind = j;
      pos = 0;
      for (int p = 0; p < split_point; p++) {
        ind_res = ind % 3;
        pos += (ind_res - 1) * (1 << p);
        tmpcount[p] = ind_res;
        ind = (ind - ind_res) / 3;
      }

      current_pos = pos;
      for (int k = 0; k < curr_n; k++) {
        tmpvec[k] = 1.0;
      }

      for (int p = 0; p < split_point; p++) {
        current_pos = current_pos - (tmpcount[p] - 1) * (1 << p);
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = tmpvec[k] * orig_coeffvec[u][p][tmpcount[p]][((k + current_pos + curr_n) % curr_n)];
        }
      }

      for (int k = 0; k < curr_n; k++) {
        fftcoeff1[u][pos + totlen1][k] += tmpvec[k];
      }
    }

    basicstep = (1 << split_point);
    tmpcount.clear();
    tmpcount.resize(curr_logn - split_point);
    all_case_count = pow(3, curr_logn - split_point);

    for (int j = 0; j < all_case_count; j++) {
      ind = j;
      pos = 0;
      for (int p = 0; p < curr_logn - split_point; p++) {
        ind_res = ind % 3;
        pos += (ind_res - 1) * (1 << p);
        tmpcount[p] = ind_res;
        ind = (ind - ind_res) / 3;
      }

      current_pos = pos;
      for (int k = 0; k < curr_n; k++) {
        tmpvec[k] = 1.0;
      }

      for (int p = 0; p < curr_logn - split_point; p++) {
        current_pos = current_pos - (tmpcount[p] - 1) * (1 << p);
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = tmpvec[k] * orig_coeffvec[u][p + split_point][tmpcount[p]][((k + basicstep * (curr_n + current_pos)) % curr_n)];
        }
      }

      for (int k = 0; k < curr_n; k++) {
        fftcoeff2[u][(pos + totlen2 + 1) % (totlen2 + 1)][k] += tmpvec[k];
      }
    }
  }
}

void Bootstrapper::geninvfftcoeff() {
  invfftcoeff1.resize(slot_vec.size());
  invfftcoeff2.resize(slot_vec.size());
  vector<complex<double>> tmpvec;
  vector<int> tmpcount;
  int totlen1, totlen2, split_point, basicstep, all_case_count, pos, current_pos, ind, ind_res, curr_logn, curr_n;
  for (int u = 0; u < slot_vec.size(); u++) {
    curr_logn = slot_vec[u];
    curr_n = (1 << curr_logn);
    split_point = ceil(curr_logn / 2.0);
    totlen1 = (1 << split_point) - 1;
    totlen2 = (1 << (curr_logn - split_point)) - 1;

    invfftcoeff1[u].resize(totlen1 + 1);
    invfftcoeff2[u].resize(2 * totlen2 + 1);

    for (int i = 0; i < totlen1 + 1; i++)
      invfftcoeff1[u][i].resize(curr_n, 0);
    for (int i = 0; i < 2 * totlen2 + 1; i++)
      invfftcoeff2[u][i].resize(2 * curr_n, 0);

    tmpvec.clear();
    tmpvec.resize(curr_n);

    tmpcount.clear();
    tmpcount.resize(split_point);
    all_case_count = pow(3, split_point);
    basicstep = (1 << (curr_logn - split_point));

    for (int j = 0; j < all_case_count; j++) {
      ind = j;
      pos = 0;
      for (int p = 0; p < split_point; p++) {
        ind_res = ind % 3;
        pos += (ind_res - 1) * (1 << (split_point - 1 - p));
        tmpcount[p] = ind_res;
        ind = (ind - ind_res) / 3;
      }
      current_pos = pos;
      for (int k = 0; k < curr_n; k++) {
        tmpvec[k] = 1.0;
      }

      for (int p = 0; p < split_point; p++) {
        current_pos = current_pos - (tmpcount[p] - 1) * (1 << (split_point - 1 - p));
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = tmpvec[k] * orig_invcoeffvec[u][p][tmpcount[p]][((k + basicstep * (curr_n + current_pos)) % curr_n)];
        }
      }

      for (int k = 0; k < curr_n; k++) {
        invfftcoeff1[u][(pos + totlen1 + 1) % (totlen1 + 1)][k] += tmpvec[k];
      }
    }
    tmpcount.clear();
    tmpcount.resize(curr_logn - split_point);
    all_case_count = pow(3, curr_logn - split_point);
    for (int j = 0; j < all_case_count; j++) {
      ind = j;
      pos = 0;
      for (int p = 0; p < curr_logn - split_point; p++) {
        ind_res = ind % 3;
        pos += (ind_res - 1) * (1 << (curr_logn - split_point - 1 - p));
        tmpcount[p] = ind_res;
        ind = (ind - ind_res) / 3;
      }
      current_pos = pos;
      for (int k = 0; k < curr_n; k++) {
        tmpvec[k] = 1.0;
      }
      for (int p = 0; p < curr_logn - split_point; p++) {
        current_pos = current_pos - (tmpcount[p] - 1) * (1 << (curr_logn - split_point - 1 - p));
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = tmpvec[k] * orig_invcoeffvec[u][p + split_point][tmpcount[p]][((k + current_pos + curr_n) % curr_n)];
        }
      }

      for (int k = 0; k < curr_n; k++) {
        invfftcoeff2[u][pos + totlen2][k] += tmpvec[k];
      }
    }

    for (int i = 0; i < totlen1 + 1; i++) {
      for (int j = 0; j < curr_n; j++)
        invfftcoeff1[u][i][j] *= 1.0 / (boundary_K * (1 << (logNh - curr_logn)));
    }

    for (int i = 0; i < 2 * totlen2 + 1; i++) {
      for (int j = 0; j < curr_n; j++) {
        invfftcoeff2[u][i][j] *= 0.5;
        invfftcoeff2[u][i][j + curr_n] = complex<double>(0, -1) * invfftcoeff2[u][i][j];
      }
    }
  }
}

void Bootstrapper::geninvfftcoeff_full() {
  invfftcoeff1.resize(slot_vec.size());
  invfftcoeff2.resize(slot_vec.size());
  vector<complex<double>> tmpvec;
  vector<int> tmpcount;
  int totlen1, totlen2, split_point, basicstep, all_case_count, pos, current_pos, ind, ind_res, curr_logn, curr_n;
  for (int u = 0; u < slot_vec.size(); u++) {
    curr_logn = slot_vec[u];
    curr_n = (1 << curr_logn);
    split_point = ceil(curr_logn / 2.0);
    totlen1 = (1 << split_point) - 1;
    totlen2 = (1 << (curr_logn - split_point)) - 1;

    invfftcoeff1.resize(totlen1 + 1);
    invfftcoeff2.resize(2 * totlen2 + 1);

    for (int i = 0; i < totlen1 + 1; i++)
      invfftcoeff1[u][i].resize(curr_n, 0);
    for (int i = 0; i < 2 * totlen2 + 1; i++)
      invfftcoeff2[u][i].resize(curr_n, 0);

    tmpvec.clear();
    tmpvec.resize(curr_n);

    tmpcount.clear();
    tmpcount.resize(split_point);
    all_case_count = pow(3, split_point);
    basicstep = (1 << (curr_logn - split_point));

    for (int j = 0; j < all_case_count; j++) {
      ind = j;
      pos = 0;
      for (int p = 0; p < split_point; p++) {
        ind_res = ind % 3;
        pos += (ind_res - 1) * (1 << (split_point - 1 - p));
        tmpcount[p] = ind_res;
        ind = (ind - ind_res) / 3;
      }
      current_pos = pos;
      for (int k = 0; k < curr_n; k++) {
        tmpvec[k] = 1.0;
      }

      for (int p = 0; p < split_point; p++) {
        current_pos = current_pos - (tmpcount[p] - 1) * (1 << (split_point - 1 - p));
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = tmpvec[k] * orig_invcoeffvec[u][p][tmpcount[p]][((k + basicstep * (curr_n + current_pos)) % curr_n)];
        }
      }

      for (int k = 0; k < curr_n; k++) {
        invfftcoeff1[u][(pos + totlen1 + 1) % (totlen1 + 1)][k] += tmpvec[k];
      }
    }
    tmpcount.clear();
    tmpcount.resize(logn - split_point);
    all_case_count = pow(3, curr_logn - split_point);
    for (int j = 0; j < all_case_count; j++) {
      ind = j;
      pos = 0;
      for (int p = 0; p < curr_logn - split_point; p++) {
        ind_res = ind % 3;
        pos += (ind_res - 1) * (1 << (curr_logn - split_point - 1 - p));
        tmpcount[p] = ind_res;
        ind = (ind - ind_res) / 3;
      }
      current_pos = pos;
      for (int k = 0; k < curr_n; k++) {
        tmpvec[k] = 1.0;
      }
      for (int p = 0; p < curr_logn - split_point; p++) {
        current_pos = current_pos - (tmpcount[p] - 1) * (1 << (curr_logn - split_point - 1 - p));
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = tmpvec[k] * orig_invcoeffvec[u][p + split_point][tmpcount[p]][((k + current_pos + curr_n) % curr_n)];
        }
      }

      for (int k = 0; k < curr_n; k++) {
        invfftcoeff2[u][pos + totlen2][k] += tmpvec[k];
      }
    }

    for (int i = 0; i < totlen1 + 1; i++) {
      for (int j = 0; j < curr_n; j++)
        invfftcoeff1[u][i][j] *= 1.0 / boundary_K;
    }

    for (int i = 0; i < 2 * totlen2 + 1; i++) {
      for (int j = 0; j < curr_n; j++)
        invfftcoeff2[u][i][j] *= 0.5;
    }
  }
}

void Bootstrapper::genfftcoeff_3() {  // not yet
  fftcoeff1.resize(slot_vec.size());
  fftcoeff2.resize(slot_vec.size());
  fftcoeff3.resize(slot_vec.size());
  vector<complex<double>> tmpvec;
  vector<int> tmpcount;
  int div_part1, div_part2, div_part3;
  int totlen1, totlen2, totlen3;
  int basicstep1, basicstep2, basicstep3;
  int all_case_count, pos, current_pos, ind, ind_res;
  int curr_logn, curr_n;

  for (int u = 0; u < slot_vec.size(); u++) {
    curr_logn = slot_vec[u];
    curr_n = (1 << curr_logn);
    if (curr_logn == logNh) {
      div_part3 = floor(curr_logn / 3.0);
      div_part2 = floor((curr_logn - div_part3) / 2.0);
      div_part1 = curr_logn - div_part3 - div_part2;

      totlen1 = (1 << div_part1) - 1;
      totlen2 = (1 << div_part2) - 1;
      totlen3 = (1 << div_part3) - 1;

      basicstep1 = 1;
      basicstep2 = (1 << div_part1);
      basicstep3 = (1 << (div_part1 + div_part2));

      fftcoeff1[u].resize(2 * totlen1 + 1);
      fftcoeff2[u].resize(2 * totlen2 + 1);
      fftcoeff3[u].resize(totlen3 + 1);

      for (int i = 0; i < 2 * totlen1 + 1; i++)
        fftcoeff1[u][i].resize(curr_n, 0);
      for (int i = 0; i < 2 * totlen2 + 1; i++)
        fftcoeff2[u][i].resize(curr_n, 0);
      for (int i = 0; i < totlen3 + 1; i++)
        fftcoeff3[u][i].resize(curr_n, 0);

      tmpvec.clear();
      tmpvec.resize(curr_n);

      tmpcount.clear();
      tmpcount.resize(div_part1);
      all_case_count = pow(3, div_part1);

      for (int j = 0; j < all_case_count; j++) {
        ind = j;
        pos = 0;
        for (int p = 0; p < div_part1; p++) {
          ind_res = ind % 3;
          pos += (ind_res - 1) * (1 << p);
          tmpcount[p] = ind_res;
          ind = (ind - ind_res) / 3;
        }

        current_pos = pos;
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = 1.0;
        }

        for (int p = 0; p < div_part1; p++) {
          current_pos = current_pos - (tmpcount[p] - 1) * (1 << p);
          for (int k = 0; k < curr_n; k++) {
            tmpvec[k] = tmpvec[k] * orig_coeffvec[u][p][tmpcount[p]][((k + basicstep1 * (curr_n + current_pos)) % curr_n)];
          }
        }

        for (int k = 0; k < curr_n; k++) {
          fftcoeff1[u][pos + totlen1][k] += tmpvec[k];
        }
      }

      tmpcount.clear();
      tmpcount.resize(div_part2);
      all_case_count = pow(3, div_part2);

      for (int j = 0; j < all_case_count; j++) {
        ind = j;
        pos = 0;
        for (int p = 0; p < div_part2; p++) {
          ind_res = ind % 3;
          pos += (ind_res - 1) * (1 << p);
          tmpcount[p] = ind_res;
          ind = (ind - ind_res) / 3;
        }

        current_pos = pos;
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = 1.0;
        }

        for (int p = 0; p < div_part2; p++) {
          current_pos = current_pos - (tmpcount[p] - 1) * (1 << p);
          for (int k = 0; k < curr_n; k++) {
            tmpvec[k] = tmpvec[k] * orig_coeffvec[u][p + div_part1][tmpcount[p]][((k + basicstep2 * (curr_n + current_pos)) % curr_n)];
          }
        }

        for (int k = 0; k < curr_n; k++) {
          fftcoeff2[u][pos + totlen2][k] += tmpvec[k];
        }
      }

      tmpcount.clear();
      tmpcount.resize(div_part3);
      all_case_count = pow(3, div_part3);

      for (int j = 0; j < all_case_count; j++) {
        ind = j;
        pos = 0;
        for (int p = 0; p < div_part3; p++) {
          ind_res = ind % 3;
          pos += (ind_res - 1) * (1 << p);
          tmpcount[p] = ind_res;
          ind = (ind - ind_res) / 3;
        }

        current_pos = pos;
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = 1.0;
        }

        for (int p = 0; p < div_part3; p++) {
          current_pos = current_pos - (tmpcount[p] - 1) * (1 << p);
          for (int k = 0; k < curr_n; k++) {
            tmpvec[k] = tmpvec[k] * orig_coeffvec[u][p + div_part1 + div_part2][tmpcount[p]][((k + basicstep3 * (curr_n + current_pos)) % curr_n)];
          }
        }

        for (int k = 0; k < curr_n; k++) {
          fftcoeff3[u][(pos + totlen3 + 1) % (totlen3 + 1)][k] += tmpvec[k];
        }
      }

    }

    else {
      div_part3 = floor(curr_logn / 3.0);
      div_part2 = floor((curr_logn - div_part3) / 2.0);
      div_part1 = curr_logn - div_part3 - div_part2;

      totlen1 = (1 << div_part1) - 1;
      totlen2 = (1 << div_part2) - 1;
      totlen3 = (1 << div_part3) - 1;

      basicstep1 = 1;
      basicstep2 = (1 << div_part1);
      basicstep3 = (1 << (div_part1 + div_part2));

      fftcoeff1[u].resize(2 * totlen1 + 1);
      fftcoeff2[u].resize(2 * totlen2 + 1);
      fftcoeff3[u].resize(2 * totlen3 + 1);

      for (int i = 0; i < 2 * totlen1 + 1; i++)
        fftcoeff1[u][i].resize(2 * curr_n, 0);
      for (int i = 0; i < 2 * totlen2 + 1; i++)
        fftcoeff2[u][i].resize(2 * curr_n, 0);
      for (int i = 0; i < 2 * totlen3 + 1; i++)
        fftcoeff3[u][i].resize(2 * curr_n, 0);

      tmpvec.clear();
      tmpvec.resize(curr_n);

      tmpcount.clear();
      tmpcount.resize(div_part1);
      all_case_count = pow(3, div_part1);

      for (int j = 0; j < all_case_count; j++) {
        ind = j;
        pos = 0;
        for (int p = 0; p < div_part1; p++) {
          ind_res = ind % 3;
          pos += (ind_res - 1) * (1 << p);
          tmpcount[p] = ind_res;
          ind = (ind - ind_res) / 3;
        }

        current_pos = pos;
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = 1.0;
        }

        for (int p = 0; p < div_part1; p++) {
          current_pos = current_pos - (tmpcount[p] - 1) * (1 << p);
          for (int k = 0; k < curr_n; k++) {
            tmpvec[k] = tmpvec[k] * orig_coeffvec[u][p][tmpcount[p]][((k + basicstep1 * (curr_n + current_pos)) % curr_n)];
          }
        }

        for (int k = 0; k < curr_n; k++) {
          fftcoeff1[u][pos + totlen1][k] += tmpvec[k];
        }
      }

      for (int i = 0; i < 2 * totlen1 + 1; i++) {
        for (int j = 0; j < curr_n; j++)
          fftcoeff1[u][i][j + curr_n] = fftcoeff1[u][i][j];
      }

      tmpcount.clear();
      tmpcount.resize(div_part2);
      all_case_count = pow(3, div_part2);

      for (int j = 0; j < all_case_count; j++) {
        ind = j;
        pos = 0;
        for (int p = 0; p < div_part2; p++) {
          ind_res = ind % 3;
          pos += (ind_res - 1) * (1 << p);
          tmpcount[p] = ind_res;
          ind = (ind - ind_res) / 3;
        }

        current_pos = pos;
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = 1.0;
        }

        for (int p = 0; p < div_part2; p++) {
          current_pos = current_pos - (tmpcount[p] - 1) * (1 << p);
          for (int k = 0; k < curr_n; k++) {
            tmpvec[k] = tmpvec[k] * orig_coeffvec[u][p + div_part1][tmpcount[p]][((k + basicstep2 * (curr_n + current_pos)) % curr_n)];
          }
        }

        for (int k = 0; k < curr_n; k++) {
          fftcoeff2[u][pos + totlen2][k] += tmpvec[k];
        }
      }

      for (int i = 0; i < 2 * totlen2 + 1; i++) {
        for (int j = 0; j < curr_n; j++)
          fftcoeff2[u][i][j + curr_n] = fftcoeff2[u][i][j];
      }

      tmpcount.clear();
      tmpcount.resize(div_part3);
      all_case_count = pow(3, div_part3);

      for (int j = 0; j < all_case_count; j++) {
        ind = j;
        pos = 0;
        for (int p = 0; p < div_part3; p++) {
          ind_res = ind % 3;
          pos += (ind_res - 1) * (1 << p);
          tmpcount[p] = ind_res;
          ind = (ind - ind_res) / 3;
        }

        current_pos = pos;
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = 1.0;
        }

        for (int p = 0; p < div_part3; p++) {
          current_pos = current_pos - (tmpcount[p] - 1) * (1 << p);
          for (int k = 0; k < curr_n; k++) {
            tmpvec[k] = tmpvec[k] * orig_coeffvec[u][p + div_part1 + div_part2][tmpcount[p]][((k + basicstep3 * (curr_n + current_pos)) % curr_n)];
          }
        }

        for (int k = 0; k < curr_n; k++) {
          fftcoeff3[u][pos + totlen3][k] += tmpvec[k];
        }
      }

      for (int i = 0; i < 2 * totlen3 + 1; i++) {
        for (int j = 0; j < curr_n; j++)
          fftcoeff3[u][i][j + curr_n] = complex<double>(0, 1) * fftcoeff3[u][i][j];
      }
    }
  }
}

// void Bootstrapper::genfftcoeff_full_3() {
// fftcoeff1.resize(slot_vec.size());
// fftcoeff2.resize(slot_vec.size());
// fftcoeff3.resize(slot_vec.size());
// vector<complex<double>> tmpvec;
// vector<int> tmpcount;
// int div_part1, div_part2, div_part3;
// int totlen1, totlen2, totlen3;
// int basicstep1, basicstep2, basicstep3;
// int all_case_count, pos, current_pos, ind, ind_res;
// int curr_logn, curr_n;

// for(int u = 0; u < slot_vec.size(); u++) {
// curr_logn = slot_vec[u];
// curr_n = (1 << curr_logn);
// div_part3 = floor(curr_logn / 3.0);
// div_part2 = floor((curr_logn - div_part3) / 2.0);
// div_part1 = curr_logn - div_part3 - div_part2;

// totlen1 = (1 << div_part1) - 1;
// totlen2 = (1 << div_part2) - 1;
// totlen3 = (1 << div_part3) - 1;

// basicstep1 = 1;
// basicstep2 = (1 << div_part1);
// basicstep3 = (1 << (div_part1 + div_part2));

// fftcoeff1[u].resize(2*totlen1 + 1);
// fftcoeff2[u].resize(2*totlen2 + 1);
// fftcoeff3[u].resize(totlen3 + 1);

// for(int i = 0; i < 2*totlen1 + 1; i++) fftcoeff1[u][i].resize(curr_n, 0);
// for(int i = 0; i < 2*totlen2 + 1; i++) fftcoeff2[u][i].resize(curr_n, 0);
// for(int i = 0; i < totlen3 + 1; i++) fftcoeff3[u][i].resize(curr_n, 0);

// tmpvec.clear();
// tmpvec.resize(curr_n);

// tmpcount.clear();
// tmpcount.resize(div_part1);
// all_case_count = pow(3, div_part1);

// for(int j = 0; j < all_case_count; j++) {
// ind = j;
// pos = 0;
// for(int p = 0; p < div_part1; p++) {
// ind_res = ind % 3;
// pos += (ind_res - 1) * (1 << p);
// tmpcount[p] = ind_res;
// ind = (ind - ind_res) / 3;
//}

// current_pos = pos;
// for(int k = 0; k < curr_n; k++) {
// tmpvec[k] = 1.0;
//}

// for(int p = 0; p < div_part1; p++) {
// current_pos = current_pos - (tmpcount[p] - 1) * (1 << p);
// for(int k = 0; k < curr_n; k++) {
// tmpvec[k] = tmpvec[k] * orig_coeffvec[u][p][tmpcount[p]][((k + basicstep1 * (curr_n + current_pos)) % curr_n)];
//}
//}

// for(int k = 0; k < curr_n; k++) {
// fftcoeff1[u][pos + totlen1][k] += tmpvec[k];
//}
//}

// tmpcount.clear();
// tmpcount.resize(div_part2);
// all_case_count = pow(3, div_part2);

// for(int j = 0; j < all_case_count; j++) {
// ind = j;
// pos = 0;
// for(int p = 0; p < div_part2; p++) {
// ind_res = ind % 3;
// pos += (ind_res - 1) * (1 << p);
// tmpcount[p] = ind_res;
// ind = (ind - ind_res) / 3;
//}

// current_pos = pos;
// for(int k = 0; k < curr_n; k++) {
// tmpvec[k] = 1.0;
//}

// for(int p = 0; p < div_part2; p++) {
// current_pos = current_pos - (tmpcount[p] - 1) * (1 << p);
// for(int k = 0; k < curr_n; k++) {
// tmpvec[k] = tmpvec[k] * orig_coeffvec[u][p + div_part1][tmpcount[p]][((k + basicstep2 * (curr_n + current_pos)) % curr_n)];
//}
//}

// for(int k = 0; k < curr_n; k++) {
// fftcoeff2[u][pos + totlen2][k] += tmpvec[k];
//}
//}

// tmpcount.clear();
// tmpcount.resize(div_part3);
// all_case_count = pow(3, div_part3);

// for(int j = 0; j < all_case_count; j++) {
// ind = j;
// pos = 0;
// for(int p = 0; p < div_part3; p++) {
// ind_res = ind % 3;
// pos += (ind_res - 1) * (1 << p);
// tmpcount[p] = ind_res;
// ind = (ind - ind_res) / 3;
//}

// current_pos = pos;
// for(int k = 0; k < curr_n; k++) {
// tmpvec[k] = 1.0;
//}

// for(int p = 0; p < div_part3; p++) {
// current_pos = current_pos - (tmpcount[p] - 1) * (1 << p);
// for(int k = 0; k < curr_n; k++) {
// tmpvec[k] = tmpvec[k] * orig_coeffvec[u][p + div_part1 + div_part2][tmpcount[p]][((k + basicstep3 * (curr_n + current_pos)) % curr_n)];
//}
//}

// for(int k = 0; k < curr_n; k++) {
// fftcoeff3[u][(pos + totlen3 + 1) % (totlen3 + 1)][k] += tmpvec[k];
//}
//}
//}
//}

void Bootstrapper::geninvfftcoeff_3() {  // not yet
  invfftcoeff1.resize(slot_vec.size());
  invfftcoeff2.resize(slot_vec.size());
  invfftcoeff3.resize(slot_vec.size());
  vector<complex<double>> tmpvec;
  vector<int> tmpcount;
  int div_part1, div_part2, div_part3;
  int totlen1, totlen2, totlen3;
  int basicstep1, basicstep2, basicstep3;
  int all_case_count, pos, current_pos, ind, ind_res;
  int curr_logn, curr_n;

  for (int u = 0; u < slot_vec.size(); u++) {
    curr_logn = slot_vec[u];
    curr_n = (1 << curr_logn);
    if (curr_logn == logNh) {
      div_part1 = floor(curr_logn / 3.0);
      div_part2 = floor((curr_logn - div_part1) / 2.0);
      div_part3 = curr_logn - div_part1 - div_part2;

      totlen1 = (1 << div_part1) - 1;
      totlen2 = (1 << div_part2) - 1;
      totlen3 = (1 << div_part3) - 1;

      basicstep1 = (1 << (curr_logn - div_part1));
      basicstep2 = (1 << (curr_logn - div_part1 - div_part2));
      basicstep3 = 1;

      invfftcoeff1[u].resize(totlen1 + 1);
      invfftcoeff2[u].resize(2 * totlen2 + 1);
      invfftcoeff3[u].resize(2 * totlen3 + 1);

      for (int i = 0; i < totlen1 + 1; i++)
        invfftcoeff1[u][i].resize(curr_n, 0);
      for (int i = 0; i < 2 * totlen2 + 1; i++)
        invfftcoeff2[u][i].resize(curr_n, 0);
      for (int i = 0; i < 2 * totlen3 + 1; i++)
        invfftcoeff3[u][i].resize(curr_n, 0);

      tmpvec.clear();
      tmpvec.resize(curr_n);

      tmpcount.clear();
      tmpcount.resize(div_part1);
      all_case_count = pow(3, div_part1);

      for (int j = 0; j < all_case_count; j++) {
        ind = j;
        pos = 0;
        for (int p = 0; p < div_part1; p++) {
          ind_res = ind % 3;
          pos += (ind_res - 1) * (1 << (div_part1 - 1 - p));
          tmpcount[p] = ind_res;
          ind = (ind - ind_res) / 3;
        }
        current_pos = pos;
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = 1.0;
        }

        for (int p = 0; p < div_part1; p++) {
          current_pos = current_pos - (tmpcount[p] - 1) * (1 << (div_part1 - 1 - p));
          for (int k = 0; k < curr_n; k++) {
            tmpvec[k] = tmpvec[k] * orig_invcoeffvec[u][p][tmpcount[p]][((k + basicstep1 * (curr_n + current_pos)) % curr_n)];
          }
        }

        for (int k = 0; k < curr_n; k++) {
          invfftcoeff1[u][(pos + totlen1 + 1) % (totlen1 + 1)][k] += tmpvec[k];
        }
      }

      tmpcount.clear();
      tmpcount.resize(div_part2);
      all_case_count = pow(3, div_part2);
      for (int j = 0; j < all_case_count; j++) {
        ind = j;
        pos = 0;
        for (int p = 0; p < div_part2; p++) {
          ind_res = ind % 3;
          pos += (ind_res - 1) * (1 << (div_part2 - 1 - p));
          tmpcount[p] = ind_res;
          ind = (ind - ind_res) / 3;
        }
        current_pos = pos;
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = 1.0;
        }
        for (int p = 0; p < div_part2; p++) {
          current_pos = current_pos - (tmpcount[p] - 1) * (1 << (div_part2 - 1 - p));
          for (int k = 0; k < curr_n; k++) {
            tmpvec[k] = tmpvec[k] * orig_invcoeffvec[u][p + div_part1][tmpcount[p]][((k + basicstep2 * (curr_n + current_pos)) % curr_n)];
          }
        }

        for (int k = 0; k < curr_n; k++) {
          invfftcoeff2[u][pos + totlen2][k] += tmpvec[k];
        }
      }

      tmpcount.clear();
      tmpcount.resize(div_part3);
      all_case_count = pow(3, div_part3);
      for (int j = 0; j < all_case_count; j++) {
        ind = j;
        pos = 0;
        for (int p = 0; p < div_part3; p++) {
          ind_res = ind % 3;
          pos += (ind_res - 1) * (1 << (div_part3 - 1 - p));
          tmpcount[p] = ind_res;
          ind = (ind - ind_res) / 3;
        }
        current_pos = pos;
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = 1.0;
        }
        for (int p = 0; p < div_part3; p++) {
          current_pos = current_pos - (tmpcount[p] - 1) * (1 << (div_part3 - 1 - p));
          for (int k = 0; k < curr_n; k++) {
            tmpvec[k] =
                tmpvec[k] * orig_invcoeffvec[u][p + div_part1 + div_part2][tmpcount[p]][((k + basicstep3 * (curr_n + current_pos)) % curr_n)];
          }
        }

        for (int k = 0; k < curr_n; k++) {
          invfftcoeff3[u][pos + totlen3][k] += tmpvec[k];
        }
      }

      for (int i = 0; i < totlen1 + 1; i++) {
        for (int j = 0; j < curr_n; j++)
          invfftcoeff1[u][i][j] *= 1.0 / boundary_K;
      }

      for (int i = 0; i < 2 * totlen3 + 1; i++) {
        for (int j = 0; j < curr_n; j++)
          invfftcoeff3[u][i][j] *= 0.5;
      }

    }

    else {
      div_part1 = floor(curr_logn / 3.0);
      div_part2 = floor((curr_logn - div_part1) / 2.0);
      div_part3 = curr_logn - div_part1 - div_part2;

      totlen1 = (1 << div_part1) - 1;
      totlen2 = (1 << div_part2) - 1;
      totlen3 = (1 << div_part3) - 1;

      basicstep1 = (1 << (curr_logn - div_part1));
      basicstep2 = (1 << (curr_logn - div_part1 - div_part2));
      basicstep3 = 1;

      invfftcoeff1[u].resize(totlen1 + 1);
      invfftcoeff2[u].resize(2 * totlen2 + 1);
      invfftcoeff3[u].resize(2 * totlen3 + 1);

      for (int i = 0; i < totlen1 + 1; i++)
        invfftcoeff1[u][i].resize(curr_n, 0);
      for (int i = 0; i < 2 * totlen2 + 1; i++)
        invfftcoeff2[u][i].resize(curr_n, 0);
      for (int i = 0; i < 2 * totlen3 + 1; i++)
        invfftcoeff3[u][i].resize(2 * curr_n, 0);

      tmpvec.clear();
      tmpvec.resize(curr_n);

      tmpcount.clear();
      tmpcount.resize(div_part1);
      all_case_count = pow(3, div_part1);

      for (int j = 0; j < all_case_count; j++) {
        ind = j;
        pos = 0;
        for (int p = 0; p < div_part1; p++) {
          ind_res = ind % 3;
          pos += (ind_res - 1) * (1 << (div_part1 - 1 - p));
          tmpcount[p] = ind_res;
          ind = (ind - ind_res) / 3;
        }
        current_pos = pos;
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = 1.0;
        }

        for (int p = 0; p < div_part1; p++) {
          current_pos = current_pos - (tmpcount[p] - 1) * (1 << (div_part1 - 1 - p));
          for (int k = 0; k < curr_n; k++) {
            tmpvec[k] = tmpvec[k] * orig_invcoeffvec[u][p][tmpcount[p]][((k + basicstep1 * (curr_n + current_pos)) % curr_n)];
          }
        }

        for (int k = 0; k < curr_n; k++) {
          invfftcoeff1[u][(pos + totlen1 + 1) % (totlen1 + 1)][k] += tmpvec[k];
        }
      }

      tmpcount.clear();
      tmpcount.resize(div_part2);
      all_case_count = pow(3, div_part2);
      for (int j = 0; j < all_case_count; j++) {
        ind = j;
        pos = 0;
        for (int p = 0; p < div_part2; p++) {
          ind_res = ind % 3;
          pos += (ind_res - 1) * (1 << (div_part2 - 1 - p));
          tmpcount[p] = ind_res;
          ind = (ind - ind_res) / 3;
        }
        current_pos = pos;
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = 1.0;
        }
        for (int p = 0; p < div_part2; p++) {
          current_pos = current_pos - (tmpcount[p] - 1) * (1 << (div_part2 - 1 - p));
          for (int k = 0; k < curr_n; k++) {
            tmpvec[k] = tmpvec[k] * orig_invcoeffvec[u][p + div_part1][tmpcount[p]][((k + basicstep2 * (curr_n + current_pos)) % curr_n)];
          }
        }

        for (int k = 0; k < curr_n; k++) {
          invfftcoeff2[u][pos + totlen2][k] += tmpvec[k];
        }
      }

      tmpcount.clear();
      tmpcount.resize(div_part3);
      all_case_count = pow(3, div_part3);
      for (int j = 0; j < all_case_count; j++) {
        ind = j;
        pos = 0;
        for (int p = 0; p < div_part3; p++) {
          ind_res = ind % 3;
          pos += (ind_res - 1) * (1 << (div_part3 - 1 - p));
          tmpcount[p] = ind_res;
          ind = (ind - ind_res) / 3;
        }
        current_pos = pos;
        for (int k = 0; k < curr_n; k++) {
          tmpvec[k] = 1.0;
        }
        for (int p = 0; p < div_part3; p++) {
          current_pos = current_pos - (tmpcount[p] - 1) * (1 << (div_part3 - 1 - p));
          for (int k = 0; k < curr_n; k++) {
            tmpvec[k] =
                tmpvec[k] * orig_invcoeffvec[u][p + div_part1 + div_part2][tmpcount[p]][((k + basicstep3 * (curr_n + current_pos)) % curr_n)];
          }
        }

        for (int k = 0; k < curr_n; k++) {
          invfftcoeff3[u][pos + totlen3][k] += tmpvec[k];
        }
      }

      for (int i = 0; i < totlen1 + 1; i++) {
        for (int j = 0; j < curr_n; j++)
          invfftcoeff1[u][i][j] *= 1.0 / (boundary_K * (1 << (logNh - curr_logn)));
      }

      for (int i = 0; i < 2 * totlen3 + 1; i++) {
        for (int j = 0; j < curr_n; j++) {
          invfftcoeff3[u][i][j] *= 0.5;
          invfftcoeff3[u][i][j + curr_n] = complex<double>(0, -1) * invfftcoeff3[u][i][j];
        }
      }
    }
  }
}

void Bootstrapper::generate_LT_coefficient() {
  genorigcoeff();
  if (logn == logNh) {
    genfftcoeff_full();
    geninvfftcoeff_full();
  } else {
    genfftcoeff();
    geninvfftcoeff();
  }
}

void Bootstrapper::generate_LT_coefficient_3() {
  genorigcoeff();
  genfftcoeff_3();
  geninvfftcoeff_3();
}

void Bootstrapper::prepare_mod_polynomial() {
  mod_reducer->generate_sin_cos_polynomial();
  mod_reducer->generate_inverse_sine_polynomial();
  // mod_reducer->write_polynomials();
}

void Bootstrapper::subsum(double scale, PhantomCiphertext &cipher) {
  int repeatcount = (1 << (logNh - logn));
  PhantomCiphertext tmp;
  int step;
  for (int i = 0; i < logNh - logn; i++) {
    step = (1 << (logn + i));
    ckks->evaluator.rotate_vector(cipher, step, *(ckks->galois_keys), tmp);
    ckks->evaluator.add_inplace_reduced_error(cipher, tmp);
  }

  PhantomPlaintext tmpplain;
  ckks->encoder.encode(1.0 / repeatcount, scale, tmpplain);
  ckks->evaluator.mod_switch_to_inplace(tmpplain, cipher.params_id());
  ckks->evaluator.multiply_plain_inplace(cipher, tmpplain);
  ckks->evaluator.rescale_to_next_inplace(cipher);
}

void Bootstrapper::bsgs_linear_transform(
    PhantomCiphertext &rtncipher, PhantomCiphertext &cipher, int totlen, int basicstep, int coeff_logn, const vector<vector<complex<double>>> &fftcoeff) {
  int gs1 = giantstep(2 * totlen + 1);
  int basicstart1 = -totlen + gs1 * floor((totlen + 0.0) / (gs1 + 0.0));
  int giantfirst1 = -floor((totlen + 0.0) / (gs1 + 0.0));
  int giantlast1 = floor((2 * totlen + 0.0) / (gs1 + 0.0)) + giantfirst1;

  vector<PhantomCiphertext> babyct(gs1, PhantomCiphertext());
  PhantomCiphertext giantct, tmpct;
  bool giantbool = false, tmpctbool = false;
  PhantomCiphertext tmptmpct;
  vector<complex<double>> rotatedcoeff;
  rotatedcoeff.reserve(Nh);

  PhantomPlaintext tmpplain;
  PhantomCiphertext addtmp1, addtmp2;

  for (int i = basicstart1; i < basicstart1 + gs1; i++) {
    if (i == 0)
      babyct[i - basicstart1] = cipher;
    else
      ckks->evaluator.rotate_vector(cipher, (Nh + i * basicstep) % Nh, *(ckks->galois_keys), babyct[i - basicstart1]);
  }

  for (int i = giantfirst1; i <= giantlast1; i++) {
    giantbool = false;
    if (i != giantlast1) {
      for (int j = basicstart1; j < basicstart1 + gs1; j++) {
        rotation(coeff_logn, Nh, (-i) * gs1 * basicstep, fftcoeff[(i * gs1 + j) + totlen], rotatedcoeff);
        ckks->evaluator.multiply_vector_reduced_error(babyct[j - basicstart1], rotatedcoeff, tmptmpct);
        if (!giantbool) {
          giantct = tmptmpct;
          giantbool = true;
        } else
          ckks->evaluator.add_inplace_reduced_error(giantct, tmptmpct);
      }
    } else {
      for (int j = basicstart1; j <= totlen - i * gs1; j++) {
        rotation(coeff_logn, Nh, (-i) * gs1 * basicstep, fftcoeff[(i * gs1 + j) + totlen], rotatedcoeff);
        ckks->evaluator.multiply_vector_reduced_error(babyct[j - basicstart1], rotatedcoeff, tmptmpct);
        if (!giantbool) {
          giantct = tmptmpct;
          giantbool = true;
        } else
          ckks->evaluator.add_inplace_reduced_error(giantct, tmptmpct);
      }
    }

    if (i != 0) {
      ckks->evaluator.rotate_vector(giantct, (Nh + i * gs1 * basicstep) % Nh, *(ckks->galois_keys), tmptmpct);
      if (!tmpctbool) {
        tmpct = tmptmpct;
        tmpctbool = true;
      } else
        ckks->evaluator.add_inplace_reduced_error(tmpct, tmptmpct);
    } else {
      if (!tmpctbool) {
        tmpct = giantct;
        tmpctbool = true;
      } else
        ckks->evaluator.add_inplace_reduced_error(tmpct, giantct);
    }
  }
  rtncipher = tmpct;
}

void Bootstrapper::rotated_bsgs_linear_transform(
    PhantomCiphertext &rtncipher, PhantomCiphertext &cipher, int totlen, int basicstep, int coeff_logn, const vector<vector<complex<double>>> &fftcoeff) {
  int gs2 = giantstep(totlen + 1);
  int giantlast2 = floor((totlen + 0.0) / (gs2 + 0.0));

  vector<PhantomCiphertext> babyct(gs2, PhantomCiphertext());
  PhantomCiphertext giantct, tmpct;
  bool giantbool = false, tmpctbool = false;
  PhantomCiphertext tmptmpct;
  vector<complex<double>> rotatedcoeff;
  rotatedcoeff.reserve(Nh);

  PhantomPlaintext tmpplain;

  PhantomCiphertext addtmp1, addtmp2;

  for (int i = 0; i < gs2; i++) {
    if (i == 0) {
      babyct[i] = cipher;
    } else {
      ckks->evaluator.rotate_vector(cipher, (Nh + i * basicstep) % Nh, *(ckks->galois_keys), babyct[i]);
    }
  }

  for (int i = 0; i <= giantlast2; i++) {
    giantbool = false;
    if (i != giantlast2) {
      for (int j = 0; j < gs2; j++) {
        rotation(coeff_logn, Nh, (-i) * gs2 * basicstep, fftcoeff[i * gs2 + j], rotatedcoeff);
        ckks->evaluator.multiply_vector_reduced_error(babyct[j], rotatedcoeff, tmptmpct);
        if (!giantbool) {
          giantct = tmptmpct;
          giantbool = true;
        } else
          ckks->evaluator.add_inplace_reduced_error(giantct, tmptmpct);
      }
    } else {
      for (int j = 0; j <= totlen - i * gs2; j++) {
        rotation(coeff_logn, Nh, (-i) * gs2 * basicstep, fftcoeff[i * gs2 + j], rotatedcoeff);
        ckks->evaluator.multiply_vector_reduced_error(babyct[j], rotatedcoeff, tmptmpct);
        if (!giantbool) {
          giantct = tmptmpct;
          giantbool = true;
        } else
          ckks->evaluator.add_inplace_reduced_error(giantct, tmptmpct);
      }
    }

    if (i != 0) {
      ckks->evaluator.rotate_vector(giantct, (Nh + i * gs2 * basicstep) % Nh, *(ckks->galois_keys), tmptmpct);
      if (!tmpctbool) {
        tmpct = tmptmpct;
        tmpctbool = true;
      } else
        ckks->evaluator.add_inplace_reduced_error(tmpct, tmptmpct);
    } else {
      if (!tmpctbool) {
        tmpct = giantct;
        tmpctbool = true;
      } else
        ckks->evaluator.add_inplace_reduced_error(tmpct, giantct);
    }
  }

  rtncipher = tmpct;
}
void Bootstrapper::bsgs_linear_transform_hoisting(
    PhantomCiphertext &rtncipher, PhantomCiphertext &cipher, int totlen, int basicstep, int coeff_logn, vector<vector<complex<double>>> fftcoeff) {
  int gs1 = giantstep(2 * totlen + 1);
  int basicstart1 = -totlen + gs1 * floor((totlen + 0.0) / (gs1 + 0.0));
  int giantfirst1 = -floor((totlen + 0.0) / (gs1 + 0.0));
  int giantlast1 = floor((2 * totlen + 0.0) / (gs1 + 0.0)) + giantfirst1;

  PhantomCiphertext *babyct = new PhantomCiphertext[gs1];
  PhantomCiphertext *giantct = 0, *tmpct = 0;
  PhantomCiphertext tmptmpct;
  vector<complex<double>> rotatedcoeff;
  rotatedcoeff.reserve(Nh);

  PhantomPlaintext tmpplain;
  PhantomCiphertext addtmp1, addtmp2;

  for (int i = basicstart1; i < basicstart1 + gs1; i++) {
    if (i == 0)
      babyct[i - basicstart1] = cipher;
    else
      ckks->evaluator.rotate_vector(cipher, (Nh + i * basicstep) % Nh, *(ckks->galois_keys), babyct[i - basicstart1]);
  }

  for (int i = giantfirst1; i <= giantlast1; i++) {
    if (!(giantct == 0))
      delete giantct;
    giantct = 0;
    if (i != giantlast1) {
      for (int j = basicstart1; j < basicstart1 + gs1; j++) {
        rotation(coeff_logn, Nh, (-i) * gs1 * basicstep, fftcoeff[(i * gs1 + j) + totlen], rotatedcoeff);
        ckks->evaluator.multiply_vector_reduced_error(babyct[j - basicstart1], rotatedcoeff, tmptmpct);
        if (giantct == 0) {
          giantct = new PhantomCiphertext();
          *giantct = tmptmpct;
        } else
          ckks->evaluator.add_inplace_reduced_error(*giantct, tmptmpct);
      }
    } else {
      for (int j = basicstart1; j <= totlen - i * gs1; j++) {
        rotation(coeff_logn, Nh, (-i) * gs1 * basicstep, fftcoeff[(i * gs1 + j) + totlen], rotatedcoeff);
        ckks->evaluator.multiply_vector_reduced_error(babyct[j - basicstart1], rotatedcoeff, tmptmpct);
        if (giantct == 0) {
          giantct = new PhantomCiphertext();
          *giantct = tmptmpct;
        } else
          ckks->evaluator.add_inplace_reduced_error(*giantct, tmptmpct);
      }
    }

    if (i != 0) {
      ckks->evaluator.rotate_vector(*giantct, (Nh + i * gs1 * basicstep) % Nh, *(ckks->galois_keys), tmptmpct);
      if (tmpct == 0) {
        tmpct = new PhantomCiphertext();
        *tmpct = tmptmpct;
      } else
        ckks->evaluator.add_inplace_reduced_error(*tmpct, tmptmpct);
    } else {
      if (tmpct == 0) {
        tmpct = new PhantomCiphertext();
        *tmpct = *giantct;
      } else
        ckks->evaluator.add_inplace_reduced_error(*tmpct, *giantct);
    }
  }

  rtncipher = *tmpct;
  delete[] babyct;
}

void Bootstrapper::rotated_bsgs_linear_transform_hoisting(
    PhantomCiphertext &rtncipher, PhantomCiphertext &cipher, int totlen, int basicstep, int coeff_logn, vector<vector<complex<double>>> fftcoeff) {
  int gs2 = giantstep(totlen + 1);
  int giantlast2 = floor((totlen + 0.0) / (gs2 + 0.0));

  PhantomCiphertext *babyct = new PhantomCiphertext[gs2];
  PhantomCiphertext *giantct = 0, *tmpct = 0;
  PhantomCiphertext tmptmpct;
  vector<complex<double>> rotatedcoeff;
  rotatedcoeff.reserve(Nh);

  PhantomPlaintext tmpplain;

  PhantomCiphertext addtmp1, addtmp2;

  if (!(tmpct == 0))
    delete tmpct;
  tmpct = 0;

  for (int i = 0; i < gs2; i++) {
    if (i == 0) {
      babyct[i] = cipher;
    } else {
      ckks->evaluator.rotate_vector(cipher, (Nh + i * basicstep) % Nh, *(ckks->galois_keys), babyct[i]);
    }
  }

  for (int i = 0; i <= giantlast2; i++) {
    if (!(giantct == 0))
      delete giantct;
    giantct = 0;
    if (i != giantlast2) {
      for (int j = 0; j < gs2; j++) {
        rotation(coeff_logn, Nh, (-i) * gs2 * basicstep, fftcoeff[i * gs2 + j], rotatedcoeff);
        ckks->evaluator.multiply_vector_reduced_error(babyct[j], rotatedcoeff, tmptmpct);
        if (giantct == 0) {
          giantct = new PhantomCiphertext();
          *giantct = tmptmpct;
        } else
          ckks->evaluator.add_inplace_reduced_error(*giantct, tmptmpct);
      }
    } else {
      for (int j = 0; j <= totlen - i * gs2; j++) {
        rotation(coeff_logn, Nh, (-i) * gs2 * basicstep, fftcoeff[i * gs2 + j], rotatedcoeff);
        ckks->evaluator.multiply_vector_reduced_error(babyct[j], rotatedcoeff, tmptmpct);
        if (giantct == 0) {
          giantct = new PhantomCiphertext();
          *giantct = tmptmpct;
        } else
          ckks->evaluator.add_inplace_reduced_error(*giantct, tmptmpct);
      }
    }

    if (i != 0) {
      ckks->evaluator.rotate_vector(*giantct, (Nh + i * gs2 * basicstep) % Nh, *(ckks->galois_keys), tmptmpct);
      if (tmpct == 0) {
        tmpct = new PhantomCiphertext;
        *tmpct = tmptmpct;
      } else
        ckks->evaluator.add_inplace_reduced_error(*tmpct, tmptmpct);
    } else {
      if (tmpct == 0) {
        tmpct = new PhantomCiphertext;
        *tmpct = *giantct;
      } else
        ckks->evaluator.add_inplace_reduced_error(*tmpct, *giantct);
    }
  }

  rtncipher = *tmpct;
  delete[] babyct;
}
void Bootstrapper::rotated_nobsgs_linear_transform(
    PhantomCiphertext &rtncipher, PhantomCiphertext &cipher, int totlen, int coeff_logn, vector<vector<complex<double>>> fftcoeff) {
  PhantomCiphertext *giantct = 0, *tmpct = 0;

  PhantomCiphertext tmptmpct;
  vector<complex<double>> rotatedcoeff;
  rotatedcoeff.reserve(Nh);

  PhantomPlaintext tmpplain;
  PhantomCiphertext addtmp1, addtmp2;

  for (int i = 0; i <= totlen; i++) {
    giantct = 0;
    rotation(coeff_logn, Nh, -i, fftcoeff[i], rotatedcoeff);
    ckks->evaluator.multiply_vector_reduced_error(cipher, rotatedcoeff, tmptmpct);
    giantct = new PhantomCiphertext();
    *giantct = tmptmpct;
    if (i != 0) {
      ckks->evaluator.rotate_vector(*giantct, i, *(ckks->galois_keys), tmptmpct);
      if (tmpct == 0) {
        tmpct = new PhantomCiphertext;
        *tmpct = tmptmpct;
      } else
        ckks->evaluator.add_inplace_reduced_error(*tmpct, tmptmpct);
    } else {
      if (tmpct == 0) {
        tmpct = new PhantomCiphertext;
        *tmpct = *giantct;
      } else
        ckks->evaluator.add_inplace_reduced_error(*tmpct, *giantct);
    }
    delete giantct;
  }
  rtncipher = *tmpct;
}

void Bootstrapper::sfl_one_depth(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  int totlen1 = (1 << logn) - 1;
  vector<vector<complex<double>>> fftcoeff1_ext(2 * totlen1 + 1);
  for (int i = 0; i < 2 * totlen1 + 1; i++) {
    fftcoeff1_ext[i].resize(2 * n);
    for (int j = 0; j < n; j++) {
      fftcoeff1_ext[i][j] = fftcoeff1[slot_index][i][j];
      fftcoeff1_ext[i][j + n] = fftcoeff1[slot_index][i][j];
    }
  }
  bsgs_linear_transform(rtncipher, cipher, totlen1, 1, logn + 1, fftcoeff1_ext);
}

void Bootstrapper::sfl_full_one_depth(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  int totlen2 = (1 << logn) - 1;
  rotated_bsgs_linear_transform(rtncipher, cipher, totlen2, 1, logn, fftcoeff2[slot_index]);
}

void Bootstrapper::sflinv_one_depth(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  int totlen1 = (1 << logn) - 1;
  rotated_bsgs_linear_transform(rtncipher, cipher, totlen1, 1, logn, invfftcoeff1[slot_index]);
}

void Bootstrapper::sflinv_one_depth_more_depth(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  int totlen1 = (1 << logn) - 1;
  rotated_nobsgs_linear_transform(rtncipher, cipher, totlen1, logn, invfftcoeff1[slot_index]);
}
void Bootstrapper::sfl(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  int split_point = floor(logn / 2.0);
  int totlen1 = (1 << split_point) - 1;
  int totlen2 = (1 << (logn - split_point)) - 1;

  PhantomCiphertext tmpct;
  bsgs_linear_transform(tmpct, cipher, totlen1, 1, logn + 1, fftcoeff1[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct);

  const auto &modulus = ckks->context->first_context_data().parms().coeff_modulus();
  auto curr_level = ckks->context->get_context_data(tmpct.params_id()).chain_depth();

  double mod_zero = (double)modulus[0].value();
  double curr_mod = (double)modulus[curr_level].value();

  vector<vector<complex<double>>> fftcoeff2_scale(2 * totlen2 + 1);
  for (int i = 0; i < 2 * totlen2 + 1; i++)
    fftcoeff2_scale[i].resize(2 * n);

  for (int i = 0; i < 2 * totlen2 + 1; i++) {
    for (int j = 0; j < 2 * n; j++) {
      fftcoeff2_scale[i][j] = fftcoeff2[slot_index][i][j] * curr_mod * mod_zero * final_scale / (tmpct.scale() * tmpct.scale() * initial_scale);
    }
  }

  int basicstep = (1 << split_point);
  bsgs_linear_transform(rtncipher, tmpct, totlen2, basicstep, logn + 1, fftcoeff2_scale);
  ckks->evaluator.rescale_to_next_inplace(rtncipher);
}

void Bootstrapper::sfl_full(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  int split_point = floor(logn / 2.0);
  int totlen1 = (1 << split_point) - 1;
  int totlen2 = (1 << (logn - split_point)) - 1;

  PhantomCiphertext tmpct;
  bsgs_linear_transform(tmpct, cipher, totlen1, 1, logn, fftcoeff1[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct);

  const auto modulus = ckks->context->first_context_data().parms().coeff_modulus();
  auto curr_level = ckks->context->get_context_data(tmpct.params_id()).chain_depth();

  double mod_zero = (double)modulus[0].value();
  double curr_mod = (double)modulus[curr_level].value();
  vector<vector<complex<double>>> fftcoeff2_scale(2 * totlen2 + 1);
  for (int i = 0; i < totlen2 + 1; i++)
    fftcoeff2_scale[i].resize(n);

  for (int i = 0; i < totlen2 + 1; i++) {
    for (int j = 0; j < n; j++) {
      fftcoeff2_scale[i][j] = fftcoeff2[slot_index][i][j] * curr_mod * mod_zero * final_scale / (tmpct.scale() * tmpct.scale() * initial_scale);
    }
  }

  int basicstep = (1 << split_point);
  rotated_bsgs_linear_transform(rtncipher, tmpct, totlen2, basicstep, logn, fftcoeff2_scale);
  ckks->evaluator.rescale_to_next_inplace(rtncipher);
}

void Bootstrapper::sflinv(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  int split_point = ceil(logn / 2.0);
  int totlen1 = (1 << split_point) - 1;
  int totlen2 = (1 << (logn - split_point)) - 1;

  PhantomCiphertext tmpct;
  int basicstep = (1 << (logn - split_point));
  rotated_bsgs_linear_transform(tmpct, cipher, totlen1, basicstep, logn, invfftcoeff1[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct);
  bsgs_linear_transform(rtncipher, tmpct, totlen2, 1, logn + 1, invfftcoeff2[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(rtncipher);
}

void Bootstrapper::sflinv_full(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  int split_point = ceil(logn / 2.0);
  int totlen1 = (1 << split_point) - 1;
  int totlen2 = (1 << (logn - split_point)) - 1;

  PhantomCiphertext tmpct;
  int basicstep = (1 << (logn - split_point));
  rotated_bsgs_linear_transform(tmpct, cipher, totlen1, basicstep, logn, invfftcoeff1[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct);
  bsgs_linear_transform(rtncipher, tmpct, totlen2, 1, logn, invfftcoeff2[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(rtncipher);
}
void Bootstrapper::sfl_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {  // not yet
  int div_part3 = floor(logn / 3.0);
  int div_part2 = floor((logn - div_part3) / 2.0);
  int div_part1 = logn - div_part3 - div_part2;

  int totlen1 = (1 << div_part1) - 1;
  int totlen2 = (1 << div_part2) - 1;
  int totlen3 = (1 << div_part3) - 1;

  int basicstep1 = 1;
  int basicstep2 = (1 << div_part1);
  int basicstep3 = (1 << (div_part1 + div_part2));

  PhantomCiphertext tmpct;
  bsgs_linear_transform(tmpct, cipher, totlen1, basicstep1, logn + 1, fftcoeff1[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct);

  PhantomCiphertext tmpct2;
  bsgs_linear_transform(tmpct2, tmpct, totlen2, basicstep2, logn + 1, fftcoeff2[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct2);

  const auto modulus = ckks->context->first_context_data().parms().coeff_modulus();
  auto curr_level = ckks->context->get_context_data(tmpct2.params_id()).chain_depth();

  double mod_zero = (double)modulus[0].value();
  double curr_mod = (double)modulus[curr_level].value();

  vector<vector<complex<double>>> fftcoeff3_scale(2 * totlen3 + 1);
  for (int i = 0; i < 2 * totlen3 + 1; i++)
    fftcoeff3_scale[i].resize(2 * n);

  for (int i = 0; i < 2 * totlen3 + 1; i++) {
    for (int j = 0; j < 2 * n; j++) {
      fftcoeff3_scale[i][j] = fftcoeff3[slot_index][i][j] * curr_mod * mod_zero * final_scale / (tmpct2.scale() * tmpct2.scale() * initial_scale);
    }
  }

  bsgs_linear_transform(rtncipher, tmpct2, totlen3, basicstep3, logn + 1, fftcoeff3_scale);
  ckks->evaluator.rescale_to_next_inplace(rtncipher);
}

void Bootstrapper::sfl_full_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  int div_part3 = floor(logn / 3.0);
  int div_part2 = floor((logn - div_part3) / 2.0);
  int div_part1 = logn - div_part3 - div_part2;

  int totlen1 = (1 << div_part1) - 1;
  int totlen2 = (1 << div_part2) - 1;
  int totlen3 = (1 << div_part3) - 1;

  int basicstep1 = 1;
  int basicstep2 = (1 << div_part1);
  int basicstep3 = (1 << (div_part1 + div_part2));

  PhantomCiphertext tmpct;
  bsgs_linear_transform(tmpct, cipher, totlen1, basicstep1, logn, fftcoeff1[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct);

  PhantomCiphertext tmpct2;
  bsgs_linear_transform(tmpct2, tmpct, totlen2, basicstep2, logn, fftcoeff2[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct2);

  const auto modulus = ckks->context->first_context_data().parms().coeff_modulus();
  auto curr_level = ckks->context->get_context_data(tmpct2.params_id()).chain_depth();

  double mod_zero = (double)modulus[0].value();
  double curr_mod = (double)modulus[curr_level].value();
  vector<vector<complex<double>>> fftcoeff3_scale(2 * totlen3 + 1);
  for (int i = 0; i < totlen2 + 1; i++)
    fftcoeff3_scale[i].resize(n);

  for (int i = 0; i < totlen2 + 1; i++) {
    for (int j = 0; j < n; j++) {
      fftcoeff3_scale[i][j] = fftcoeff3[slot_index][i][j] * curr_mod * mod_zero * final_scale / (tmpct2.scale() * tmpct2.scale() * initial_scale);
    }
  }

  rotated_bsgs_linear_transform(rtncipher, tmpct2, totlen3, basicstep3, logn, fftcoeff3_scale);
  ckks->evaluator.rescale_to_next_inplace(rtncipher);
}
void Bootstrapper::sfl_half_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {  // not yet
  int div_part3 = floor(logn / 3.0);
  int div_part2 = floor((logn - div_part3) / 2.0);
  int div_part1 = logn - div_part3 - div_part2;

  int totlen1 = (1 << div_part1) - 1;
  int totlen2 = (1 << div_part2) - 1;
  int totlen3 = (1 << div_part3) - 1;

  int basicstep1 = 1;
  int basicstep2 = (1 << div_part1);
  int basicstep3 = (1 << (div_part1 + div_part2));

  PhantomCiphertext tmpct;
  bsgs_linear_transform(tmpct, cipher, totlen1, basicstep1, logn + 1, fftcoeff1[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct);

  PhantomCiphertext tmpct2;
  bsgs_linear_transform(tmpct2, tmpct, totlen2, basicstep2, logn + 1, fftcoeff2[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct2);

  const auto modulus = ckks->context->first_context_data().parms().coeff_modulus();
  auto curr_level = ckks->context->get_context_data(tmpct2.params_id()).chain_depth();

  double mod_zero = (double)modulus[0].value();
  double curr_mod = (double)modulus[curr_level].value();
  vector<vector<complex<double>>> fftcoeff3_scale(2 * totlen3 + 1);
  for (int i = 0; i < 2 * totlen3 + 1; i++)
    fftcoeff3_scale[i].resize(2 * n);

  for (int i = 0; i < 2 * totlen3 + 1; i++) {
    for (int j = 0; j < 2 * n; j++) {
      fftcoeff3_scale[i][j] = fftcoeff3[slot_index][i][j] * curr_mod * mod_zero * final_scale / (2 * tmpct2.scale() * tmpct2.scale() * initial_scale);
    }
  }

  bsgs_linear_transform(rtncipher, tmpct2, totlen3, basicstep3, logn + 1, fftcoeff3_scale);
  ckks->evaluator.rescale_to_next_inplace(rtncipher);
}

void Bootstrapper::sfl_full_half_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  int div_part3 = floor(logn / 3.0);
  int div_part2 = floor((logn - div_part3) / 2.0);
  int div_part1 = logn - div_part3 - div_part2;

  int totlen1 = (1 << div_part1) - 1;
  int totlen2 = (1 << div_part2) - 1;
  int totlen3 = (1 << div_part3) - 1;

  int basicstep1 = 1;
  int basicstep2 = (1 << div_part1);
  int basicstep3 = (1 << (div_part1 + div_part2));

  PhantomCiphertext tmpct;
  bsgs_linear_transform(tmpct, cipher, totlen1, basicstep1, logn, fftcoeff1[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct);

  PhantomCiphertext tmpct2;
  bsgs_linear_transform(tmpct2, tmpct, totlen2, basicstep2, logn, fftcoeff2[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct2);

  const auto modulus = ckks->context->first_context_data().parms().coeff_modulus();
  auto curr_level = ckks->context->get_context_data(tmpct2.params_id()).chain_depth();

  double mod_zero = (double)modulus[0].value();
  double curr_mod = (double)modulus[curr_level].value();
  vector<vector<complex<double>>> fftcoeff3_scale(2 * totlen3 + 1);
  for (int i = 0; i < totlen2 + 1; i++)
    fftcoeff3_scale[i].resize(n);

  for (int i = 0; i < totlen2 + 1; i++) {
    for (int j = 0; j < n; j++) {
      fftcoeff3_scale[i][j] = fftcoeff3[slot_index][i][j] * curr_mod * mod_zero * final_scale / (2 * tmpct2.scale() * tmpct2.scale() * initial_scale);
    }
  }

  rotated_bsgs_linear_transform(rtncipher, tmpct2, totlen3, basicstep3, logn, fftcoeff3_scale);
  ckks->evaluator.rescale_to_next_inplace(rtncipher);
}

void Bootstrapper::sflinv_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {  // not yet
  int div_part1 = floor(logn / 3.0);
  int div_part2 = floor((logn - div_part1) / 2.0);
  int div_part3 = logn - div_part1 - div_part2;

  int totlen1 = (1 << div_part1) - 1;
  int totlen2 = (1 << div_part2) - 1;
  int totlen3 = (1 << div_part3) - 1;

  int basicstep1 = (1 << (logn - div_part1));
  int basicstep2 = (1 << (logn - div_part1 - div_part2));
  int basicstep3 = 1;

  PhantomCiphertext tmpct;
  rotated_bsgs_linear_transform(tmpct, cipher, totlen1, basicstep1, logn, invfftcoeff1[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct);

  PhantomCiphertext tmpct2;
  bsgs_linear_transform(tmpct2, tmpct, totlen2, basicstep2, logn, invfftcoeff2[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct2);

  bsgs_linear_transform(rtncipher, tmpct2, totlen3, basicstep3, logn + 1, invfftcoeff3[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(rtncipher);
}

void Bootstrapper::sflinv_full_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  int div_part1 = floor(logn / 3.0);
  int div_part2 = floor((logn - div_part1) / 2.0);
  int div_part3 = logn - div_part1 - div_part2;

  int totlen1 = (1 << div_part1) - 1;
  int totlen2 = (1 << div_part2) - 1;
  int totlen3 = (1 << div_part3) - 1;

  int basicstep1 = (1 << (logn - div_part1));
  int basicstep2 = (1 << (logn - div_part1 - div_part2));
  int basicstep3 = 1;

  PhantomCiphertext tmpct;
  rotated_bsgs_linear_transform(tmpct, cipher, totlen1, basicstep1, logn, invfftcoeff1[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct);
  PhantomCiphertext tmpct2;
  bsgs_linear_transform(tmpct2, tmpct, totlen2, basicstep2, logn, invfftcoeff2[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(tmpct2);
  bsgs_linear_transform(rtncipher, tmpct2, totlen3, basicstep3, logn, invfftcoeff3[slot_index]);
  ckks->evaluator.rescale_to_next_inplace(rtncipher);
}

void Bootstrapper::sfl_hoisting(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  int split_point = floor(logn / 2.0);
  int totlen1 = (1 << split_point) - 1;
  int totlen2 = (1 << (logn - split_point)) - 1;

  PhantomCiphertext tmpct;
  bsgs_linear_transform(tmpct, cipher, totlen1, 1, logn + 1, fftcoeff1[slot_index]);

  int basicstep = (1 << split_point);
  bsgs_linear_transform(rtncipher, tmpct, totlen2, basicstep, logn + 1, fftcoeff2[slot_index]);
}

void Bootstrapper::sfl_full_hoisting(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  int split_point = floor(logn / 2.0);
  int totlen1 = (1 << split_point) - 1;
  int totlen2 = (1 << (logn - split_point)) - 1;

  PhantomCiphertext tmpct;
  bsgs_linear_transform(tmpct, cipher, totlen1, 1, logn, fftcoeff1[slot_index]);

  int basicstep = (1 << split_point);
  rotated_bsgs_linear_transform(rtncipher, tmpct, totlen2, basicstep, logn, fftcoeff2[slot_index]);
}

void Bootstrapper::sflinv_hoisting(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  int split_point = ceil(logn / 2.0);
  int totlen1 = (1 << split_point) - 1;
  int totlen2 = (1 << (logn - split_point)) - 1;

  PhantomCiphertext tmpct;
  int basicstep = (1 << (logn - split_point));
  rotated_bsgs_linear_transform(tmpct, cipher, totlen1, basicstep, logn, invfftcoeff1[slot_index]);
  bsgs_linear_transform(rtncipher, tmpct, totlen2, 1, logn + 1, invfftcoeff2[slot_index]);
}

void Bootstrapper::sflinv_full_hoisting(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  int split_point = ceil(logn / 2.0);
  int totlen1 = (1 << split_point) - 1;
  int totlen2 = (1 << (logn - split_point)) - 1;

  PhantomCiphertext tmpct;
  int basicstep = (1 << (logn - split_point));
  rotated_bsgs_linear_transform(tmpct, cipher, totlen1, basicstep, logn, invfftcoeff1[slot_index]);
  bsgs_linear_transform(rtncipher, tmpct, totlen2, 1, logn, invfftcoeff2[slot_index]);
}

void Bootstrapper::coefftoslot(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3;
  sflinv(tmpct1, cipher);
  ckks->evaluator.complex_conjugate(tmpct1, *(ckks->galois_keys), tmpct2);
  ckks->evaluator.add_reduced_error(tmpct1, tmpct2, rtncipher);
}

void Bootstrapper::slottocoeff(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3;
  sfl(tmpct1, cipher);
  ckks->evaluator.rotate_vector(tmpct1, n, *(ckks->galois_keys), tmpct2);
  ckks->evaluator.add_reduced_error(tmpct1, tmpct2, rtncipher);
}

void Bootstrapper::coefftoslot_full(PhantomCiphertext &rtncipher1, PhantomCiphertext &rtncipher2, PhantomCiphertext &cipher) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3, tmpct4, tmpct5;
  sflinv_full(tmpct1, cipher);
  complex<double> iunit(0.0, 1.0);
  vector<complex<double>> tmpvec(Nh, 0);
  for (int i = 0; i < Nh; i++) {
    tmpvec[i] -= iunit;
  }
  PhantomPlaintext tmpplain;
  ckks->encoder.encode(tmpvec, 1.0, tmpplain);
  ckks->evaluator.mod_switch_to_inplace(tmpplain, tmpct1.params_id());
  ckks->evaluator.multiply_plain(tmpct1, tmpplain, tmpct2);

  ckks->evaluator.complex_conjugate(tmpct2, *(ckks->galois_keys), tmpct3);
  ckks->evaluator.complex_conjugate(tmpct1, *(ckks->galois_keys), tmpct4);
  ckks->evaluator.add_reduced_error(tmpct1, tmpct4, rtncipher1);
  ckks->evaluator.add_reduced_error(tmpct2, tmpct3, rtncipher2);
}
void Bootstrapper::slottocoeff_full(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher1, PhantomCiphertext &cipher2) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3;
  complex<double> iunit(0.0, 1.0);

  vector<complex<double>> tmpvec(Nh, 0);
  for (int i = 0; i < Nh; i++) {
    tmpvec[i] += iunit;
  }

  PhantomPlaintext tmpplain;
  ckks->encoder.encode(tmpvec, 1.0, tmpplain);

  ckks->evaluator.mod_switch_to_inplace(tmpplain, cipher2.params_id());
  ckks->evaluator.multiply_plain(cipher2, tmpplain, tmpct1);
  ckks->evaluator.add_reduced_error(cipher1, tmpct1, tmpct3);

  sfl_full(rtncipher, tmpct3);
}

// NOTE:
void Bootstrapper::coefftoslot_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3;
  sflinv_3(tmpct1, cipher);
  ckks->evaluator.complex_conjugate(tmpct1, *(ckks->galois_keys), tmpct2);
  ckks->evaluator.add_reduced_error(tmpct1, tmpct2, rtncipher);
}

// NOTE:
void Bootstrapper::slottocoeff_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3;
  sfl_3(tmpct1, cipher);
  ckks->evaluator.rotate_vector(tmpct1, n, *(ckks->galois_keys), tmpct2);
  ckks->evaluator.add_reduced_error(tmpct1, tmpct2, rtncipher);
}

void Bootstrapper::slottocoeff_half_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3;
  sfl_half_3(tmpct1, cipher);
  ckks->evaluator.rotate_vector(tmpct1, n, *(ckks->galois_keys), tmpct2);
  ckks->evaluator.add_reduced_error(tmpct1, tmpct2, rtncipher);
}

void Bootstrapper::coefftoslot_full_3(PhantomCiphertext &rtncipher1, PhantomCiphertext &rtncipher2, PhantomCiphertext &cipher) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3, tmpct4, tmpct5;
  sflinv_full_3(tmpct1, cipher);
  complex<double> iunit(0.0, 1.0);
  vector<complex<double>> tmpvec(Nh, 0);
  for (int i = 0; i < Nh; i++) {
    tmpvec[i] -= iunit;
  }
  PhantomPlaintext tmpplain;
  ckks->encoder.encode(tmpvec, 1.0, tmpplain);
  ckks->evaluator.mod_switch_to_inplace(tmpplain, tmpct1.params_id());
  ckks->evaluator.multiply_plain(tmpct1, tmpplain, tmpct2);

  ckks->evaluator.complex_conjugate(tmpct2, *(ckks->galois_keys), tmpct3);
  ckks->evaluator.complex_conjugate(tmpct1, *(ckks->galois_keys), tmpct4);

  ckks->evaluator.add_reduced_error(tmpct1, tmpct4, rtncipher1);
  ckks->evaluator.add_reduced_error(tmpct2, tmpct3, rtncipher2);
}
void Bootstrapper::slottocoeff_full_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher1, PhantomCiphertext &cipher2) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3;
  complex<double> iunit(0.0, 1.0);

  vector<complex<double>> tmpvec(Nh, 0);
  for (int i = 0; i < Nh; i++) {
    tmpvec[i] += iunit;
  }

  PhantomPlaintext tmpplain;
  ckks->encoder.encode(tmpvec, 1.0, tmpplain);

  ckks->evaluator.mod_switch_to_inplace(tmpplain, cipher2.params_id());
  ckks->evaluator.multiply_plain(cipher2, tmpplain, tmpct1);
  ckks->evaluator.add_reduced_error(cipher1, tmpct1, tmpct3);

  sfl_full_3(rtncipher, tmpct3);
}
void Bootstrapper::slottocoeff_full_half_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher1, PhantomCiphertext &cipher2) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3;
  complex<double> iunit(0.0, 1.0);

  vector<complex<double>> tmpvec(Nh, 0);
  for (int i = 0; i < Nh; i++) {
    tmpvec[i] += iunit;
  }

  PhantomPlaintext tmpplain;
  ckks->encoder.encode(tmpvec, 1.0, tmpplain);

  ckks->evaluator.mod_switch_to_inplace(tmpplain, cipher2.params_id());
  ckks->evaluator.multiply_plain(cipher2, tmpplain, tmpct1);
  ckks->evaluator.add_reduced_error(cipher1, tmpct1, tmpct3);

  sfl_full_half_3(rtncipher, tmpct3);
}
void Bootstrapper::coefftoslot_hoisting(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3;
  sflinv_hoisting(tmpct1, cipher);
  ckks->evaluator.complex_conjugate(tmpct1, *(ckks->galois_keys), tmpct2);
  ckks->evaluator.add_reduced_error(tmpct1, tmpct2, rtncipher);
}

void Bootstrapper::slottocoeff_hoisting(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3;
  sfl_hoisting(tmpct1, cipher);
  ckks->evaluator.rotate_vector(tmpct1, n, *(ckks->galois_keys), tmpct2);
  ckks->evaluator.add_reduced_error(tmpct1, tmpct2, rtncipher);
}

void Bootstrapper::coefftoslot_full_hoisting(PhantomCiphertext &rtncipher1, PhantomCiphertext &rtncipher2, PhantomCiphertext &cipher) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3, tmpct4, tmpct5;
  sflinv_full_hoisting(tmpct1, cipher);
  complex<double> iunit(0.0, 1.0);
  vector<complex<double>> tmpvec(Nh, 0);
  for (int i = 0; i < Nh; i++) {
    tmpvec[i] -= iunit;
  }
  PhantomPlaintext tmpplain;
  ckks->encoder.encode(tmpvec, 1.0, tmpplain);
  ckks->evaluator.mod_switch_to_inplace(tmpplain, tmpct1.params_id());
  ckks->evaluator.multiply_plain(tmpct1, tmpplain, tmpct2);

  ckks->evaluator.complex_conjugate(tmpct2, *(ckks->galois_keys), tmpct3);
  ckks->evaluator.complex_conjugate(tmpct1, *(ckks->galois_keys), tmpct4);
  ckks->evaluator.add_reduced_error(tmpct1, tmpct4, rtncipher1);
  ckks->evaluator.add_reduced_error(tmpct2, tmpct3, rtncipher2);
}
void Bootstrapper::slottocoeff_full_hoisting(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher1, PhantomCiphertext &cipher2) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3;
  complex<double> iunit(0.0, 1.0);

  vector<complex<double>> tmpvec(Nh, 0);
  for (int i = 0; i < Nh; i++) {
    tmpvec[i] += iunit;
  }

  PhantomPlaintext tmpplain;
  ckks->encoder.encode(tmpvec, 1.0, tmpplain);

  ckks->evaluator.mod_switch_to_inplace(tmpplain, cipher2.params_id());
  ckks->evaluator.multiply_plain(cipher2, tmpplain, tmpct1);
  ckks->evaluator.add_reduced_error(cipher1, tmpct1, tmpct3);

  sfl_full_hoisting(rtncipher, tmpct3);
}

void Bootstrapper::coefftoslot_one_depth(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3;
  sflinv_one_depth(tmpct1, cipher);
  complex<double> iunit(0.0, 1.0);
  vector<complex<double>> tmpvec(Nh, 0);
  int repeatcount = Nh / (2 * n);

  for (int i = 0; i < repeatcount; i++) {
    for (int j = 0; j < n; j++) {
      tmpvec[i * (2 * n) + j] = 0.5;
      tmpvec[i * (2 * n) + j + n] -= 0.5 * iunit;
    }
  }

  ckks->evaluator.multiply_vector_reduced_error(tmpct1, tmpvec, tmpct2);
  ckks->evaluator.complex_conjugate(tmpct2, *(ckks->galois_keys), tmpct1);
  ckks->evaluator.add_reduced_error(tmpct1, tmpct2, rtncipher);
}
void Bootstrapper::slottocoeff_one_depth(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3;
  complex<double> iunit(0.0, 1.0);

  vector<complex<double>> tmpvec(Nh, 0);
  int repeatcount = Nh / (2 * n);
  for (int i = 0; i < repeatcount; i++) {
    for (int j = 0; j < n; j++) {
      tmpvec[i * (2 * n) + j] = 1.0;
      tmpvec[i * (2 * n) + j + n] += iunit;
    }
  }
  ckks->evaluator.multiply_vector_reduced_error(cipher, tmpvec, tmpct3);
  ckks->evaluator.rotate_vector(tmpct3, n, *(ckks->galois_keys), tmpct2);
  ckks->evaluator.add_reduced_error(tmpct3, tmpct2, tmpct1);

  sfl_one_depth(rtncipher, tmpct1);
}

void Bootstrapper::coefftoslot_full_one_depth(PhantomCiphertext &rtncipher1, PhantomCiphertext &rtncipher2, PhantomCiphertext &cipher) {}

void Bootstrapper::slottocoeff_full_one_depth(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher1, PhantomCiphertext &cipher2) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3;
  complex<double> iunit(0.0, 1.0);

  vector<complex<double>> tmpvec(Nh, 0);
  int repeatcount = Nh / (2 * n);
  for (int i = 0; i < repeatcount; i++) {
    for (int j = 0; j < n; j++) {
      tmpvec[i * (2 * n) + j] += iunit;
      tmpvec[i * (2 * n) + j + n] += iunit;
    }
  }
  PhantomPlaintext tmpplain;
  ckks->encoder.encode(tmpvec, 1.0, tmpplain);
  ckks->evaluator.mod_switch_to_inplace(tmpplain, cipher2.params_id());
  ckks->evaluator.multiply_plain(cipher2, tmpplain, tmpct1);
  ckks->evaluator.add_reduced_error(cipher1, tmpct1, tmpct3);

  sfl_full_one_depth(rtncipher, tmpct3);
}

void Bootstrapper::coefftoslot_full_mul_first(PhantomCiphertext &rtncipher1, PhantomCiphertext &rtncipher2, PhantomCiphertext &cipher) {
  PhantomCiphertext tmpct1, tmpct2, tmpct3, tmpct4, tmpct5;
  sflinv_one_depth_more_depth(tmpct1, cipher);
  complex<double> iunit(0.0, 1.0);
  vector<complex<double>> tmpvec(Nh, 0);
  int repeatcount = Nh / (2 * n);

  for (int i = 0; i < repeatcount; i++) {
    for (int j = 0; j < n; j++) {
      tmpvec[i * (2 * n) + j] = -0.5 * iunit;
      tmpvec[i * (2 * n) + j + n] = -0.5 * iunit;
    }
  }

  ckks->evaluator.multiply_vector_reduced_error(tmpct1, tmpvec, tmpct2);
  vector<complex<double>> tmpvec_1(Nh, 0);

  for (int i = 0; i < repeatcount; i++) {
    for (int j = 0; j < n; j++) {
      tmpvec_1[i * (2 * n) + j] = 0.5;
      tmpvec_1[i * (2 * n) + j + n] = 0.5;
    }
  }

  ckks->evaluator.multiply_vector_reduced_error(tmpct1, tmpvec_1, tmpct3);
  ckks->evaluator.complex_conjugate(tmpct2, *(ckks->galois_keys), tmpct4);
  ckks->evaluator.complex_conjugate(tmpct3, *(ckks->galois_keys), tmpct5);
  ckks->evaluator.add_reduced_error(tmpct3, tmpct5, rtncipher1);
  ckks->evaluator.add_reduced_error(tmpct2, tmpct4, rtncipher2);
}

__global__ void kernel_modraise_inplace(uint64_t *poly_dest, const uint64_t *poly_src_zero, const uint64_t modulus,
                                        const uint64_t minus_q0_modulus, size_t N, uint64_t q0) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < N) {
    auto q = modulus;
    poly_dest[i] = poly_src_zero[i] % q;
    if (poly_src_zero[i] > (q0 >> 1)) {
      poly_dest[i] += minus_q0_modulus;
      poly_dest[i] -= (poly_dest[i] >= q) ? q : 0;
    }
  }
}

void Bootstrapper::modraise_inplace(PhantomCiphertext &cipher) {
  if (cipher.size() != 2) {
    throw invalid_argument("Ciphertexts of size 2 are supported only!");
  }

  if (cipher.coeff_modulus_size() != 1) {
    throw invalid_argument("Ciphertexts in the lowest level are supported only!");
  }

  if (cipher.is_ntt_form()) {
    ckks->evaluator.transform_from_ntt_inplace(cipher);
  }

  const auto &stream = phantom::util::global_variables::default_stream->get_stream();
  auto N = cipher.poly_modulus_degree();
  const auto &context_data = ckks->context->first_context_data();
  const auto &modulus = context_data.parms().coeff_modulus();
  auto mod_count = modulus.size();
  auto rns_coeff_count = mod_count * N;

  uint64_t q0 = modulus[0].value();
  vector<uint64_t> minus_q0(mod_count);
  minus_q0[0] = 0;

  for (size_t l = 1; l < mod_count; l++) {
    minus_q0[l] = modulus[l].value() - q0 % modulus[l].value();
  }

  // Resize to the full level.
  cipher.resize(*ckks->context, context_data.chain_index(), stream);

  auto ciphertext_size = cipher.size();
  PhantomCiphertext cipher_copy = cipher;

  uint64_t gridDimGlb = N / blockDimGlb.x;
  for (size_t poly_idx = 0; poly_idx < ciphertext_size; poly_idx++) {
    const auto rns_poly_src = cipher_copy.data() + poly_idx * rns_coeff_count;
    const auto rns_poly_dest = cipher.data() + poly_idx * rns_coeff_count;

    for (size_t j = 0; j < mod_count; j++) {
      const auto poly_dest = rns_poly_dest + j * N;
      kernel_modraise_inplace<<<gridDimGlb, blockDimGlb, 0, stream>>>(
          poly_dest, rns_poly_src, modulus[j].value(), minus_q0[j], N, q0);
    }
  }

  ckks->evaluator.transform_to_ntt_inplace(cipher);
}

void Bootstrapper::bootstrap_sparse(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  std::cout << "Modulus Raising..." << endl;
  modraise_inplace(cipher);

  const auto modulus = ckks->context->first_context_data().parms().coeff_modulus();
  cipher.scale() = ((double)modulus[0].value());

  std::cout << "Subsum..." << endl;
  PhantomCiphertext rot;
  for (long i = logn; i < logNh; ++i) {
    ckks->evaluator.rotate_vector(cipher, (1 << i), *(ckks->galois_keys), rot);
    // ckks->evaluator.add_inplace_original(cipher, rot);
    ckks->evaluator.add_inplace(cipher, rot);
  }

  PhantomCiphertext rtn;
  if (logn == 0) {
    vector<complex<double>> cts_vec(Nh, 0);
    for (int i = 0; i < Nh; i++) {
      if (i % 2 == 0)
        cts_vec[i] = 1.0 / (2.0 * boundary_K * (1 << logNh));
      else
        cts_vec[i] = -complex<double>(0, 1.0) / (2.0 * boundary_K * (1 << logNh));
    }
    ckks->evaluator.multiply_vector_reduced_error(cipher, cts_vec, rtn);
    ckks->evaluator.rescale_to_next_inplace(rtn);

    PhantomCiphertext conjrtn;
    ckks->evaluator.complex_conjugate(rtn, *(ckks->galois_keys), conjrtn);
    ckks->evaluator.add_inplace_reduced_error(rtn, conjrtn);
  }

  else {
    std::cout << "Coefftoslot..." << endl;
    coefftoslot(rtn, cipher);
  }

  std::cout << "Modular reduction..." << endl;
  PhantomCiphertext modrtn;
  mod_reducer->modular_reduction(modrtn, rtn);
  std::cout << "mod end" << endl;

  if (logn == 0) {
    const auto modulus = ckks->context->first_context_data().parms().coeff_modulus();
    auto curr_level = ckks->context->get_context_data(modrtn.params_id()).chain_depth();

    double mod_zero = (double)modulus[0].value();
    double curr_mod = (double)modulus[curr_level].value();

    double scale_adj = curr_mod * mod_zero * final_scale / (modrtn.scale() * modrtn.scale() * initial_scale);

    vector<complex<double>> stc_vec(Nh, 0);
    for (int i = 0; i < Nh; i++) {
      if (i % 2 == 0)
        stc_vec[i] = scale_adj;
      else
        stc_vec[i] = complex<double>(0, 1.0) * scale_adj;
    }
    ckks->evaluator.multiply_vector_reduced_error(modrtn, stc_vec, rtncipher);
    ckks->evaluator.rescale_to_next_inplace(rtncipher);

    PhantomCiphertext rotrtncipher;
    ckks->evaluator.rotate_vector(rtncipher, 1, *(ckks->galois_keys), rotrtncipher);
    ckks->evaluator.add_inplace_reduced_error(rtncipher, rotrtncipher);
  }

  else {
    std::cout << "Slottocoeff..." << endl;
    slottocoeff(rtncipher, modrtn);
  }
  rtncipher.scale() = final_scale;
}

void Bootstrapper::bootstrap_full(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  // // for debugging
  // chrono::high_resolution_clock::time_point start, end;
  // chrono::microseconds diff;
  // start = chrono::high_resolution_clock::now();

  std::cout << "Modulus Raising..." << endl;
  modraise_inplace(cipher);

  // // for debugging
  // end = chrono::high_resolution_clock::now();
  // diff = chrono::duration_cast<chrono::milliseconds>(end - start);
  // std::cout << "mod raise time : " << diff.count() / 1000 << " ms" << endl;

  const auto modulus = ckks->context->first_context_data().parms().coeff_modulus();
  cipher.scale() = ((double)modulus[0].value());

  // // for debugging
  // start = chrono::high_resolution_clock::now();

  std::cout << "Coefftoslot...2" << endl;
  PhantomCiphertext rtn1, rtn2;
  coefftoslot_full(rtn1, rtn2, cipher);

  // // for debugging
  // end = chrono::high_resolution_clock::now();
  // diff = chrono::duration_cast<chrono::milliseconds>(end - start);
  // std::cout << "coefftoslot time : " << diff.count() / 1000 << " ms" << endl;
  // start = chrono::high_resolution_clock::now();

  std::cout << "Modular reduction..." << endl;
  PhantomCiphertext modrtn1, modrtn2;
  mod_reducer->modular_reduction(modrtn1, rtn1);
  mod_reducer->modular_reduction(modrtn2, rtn2);

  // // for debugging
  // end = chrono::high_resolution_clock::now();
  // diff = chrono::duration_cast<chrono::milliseconds>(end - start);
  // std::cout << "mod reduction time : " << diff.count() / 1000 << " ms" << endl;
  // start = chrono::high_resolution_clock::now();

  std::cout << "Slottocoeff..." << endl;
  slottocoeff_full(rtncipher, modrtn1, modrtn2);

  // // for debugging
  // end = chrono::high_resolution_clock::now();
  // diff = chrono::duration_cast<chrono::milliseconds>(end - start);
  // std::cout << "slottocoeff time : " << diff.count() / 1000 << " ms" << endl;

  rtncipher.scale() = final_scale;
}

void Bootstrapper::bootstrap_sparse_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  // ModRaise
  auto timer = Timer();

  std::cout << "Modulus Raising... ";
  modraise_inplace(cipher);

  timer.stop();
  std::cout << timer.duration<milliseconds>() << "ms\n";
  // ckks->print_decrypted_ct(cipher, 10);

  const auto modulus = ckks->context->first_context_data().parms().coeff_modulus();
  cipher.scale() = ((double)modulus[0].value());

  timer.start();

  std::cout << "Subsum... ";
  PhantomCiphertext rot;
  for (auto i = logn; i < logNh; i++) {
    ckks->evaluator.rotate_vector(cipher, (1 << i), *(ckks->galois_keys), rot);
    ckks->evaluator.add_inplace(cipher, rot);
  }

  timer.stop();
  std::cout << timer.duration<milliseconds>() << "ms\n";
  // ckks->print_decrypted_ct(cipher, 10);

  PhantomCiphertext rtn;
  if (logn == 0) {
    vector<complex<double>> cts_vec(Nh, 0);
    for (int i = 0; i < Nh; i++) {
      if (i % 2 == 0)
        cts_vec[i] = 1.0 / (2.0 * boundary_K * (1 << logNh));
      else
        cts_vec[i] = -complex<double>(0, 1.0) / (2.0 * boundary_K * (1 << logNh));
    }
    ckks->evaluator.multiply_vector_reduced_error(cipher, cts_vec, rtn);
    ckks->evaluator.rescale_to_next_inplace(rtn);

    PhantomCiphertext conjrtn;
    ckks->evaluator.complex_conjugate(rtn, *(ckks->galois_keys), conjrtn);
    ckks->evaluator.add_inplace_reduced_error(rtn, conjrtn);
  }

  // Coefficient to Slots
  else {
    timer.start();

    std::cout << "Coeff-to-slot... ";
    coefftoslot_3(rtn, cipher);

    timer.stop();
    std::cout << timer.duration<milliseconds>() << "ms\n";
    // ckks->print_decrypted_ct(rtn, 10);
  }

  // Modular Reduction
  timer.start();

  std::cout << "Modular reduction... ";
  PhantomCiphertext modrtn;
  mod_reducer->modular_reduction(modrtn, rtn);

  timer.stop();
  std::cout << timer.duration<milliseconds>() << "ms\n";
  // ckks->print_decrypted_ct(modrtn, 10);

  if (logn == 0) {
    const auto modulus = ckks->context->first_context_data().parms().coeff_modulus();
    auto curr_level = ckks->context->get_context_data(modrtn.params_id()).chain_depth();

    double mod_zero = (double)modulus[0].value();
    double curr_mod = (double)modulus[curr_level].value();

    double scale_adj = curr_mod * mod_zero * final_scale / (modrtn.scale() * modrtn.scale() * initial_scale);

    vector<complex<double>> stc_vec(Nh, 0);
    for (int i = 0; i < Nh; i++) {
      if (i % 2 == 0)
        stc_vec[i] = scale_adj;
      else
        stc_vec[i] = complex<double>(0, 1.0) * scale_adj;
    }
    ckks->evaluator.multiply_vector_reduced_error(modrtn, stc_vec, rtncipher);
    ckks->evaluator.rescale_to_next_inplace(rtncipher);

    PhantomCiphertext rotrtncipher;
    ckks->evaluator.rotate_vector(rtncipher, 1, *(ckks->galois_keys), rotrtncipher);
    ckks->evaluator.add_inplace_reduced_error(rtncipher, rotrtncipher);
  }

  // Slots to Coefficients
  else {
    timer.start();

    std::cout << "Slot-to-coeff... ";
    slottocoeff_3(rtncipher, modrtn);

    timer.stop();
    std::cout << timer.duration<milliseconds>() << "ms\n";
    // ckks->print_decrypted_ct(rtncipher, 10);
  }

  rtncipher.scale() = final_scale;
}

void Bootstrapper::bootstrap_full_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  std::cout << "Modulus Raising..." << endl;
  modraise_inplace(cipher);

  const auto modulus = ckks->context->first_context_data().parms().coeff_modulus();
  cipher.scale() = ((double)modulus[0].value());

  std::cout << "Coefftoslot...4" << endl;
  PhantomCiphertext rtn1, rtn2;
  coefftoslot_full_3(rtn1, rtn2, cipher);

  std::cout << "Modular reduction..." << endl;
  PhantomCiphertext modrtn1, modrtn2;
  mod_reducer->modular_reduction(modrtn1, rtn1);
  mod_reducer->modular_reduction(modrtn2, rtn2);

  std::cout << "Slottocoeff..." << endl;
  slottocoeff_full_3(rtncipher, modrtn1, modrtn2);

  rtncipher.scale() = final_scale;
}

void Bootstrapper::bootstrap_sparse_real_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  std::cout << "Modulus Raising..." << endl;
  modraise_inplace(cipher);

  const auto modulus = ckks->context->first_context_data().parms().coeff_modulus();
  cipher.scale() = ((double)modulus[0].value());

  std::cout << "Subsum..." << endl;
  PhantomCiphertext rot;
  for (long i = logn; i < logNh; ++i) {
    ckks->evaluator.rotate_vector(cipher, (1 << i), *(ckks->galois_keys), rot);
    // ckks->evaluator.add_inplace_original(cipher, rot);
    ckks->evaluator.add_inplace(cipher, rot);
  }

  PhantomCiphertext rtn;
  if (logn == 0) {
    vector<complex<double>> cts_vec(Nh, 0);
    for (int i = 0; i < Nh; i++) {
      if (i % 2 == 0)
        cts_vec[i] = 1.0 / (2.0 * boundary_K * (1 << logNh));
      else
        cts_vec[i] = -complex<double>(0, 1.0) / (2.0 * boundary_K * (1 << logNh));
    }
    ckks->evaluator.multiply_vector_reduced_error(cipher, cts_vec, rtn);
    ckks->evaluator.rescale_to_next_inplace(rtn);

    PhantomCiphertext conjrtn;
    ckks->evaluator.complex_conjugate(rtn, *(ckks->galois_keys), conjrtn);
    ckks->evaluator.add_inplace_reduced_error(rtn, conjrtn);
  }

  else {
    std::cout << "Coefftoslot...5" << endl;
    coefftoslot_3(rtn, cipher);
  }

  std::cout << "Modular reduction..." << endl;
  PhantomCiphertext modrtn;
  mod_reducer->modular_reduction(modrtn, rtn);

  if (logn == 0) {
    const auto modulus = ckks->context->first_context_data().parms().coeff_modulus();
    auto curr_level = ckks->context->get_context_data(modrtn.params_id()).chain_depth();

    double mod_zero = (double)modulus[0].value();
    double curr_mod = (double)modulus[curr_level].value();

    double scale_adj = curr_mod * mod_zero * final_scale / (modrtn.scale() * modrtn.scale() * initial_scale);

    vector<complex<double>> stc_vec(Nh, 0);
    for (int i = 0; i < Nh; i++) {
      if (i % 2 == 0)
        stc_vec[i] = scale_adj;
      else
        stc_vec[i] = complex<double>(0, 1.0) * scale_adj;
    }
    ckks->evaluator.multiply_vector_reduced_error(modrtn, stc_vec, rtncipher);
    ckks->evaluator.rescale_to_next_inplace(rtncipher);

    PhantomCiphertext rotrtncipher;
    ckks->evaluator.rotate_vector(rtncipher, 1, *(ckks->galois_keys), rotrtncipher);
    ckks->evaluator.add_inplace_reduced_error(rtncipher, rotrtncipher);
  }

  else {
    std::cout << "Slottocoeff..." << endl;
    slottocoeff_half_3(rtncipher, modrtn);
  }
  rtncipher.scale() = final_scale;
  PhantomCiphertext conjct;
  ckks->evaluator.complex_conjugate(rtncipher, *(ckks->galois_keys), conjct);
  ckks->evaluator.add_inplace_reduced_error(rtncipher, conjct);
}

void Bootstrapper::bootstrap_full_real_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  std::cout << "Modulus Raising..." << endl;
  modraise_inplace(cipher);

  const auto &modulus = ckks->context->first_context_data().parms().coeff_modulus();
  cipher.scale() = ((double)modulus[0].value());

  std::cout << "Coefftoslot...6" << endl;
  PhantomCiphertext rtn1, rtn2;
  coefftoslot_full_3(rtn1, rtn2, cipher);

  std::cout << "Modular reduction..." << endl;
  PhantomCiphertext modrtn1, modrtn2;
  mod_reducer->modular_reduction(modrtn1, rtn1);
  mod_reducer->modular_reduction(modrtn2, rtn2);

  std::cout << "Slottocoeff..." << endl;
  slottocoeff_full_half_3(rtncipher, modrtn1, modrtn2);

  rtncipher.scale() = final_scale;
  PhantomCiphertext conjct;
  ckks->evaluator.complex_conjugate(rtncipher, *(ckks->galois_keys), conjct);
  ckks->evaluator.add_inplace_reduced_error(rtncipher, conjct);
}

void Bootstrapper::bootstrap_sparse_hoisting(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  std::cout << "Modulus Raising..." << endl;
  modraise_inplace(cipher);

  const auto &modulus = ckks->context->first_context_data().parms().coeff_modulus();
  cipher.scale() = ((double)modulus[0].value());

  std::cout << "Subsum..." << endl;
  PhantomCiphertext rot;
  for (long i = logn; i < logNh; ++i) {
    ckks->evaluator.rotate_vector(cipher, (1 << i), *(ckks->galois_keys), rot);
    // ckks->evaluator.add_inplace_original(cipher, rot);
    ckks->evaluator.add_inplace(cipher, rot);
  }

  PhantomPlaintext tmpplain;

  ckks->encoder.encode(1.0 / (boundary_K * (1 << (logNh - logn))), cipher.scale(), tmpplain);  // scale of tmpplain = Delta
  ckks->evaluator.mod_switch_to_inplace(tmpplain, cipher.params_id());
  ckks->evaluator.multiply_plain_inplace(cipher, tmpplain);  // scale of cipher = Delta^2

  // Lazy rescaling from now on: No rescaling!!
  // NO NEED TO USE "ORIGINAL" METHODS OF EVALUATOR FROM NOW ON, AS THE SCALING FACTOR IS Delta^2

  std::cout << "Coefftoslot...7" << endl;
  PhantomCiphertext rtn;
  coefftoslot_hoisting(rtn, cipher);

  std::cout << "Modular reduction..." << endl;
  PhantomCiphertext modrtn;
  mod_reducer->modular_reduction(modrtn, rtn);

  std::cout << "Slottocoeff..." << endl;
  slottocoeff_hoisting(rtncipher, modrtn);
}
void Bootstrapper::bootstrap_full_hoisting(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  std::cout << "Modulus Raising..." << endl;
  modraise_inplace(cipher);

  const auto &modulus = ckks->context->first_context_data().parms().coeff_modulus();
  cipher.scale() = ((double)modulus[0].value());

  PhantomPlaintext tmpplain;

  ckks->encoder.encode(1.0 / boundary_K, cipher.scale(), tmpplain);  // scale of tmpplain = Delta
  ckks->evaluator.mod_switch_to_inplace(tmpplain, cipher.params_id());
  ckks->evaluator.multiply_plain_inplace(cipher, tmpplain);  // scale of cipher = Delta^2

  // Lazy rescaling from now on: No rescaling!!
  // NO NEED TO USE "ORIGINAL" METHODS OF EVALUATOR FROM NOW ON, AS THE SCALING FACTOR IS Delta^2

  std::cout << "Coefftoslot...8" << endl;
  PhantomCiphertext rtn1, rtn2;
  coefftoslot_full_hoisting(rtn1, rtn2, cipher);

  std::cout << "Modular reduction..." << endl;
  PhantomCiphertext modrtn1, modrtn2;
  mod_reducer->modular_reduction(modrtn1, rtn1);
  mod_reducer->modular_reduction(modrtn2, rtn2);

  std::cout << "Slottocoeff..." << endl;
  slottocoeff_full_hoisting(rtncipher, modrtn1, modrtn2);
}

void Bootstrapper::bootstrap_one_depth(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  std::cout << "Modulus Raising..." << endl;
  modraise_inplace(cipher);

  const auto &modulus = ckks->context->first_context_data().parms().coeff_modulus();
  cipher.scale() = ((double)modulus[0].value());

  std::cout << "Subsum..." << endl;
  PhantomCiphertext rot;
  for (long i = logn; i < logNh; ++i) {
    ckks->evaluator.rotate_vector(cipher, (1 << i), *(ckks->galois_keys), rot);
    // ckks->evaluator.add_inplace_original(cipher, rot);
    ckks->evaluator.add_inplace(cipher, rot);
  }

  PhantomPlaintext tmpplain;

  ckks->encoder.encode(1.0 / (boundary_K * (1 << (logNh - logn))), cipher.scale(), tmpplain);  // scale of tmpplain = Delta
  ckks->evaluator.mod_switch_to_inplace(tmpplain, cipher.params_id());
  ckks->evaluator.multiply_plain_inplace(cipher, tmpplain);  // scale of cipher = Delta^2

  // Lazy rescaling from now on: No rescaling!!
  // NO NEED TO USE "ORIGINAL" METHODS OF EVALUATOR FROM NOW ON, AS THE SCALING FACTOR IS Delta^2

  std::cout << "Coefftoslot...9" << endl;
  PhantomCiphertext rtn;
  coefftoslot_one_depth(rtn, cipher);

  std::cout << "Modular reduction..." << endl;
  PhantomCiphertext modrtn;
  mod_reducer->modular_reduction(modrtn, rtn);

  std::cout << "Slottocoeff..." << endl;
  slottocoeff_one_depth(rtncipher, modrtn);
}

void Bootstrapper::bootstrap_more_depth(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  std::cout << "Modulus Raising..." << endl;
  modraise_inplace(cipher);

  const auto &modulus = ckks->context->first_context_data().parms().coeff_modulus();
  cipher.scale() = ((double)modulus[0].value());

  PhantomPlaintext tmpplain;

  ckks->encoder.encode(1.0 / boundary_K, cipher.scale(), tmpplain);  // scale of tmpplain = Delta
  ckks->evaluator.mod_switch_to_inplace(tmpplain, cipher.params_id());
  ckks->evaluator.multiply_plain_inplace(cipher, tmpplain);  // scale of cipher = Delta^2

  // Lazy rescaling from now on: No rescaling!!
  // NO NEED TO USE "ORIGINAL" METHODS OF EVALUATOR FROM NOW ON, AS THE SCALING FACTOR IS Delta^2

  std::cout << "Coefftoslot...10" << endl;
  PhantomCiphertext rtn1, rtn2;
  coefftoslot_full_mul_first(rtn1, rtn2, cipher);

  std::cout << "Modular reduction..." << endl;
  PhantomCiphertext modrtn1, modrtn2;
  mod_reducer->modular_reduction(modrtn1, rtn1);
  mod_reducer->modular_reduction(modrtn2, rtn2);

  std::cout << "Slottocoeff..." << endl;
  slottocoeff_full_one_depth(rtncipher, modrtn1, modrtn2);
}

void Bootstrapper::bootstrap(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  initial_scale = cipher.scale();
  if (logn == logNh)
    bootstrap_full(rtncipher, cipher);
  else
    bootstrap_sparse(rtncipher, cipher);
}

void Bootstrapper::bootstrap_inplace(PhantomCiphertext &cipher) {
  PhantomCiphertext rtncipher;
  bootstrap(rtncipher, cipher);
  cipher = rtncipher;
}

void Bootstrapper::bootstrap_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  if (rtncipher.chain_index()) {
    throw invalid_argument("Return cipher should initially be a new ciphertext.");
  }

  initial_scale = cipher.scale();
  if (logn == logNh)
    bootstrap_full_3(rtncipher, cipher);
  else
    bootstrap_sparse_3(rtncipher, cipher);

  cudaStreamSynchronize(cipher.data_ptr().get_stream());
}

void Bootstrapper::bootstrap_inplace_3(PhantomCiphertext &cipher) {
  PhantomCiphertext rtncipher;
  bootstrap_3(rtncipher, cipher);
  cipher = rtncipher;
}

void Bootstrapper::bootstrap_real_3(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  initial_scale = cipher.scale();
  if (logn == logNh)
    bootstrap_full_real_3(rtncipher, cipher);
  else
    bootstrap_sparse_real_3(rtncipher, cipher);
}

void Bootstrapper::bootstrap_inplace_real_3(PhantomCiphertext &cipher) {
  PhantomCiphertext rtncipher;
  bootstrap_real_3(rtncipher, cipher);
  cipher = rtncipher;
}

void Bootstrapper::bootstrap_hoisting(PhantomCiphertext &rtncipher, PhantomCiphertext &cipher) {
  if (logn == logNh)
    bootstrap_full_hoisting(rtncipher, cipher);
  else
    bootstrap_sparse_hoisting(rtncipher, cipher);
}

void Bootstrapper::bootstrap_inplace_hoisting(PhantomCiphertext &cipher) {
  PhantomCiphertext rtncipher;
  bootstrap_hoisting(rtncipher, cipher);
  cipher = rtncipher;
}
