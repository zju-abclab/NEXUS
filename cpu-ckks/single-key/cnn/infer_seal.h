#pragma	once
#include "seal/seal.h"
#include "SEALcomp.h"
#include "MinicompFunc.h"
#include "func.h"
#include "PolyUpdate.h"
#include "program.h"
#include "Bootstrapper.h"
#include "cnn_seal.h"
#include <omp.h>
#include <NTL/RR.h>
#include <fstream>
#include <vector>
#include <chrono>

// import parameters
void import_parameters_cifar10(vector<double> &linear_weight, vector<double> &linear_bias, vector<vector<double>> &conv_weight, vector<vector<double>> &bn_bias, vector<vector<double>> &bn_running_mean, vector<vector<double>> &bn_running_var, vector<vector<double>> &bn_weight, size_t layer_num, size_t end_num);
void import_parameters_cifar100(vector<double> &linear_weight, vector<double> &linear_bias, vector<vector<double>> &conv_weight, vector<vector<double>> &bn_bias, vector<vector<double>> &bn_running_mean, vector<vector<double>> &bn_running_var, vector<vector<double>> &bn_weight, vector<vector<double>> &shortcut_weight, vector<vector<double>> &shortcut_bn_bias, vector<vector<double>> &shortcut_bn_mean, vector<vector<double>> &shortcut_bn_var, vector<vector<double>> &shortcut_bn_weight, size_t layer_num, size_t end_num);

// cifar10, cifar100 integrated
void ResNet_cifar10_seal_sparse(size_t layer_num, size_t start_image_id, size_t end_image_id);
// void ResNet_cifar100_seal_sparse(size_t layer_num, size_t start_image_id, size_t end_image_id);