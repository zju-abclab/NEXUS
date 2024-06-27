#include "infer_seal.h"

void import_parameters_cifar10(
    vector<double> &linear_weight,
    vector<double> &linear_bias,
    vector<vector<double>> &conv_weight,
    vector<vector<double>> &bn_bias,
    vector<vector<double>> &bn_running_mean,
    vector<vector<double>> &bn_running_var,
    vector<vector<double>> &bn_weight,
    size_t layer_num,
    size_t end_num)
{
    string dir;
    if (layer_num != 20 && layer_num != 32 && layer_num != 44 && layer_num != 56 && layer_num != 110)
        throw std::invalid_argument("layer number is not valid");
    if (layer_num == 20)
        dir = "resnet20_new";
    else if (layer_num == 32)
        dir = "resnet32_new";
    else if (layer_num == 44)
        dir = "resnet44_new";
    else if (layer_num == 56)
        dir = "resnet56_new";
    else if (layer_num == 110)
        dir = "resnet110_new";

    ifstream in;
    double val;
    size_t num_c = 0, num_b = 0, num_m = 0, num_v = 0, num_w = 0;

    conv_weight.clear();
    conv_weight.resize(layer_num - 1);
    bn_bias.clear();
    bn_bias.resize(layer_num - 1);
    bn_running_mean.clear();
    bn_running_mean.resize(layer_num - 1);
    bn_running_var.clear();
    bn_running_var.resize(layer_num - 1);
    bn_weight.clear();
    bn_weight.resize(layer_num - 1);

    int fh = 3, fw = 3;
    int ci = 0, co = 0;

    // convolution parameters
    ci = 3, co = 16;
    in.open("../../pretrained_parameters/" + dir + "/conv1_weight.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < fh * fw * ci * co; i++) {
        in >> val;
        conv_weight[num_c].emplace_back(val);
    }
    in.close();
    num_c++;

    // convolution parameters
    for (int j = 1; j <= 3; j++) {
        for (int k = 0; k <= end_num; k++) {
            // co setting
            if (j == 1)
                co = 16;
            else if (j == 2)
                co = 32;
            else if (j == 3)
                co = 64;

            // ci setting
            if (j == 1 || (j == 2 && k == 0))
                ci = 16;
            else if ((j == 2 && k != 0) || (j == 3 && k == 0))
                ci = 32;
            else
                ci = 64;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_conv1_weight.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < fh * fw * ci * co; i++) {
                in >> val;
                conv_weight[num_c].emplace_back(val);
            }
            in.close();
            num_c++;

            // ci setting
            if (j == 1)
                ci = 16;
            else if (j == 2)
                ci = 32;
            else if (j == 3)
                ci = 64;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_conv2_weight.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < fh * fw * ci * co; i++) {
                in >> val;
                conv_weight[num_c].emplace_back(val);
            }
            in.close();
            num_c++;
        }
    }

    // batch_normalization parameters
    ci = 16;
    in.open("../../pretrained_parameters/" + dir + "/bn1_bias.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        bn_bias[num_b].emplace_back(val);
    }
    in.close();
    num_b++;
    in.open("../../pretrained_parameters/" + dir + "/bn1_running_mean.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        bn_running_mean[num_m].emplace_back(val);
    }
    in.close();
    num_m++;
    in.open("../../pretrained_parameters/" + dir + "/bn1_running_var.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        bn_running_var[num_v].emplace_back(val);
    }
    in.close();
    num_v++;
    in.open("../../pretrained_parameters/" + dir + "/bn1_weight.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        bn_weight[num_w].emplace_back(val);
    }
    in.close();
    num_w++;

    // batch_normalization parameters
    for (int j = 1; j <= 3; j++) {
        int ci;
        if (j == 1)
            ci = 16;
        else if (j == 2)
            ci = 32;
        else if (j == 3)
            ci = 64;

        for (int k = 0; k <= end_num; k++) {
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn1_bias.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_bias[num_b].emplace_back(val);
            }
            in.close();
            num_b++;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn1_running_mean.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_running_mean[num_m].emplace_back(val);
            }
            in.close();
            num_m++;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn1_running_var.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_running_var[num_v].emplace_back(val);
            }
            in.close();
            num_v++;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn1_weight.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_weight[num_w].emplace_back(val);
            }
            in.close();
            num_w++;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn2_bias.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_bias[num_b].emplace_back(val);
            }
            in.close();
            num_b++;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn2_running_mean.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_running_mean[num_m].emplace_back(val);
            }
            in.close();
            num_m++;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn2_running_var.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_running_var[num_v].emplace_back(val);
            }
            in.close();
            num_v++;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn2_weight.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_weight[num_w].emplace_back(val);
            }
            in.close();
            num_w++;
        }
    }

    // FC
    in.open("../../pretrained_parameters/" + dir + "/linear_weight.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < 10 * 64; i++) {
        in >> val;
        linear_weight.emplace_back(val);
    }
    in.close();
    in.open("../../pretrained_parameters/" + dir + "/linear_bias.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < 10; i++) {
        in >> val;
        linear_bias.emplace_back(val);
    }
    in.close();
}
void import_parameters_cifar100(
    vector<double> &linear_weight,
    vector<double> &linear_bias,
    vector<vector<double>> &conv_weight,
    vector<vector<double>> &bn_bias,
    vector<vector<double>> &bn_running_mean,
    vector<vector<double>> &bn_running_var,
    vector<vector<double>> &bn_weight,
    vector<vector<double>> &shortcut_weight,
    vector<vector<double>> &shortcut_bn_bias,
    vector<vector<double>> &shortcut_bn_mean,
    vector<vector<double>> &shortcut_bn_var,
    vector<vector<double>> &shortcut_bn_weight,
    size_t layer_num,
    size_t end_num)
{
    string dir;
    if (layer_num != 32)
        throw std::invalid_argument("layer number is not valid");
    dir = "resnet32_cifar100";

    ifstream in;
    double val;
    size_t num_c = 0, num_b = 0, num_m = 0, num_v = 0, num_w = 0;

    conv_weight.clear();
    conv_weight.resize(layer_num - 1);
    bn_bias.clear();
    bn_bias.resize(layer_num - 1);
    bn_running_mean.clear();
    bn_running_mean.resize(layer_num - 1);
    bn_running_var.clear();
    bn_running_var.resize(layer_num - 1);
    bn_weight.clear();
    bn_weight.resize(layer_num - 1);

    shortcut_weight.clear();
    shortcut_weight.resize(2);
    shortcut_bn_bias.clear();
    shortcut_bn_bias.resize(2);
    shortcut_bn_mean.clear();
    shortcut_bn_mean.resize(2);
    shortcut_bn_var.clear();
    shortcut_bn_var.resize(2);
    shortcut_bn_weight.clear();
    shortcut_bn_weight.resize(2);

    int fh = 3, fw = 3;
    int ci = 0, co = 0;

    // convolution parameters
    ci = 3, co = 16;
    in.open("../../pretrained_parameters/" + dir + "/conv1_weight.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < fh * fw * ci * co; i++) {
        in >> val;
        conv_weight[num_c].emplace_back(val);
    }
    in.close();
    num_c++;

    // convolution parameters
    for (int j = 1; j <= 3; j++) {
        for (int k = 0; k <= end_num; k++) {
            // co setting
            if (j == 1)
                co = 16;
            else if (j == 2)
                co = 32;
            else if (j == 3)
                co = 64;

            // ci setting
            if (j == 1 || (j == 2 && k == 0))
                ci = 16;
            else if ((j == 2 && k != 0) || (j == 3 && k == 0))
                ci = 32;
            else
                ci = 64;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_conv1_weight.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < fh * fw * ci * co; i++) {
                in >> val;
                conv_weight[num_c].emplace_back(val);
            }
            in.close();
            num_c++;

            // ci setting
            if (j == 1)
                ci = 16;
            else if (j == 2)
                ci = 32;
            else if (j == 3)
                ci = 64;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_conv2_weight.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < fh * fw * ci * co; i++) {
                in >> val;
                conv_weight[num_c].emplace_back(val);
            }
            in.close();
            num_c++;
        }
    }

    // shortcut convolution parameters
    fh = 1, fw = 1;
    ci = 16, co = 32;
    in.open("../../pretrained_parameters/" + dir + "/layer2_0_shortcut_0_weight.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < fh * fw * ci * co; i++) {
        in >> val;
        shortcut_weight[0].emplace_back(val);
    }
    in.close();
    ci = 32, co = 64;
    in.open("../../pretrained_parameters/" + dir + "/layer3_0_shortcut_0_weight.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < fh * fw * ci * co; i++) {
        in >> val;
        shortcut_weight[1].emplace_back(val);
    }
    in.close();

    // batch_normalization parameters
    ci = 16;
    in.open("../../pretrained_parameters/" + dir + "/bn1_bias.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        bn_bias[num_b].emplace_back(val);
    }
    in.close();
    num_b++;
    in.open("../../pretrained_parameters/" + dir + "/bn1_running_mean.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        bn_running_mean[num_m].emplace_back(val);
    }
    in.close();
    num_m++;
    in.open("../../pretrained_parameters/" + dir + "/bn1_running_var.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        bn_running_var[num_v].emplace_back(val);
    }
    in.close();
    num_v++;
    in.open("../../pretrained_parameters/" + dir + "/bn1_weight.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        bn_weight[num_w].emplace_back(val);
    }
    in.close();
    num_w++;

    // batch_normalization parameters
    for (int j = 1; j <= 3; j++) {
        int ci;
        if (j == 1)
            ci = 16;
        else if (j == 2)
            ci = 32;
        else if (j == 3)
            ci = 64;

        for (int k = 0; k <= end_num; k++) {
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn1_bias.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_bias[num_b].emplace_back(val);
            }
            in.close();
            num_b++;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn1_running_mean.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_running_mean[num_m].emplace_back(val);
            }
            in.close();
            num_m++;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn1_running_var.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_running_var[num_v].emplace_back(val);
            }
            in.close();
            num_v++;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn1_weight.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_weight[num_w].emplace_back(val);
            }
            in.close();
            num_w++;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn2_bias.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_bias[num_b].emplace_back(val);
            }
            in.close();
            num_b++;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn2_running_mean.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_running_mean[num_m].emplace_back(val);
            }
            in.close();
            num_m++;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn2_running_var.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_running_var[num_v].emplace_back(val);
            }
            in.close();
            num_v++;
            in.open("../../pretrained_parameters/" + dir + "/layer" + to_string(j) + "_" + to_string(k) + "_bn2_weight.txt");
            if (!in.is_open())
                throw std::runtime_error("file is not open");
            for (long i = 0; i < ci; i++) {
                in >> val;
                bn_weight[num_w].emplace_back(val);
            }
            in.close();
            num_w++;
        }
    }

    // shortcut batch normalization parameters
    ci = 32;
    in.open("../../pretrained_parameters/" + dir + "/layer2_0_shortcut_1_bias.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        shortcut_bn_bias[0].emplace_back(val);
    }
    in.close(); // layer 1
    in.open("../../pretrained_parameters/" + dir + "/layer2_0_shortcut_1_running_mean.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        shortcut_bn_mean[0].emplace_back(val);
    }
    in.close();
    in.open("../../pretrained_parameters/" + dir + "/layer2_0_shortcut_1_running_var.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        shortcut_bn_var[0].emplace_back(val);
    }
    in.close();
    in.open("../../pretrained_parameters/" + dir + "/layer2_0_shortcut_1_weight.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        shortcut_bn_weight[0].emplace_back(val);
    }
    in.close();

    ci = 64;
    in.open("../../pretrained_parameters/" + dir + "/layer3_0_shortcut_1_bias.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        shortcut_bn_bias[1].emplace_back(val);
    }
    in.close(); // layer 1
    in.open("../../pretrained_parameters/" + dir + "/layer3_0_shortcut_1_running_mean.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        shortcut_bn_mean[1].emplace_back(val);
    }
    in.close();
    in.open("../../pretrained_parameters/" + dir + "/layer3_0_shortcut_1_running_var.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        shortcut_bn_var[1].emplace_back(val);
    }
    in.close();
    in.open("../../pretrained_parameters/" + dir + "/layer3_0_shortcut_1_weight.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < ci; i++) {
        in >> val;
        shortcut_bn_weight[1].emplace_back(val);
    }
    in.close();

    // FC
    in.open("../../pretrained_parameters/" + dir + "/linear_weight.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < 100 * 64; i++) {
        in >> val;
        linear_weight.emplace_back(val);
    }
    in.close();
    in.open("../../pretrained_parameters/" + dir + "/linear_bias.txt");
    if (!in.is_open())
        throw std::runtime_error("file is not open");
    for (long i = 0; i < 100; i++) {
        in >> val;
        linear_bias.emplace_back(val);
    }
    in.close();
}

void ResNet_cifar10_seal_sparse(size_t layer_num, size_t start_image_id, size_t end_image_id)
{
    // approximation boundary setting
    double B = 40.0; // approximation boundary

    // approx ReLU setting
    long alpha = 13;                  // precision parameter alpha
    long comp_no = 3;                 // number of compositions
    vector<int> deg = { 15, 15, 27 }; // degrees of component polynomials
    // double eta = pow(2.0,-15);		// margin
    double scaled_val = 1.7; // scaled_val: the last scaled value
    // double max_factor = 16;		// max_factor = 1 for comparison operation. max_factor > 1 for max or ReLU function
    vector<Tree> tree; // structure of polynomial evaluation
    evaltype ev_type = evaltype::oddbaby;
    // RR::SetOutputPrecision(25);

    // generate tree
    for (int i = 0; i < comp_no; i++) {
        Tree tr;
        if (ev_type == evaltype::oddbaby)
            upgrade_oddbaby(deg[i], tr);
        else if (ev_type == evaltype::baby)
            upgrade_baby(deg[i], tr);
        else
            std::invalid_argument("evaluation type is not correct");
        tree.emplace_back(tr);
        // tr.print();
    }

    // all threads output files
    ofstream out_share;
    if (layer_num == 20)
        out_share.open("../../result/resnet20_cifar10_label_" + to_string(start_image_id) + "_" + to_string(end_image_id));
    else if (layer_num == 32)
        out_share.open("../../result/resnet32_cifar10_label_" + to_string(start_image_id) + "_" + to_string(end_image_id));
    else if (layer_num == 44)
        out_share.open("../../result/resnet44_cifar10_label_" + to_string(start_image_id) + "_" + to_string(end_image_id));
    else if (layer_num == 56)
        out_share.open("../../result/resnet56_cifar10_label_" + to_string(start_image_id) + "_" + to_string(end_image_id));
    else if (layer_num == 110)
        out_share.open("../../result/resnet110_cifar10_label_" + to_string(start_image_id) + "_" + to_string(end_image_id));
    else
        throw std::invalid_argument("layer_num is not correct");

    // SEAL and bootstrapping setting
    long boundary_K = 25;
    long boot_deg = 59;
    long scale_factor = 2;
    long inverse_deg = 1;
    long logN = 16;
    long loge = 10;
    long logn = 15;   // full slots
    long logn_1 = 14; // sparse slots
    long logn_2 = 13;
    long logn_3 = 12;
    int logp = 46;
    int logq = 51;
    int log_special_prime = 51;
    int log_integer_part = logq - logp - loge + 5;
    int remaining_level = 16; // Calculation required
    int boot_level = 14;      //
    int total_level = remaining_level + boot_level;

    vector<int> coeff_bit_vec;
    coeff_bit_vec.push_back(logq);
    for (int i = 0; i < remaining_level; i++)
        coeff_bit_vec.push_back(logp);
    for (int i = 0; i < boot_level; i++)
        coeff_bit_vec.push_back(logq);
    coeff_bit_vec.push_back(log_special_prime);

    cout << "Setting Parameters" << endl;
    EncryptionParameters parms(scheme_type::ckks);
    size_t poly_modulus_degree = (size_t)(1 << logN);
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_bit_vec));

    // added
    size_t secret_key_hamming_weight = 192;
    parms.set_secret_key_hamming_weight(secret_key_hamming_weight);
    // parms.set_sparse_slots(1 << logn_1);
    double scale = pow(2.0, logp);

    SEALContext context(parms);
    // KeyGenerator keygen(context, 192);
    KeyGenerator keygen(context);
    PublicKey public_key;
    keygen.create_public_key(public_key);
    auto secret_key = keygen.secret_key();
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);
    GaloisKeys gal_keys;

    CKKSEncoder encoder(context);
    Encryptor encryptor(context, public_key);
    Evaluator evaluator(context, encoder);
    Decryptor decryptor(context, secret_key);
    // ScaleInvEvaluator scale_evaluator(context, encoder, relin_keys);

    Bootstrapper bootstrapper_1(
        loge,
        logn_1,
        logN - 1,
        total_level,
        scale,
        boundary_K,
        boot_deg,
        scale_factor,
        inverse_deg,
        context,
        keygen,
        encoder,
        encryptor,
        decryptor,
        evaluator,
        relin_keys,
        gal_keys);
    Bootstrapper bootstrapper_2(
        loge,
        logn_2,
        logN - 1,
        total_level,
        scale,
        boundary_K,
        boot_deg,
        scale_factor,
        inverse_deg,
        context,
        keygen,
        encoder,
        encryptor,
        decryptor,
        evaluator,
        relin_keys,
        gal_keys);
    Bootstrapper bootstrapper_3(
        loge,
        logn_3,
        logN - 1,
        total_level,
        scale,
        boundary_K,
        boot_deg,
        scale_factor,
        inverse_deg,
        context,
        keygen,
        encoder,
        encryptor,
        decryptor,
        evaluator,
        relin_keys,
        gal_keys);

    //	additional rotation kinds for CNN
    vector<int> rotation_kinds = { 0,
                                   1,
                                   2,
                                   3,
                                   4,
                                   5,
                                   6,
                                   7,
                                   8,
                                   9,
                                   10,
                                   11,
                                   12,
                                   13,
                                   14,
                                   15,
                                   16,
                                   17,
                                   18,
                                   19,
                                   20,
                                   21,
                                   22,
                                   23,
                                   24,
                                   25,
                                   26,
                                   27,
                                   28,
                                   29,
                                   30,
                                   31,
                                   32,
                                   33
                                   // ,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55
                                   ,
                                   56
                                   // ,57,58,59,60,61
                                   ,
                                   62,
                                   63,
                                   64,
                                   66,
                                   84,
                                   124,
                                   128,
                                   132,
                                   256,
                                   512,
                                   959,
                                   960,
                                   990,
                                   991,
                                   1008,
                                   1023,
                                   1024,
                                   1036,
                                   1064,
                                   1092,
                                   1952,
                                   1982,
                                   1983,
                                   2016,
                                   2044,
                                   2047,
                                   2048,
                                   2072,
                                   2078,
                                   2100,
                                   3007,
                                   3024,
                                   3040,
                                   3052,
                                   3070,
                                   3071,
                                   3072,
                                   3080,
                                   3108,
                                   4031,
                                   4032,
                                   4062,
                                   4063,
                                   4095,
                                   4096,
                                   5023,
                                   5024,
                                   5054,
                                   5055,
                                   5087,
                                   5118,
                                   5119,
                                   5120,
                                   6047,
                                   6078,
                                   6079,
                                   6111,
                                   6112,
                                   6142,
                                   6143,
                                   6144,
                                   7071,
                                   7102,
                                   7103,
                                   7135,
                                   7166,
                                   7167,
                                   7168,
                                   8095,
                                   8126,
                                   8127,
                                   8159,
                                   8190,
                                   8191,
                                   8192,
                                   9149,
                                   9183,
                                   9184,
                                   9213,
                                   9215,
                                   9216,
                                   10173,
                                   10207,
                                   10208,
                                   10237,
                                   10239,
                                   10240,
                                   11197,
                                   11231,
                                   11232,
                                   11261,
                                   11263,
                                   11264,
                                   12221,
                                   12255,
                                   12256,
                                   12285,
                                   12287,
                                   12288,
                                   13214,
                                   13216,
                                   13246,
                                   13278,
                                   13279,
                                   13280,
                                   13310,
                                   13311,
                                   13312,
                                   14238,
                                   14240,
                                   14270,
                                   14302,
                                   14303,
                                   14304,
                                   14334,
                                   14335,
                                   15262,
                                   15264,
                                   15294,
                                   15326,
                                   15327,
                                   15328,
                                   15358,
                                   15359,
                                   15360,
                                   16286,
                                   16288,
                                   16318,
                                   16350,
                                   16351,
                                   16352,
                                   16382,
                                   16383,
                                   16384,
                                   17311,
                                   17375,
                                   18335,
                                   18399,
                                   18432,
                                   19359,
                                   19423,
                                   20383,
                                   20447,
                                   20480,
                                   21405,
                                   21406,
                                   21437,
                                   21469,
                                   21470,
                                   21471,
                                   21501,
                                   21504,
                                   22429,
                                   22430,
                                   22461,
                                   22493,
                                   22494,
                                   22495,
                                   22525,
                                   22528,
                                   23453,
                                   23454,
                                   23485,
                                   23517,
                                   23518,
                                   23519,
                                   23549,
                                   24477,
                                   24478,
                                   24509,
                                   24541,
                                   24542,
                                   24543,
                                   24573,
                                   24576,
                                   25501,
                                   25565,
                                   25568,
                                   25600,
                                   26525,
                                   26589,
                                   26592,
                                   26624,
                                   27549,
                                   27613,
                                   27616,
                                   27648,
                                   28573,
                                   28637,
                                   28640,
                                   28672,
                                   29600,
                                   29632,
                                   29664,
                                   29696,
                                   30624,
                                   30656,
                                   30688,
                                   30720,
                                   31648,
                                   31680,
                                   31712,
                                   31743,
                                   31744,
                                   31774,
                                   32636,
                                   32640,
                                   32644,
                                   32672,
                                   32702,
                                   32704,
                                   32706,
                                   32735,
                                   32736,
                                   32737,
                                   32759,
                                   32760,
                                   32761,
                                   32762,
                                   32763,
                                   32764,
                                   32765,
                                   32766,
                                   32767 };

    // bootstrapping preprocessing
    cout << "Generating Optimal Minimax Polynomials..." << endl;
    bootstrapper_1.prepare_mod_polynomial();
    bootstrapper_2.prepare_mod_polynomial();
    bootstrapper_3.prepare_mod_polynomial();

    cout << "Adding Bootstrapping Keys..." << endl;
    vector<int> gal_steps_vector;
    gal_steps_vector.push_back(0);
    for (int i = 0; i < logN - 1; i++)
        gal_steps_vector.push_back((1 << i));
    for (auto rot : rotation_kinds) {
        if (find(gal_steps_vector.begin(), gal_steps_vector.end(), rot) == gal_steps_vector.end())
            gal_steps_vector.push_back(rot);
    }
    bootstrapper_1.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
    bootstrapper_2.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
    bootstrapper_3.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
    keygen.create_galois_keys(gal_steps_vector, gal_keys);

    bootstrapper_1.slot_vec.push_back(logn_1);
    bootstrapper_2.slot_vec.push_back(logn_2);
    bootstrapper_3.slot_vec.push_back(logn_3);

    cout << "Generating Linear Transformation Coefficients..." << endl;
    bootstrapper_1.generate_LT_coefficient_3();
    bootstrapper_2.generate_LT_coefficient_3();
    bootstrapper_3.generate_LT_coefficient_3();

    // time setting
    chrono::high_resolution_clock::time_point all_time_start, all_time_end;
    chrono::microseconds all_time_diff;
    all_time_start = chrono::high_resolution_clock::now();

    // end number
    int end_num = 0;
    if (layer_num == 20)
        end_num = 2; // 0 ~ 2
    else if (layer_num == 32)
        end_num = 4; // 0 ~ 4
    else if (layer_num == 44)
        end_num = 6; // 0 ~ 6
    else if (layer_num == 56)
        end_num = 8; // 0 ~ 8
    else if (layer_num == 110)
        end_num = 17; // 0 ~ 17
    else
        throw std::invalid_argument("layer_num is not correct");

#pragma omp parallel for num_threads(50)
    for (size_t image_id = start_image_id; image_id <= end_image_id; image_id++) {
        // each thread output result file
        ofstream output;
        if (layer_num == 20)
            output.open("../../result/resnet20_cifar10_image" + to_string(image_id) + ".txt");
        else if (layer_num == 32)
            output.open("../../result/resnet32_cifar10_image" + to_string(image_id) + ".txt");
        else if (layer_num == 44)
            output.open("../../result/resnet44_cifar10_image" + to_string(image_id) + ".txt");
        else if (layer_num == 56)
            output.open("../../result/resnet56_cifar10_image" + to_string(image_id) + ".txt");
        else if (layer_num == 110)
            output.open("../../result/resnet110_cifar10_image" + to_string(image_id) + ".txt");
        else
            throw std::invalid_argument("layer_num is not correct");
        string dir = "resnet" + to_string(layer_num) + "_new";

        // ciphertext pool generation
        vector<Ciphertext> cipher_pool(14);

        // time setting
        chrono::high_resolution_clock::time_point time_start, time_end, total_time_start, total_time_end;
        chrono::microseconds time_diff, total_time_diff;

        // variables
        TensorCipher cnn, temp;

        // deep learning parameters and import
        int co = 0, st = 0, fh = 3, fw = 3;
        long init_p = 8, n = 1 << logn;
        int stage = 0;
        double epsilon = 0.00001;
        vector<double> image, linear_weight, linear_bias;
        vector<vector<double>> conv_weight, bn_bias, bn_running_mean, bn_running_var, bn_weight;
        import_parameters_cifar10(linear_weight, linear_bias, conv_weight, bn_bias, bn_running_mean, bn_running_var, bn_weight, layer_num, end_num);

        // pack images compactly
        ifstream in;
        double val;
        in.open("../../../testFile/test_values.txt");
        for (long i = 0; i < 1 << logn; i++)
            image.emplace_back(0);
        for (long i = 0; i < 32 * 32 * 3 * image_id; i++) {
            in >> val;
        }
        for (long i = 0; i < 32 * 32 * 3; i++) {
            in >> val;
            image[i] = val;
        }
        in.close();
        for (long i = n / init_p; i < n; i++)
            image[i] = image[i % (n / init_p)];
        for (long i = 0; i < n; i++)
            image[i] /= B; // for boundary [-1,1]

        ifstream in_label;
        int image_label;
        in_label.open("../../../testFile/test_label.txt");
        for (long i = 0; i < image_id; i++) {
            in_label >> image_label;
        }
        in_label >> image_label;

        // generate CIFAR-10 image
        cnn = TensorCipher(logn, 1, 32, 32, 3, 3, init_p, image, encryptor, encoder, logq);
        // decrypt_and_print(cnn.cipher(), decryptor, encoder, 1<<logn, 256, 2); cnn.print_parms();
        cout << "remaining level : " << context.get_context_data(cnn.cipher().parms_id())->chain_index() << endl;
        cout << "scale: " << cnn.cipher().scale() << endl;
        total_time_start = chrono::high_resolution_clock::now();

        // modulus down
        Ciphertext ctxt;
        ctxt = cnn.cipher();
        for (int i = 0; i < boot_level - 3; i++)
            evaluator.mod_switch_to_next_inplace(ctxt);
        cnn.set_ciphertext(ctxt);

        // layer 0
        cout << "layer 0" << endl;
        output << "layer 0" << endl;
        multiplexed_parallel_convolution_print(
            cnn,
            cnn,
            16,
            1,
            fh,
            fw,
            conv_weight[stage],
            bn_running_var[stage],
            bn_weight[stage],
            epsilon,
            encoder,
            encryptor,
            evaluator,
            gal_keys,
            cipher_pool,
            output,
            decryptor,
            context,
            stage);

        // scaling factor ~2^51 -> 2^46
        const auto &modulus = iter(context.first_context_data()->parms().coeff_modulus());
        ctxt = cnn.cipher();
        size_t cur_level = ctxt.coeff_modulus_size();
        Plaintext scaler;
        double scale_change = pow(2.0, 46) * ((double)modulus[cur_level - 1].value()) / ctxt.scale();
        encoder.encode(1, scale_change, scaler);
        evaluator.mod_switch_to_inplace(scaler, ctxt.parms_id());
        evaluator.multiply_plain_inplace(ctxt, scaler);
        evaluator.rescale_to_next_inplace(ctxt);
        ctxt.scale() = pow(2.0, 46);
        cnn.set_ciphertext(ctxt);

        multiplexed_parallel_batch_norm_seal_print(
            cnn,
            cnn,
            bn_bias[stage],
            bn_running_mean[stage],
            bn_running_var[stage],
            bn_weight[stage],
            epsilon,
            encoder,
            encryptor,
            evaluator,
            B,
            output,
            decryptor,
            context,
            stage);
        approx_ReLU_seal_print(
            cnn,
            cnn,
            comp_no,
            deg,
            alpha,
            tree,
            scaled_val,
            logp,
            encryptor,
            evaluator,
            decryptor,
            encoder,
            public_key,
            secret_key,
            relin_keys,
            B,
            output,
            context,
            gal_keys,
            stage);

        for (int j = 0; j < 3; j++) // layer 1_x, 2_x, 3_x
        {
            if (j == 0)
                co = 16;
            else if (j == 1)
                co = 32;
            else if (j == 2)
                co = 64;

            // // sparse slot
            // if(j==0) {
            // 	parms.set_sparse_slots(1<<logn_1);
            // 	encoder.set_sparse_slots(1<<logn_1);
            // } else if(j==1) {
            // 	parms.set_sparse_slots(1<<logn_2);
            // 	encoder.set_sparse_slots(1<<logn_2);
            // } else if(j==2) {
            // 	parms.set_sparse_slots(1<<logn_3);
            // 	encoder.set_sparse_slots(1<<logn_3);
            // }

            for (int k = 0; k <= end_num; k++) // 0 ~ 2/4/6/8/17
            {
                stage = 2 * ((end_num + 1) * j + k) + 1;
                cout << "layer " << stage << endl;
                output << "layer " << stage << endl;
                temp = cnn;
                if (j >= 1 && k == 0)
                    st = 2;
                else
                    st = 1;
                multiplexed_parallel_convolution_print(
                    cnn,
                    cnn,
                    co,
                    st,
                    fh,
                    fw,
                    conv_weight[stage],
                    bn_running_var[stage],
                    bn_weight[stage],
                    epsilon,
                    encoder,
                    encryptor,
                    evaluator,
                    gal_keys,
                    cipher_pool,
                    output,
                    decryptor,
                    context,
                    stage);
                multiplexed_parallel_batch_norm_seal_print(
                    cnn,
                    cnn,
                    bn_bias[stage],
                    bn_running_mean[stage],
                    bn_running_var[stage],
                    bn_weight[stage],
                    epsilon,
                    encoder,
                    encryptor,
                    evaluator,
                    B,
                    output,
                    decryptor,
                    context,
                    stage);
                if (j == 0)
                    bootstrap_print(cnn, cnn, bootstrapper_1, output, decryptor, encoder, context, stage);
                else if (j == 1)
                    bootstrap_print(cnn, cnn, bootstrapper_2, output, decryptor, encoder, context, stage);
                else if (j == 2)
                    bootstrap_print(cnn, cnn, bootstrapper_3, output, decryptor, encoder, context, stage);
                approx_ReLU_seal_print(
                    cnn,
                    cnn,
                    comp_no,
                    deg,
                    alpha,
                    tree,
                    scaled_val,
                    logp,
                    encryptor,
                    evaluator,
                    decryptor,
                    encoder,
                    public_key,
                    secret_key,
                    relin_keys,
                    B,
                    output,
                    context,
                    gal_keys,
                    stage);

                stage = 2 * ((end_num + 1) * j + k) + 2;
                cout << "layer " << stage << endl;
                output << "layer " << stage << endl;
                st = 1;

                multiplexed_parallel_convolution_print(
                    cnn,
                    cnn,
                    co,
                    st,
                    fh,
                    fw,
                    conv_weight[stage],
                    bn_running_var[stage],
                    bn_weight[stage],
                    epsilon,
                    encoder,
                    encryptor,
                    evaluator,
                    gal_keys,
                    cipher_pool,
                    output,
                    decryptor,
                    context,
                    stage);
                multiplexed_parallel_batch_norm_seal_print(
                    cnn,
                    cnn,
                    bn_bias[stage],
                    bn_running_mean[stage],
                    bn_running_var[stage],
                    bn_weight[stage],
                    epsilon,
                    encoder,
                    encryptor,
                    evaluator,
                    B,
                    output,
                    decryptor,
                    context,
                    stage);
                if (j >= 1 && k == 0)
                    multiplexed_parallel_downsampling_seal_print(temp, temp, evaluator, decryptor, encoder, context, gal_keys, output);
                cipher_add_seal_print(temp, cnn, cnn, evaluator, output, decryptor, encoder, context);
                if (j == 0)
                    bootstrap_print(cnn, cnn, bootstrapper_1, output, decryptor, encoder, context, stage);
                else if (j == 1)
                    bootstrap_print(cnn, cnn, bootstrapper_2, output, decryptor, encoder, context, stage);
                else if (j == 2)
                    bootstrap_print(cnn, cnn, bootstrapper_3, output, decryptor, encoder, context, stage);
                approx_ReLU_seal_print(
                    cnn,
                    cnn,
                    comp_no,
                    deg,
                    alpha,
                    tree,
                    scaled_val,
                    logp,
                    encryptor,
                    evaluator,
                    decryptor,
                    encoder,
                    public_key,
                    secret_key,
                    relin_keys,
                    B,
                    output,
                    context,
                    gal_keys,
                    stage);
            }
        }
        cout << "layer " << layer_num - 1 << endl;
        output << "layer " << layer_num - 1 << endl;
        averagepooling_seal_scale_print(cnn, cnn, evaluator, gal_keys, B, output, decryptor, encoder, context);
        fully_connected_seal_print(cnn, cnn, linear_weight, linear_bias, 10, 64, evaluator, gal_keys, output, decryptor, encoder, context);

        total_time_end = chrono::high_resolution_clock::now();
        total_time_diff = chrono::duration_cast<chrono::milliseconds>(total_time_end - total_time_start);

        // final text file print
        Plaintext plain;
        decryptor.decrypt(cnn.cipher(), plain);
        vector<complex<double>> rtn_vec;
        // encoder.decode(plain, rtn_vec, 1<<logn);
        encoder.decode(plain, rtn_vec);
        cout << "( ";
        output << "( ";
        for (size_t i = 0; i < 9; i++) {
            cout << rtn_vec[i] << ", ";
            output << rtn_vec[i] << ", ";
        }
        cout << rtn_vec[9] << ")" << endl;
        output << rtn_vec[9] << ")" << endl;
        cout << "total time : " << total_time_diff.count() / 1000 << " ms" << endl;
        output << "total time : " << total_time_diff.count() / 1000 << " ms" << endl;

        size_t label = 0;
        double max_score = -100.0;
        for (size_t i = 0; i < 10; i++) {
            if (max_score < rtn_vec[i].real()) {
                label = i;
                max_score = rtn_vec[i].real();
            }
        }
        cout << "image label: " << image_label << endl;
        cout << "inferred label: " << label << endl;
        cout << "max score: " << max_score << endl;
        output << "image label: " << image_label << endl;
        output << "inferred label: " << label << endl;
        output << "max score: " << max_score << endl;
        out_share << "image_id: " << image_id << ", "
                  << "image label: " << image_label << ", inferred label: " << label << endl;
    }

    all_time_end = chrono::high_resolution_clock::now();
    all_time_diff = chrono::duration_cast<chrono::milliseconds>(all_time_end - all_time_start);
    cout << "all threads time : " << all_time_diff.count() / 1000 << " ms" << endl;
    out_share << endl << "all threads time : " << all_time_diff.count() / 1000 << " ms" << endl;
}
// void ResNet_cifar100_seal_sparse(size_t layer_num, size_t start_image_id, size_t end_image_id)
// {
// 	// approximation boundary setting
// 	double B = 65.0;	// approximation boundary

// 	// approx ReLU setting
// 	long alpha = 13;			// precision parameter alpha
// 	long comp_no = 3;		// number of compositions
// 	vector<int> deg = {15,15,27};		// degrees of component polynomials
// 	// double eta = pow(2.0,-15);		// margin
// 	double scaled_val = 1.7;		// scaled_val: the last scaled value
// 	// double max_factor = 16;		// max_factor = 1 for comparison operation. max_factor > 1 for max or ReLU function
// 	vector<Tree> tree;		// structure of polynomial evaluation
// 	evaltype ev_type = evaltype::oddbaby;
// 	// RR::SetOutputPrecision(25);

// 	// generate tree
// 	for(int i=0; i<comp_no; i++)
// 	{
// 		Tree tr;
// 		if(ev_type == evaltype::oddbaby) upgrade_oddbaby(deg[i], tr);
// 		else if(ev_type == evaltype::baby) upgrade_baby(deg[i], tr);
// 		else std::invalid_argument("evaluation type is not correct");
// 		tree.emplace_back(tr);
// 		// tr.print();
// 	}

// 	// all threads output files
// 	ofstream out_share;
// 	if(layer_num == 32) out_share.open("../../result/resnet32_cifar100_label_" + to_string(start_image_id) + "_" + to_string(end_image_id));
// 	else throw std::invalid_argument("layer number is not correct");

// 	// SEAL and bootstrapping setting
// 	long boundary_K = 25;
// 	long boot_deg = 59;
//     long scale_factor = 2;
//     long inverse_deg = 1;
// 	long logN = 16;
// 	long loge = 10;
// 	long logn = 15;		// full slots
// 	long logn_1 = 14;	// sparse slots
// 	long logn_2 = 13;
// 	long logn_3 = 12;
// 	int logp = 46;
// 	int logq = 51;
// 	int log_special_prime = 51;
//     int log_integer_part = logq - logp - loge + 5;
// 	int remaining_level = 16; // Calculation required
// 	int boot_level = 14; //
// 	int total_level = remaining_level + boot_level;

// 	vector<int> coeff_bit_vec;
// 	coeff_bit_vec.push_back(logq);
// 	for (int i = 0; i < remaining_level; i++) coeff_bit_vec.push_back(logp);
// 	for (int i = 0; i < boot_level; i++) coeff_bit_vec.push_back(logq);
// 	coeff_bit_vec.push_back(log_special_prime);

// 	cout << "Setting Parameters" << endl;
// 	EncryptionParameters parms(scheme_type::ckks);
// 	size_t poly_modulus_degree = (size_t)(1 << logN);
// 	parms.set_poly_modulus_degree(poly_modulus_degree);
// 	parms.set_coeff_modulus(CoeffModulus::Create(poly_modulus_degree, coeff_bit_vec));
// 	double scale = pow(2.0, logp);

// 	SEALContext context(parms);
// 	KeyGenerator keygen(context, 192);
//     PublicKey public_key;
// 	keygen.create_public_key(public_key);
// 	auto secret_key = keygen.secret_key();
//     RelinKeys relin_keys;
// 	keygen.create_relin_keys(relin_keys);
// 	GaloisKeys gal_keys;

// 	Encryptor encryptor(context, public_key);
// 	Evaluator evaluator(context);
// 	Decryptor decryptor(context, secret_key);
// 	CKKSEncoder encoder(context);
// 	ScaleInvEvaluator scale_evaluator(context, encoder, relin_keys);

// 	Bootstrapper bootstrapper_1(loge, logn_1, logN - 1, total_level, scale, boundary_K, boot_deg, scale_factor, inverse_deg, context, keygen, encoder,
// encryptor, decryptor, scale_evaluator, relin_keys, gal_keys); 	Bootstrapper bootstrapper_2(loge, logn_2, logN - 1, total_level, scale, boundary_K,
// boot_deg, scale_factor, inverse_deg, context, keygen, encoder, encryptor, decryptor, scale_evaluator, relin_keys, gal_keys); 	Bootstrapper
// bootstrapper_3(loge, logn_3, logN - 1, total_level, scale, boundary_K, boot_deg, scale_factor, inverse_deg, context, keygen, encoder, encryptor, decryptor,
// scale_evaluator, relin_keys, gal_keys);

// //	additional rotation kinds for CNN
// 	vector<int> rotation_kinds = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33
// 		// ,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55
// 		,56
// 		// ,57,58,59,60,61
// 		,62,63,64,66,84,124,128,132,256,512,959,960,990,991,1008
// 		,1023,1024,1036,1064,1092,1952,1982,1983,2016,2044,2047,2048,2072,2078,2100,3007,3024,3040,3052,3070,3071,3072,3080,3108,4031
// 		,4032,4062,4063,4095,4096,5023,5024,5054,5055,5087,5118,5119,5120,6047,6078,6079,6111,6112,6142,6143,6144,7071,7102,7103,7135
// 		,7166,7167,7168,8095,8126,8127,8159,8190,8191,8192,9149,9183,9184,9213,9215,9216,10173,10207,10208,10237,10239,10240,11197,11231
// 		,11232,11261,11263,11264,12221,12255,12256,12285,12287,12288,13214,13216,13246,13278,13279,13280,13310,13311,13312,14238,14240
// 		,14270,14302,14303,14304,14334,14335,15262,15264,15294,15326,15327,15328,15358,15359,15360,16286,16288,16318,16350,16351,16352
// 		,16382,16383,16384,17311,17375,18335,18399,18432,19359,19423,20383,20447,20480,21405,21406,21437,21469,21470,21471,21501,21504
// 		,22429,22430,22461,22493,22494,22495,22525,22528,23453,23454,23485,23517,23518,23519,23549,24477,24478,24509,24541,24542,24543
// 		,24573,24576,25501,25565,25568,25600,26525,26589,26592,26624,27549,27613,27616,27648,28573,28637,28640,28672,29600,29632,29664
// 		,29696,30624,30656,30688,30720,31648,31680,31712,31743,31744,31774,32636,32640,32644,32672,32702,32704,32706,32735
// 		,32736,32737,32759,32760,32761,32762,32763,32764,32765,32766,32767
// 	};

// 	// bootstrapping preprocessing
// 	cout << "Generating Optimal Minimax Polynomials..." << endl;
// 	bootstrapper_1.prepare_mod_polynomial();
// 	bootstrapper_2.prepare_mod_polynomial();
// 	bootstrapper_3.prepare_mod_polynomial();

// 	cout << "Adding Bootstrapping Keys..." << endl;
// 	vector<int> gal_steps_vector;
// 	gal_steps_vector.push_back(0);
// 	for(int i=0; i<logN-1; i++) gal_steps_vector.push_back((1 << i));
// 	for(auto rot: rotation_kinds)
// 	{
// 		if(find(gal_steps_vector.begin(), gal_steps_vector.end(), rot) == gal_steps_vector.end()) gal_steps_vector.push_back(rot);
// 	}
// 	bootstrapper_1.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
// 	bootstrapper_2.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
// 	bootstrapper_3.addLeftRotKeys_Linear_to_vector_3(gal_steps_vector);
// 	keygen.create_galois_keys(gal_steps_vector, gal_keys);

// 	bootstrapper_1.slot_vec.push_back(logn_1);
// 	bootstrapper_2.slot_vec.push_back(logn_2);
// 	bootstrapper_3.slot_vec.push_back(logn_3);

// 	cout << "Generating Linear Transformation Coefficients..." << endl;
// 	bootstrapper_1.generate_LT_coefficient_3();
// 	bootstrapper_2.generate_LT_coefficient_3();
// 	bootstrapper_3.generate_LT_coefficient_3();

// 	// time setting
// 	chrono::high_resolution_clock::time_point all_time_start, all_time_end;
// 	chrono::microseconds all_time_diff;
// 	all_time_start = chrono::high_resolution_clock::now();

// 	// end number
// 	int end_num = 4;

// 	#pragma omp parallel for num_threads(50)
// 	for(size_t image_id = start_image_id; image_id <=end_image_id; image_id++)
// 	{
// 		// each thread output result file
// 		ofstream output;
// 		output.open("../../result/resnet32_cifar100_image" + to_string(image_id) + ".txt");
// 		string dir = "resnet32_cifar100";		// parameter directory name

// 		// ciphertext pool generation
// 		vector<Ciphertext> cipher_pool(14);

// 		// time setting
// 		chrono::high_resolution_clock::time_point time_start, time_end, total_time_start, total_time_end;
// 		chrono::microseconds time_diff, total_time_diff;

// 		// variables
// 		TensorCipher cnn, temp;

// 		// deep learning parameters and import
// 		int co = 0, st = 0, fh = 3, fw = 3;
// 		long init_p = 8, n = 1<<logn;
// 		int stage = 0;
// 		double epsilon = 0.00001;
// 		vector<double> image, linear_weight, linear_bias;
// 		vector<vector<double>> conv_weight, bn_bias, bn_running_mean, bn_running_var, bn_weight, shortcut_weight, shortcut_bn_bias, shortcut_bn_mean,
// shortcut_bn_var, shortcut_bn_weight; 		import_parameters_cifar100(linear_weight, linear_bias, conv_weight, bn_bias, bn_running_mean, bn_running_var,
// bn_weight, shortcut_weight, shortcut_bn_bias, shortcut_bn_mean, shortcut_bn_var, shortcut_bn_weight, layer_num, end_num);

// 		// pack images compactly
// 		ifstream in;
// 		double val;
// 		in.open("../../../testFile_cifar100/test_values.txt");
// 		for(long i=0; i<1<<logn; i++) image.emplace_back(0);
// 		for(long i=0; i<32*32*3 *image_id; i++) {in>>val;}
// 		for(long i=0; i<32*32*3; i++) {in>>val; image[i]=val;}  in.close();
// 		for(long i=n/init_p; i<n; i++) image[i] = image[i%(n/init_p)];
// 		for(long i=0; i<n; i++) image[i] /= B;		// for boundary [-1,1]

// 		ifstream in_label;
// 		int image_label;
// 		in_label.open("../../../testFile_cifar100/test_label.txt");
// 		for(long i=0; i<image_id; i++) {in_label>>image_label;}
// 		in_label >> image_label;

// 		// generate CIFAR-100 image
// 		cnn = TensorCipher(logn, 1, 32, 32, 3, 3, init_p, image, encryptor, encoder, logq);
// 		// decrypt_and_print(cnn.cipher(), decryptor, encoder, 1<<logn, 256, 2); cnn.print_parms();
// 		cout << "remaining level : " << context.get_context_data(cnn.cipher().parms_id())->chain_index() << endl;
// 		cout << "scale: " << cnn.cipher().scale() << endl;
// 		total_time_start = chrono::high_resolution_clock::now();

// 		// modulus down
// 		Ciphertext ctxt;
// 		ctxt = cnn.cipher();
// 		for(int i=0; i<boot_level-3; i++) evaluator.mod_switch_to_next_inplace(ctxt);
// 		cnn.set_ciphertext(ctxt);

// 		// layer 0
// 		cout << "layer 0" << endl;
// 		output << "layer 0" << endl;
// 		multiplexed_parallel_convolution_print(cnn, cnn, 16, 1, fh, fw, conv_weight[stage], bn_running_var[stage], bn_weight[stage], epsilon, encoder,
// encryptor, scale_evaluator, gal_keys, cipher_pool, output, decryptor, context, stage);

// 		// scaling factor ~2^51 -> 2^46
// 		const auto &modulus = iter(context.first_context_data()->parms().coeff_modulus());
// 		ctxt = cnn.cipher();
// 		size_t cur_level = ctxt.coeff_modulus_size();
// 		Plaintext scaler;
// 		double scale_change = pow(2.0,46) * ((double)modulus[cur_level-1].value()) / ctxt.scale();
// 		encoder.encode(1, scale_change, scaler);
// 		evaluator.mod_switch_to_inplace(scaler, ctxt.parms_id());
// 		evaluator.multiply_plain_inplace(ctxt, scaler);
// 		evaluator.rescale_to_next_inplace(ctxt);
// 		ctxt.scale() = pow(2.0,46);
// 		cnn.set_ciphertext(ctxt);

// 		multiplexed_parallel_batch_norm_seal_print(cnn, cnn, bn_bias[stage], bn_running_mean[stage], bn_running_var[stage], bn_weight[stage], epsilon, encoder,
// encryptor, scale_evaluator, B, output, decryptor, context, stage); 		approx_ReLU_seal_print(cnn, cnn, comp_no, deg, alpha, tree, scaled_val, logp,
// encryptor, evaluator, scale_evaluator, decryptor, encoder, public_key, secret_key, relin_keys, B, output, context, gal_keys, stage);

// 		for(int j=0; j<3; j++)		// layer 1_x, 2_x, 3_x
// 		{
// 			if(j==0) co = 16;
// 			else if(j==1) co = 32;
// 			else if(j==2) co = 64;

// 			for(int k=0; k<=end_num; k++)	// 0 ~ 4
// 			{
// 				stage = 2*((end_num+1)*j+k)+1;
// 				cout << "layer " << stage << endl;
// 				output << "layer " << stage << endl;
// 				temp = cnn;
// 				if(j>=1 && k==0) st = 2;
// 				else st = 1;
// 				fh = 3, fw = 3;
// 				multiplexed_parallel_convolution_print(cnn, cnn, co, st, fh, fw, conv_weight[stage], bn_running_var[stage], bn_weight[stage], epsilon, encoder,
// encryptor, scale_evaluator, gal_keys, cipher_pool, output, decryptor, context, stage); 				multiplexed_parallel_batch_norm_seal_print(cnn, cnn,
// bn_bias[stage], bn_running_mean[stage], bn_running_var[stage], bn_weight[stage], epsilon, encoder, encryptor, scale_evaluator, B, output, decryptor, context,
// stage); 				if(j==0) bootstrap_print(cnn, cnn, bootstrapper_1, output, decryptor, encoder, context, stage); 				else if(j==1)
// bootstrap_print(cnn, cnn, bootstrapper_2, output, decryptor, encoder, context, stage); 				else if(j==2) bootstrap_print(cnn, cnn, bootstrapper_3,
// output, decryptor,
// encoder, context, stage); 				approx_ReLU_seal_print(cnn, cnn, comp_no, deg, alpha, tree, scaled_val, logp, encryptor, evaluator, scale_evaluator,
// decryptor, encoder, public_key, secret_key, relin_keys, B, output, context, gal_keys, stage);

// 				stage = 2*((end_num+1)*j+k)+2;
// 				cout << "layer " << stage << endl;
// 				output << "layer " << stage << endl;
// 				st = 1, fh = 3, fw = 3;
// 				multiplexed_parallel_convolution_print(cnn, cnn, co, st, fh, fw, conv_weight[stage], bn_running_var[stage], bn_weight[stage], epsilon, encoder,
// encryptor, scale_evaluator, gal_keys, cipher_pool, output, decryptor, context, stage); 				multiplexed_parallel_batch_norm_seal_print(cnn, cnn,
// bn_bias[stage], bn_running_mean[stage], bn_running_var[stage], bn_weight[stage], epsilon, encoder, encryptor, scale_evaluator, B, output, decryptor, context,
// stage); 				if(j>=1
// && k==0)
// 				{
// 					st = 2, fh = 1, fw = 1;
// 					multiplexed_parallel_convolution_print(temp, temp, co, st, fh, fw, shortcut_weight[j-1], shortcut_bn_var[j-1], shortcut_bn_weight[j-1],
// epsilon, encoder, encryptor, scale_evaluator, gal_keys, cipher_pool, output, decryptor, context, stage);
// multiplexed_parallel_batch_norm_seal_print(temp, temp, shortcut_bn_bias[j-1], shortcut_bn_mean[j-1], shortcut_bn_var[j-1], shortcut_bn_weight[j-1], epsilon,
// encoder, encryptor, scale_evaluator, B, output, decryptor, context, stage);
// 				}
// 				cipher_add_seal_print(temp, cnn, cnn, scale_evaluator, output, decryptor, encoder, context);
// 				if(j==0) bootstrap_print(cnn, cnn, bootstrapper_1, output, decryptor, encoder, context, stage);
// 				else if(j==1) bootstrap_print(cnn, cnn, bootstrapper_2, output, decryptor, encoder, context, stage);
// 				else if(j==2) bootstrap_print(cnn, cnn, bootstrapper_3, output, decryptor, encoder, context, stage);
// 				approx_ReLU_seal_print(cnn, cnn, comp_no, deg, alpha, tree, scaled_val, logp, encryptor, evaluator, scale_evaluator, decryptor, encoder,
// public_key, secret_key, relin_keys, B, output, context, gal_keys, stage);
// 			}
// 		}
// 		cout << "layer " << layer_num - 1 << endl;
// 		output << "layer " << layer_num - 1 << endl;
// 		averagepooling_seal_scale_print(cnn, cnn, scale_evaluator, gal_keys, B, output, decryptor, encoder, context);
// 		fully_connected_seal_print(cnn, cnn, linear_weight, linear_bias, 100, 64, scale_evaluator, gal_keys, output, decryptor, encoder, context);

// 		total_time_end = chrono::high_resolution_clock::now();
// 		total_time_diff = chrono::duration_cast<chrono::milliseconds>(total_time_end - total_time_start);

// 		// final text file print
// 		Plaintext plain;
// 		decryptor.decrypt(cnn.cipher(), plain);
// 		vector<complex<double>> rtn_vec;
// 		encoder.decode(plain, rtn_vec, 1<<logn);
// 		cout << "final decrypted values: " << endl;
// 		output << "final decrypted values: " << endl;
// 		cout << "( ";
// 		output << "( ";
// 		for (size_t i = 0; i < 99; i++) {
// 			cout << rtn_vec[i] << ", ";
// 			output << rtn_vec[i] << ", ";
// 		}
// 		cout << rtn_vec[99] << ")" << endl;
// 		output << rtn_vec[99] << ")" << endl;
// 		cout << "total time : " << total_time_diff.count() / 1000 << " ms" << endl;
// 		output << "total time : " << total_time_diff.count() / 1000 << " ms" << endl;

// 		size_t label = 0;
// 		double max_score = -100.0;
// 		for(size_t i=0; i<100; i++)
// 		{
// 			if(max_score < rtn_vec[i].real())
// 			{
// 				label = i;
// 				max_score = rtn_vec[i].real();
// 			}
// 		}
// 		cout << "image label: " << image_label << endl;
// 		cout << "inferred label: " << label << endl;
// 		cout << "max score: " << max_score << endl;
// 		output << "image label: " << image_label << endl;
// 		output << "inferred label: " << label << endl;
// 		output << "max score: " << max_score << endl;
// 		out_share << "image_id: " << image_id << ", " << "image label: " << image_label << ", inferred label: " << label << endl;

// 	}

// 	all_time_end = chrono::high_resolution_clock::now();
// 	all_time_diff = chrono::duration_cast<chrono::milliseconds>(all_time_end - all_time_start);
// 	cout << "all threads time : " << all_time_diff.count() / 1000 << " ms" << endl;
// 	out_share << endl << "all threads time : " << all_time_diff.count() / 1000 << " ms" << endl;

// }
