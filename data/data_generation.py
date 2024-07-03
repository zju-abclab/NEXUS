import numpy as np
import matplotlib.pyplot as plt


# Define the GELU function
def gelu(x):
    return x * 0.5 * (1 + np.tanh(np.sqrt(2 / np.pi) * (x + 0.044715 * x**3)))


def layernorm(x):
    return np.sqrt(len(x)) * x / np.sqrt(np.sum(x**2))


def softmax(x):
    return np.exp(x) / np.sum(np.exp(x))


# GELU
interval = [-5, 5]
x = np.linspace(interval[0], interval[1], 32768)
y = gelu(x)
np.savetxt("data/input/gelu_input_32768.txt", x, fmt="%f", delimiter=" ")
np.savetxt("data/calibration/gelu_calibration_32768.txt", y, fmt="%f", delimiter=" ")

# LayerNorm
interval = [-10, 10]
x = np.linspace(interval[0], interval[1], 768)
y = layernorm(x)
np.savetxt("data/input/layernorm_input_16_768.txt", x, fmt="%f", delimiter=" ")
np.savetxt("data/calibration/layernorm_calibration_16_768.txt", y, fmt="%f", delimiter=" ")

# Softmax
interval = [-10, 0]
x = np.linspace(interval[0], interval[1], 128)
y = softmax(x)
np.savetxt("data/input/softmax_input_128_128.txt", x, fmt="%f", delimiter=" ")
np.savetxt("data/calibration/softmax_calibration_128_128.txt", y, fmt="%f", delimiter=" ")


def generate_matrices():
    np.random.seed(42)
    matrix_4096x768 = np.random.uniform(-1, 1, (128 * 32, 768))
    matrix_768x64 = np.random.uniform(-1, 1, (768, 64))

    return matrix_4096x768, matrix_768x64


def save_matrix_to_txt(matrix, filename):
    np.savetxt(filename, matrix, fmt="%.4f", delimiter=" ")


def multiply_matrices(matrix1, matrix2):
    return np.dot(matrix1, matrix2)



matrix_4096x768, matrix_768x64 = generate_matrices()
result_matrix = multiply_matrices(matrix_4096x768, matrix_768x64)


save_matrix_to_txt(matrix_4096x768, "data/input/matrixmul_input_m_128_n_768_k_64_batch_128.txt")
save_matrix_to_txt(matrix_768x64, "data/input/matrix_input_n_768_k_64.txt")
save_matrix_to_txt(result_matrix, "data/calibration/matrix_output_m_128_k_64_batch_128.txt")
