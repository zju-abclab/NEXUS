import numpy as np
import matplotlib.pyplot as plt


# Define the GELU function
def gelu(x):
    return x * 0.5 * (1 + np.tanh(np.sqrt(2 / np.pi) * (x + 0.044715 * x**3)))


def layernorm(x):
    return np.sqrt(len(x)) * x / np.sqrt(np.sum(x**2))


def softmax(x):
    return np.exp(x) / np.sum(np.exp(x))

def argmax(x):
    max_idx = np.argmax(x)
    return np.eye(len(x))[max_idx]


# GELU
interval = [-5, 5]
x = np.linspace(interval[0], interval[1], 32768)
y = gelu(x)
np.savetxt("input/gelu_input_32768.txt", x, fmt="%f", delimiter=" ")
np.savetxt("calibration/gelu_calibration_32768.txt", y, fmt="%f", delimiter=" ")

# LayerNorm
interval = [-10, 10]
x = np.linspace(interval[0], interval[1], 768)
y = layernorm(x)
np.savetxt("input/layernorm_input_16_768.txt", x, fmt="%f", delimiter=" ")
np.savetxt("calibration/layernorm_calibration_16_768.txt", y, fmt="%f", delimiter=" ")

# Softmax
interval = [-10, 0]
x = np.linspace(interval[0], interval[1], 128)
y = softmax(x)
np.savetxt("input/softmax_input_128_128.txt", x, fmt="%f", delimiter=" ")
np.savetxt("calibration/softmax_calibration_128_128.txt", y, fmt="%f", delimiter=" ")

# Argmax
interval = [-0.5, 0.499999]
x = np.linspace(interval[0], interval[1], 8)
x[3] = 0.5
y = argmax(x)
np.savetxt("input/argmax_input_8.txt", x, fmt="%f", delimiter=" ")
np.savetxt("calibration/argmax_calibration_8.txt", y, fmt="%d", delimiter=" ")

# Matrix Multiplication
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


save_matrix_to_txt(matrix_4096x768, "input/matrixmul_input_m_128_n_768_k_64_batch_128.txt")
save_matrix_to_txt(matrix_768x64, "input/matrix_input_n_768_k_64.txt")
save_matrix_to_txt(result_matrix, "calibration/matrix_output_m_128_k_64_batch_128.txt")
