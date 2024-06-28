import numpy as np
import matplotlib.pyplot as plt

# Define the GELU function
def gelu(x):
    return x * 0.5 * (1 + np.tanh(np.sqrt(2 / np.pi) * (x + 0.044715 * x**3)))

def layernorm(x): 
    return np.sqrt(len(x)) * x / np.sqrt(np.sum(x**2) )

def softmax(x): 
    return np.exp(x) / np.sum(np.exp(x))


# GELU
interval = [-5, 5]
x = np.linspace(interval[0], interval[1], 32768)
y = gelu(x)
np.savetxt("data/input/gelu_input_32768.txt", x, fmt='%f',delimiter=' ')
np.savetxt("data/calibration/gelu_calibration_32768.txt", y, fmt='%f',delimiter=' ')

# LayerNorm
interval = [-10, 10]
x = np.linspace(interval[0], interval[1], 768)
y = layernorm(x)
np.savetxt("data/input/layernorm_input_16_768.txt", x, fmt='%f',delimiter=' ')
np.savetxt("data/calibration/layernorm_calibration_16_768.txt", y, fmt='%f',delimiter=' ')

# Softmax
interval = [-10, 0]
x = np.linspace(interval[0], interval[1], 128)
y = softmax(x)
np.savetxt("data/input/softmax_input_128_128.txt", x, fmt='%f',delimiter=' ')
np.savetxt("data/calibration/softmax_calibration_128_128.txt", y, fmt='%f',delimiter=' ')

