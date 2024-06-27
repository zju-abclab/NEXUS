import numpy as np
import matplotlib.pyplot as plt

# Define the GELU function
def gelu(x):
    return x * 0.5 * (1 + np.tanh(np.sqrt(2 / np.pi) * (x + 0.044715 * x**3)))

# def layernorm(x):


# Define the interval over which to approximate
interval = [-5, 5]

# Generate points within the interval
x = np.linspace(interval[0], interval[1], 32768)
y = gelu(x)

np.savetxt("data/input/gelu_input_32768.txt", x, fmt='%f',delimiter=' ')
np.savetxt("data/calibration/gelu_calibration_32768.txt", y, fmt='%f',delimiter=' ')


x = np.linspace(interval[0], interval[1], 768)

