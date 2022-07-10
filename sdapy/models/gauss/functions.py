import numpy as np

def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def double_gauss(x, H, A1, x01, sigma1, A2, x02, sigma2): 
    return H + A1 * np.exp(-(x - x01) ** 2 / (2 * sigma1 ** 2)) \
        + A2 * np.exp(-(x - x02) ** 2 / (2 * sigma2 ** 2))
