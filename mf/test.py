import numpy as np

if __name__ == "__main__":
    B = 10
    A = 5
    a = 0.25 * np.arcsin(B/(np.sqrt(A**2 + B**2)))
    a2 = 0.25 * np.arccos(-A/(np.sqrt(A**2 + B**2)))

    print(a, a2)
    print(np.pi/4)
