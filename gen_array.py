import numpy as np

NX = 100
NY = 100

xlim = (-10, 10)
ylim = (-10, 10)
x = np.linspace(xlim[0], xlim[1], NX, dtype=np.float32)
y = np.linspace(ylim[0], ylim[1], NY, dtype=np.float32)
X, Y = np.meshgrid(x, y)

R = np.sqrt(X**2 + Y**2)
sig = np.exp(-R**2)

# Sanity check - see fortran.
print(sig)

print(sig[5, 5])

sig.T.tofile('data/sig.bin')
X.T.tofile('data/X.bin')
Y.T.tofile('data/Y.bin')
