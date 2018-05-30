"""Calculates 2d Fourier Series numerically

Uses sines/cosines to calculate a given number of term of Fourier Series for any input signal. """
import sys

import numpy as np
import pylab as plt

NX = 100
NY = 100


def rmse(a1, a2):
    return np.sqrt(np.sum((a1 - a2)**2) / a1.size)


def calc_2d_fs(X, Y, sig, N):
    """Calculates 2D Fourier series to N terms over domain defined by X, Y

    Uses sines/cosines, and takes calculations from here:
    https://math.stackexchange.com/a/1695296/109131
    """
    dx = X[0, 1] - X[0, 0]
    dy = Y[1, 0] - Y[0, 0]
    dA = dx * dy
    Lx = X[0, -1] - X[0, 0]
    Ly = Y[-1, 0] - Y[0, 0]

    # N.B. closures over X/Y/Lx/Ly for convenience. i.e. they know about domain size.
    def sx(n):
        return np.sin(2 * np.pi * n * X / Lx)
    def sy(n):
        return np.sin(2 * np.pi * n * Y / Ly)
    def cx(n):
        return np.cos(2 * np.pi * n * X / Lx)
    def cy(n):
        return np.cos(2 * np.pi * n * Y / Ly)

    alpha = np.zeros((N, N))
    beta = np.zeros((N, N))
    gamma = np.zeros((N, N))
    delta = np.zeros((N, N))

    # Calculate components:
    for i in range(N):
        print('comp: {}/{}'.format(i + 1, N))
        for j in range(N):
            if i == 0 and j == 0:
                kappa = 1
            elif i == 0 or j == 0:
                kappa = 2
            else:
                kappa = 4
            alpha[i, j] = kappa * (sig * cx(i) * cy(j)).sum() * dA / (Lx * Ly)
            beta[i, j] = kappa * (sig * cx(i) * sy(j)).sum() * dA / (Lx * Ly)
            gamma[i, j] = kappa * (sig * sx(i) * cy(j)).sum() * dA / (Lx * Ly)
            delta[i, j] = kappa * (sig * sx(i) * sy(j)).sum() * dA / (Lx * Ly)
            
    sig_fs = np.zeros_like(sig)

    # Calculate Fourier Series:
    for i in range(N):
        print('fs: {}/{}'.format(i + 1, N))
        for j in range(N):
            sig_fs += (alpha[i, j] * cx(i) * cy(j) +
                       beta[i, j] * cx(i) * sy(j) + 
                       gamma[i, j] * sx(i) * cy(j) + 
                       delta[i, j] * sx(i) * sy(j))

    return alpha, beta, gamma, delta, sig_fs
                       

def plot_results(signame, X, Y, sig, sig_fs, N, alpha, beta, gamma, delta):
    extent = (np.min(X), np.max(X), np.min(Y), np.max(Y))
    plt.figure(1)
    plt.clf()
    plt.title('input signal: {}'.format(signame))
    plt.imshow(sig, extent=extent, interpolation='nearest')
    plt.colorbar()
    plt.pause(0.001)

    plt.figure(2)
    plt.clf()
    plt.title('fourier series ({} terms)'.format(N))
    plt.imshow(sig_fs, extent=extent, interpolation='nearest')
    plt.colorbar()
    plt.pause(0.001)

    plt.figure(3)
    plt.clf()
    plt.title('signal - fourier series ({} terms), RMSE: {:.5f}'.format(N, rmse(sig, sig_fs)))
    plt.imshow(sig - sig_fs, extent=extent, interpolation='nearest')
    plt.colorbar()
    plt.pause(0.001)

    plt.figure(4)
    plt.clf()
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row', num=4)
    # f.set_title('compontents')

    min_all = min(map(np.min, [alpha, beta, gamma, delta]))
    max_all = max(map(np.max, [alpha, beta, gamma, delta]))

    ax1.set_title('$\\alpha$ [cos(nx).cos(my)]')
    ax1.imshow(alpha, vmin=min_all, vmax=max_all, interpolation='nearest')

    ax2.set_title('$\\beta$ [cos(nx).sin(my)]')
    ax2.imshow(beta, vmin=min_all, vmax=max_all, interpolation='nearest')

    ax3.set_title('$\\gamma$ [cos(nx).sin(my)]')
    ax3.imshow(gamma, vmin=min_all, vmax=max_all, interpolation='nearest')

    ax4.set_title('$\\delta$ [sin(nx).sin(my)]')
    im = ax4.imshow(delta, vmin=min_all, vmax=max_all, interpolation='nearest')

    f.subplots_adjust(right=0.8)
    cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
    f.colorbar(im, cax=cbar_ax)

    plt.pause(0.001)


def main(signame, N, xlim=(-10, 10), ylim=(-10, 10)):
    print('Run for {}, N={}'.format(signame, N))
    x = np.linspace(xlim[0], xlim[1], NX)
    y = np.linspace(ylim[0], ylim[1], NY)

    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X**2 + Y**2)
    R_off = np.sqrt((X - 3)**2 + (Y - 4)**2)

    if signame == 'gauss':
        sig = np.exp(-R**2)
    if signame == 'offset_gauss':
        sig = np.exp(-R_off**2)
    elif signame == 'cylinder':
        sig = (R < 3).astype(float)
    elif signame == 'square':
        sig = ((X < 0) & (Y < 0)).astype(float) + ((X > 0) & (Y > 0)).astype(float) * -1
    elif signame == 'saw_x':
        sig = X
    elif signame == 'saw_xy':
        sig = X + Y
    elif signame == 'sin_xy':
        sig = np.sin(X + Y)
    elif signame == 'sinc':
        # Stop div by zero errors.
        sig = np.sin(R + 0.00001) / R + 0.00001
    else:
        raise Exception('Unkown signame: {}'.format(signame))

    alpha, beta, gamma, delta, sig_fs = calc_2d_fs(X, Y, sig, N=N)
    plot_results(signame, X, Y, sig, sig_fs, N, alpha, beta, gamma, delta)
    return alpha, beta, gamma, delta, sig_fs


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print('Usage: python {} <runtype> <signame> <n_terms>'.format(sys.argv[0]))
        sys.exit(1)

    runtype = sys.argv[1]
    signame = sys.argv[2]
    N = int(sys.argv[3])
    if runtype == 'single':
        alpha, beta, gamma, delta, sig_fs = main(signame, N)
        # zero out low values:
        for arr in [alpha, beta, gamma, delta]:
            arr[np.abs(arr) < 1e-12] = 0
    elif runtype == 'up_to_N':
        r = ''
        for i in range(1, N):
            alpha, beta, gamma, delta, sig_fs = main(signame, i)
            if r == 'c':
                plt.pause(0.001)
            else:
                r = input('q to quit, c to cont: ')
                if r == 'q':
                    break
    else:
        raise Exception('Unkown runtype: {}'.format(runtype))
    plt.show()
