"""Calculates 2d Fourier Series numerically

Uses sines/cosines to calculate a given number of term of Fourier Series for any input signal. """
import sys

import numpy as np
import pylab as plt
from scipy.io import FortranFile

NX = 100
NY = 100


def rmse(a1, a2):
    return np.sqrt(np.sum((a1 - a2)**2) / a1.size)


def calc_domain_details(X, Y):
    dx = X[0, 1] - X[0, 0]
    dy = Y[1, 0] - Y[0, 0]
    dA = dx * dy
    Lx = X[0, -1] - X[0, 0]
    Ly = Y[-1, 0] - Y[0, 0]
    return dx, dy, dA, Lx, Ly

def calc_2d_fs(X, Y, sig, N):
    """Calculates 2D Fourier series to N terms over domain defined by X, Y

    Uses sines/cosines, and takes calculations from here:
    https://math.stackexchange.com/a/1695296/109131
    """
    dx, dy, dA, Lx, Ly = calc_domain_details(X, Y)

    alpha = np.zeros((N, N))
    beta = np.zeros((N, N))
    gamma = np.zeros((N, N))
    delta = np.zeros((N, N))

    sig_fs = np.zeros_like(sig)

    # Calculate components:
    # Calculate Fourier Series:
    for i in range(N):
        print('{}/{}'.format(i + 1, N))
        for j in range(N):

            snx = np.sin(2 * np.pi * i * X / Lx)
            sny = np.sin(2 * np.pi * i * Y / Ly)
            cnx = np.cos(2 * np.pi * j * X / Lx)
            cny = np.cos(2 * np.pi * j * Y / Ly)

            if i == 0 and j == 0:
                kappa = 1
            elif i == 0 or j == 0:
                kappa = 2
            else:
                kappa = 4
            alpha[i, j] = kappa * (sig * cnx * cny).sum() * dA / (Lx * Ly)
            beta[i, j]  = kappa * (sig * cnx * sny).sum() * dA / (Lx * Ly)
            gamma[i, j] = kappa * (sig * snx * cny).sum() * dA / (Lx * Ly)
            delta[i, j] = kappa * (sig * snx * sny).sum() * dA / (Lx * Ly)

            sig_fs += (alpha[i, j] * cnx * cny +
                       beta[i, j]  * cnx * sny + 
                       gamma[i, j] * snx * cny + 
                       delta[i, j] * snx * sny)

    return alpha, beta, gamma, delta, sig_fs
                       

def plot_results(signame, X, Y, sig, sig_fs, N, alpha=None, beta=None, gamma=None, delta=None):
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


    if alpha is not None:
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


def load_sig(signame, xlim=(-10, 10), ylim=(-10, 10)):
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
    return sig, X, Y

def main(signame, N):
    print('Run for {}, N={}'.format(signame, N))

    sig, X, Y = load_sig(signame)

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
    elif runtype == 'save':
        sig, X, Y = load_sig(signame)
        dx, dy, dA, Lx, Ly = calc_domain_details(X, Y)
        print(f'dx = {dx}, dy = {dy}, dA = {dA}, Lx = {Lx}, Ly = {Ly}')
        sig.T.astype(np.float64).tofile('data/sig.bin')
        X.T.astype(np.float64).tofile('data/X.bin')
        Y.T.astype(np.float64).tofile('data/Y.bin')
    elif runtype == 'plot_saved':
        ff = FortranFile('data/sig_fs.bin')
        sig_fs = []
        while True:
            try:
                sig_fs.append(ff.read_reals(dtype=np.float64))
            except:
                break
        sig_fs = np.array(sig_fs)
        sig, X, Y = load_sig(signame)
        plot_results(signame, X, Y, sig, sig_fs, N)
    else:
        raise Exception('Unkown runtype: {}'.format(runtype))
    plt.show()
