#!/usr/bin/env python3
"""
Сверка rk.dat / ae.dat / ai.dat с точным решением гармонического осциллятора.
"""
import numpy as np
import sys


def read_ode(path):
    d = np.loadtxt(path)
    return d[:, 0], d[:, 1:]


def exact(t):
    return np.cos(t), -np.sin(t)


def main():
    ok = True
    results = {}
    for fname in ['rk.dat', 'ae.dat', 'ai.dat']:
        t, X = read_ode(fname)
        xe, ve = exact(t)
        err_x = np.max(np.abs(X[:, 0] - xe))
        err_v = np.max(np.abs(X[:, 1] - ve))
        print(f'{fname}: N={len(t):4d}')
        print(f'\tmax|x-cos t|={err_x:.3e},  max|v+sin t|={err_v:.3e}')
        results[fname] = (t, X)
        if err_x > 1e-4 or err_v > 1e-4:
            ok = False

    t_rk, X_rk = results['rk.dat']
    t_ae, X_ae = results['ae.dat']
    t_ai, X_ai = results['ai.dat']
    if len(t_rk) == len(t_ae) == len(t_ai):
        print(f'|AE - RK4|_max = {np.max(np.abs(X_ae - X_rk)):.3e}')
        print(f'|AI - RK4|_max = {np.max(np.abs(X_ai - X_rk)):.3e}')

    print('OK' if ok else 'FAIL')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
