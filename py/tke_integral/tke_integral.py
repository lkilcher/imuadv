from __future__ import division
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

a, b, pow, k, tke, eps, alpha, beta, beta2 = sp.var('a b pow k tke eps alpha beta beta2',
                                                    real=True, positive=True)


def sinc(x):
    return sp.sin(x) / x

spec = a / (1 + b * k ** pow)

intg = sp.integrate(spec, (k, 0, sp.oo), conds='separate')[0]

# beta = 3 pi / 5 = pi / pow
tmp = sp.solve(intg.replace(b, a / (alpha * eps ** (pow - 1))) - tke,
               a)[0]
# print sp.N(tmp.replace(pow, 5 / 3))
# print tmp.replace(sp.pi, beta * pow).replace(pow, 5 / 3)


def spec_epsilon(epsilon, tke, k, alpha=0.5, ):
    beta = np.pi * 3 / 5
    a = (np.sin(beta) / beta * tke) ** (5 / 2) / (alpha ** (3 / 2) * epsilon)
    b = a / (alpha * epsilon ** (2 / 3))
    return a / (1 + b * k ** (5 / 3))

line = {'f': np.logspace(-9, 5, 1000)}
line['U'] = 1.0
line['epsilon'] = 1e-4
line['tke'] = 0.03
line['k'] = line['f'] * 2 * np.pi * line['U']
line['spec'] = spec_epsilon(
    line['epsilon'], line['tke'], line['k'])

plt.figure(1)
plt.clf()
ax = plt.gca()
ax.loglog(line['k'], line['spec'])
ax.set_xlim([1e-3, 10])
ax.set_ylim([1e-5, 10])

plt.figure(2)
plt.clf()
ax = plt.gca()
ax.plot(line['k'], line['spec'])
ax.set_xlim([1e-3, 10])
ax.set_ylim([1e-5, 10])

assert np.abs(np.trapz(line['spec'], line['k']) - line['tke']) / line['tke'] < 0.01, 'The spec_epsilon routine is not working correctly.'
