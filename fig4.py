import numpy as np
import matplotlib.pyplot as plt
from nuclear import reaction_rate

label=[
    r'n+p$\to$d',
    r'd+p$\to$He$^3$',
    r'd+d$\to$He$^3$+n',
    r'd+d$\to$t+p',
    r'He$^3$+n$\to$t+p',
    r't+d$\to$He$^4$+n',
    r'He$^3$+d$\to$He$^4$+p',
    r'He$^4$+He$^3\to$Be$^7$',
    r'He$^4$+t$\to$Li$^7$',
    r'Be$^7$+n$\to$Li$^7$+p',
    r'Li$^7$+p$\to$He$^4$+He$^4$'
]

T = np.geomspace(1e-2, 1e1, 128)
R = reaction_rate(T)
plt.axis([1e1, 2e-4, 1e-24, 1e-14])
plt.loglog(T, R.T)
plt.xlabel('T = temperature / MeV')
plt.ylabel(r'$\langle\sigma v\rangle$ = nuclear reaction rate / cm$^3$ sec$^{-1}$')
plt.legend(label)
plt.show()
