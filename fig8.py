import numpy as np
import matplotlib.pyplot as plt
from BBN import BBN,initialize

index = ['helium4']
label = [r'$N_{\nu}$ = 2', '3', '4']

eta = np.geomspace(3e-11, 1e-8, 50)
X = []
for N_nu in [2,3,4]:# number of neutrino generation
    initialize(N_nu=N_nu)
    X1 = []
    for e in eta:
        T,X2 = BBN(e, index, 2, rtol=1e-6, atol=1e-9)
        X1.append(X2[0,-1])
        print(e, X1[-1])

    X.append(X1)

plt.axis([eta[0], eta[-1], 0.17, 0.27])
plt.semilogx(eta, np.asarray(X).T)
plt.xlabel(r'$\eta$ = baryon to photon ratio')
plt.ylabel(r'$X_4$ = mass fraction of He$^4$')
plt.legend(label, markerfirst=False)
plt.show()
