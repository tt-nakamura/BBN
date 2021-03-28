import numpy as np
import matplotlib.pyplot as plt
from BBN import BBN

index = ['proton', 'deutron', 'tritium', 'helium3', 'helium4']
label = ['p', 'd', 't', r'He$^3$', r'He$^4$']

eta = np.geomspace(1e-11, 1e-8, 50)
X = []
for e in eta:
    T,X1 = BBN(e, index, 2, rtol=1e-5)
    X.append(X1[:,-1])
    print(X[-1])

plt.axis([eta[0], eta[-1], 1e-6, 1])
plt.loglog(eta, X)
plt.xlabel(r'$\eta$ = baryon to photon ratio')
plt.ylabel('X = mass fraction')
plt.legend(label)
plt.show()
