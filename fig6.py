import matplotlib.pyplot as plt
from BBN import BBN

index = ['neutron', 'proton', 'helium4']
label = ['n', 'p', r'He$^4$']

T,X = BBN(5e-10, index)
plt.axis([1e1, 1e-2, 8e-2, 1])
plt.loglog(T, X.T)
plt.xlabel('T = temperature / MeV')
plt.ylabel('X = mass fraction')
plt.legend(label)
plt.show()
