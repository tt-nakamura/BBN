import matplotlib.pyplot as plt
from BBN import BBN

index = ['neutron', 'proton', 'deutron', 'tritium',
         'helium3', 'helium4', 'lithium7', 'beryllium7']
label = ['n', 'p', 'd', 't', r'He$^3$',
         r'He$^4$', r'Li$^7$', r'Be$^7$']

T,X = BBN(5e-10, index, atol=1e-13)
plt.axis([2, 1e-2, 1e-13, 2])
plt.loglog(T, X.T)
plt.xlabel('T = temperature / MeV')
plt.ylabel('X = mass fraction')
plt.legend(label)
plt.show()
