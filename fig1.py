import matplotlib.pyplot as plt
from expansion import expansion

T0,T1 = 100,0.01 # initial and final temperature / MeV
T,T_nu,t = expansion(T0,T1)
plt.axis([1e-2,1e4,1e-2,1e1])
plt.loglog(t,T, label='photon and electron')
plt.loglog(t,T_nu, label='neutrino')
plt.xlabel('t = time / sec')
plt.ylabel('T = temperature / MeV')
plt.legend()
plt.show()
