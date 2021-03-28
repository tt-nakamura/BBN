import matplotlib.pyplot as plt
from expansion import expansion, expansion_eq, weak_rate

T0,T1 = 100,0.01 # initial and final temperature / MeV
T,T_nu,t = expansion(T0,T1)
p_n,n_p = weak_rate(T, T_nu)
dT_nu_dT, dt_dT = expansion_eq(T,[T_nu,t])
H = -dT_nu_dT/T_nu/dt_dT # cosmic expansion rate
plt.axis([1e2,1e-2,1e-10,1e10])
plt.loglog(T, p_n, label='proton to neutron')
plt.loglog(T, n_p, label='neutron to proton')
plt.loglog(T, H, label='cosmic expansion rate')
plt.xlabel('T = temperature / MeV')
plt.ylabel(r'$\Gamma$ = weak reaction rate / sec$^{-1}$')
plt.legend()
plt.show()
