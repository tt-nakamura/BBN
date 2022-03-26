#include<fstream>
#include<cmath>
#include "BBN.h"

main() {
    int i,n(256);
    std::ofstream f("fig2.txt");
    double T0(100), T1(1e-2), dT(pow(T1/T0, 1./n));
    double p_n, n_p, T, T_nu, t, H;
    expansion_init(T0);
    for(i=0; i<=n; i++) {
        T = T0*pow(dT,i);
        expansion(t, T_nu, T);
        weak_rate(p_n, n_p, T, T_nu);
        H = expansion_rate(T, T_nu);
        f << T << ' ';
        f << p_n << ' ';
        f << n_p << ' ';
        f << H << '\n';
	}
}