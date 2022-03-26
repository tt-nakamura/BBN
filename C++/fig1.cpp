#include<cmath>
#include<fstream>
#include "BBN.h"

main() {
    std::ofstream f("fig1.txt");
    int i,n(256);
    double T0(1e2), T1(1e-2), dT(pow(T1/T0, 1./n));
    double T, T_nu, t;
    expansion_init(T0);
    for(i=0; i<=n; i++) {
        T = T0*pow(dT, i);
        expansion(t, T_nu, T);
        f << t << ' ';
        f << T << ' ';
        f << T_nu << '\n';
    }
}