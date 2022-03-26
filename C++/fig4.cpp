#include<cmath>
#include<fstream>
#include "BBN.h"

main() {
    std::ofstream f("fig4.txt");
    int i,j,n(100);
    double T1(10), T2(0.01), dT(pow(T2/T1,1./n));
    double T;
    Vec_DP r(BBN::N_reaction);
    for(int i=0; i<=n; i++) {
        T = T1*pow(dT,i);
        BBN::reaction_rate(r,T);
        f << T;
        for(j=0; j<BBN::N_reaction; j++) f << ' ' << r[j];
        f << '\n';
    }
}