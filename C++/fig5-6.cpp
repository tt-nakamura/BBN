#include<cmath>
#include<fstream>
#include "BBN.h"

main() {
    std::ofstream f("fig5-6.txt");
    int i,j,n(256);
    double eta(5e-10), T0(10), T1(0.01);
    double T, dT(pow(T1/T0, 1./n));
    BBN::init(eta, T0, T1);
    for(i=0; i<=n; i++) {
        T = T0*pow(dT,i);
        BBN::set_temperature(T);
        f << T;
        for(j=0; j<BBN::N_element; j++)
            f << ' ' << BBN::mass_fraction(j);
        f << '\n';
    }
}