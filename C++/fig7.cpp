#include<cmath>
#include<fstream>
#include "BBN.h"

main() {
    std::ofstream f("fig7.txt");
    int i,j,n(100);
    double eta0(1e-11), eta1(1e-8), T0(10), T1(0.01);
    double eta, de(pow(eta1/eta0, 1./n));
    for(i=0; i<=n; i++) {
        eta = eta0*pow(de,i);
        BBN::init(eta, T0, T1);
        BBN::set_temperature(T1);
        f << eta;
        for(j=0; j<BBN::N_element; j++)
            f << ' ' << BBN::mass_fraction(j);
        f << '\n';
        std::cout << eta << '\n';
    }
}