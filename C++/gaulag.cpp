// Gauss-Laguerre quadrature
// W. H. Press, et al, "Numerical Recipes" section 4.5

#include <cmath>
#include "nr.h"
using namespace std;

DP gammln(const DP xx)
{
    int j;
    DP x,y,tmp,ser;
    static const DP cof[6]={76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
        -0.5395239384953e-5};
    
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<6;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

void gaulag(Vec_O_DP &x, Vec_O_DP &w, const DP alf)
{
    const int MAXIT=10;
    const DP EPS=1.0e-14;
    int i,its,j;
    DP ai,p1,p2,p3,pp,z,z1;
    
    int n=x.size();
    for (i=0;i<n;i++) {
        if (i == 0) {
            z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
        } else if (i == 1) {
            z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
        } else {
            ai=i-1;
            z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
                  (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
        }
        for (its=0;its<MAXIT;its++) {
            p1=1.0;
            p2=0.0;
            for (j=0;j<n;j++) {
                p3=p2;
                p2=p1;
                p1=((2*j+1+alf-z)*p2-(j+alf)*p3)/(j+1);
            }
            pp=(n*p1-(n+alf)*p2)/z;
            z1=z;
            z=z1-p1/pp;
            if (fabs(z-z1) <= EPS) break;
        }
        if (its >= MAXIT) nrerror("too many iterations in gaulag");
        x[i]=z;
        w[i] = -exp(gammln(alf+n)-gammln(DP(n)))/(pp*n*p2);
    }
}