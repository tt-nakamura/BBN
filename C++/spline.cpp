// cubic spline interpolation
// W. H. Press, et al, "Numerical Recipes" section 3.3

#include "nr.h"

void spline(Vec_I_DP &x, Vec_I_DP &y, const DP yp1, const DP ypn,
    Vec_O_DP &y2)
{
    int i,k;
    DP p,qn,sig,un;

    int n=y2.size();
    Vec_DP u(n-1);
    if (yp1 > 0.99e30)
        y2[0]=u[0]=0.0;
    else {
        y2[0] = -0.5;
        u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }
    for (i=1;i<n-1;i++) {
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p=sig*y2[i-1]+2.0;
        y2[i]=(sig-1.0)/p;
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if (ypn > 0.99e30)
        qn=un=0.0;
    else {
        qn=0.5;
        un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }
    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
    for (k=n-2;k>=0;k--)
        y2[k]=y2[k]*y2[k+1]+u[k];
}

void splint(Vec_I_DP &xa, Vec_I_DP &ya, Vec_I_DP &y2a, const DP x, DP &y)
{
    int k;
    DP h,b,a;

    int n=xa.size();
    int klo=0;
    int khi=n-1;
    while (khi-klo > 1) {
        k=(khi+klo) >> 1;
        if (xa[k] > x) khi=k;
        else klo=k;
    }
    h=xa[khi]-xa[klo];
    if (h == 0.0) nrerror("Bad xa input to routine splint");
    a=(xa[khi]-x)/h;
    b=(x-xa[klo])/h;
    y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]
        +(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

double splint(Vec_I_DP &xa, Vec_I_DP &ya, Vec_I_DP &y2a, const DP x)
{
    double y;
    splint(xa, ya, y2a, x, y);
    return y;
}