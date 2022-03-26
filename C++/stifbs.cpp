// integration of stiff differential equation
// W. H. Press, et al, "Numerical Recipes" section 16.6

#include <cmath>
#include "nr.h"
using namespace std;

void ludcmp(Mat_DP &a, Vec_INT &indx, double &d);
void lubksb(const Mat_DP &a, const Vec_INT &indx, Vec_DP &b);

void simpr(Vec_I_DP &y, Vec_I_DP &dydx, Vec_I_DP &dfdx, Mat_I_DP &dfdy,
    const DP xs, const DP htot, const int nstep, Vec_O_DP &yout,
    void derivs(const DP, Vec_I_DP &, Vec_O_DP &))
{
    int i,j,nn;
    DP d,h,x;

    int n=y.size();
    Mat_DP a(n,n);
    Vec_INT indx(n);
    Vec_DP del(n),ytemp(n);
    h=htot/nstep;
    for (i=0;i<n;i++) {
        for (j=0;j<n;j++) a[i][j] = -h*dfdy[i][j];
        ++a[i][i];
    }
    ludcmp(a,indx,d);
    for (i=0;i<n;i++)
        yout[i]=h*(dydx[i]+h*dfdx[i]);
    lubksb(a,indx,yout);
    for (i=0;i<n;i++)
        ytemp[i]=y[i]+(del[i]=yout[i]);
    x=xs+h;
    derivs(x,ytemp,yout);
    for (nn=2;nn<=nstep;nn++) {
        for (i=0;i<n;i++)
            yout[i]=h*yout[i]-del[i];
        lubksb(a,indx,yout);
        for (i=0;i<n;i++) ytemp[i] += (del[i] += 2.0*yout[i]);
        x += h;
        derivs(x,ytemp,yout);
    }
    for (i=0;i<n;i++)
        yout[i]=h*yout[i]-del[i];
    lubksb(a,indx,yout);
    for (i=0;i<n;i++)
        yout[i] += ytemp[i];
}

Vec_DP *x_p;
Mat_DP *d_p;

void pzextr(const int iest, const DP xest, Vec_I_DP &yest, Vec_O_DP &yz,
    Vec_O_DP &dy)
{
    int j,k1;
    DP q,f2,f1,delta;

    int nv=yz.size();
    Vec_DP c(nv);
    Vec_DP &x=*x_p;
    Mat_DP &d=*d_p;
    x[iest]=xest;
    for (j=0;j<nv;j++) dy[j]=yz[j]=yest[j];
    if (iest == 0) {
        for (j=0;j<nv;j++) d[j][0]=yest[j];
    } else {
        for (j=0;j<nv;j++) c[j]=yest[j];
        for (k1=0;k1<iest;k1++) {
            delta=1.0/(x[iest-k1-1]-xest);
            f1=xest*delta;
            f2=x[iest-k1-1]*delta;
            for (j=0;j<nv;j++) {
                q=d[j][k1];
                d[j][k1]=dy[j];
                delta=c[j]-q;
                dy[j]=f1*delta;
                c[j]=f2*delta;
                yz[j] += dy[j];
            }
        }
        for (j=0;j<nv;j++) d[j][iest]=dy[j];
    }
}

void (*jacobn_s)(double, const Vec_DP&, Vec_DP&, Mat_DP&);

void stifbs(Vec_IO_DP &y, Vec_IO_DP &dydx, DP &xx, const DP htry,
    const DP eps, Vec_I_DP &yscal, DP &hdid, DP &hnext,
    void derivs(const DP, Vec_I_DP &, Vec_O_DP &))
{
    const int KMAXX=7,IMAXX=KMAXX+1;
    const DP SAFE1=0.25,SAFE2=0.7,REDMAX=1.0e-5,REDMIN=0.7;
    const DP TINY=1.0e-30,SCALMX=0.1;
    bool exitflag=false;
    int i,iq,k,kk,km,reduct;
    static int first=1,kmax,kopt,nvold = -1;
    DP eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
    static DP epsold = -1.0,xnew;
    static Vec_DP a(IMAXX);
    static Mat_DP alf(KMAXX,KMAXX);
    static int nseq_d[IMAXX]={2,6,10,14,22,34,50,70};
    Vec_INT nseq(nseq_d,IMAXX);

    int nv=y.size();
    d_p=new Mat_DP(nv,KMAXX);
    x_p=new Vec_DP(KMAXX);
    Vec_DP dfdx(nv),err(KMAXX),yerr(nv),ysav(nv),yseq(nv);
    Mat_DP dfdy(nv,nv);
    if (eps != epsold || nv != nvold) {
        hnext = xnew = -1.0e29;
        eps1=SAFE1*eps;
        a[0]=nseq[0]+1;
        for (k=0;k<KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
        for (iq=1;iq<KMAXX;iq++) {
            for (k=0;k<iq;k++)
                alf[k][iq]=pow(eps1,(a[k+1]-a[iq+1])/
                    ((a[iq+1]-a[0]+1.0)*(2*k+3)));
        }
        epsold=eps;
        nvold=nv;
        a[0] += nv;
        for (k=0;k<KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
        for (kopt=1;kopt<KMAXX-1;kopt++)
            if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
        kmax=kopt;
    }
    h=htry;
    for (i=0;i<nv;i++) ysav[i]=y[i];
    jacobn_s(xx,y,dfdx,dfdy);
    if (xx != xnew || h != hnext) {
        first=1;
        kopt=kmax;
    }
    reduct=0;
    for (;;) {
        for (k=0;k<=kmax;k++) {
            xnew=xx+h;
//            if (xnew == xx) nrerror("step size underflow in stifbs");
            simpr(ysav,dydx,dfdx,dfdy,xx,h,nseq[k],yseq,derivs);
            xest=SQR(h/nseq[k]);
            pzextr(k,xest,yseq,y,yerr);
            if (k != 0) {
                errmax=TINY;
                for (i=0;i<nv;i++) errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));
                errmax /= eps;
                km=k-1;
                err[km]=pow(errmax/SAFE1,1.0/(2*km+3));
            }
            if (k != 0 && (k >= kopt-1 || first)) {
                if (errmax < 1.0) {
                    exitflag=true;
                    break;
                }
                if (k == kmax || k == kopt+1) {
                    red=SAFE2/err[km];
                    break;
                }
                else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
                    red=1.0/err[km];
                    break;
                }
                else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
                    red=alf[km][kmax-1]*SAFE2/err[km];
                    break;
                }
                else if (alf[km][kopt] < err[km]) {
                    red=alf[km][kopt-1]/err[km];
                    break;
                }
            }
        }
        if (exitflag) break;
        red=MIN(red,REDMIN);
        red=MAX(red,REDMAX);
        h *= red;
        reduct=1;
    }
    xx=xnew;
    hdid=h;
    first=0;
    wrkmin=1.0e35;
    for (kk=0;kk<=km;kk++) {
        fact=MAX(err[kk],SCALMX);
        work=fact*a[kk+1];
        if (work < wrkmin) {
            scale=fact;
            wrkmin=work;
            kopt=kk+1;
        }
    }
    hnext=h/scale;
    if (kopt >= k && kopt != kmax && !reduct) {
        fact=MAX(scale/alf[kopt-1][kopt],SCALMX);
        if (a[kopt+1]*fact <= wrkmin) {
            hnext=h/fact;
            kopt++;
        }
    }
    delete d_p;
    delete x_p;
}

void odeint(Vec_IO_DP &ystart, const DP x1, const DP x2, const DP eps,
            const DP h1, const DP hmin, int &nok, int &nbad,
            void derivs(const DP, Vec_I_DP &, Vec_O_DP &),
            void rkqs(Vec_IO_DP &, Vec_IO_DP &, DP &, const DP, const DP,
                      Vec_I_DP &, DP &, DP &,
                      void (*)(const DP, Vec_I_DP &, Vec_O_DP &)));

void odeint(Vec_DP& y, void f(double, const Vec_DP& y, Vec_DP& f),
                 void jac(double, const Vec_DP& fx, Vec_DP&, Mat_DP& fy),
                 double a, double b, double eps)
// solve initial value problem for stiff differential equation
// input:
//   y = initial value of dependent variables y at x=a
//   f = right hand side of differential equaiton dy/dx = f(x,y)
//   jac = jacobian dfdx and dfdy (used in stiff solver)
//   a,b = span of independent variable x for integration
//   eps = error tolerance (Numerical Recipes, section 16.2)
// output:
//   y = final value of dependent variables y at x=b
{
    if(a==b) return;
    int dum;
    extern int kmax;
    kmax = 0;
    jacobn_s = jac;
    odeint(y,a,b, eps, (b-a)*eps, 0,dum,dum,f,stifbs);
}