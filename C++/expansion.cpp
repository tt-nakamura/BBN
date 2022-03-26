// expansion of the early universe
// reference:
//   R. A. Alpher, J. W. Follin and R. C. Herman
//     Physical Review 92 (1953) 1347
//   E. W. Kolb and M. S. Turner
//     "The Early Universe" chapter 3

#include<cmath>
#include "BBN.h"

static int N(64);// number of nodes for quadrature
static Vec_DP y(2);// dependent variable of diff_eq
static Vec_DP node(N), weight(0.,N), ex(N);
static double temperature;// current temperature
static double a_nu;// neutrino radiation const

void expansion_init(double T0, double N_nu)
// T0 = initial temperature / MeV
// N_nu = number of neutrino generation
{
    int i;
    temperature = T0;
    a_nu = a_rad*0.875*N_nu;// neutrino radiation const
    y[0] = 0; // time
    y[1] = T0;// neutrino temperature
    if(weight[0]) return;
    gaulag(node, weight, 0);// execute only once
    for(i=0; i<N; i++) ex[i] = exp(-node[i]);
}

void expansion_eq(double T, const Vec_DP& y, Vec_DP& f)
// input: T = temperature / MeV
//        y[0] = time since T=T0 / sec
//        y[1] = neutrino temperature / MeV
// output: f[0] = dt/dT
//         f[1] = d(T_nu)/dT
{
    static double GP83(sqrt(8*PI*Grav/3)/hbar);
    int i;
    double a,b,x2,z,u,v,T4;
    double E_nu,E_r,E(0),P(0),C(0),H;
    const double& T_nu(y[1]);
    const Vec_DP& x(node);

    T4 = pow(T,4);
    E_r = a_rad*T4;// photon energy density
    E_nu = a_nu*pow(T_nu, 4);// neutrino enegy density

    a = pow(me/T, 2);
    b = T4*2/PI/PI;
    for(i=0; i<N; i++) {// quadrature
        x2 = x[i]*x[i];
        v = sqrt(x2 + a);
        z = exp(a/(v + x[i]));
        u = x2/(z + ex[i])*weight[i];
        E += v*u;
        P += x2/v/3*u;
        C += v*v*z/(z + ex[i])*u;
    }
    E *= b;// electron energy density / MeV^4/(hbar*c)^3
    P *= b;// electron pressure / MeV^4/(hbar*c)^3
    C *= b;// dE/dlnT = (electron specific heat)*T

    E += E_r;// photon energy density
    P += E_r/3;// photon pressure
    C += 4*E_r;// photon specific heat
    H = GP83*sqrt(E + E_nu);// Hubble expansion rate
    C /= 3*T*(E+P); // dln(T_nu)/dT = -dln(a)/dT
    f[0] = -C/H;
    f[1] = C*T_nu;
}

void expansion(double& t, double& T_nu, double T)
// solve ivp until temperature==T
// input: T = temperature / MeV
// output: t = time since T=T0 / sec
//         T_nu = neutrino temperature / MeV
{
    odeint(y, expansion_eq, temperature, T, 1e-9);
    temperature = T;
    t = y[0];
    T_nu = y[1];
}

double expansion_rate(double T, double T_nu)
// input: T = temperatures / MeV
//        T_nu = neutrino temperature / MeV
// return: Hubble expansion rate dln(a)/dt / sec^-1
{
    Vec_DP y(2),f(2);
    y[1] = T_nu;
    expansion_eq(T,y,f);
    return -f[1]/f[0]/T_nu;
}

void weak_rate(double& p_n, double& n_p, double T, double T_nu, double tau)
// input: T = temperature / MeV
//        T_nu = neutrino temperature / MeV
//        tau = neutron lifetime / sec
// output: p_n = proton to neutron conversion rate / sec^-1
//         n_p = neutron to proton conversion rate / sec^-1
{
    int i;
    double a,b,e,q,z,v,v1,v2,z1,z2,u,kappa;
    const Vec_DP& x(node);

    a = me/T;
    b = T/T_nu;
    q = (mn-mp)/T;
    e = exp(-a);
    p_n = n_p = 0;
    for(i=0; i<N; i++) {// quadrature
        v = x[i] + a;
        z = exp(-v);
        z1 = exp((v+q)*b - v);
        z2 = exp((v-q)*b - v);
        v1 = pow(v+q, 2)/(z1 + z);
        v2 = pow(v-q, 2)/(z2 + z);
        u = v*sqrt(x[i]*(x[i] + 2*a))*e/(1+z)*weight[i];
        p_n += (v1 + v2*z2)*u;
        n_p += (v2 + v1*z1)*u;
    }
    v = 1./1.6361/tau/pow(a,5);
    p_n *= v;
    n_p *= v;
}

static int M(256);// number of grid points for interpolation
Vec_DP BBN::x0(M), BBN::x1(M);// interplation variables
Vec_DP BBN::y0(M), BBN::dy0(M);
Vec_DP BBN::y1(M), BBN::dy1(M);
Vec_DP BBN::y2(M), BBN::dy2(M);
Vec_DP BBN::y3(M), BBN::dy3(M);
Vec_DP BBN::y4(M), BBN::dy4(M);

void BBN::interp_init(double T_init, double T_final, double N_nu)
// initialize interpolation variables
// T_init = max temperature / MeV
// T_final = min temperature / MeV
// N_nu = number of neutrino generation
{
    static double T0(100);// temperature at time=0 / MeV
    static double T_init_(0), T_final_, N_nu_;// previous inputs
    int i,j;
    double T, dT, t, T_nu, p_n, n_p;
    // check if parameters change
    if(T_init == T_init_ && T_final == T_final_ && N_nu == N_nu_)
        return;
    T_init_ = T_init; T_final_ = T_final; N_nu_ = N_nu;

    expansion_init(T0, N_nu);
    dT = pow(T_final/T_init, 1./(M-1));
    for(i=0, j=M-1; i<M; i++, j--) {
        T = T_init*pow(dT,i);
        expansion(t, T_nu, T);
        weak_rate(p_n, n_p, T, T_nu);
        x0[i] = t;// time (incresing order)
        x1[j] = T;// temperature (incresing order)
        y0[i] = T;// temperature (decresing order)
        y1[j] = t;// time (decresing order)
        y2[j] = T_nu/T;
        y3[j] = p_n;
        y4[j] = n_p;
    }
    // spline interpolation
    spline(x0,y0,1e30,1e30,dy0);
    spline(x1,y1,1e30,1e30,dy1);
    spline(x1,y2,1e30,1e30,dy2);
    spline(x1,y3,1e30,1e30,dy3);
    spline(x1,y4,1e30,1e30,dy4);
}