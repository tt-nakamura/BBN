// Big-Bang Nucleosynthesis
// reference:
//   P. J. E. Peebles
//     The Astrophysical Journal 146 (1966) 542
//   E. W. Kolb and M. S. Turner
//     "The Early Universe" chapter 4

#include<cmath>
#include "BBN.h"

double BBN::time;// time since T=T0 / sec
double BBN::n0;// number density of nucleons at T_nu=1MeV
double BBN::weak;// weak interaction strength (1 if tau=tau_n)

void BBN::init(double eta, double T_init, double T_final, double N_nu, double tau)
// eta = baryon to photon ratio (Kolb & Turner eq.3.104)
// T_init = initial temperature / MeV
// T_final = final temperature / MeV
// N_nu = number of neutrino generation
// tau = neutron lifetime / sec
{
    static double Q(mn-mp);// mass difference of neutron and proton
    // number density of nucleons at T_nu=1MeV
    n0 = 11./4.*eta*2*zeta3/PI/PI/pow(hbar*c,3);
    weak = tau_n/tau;// weak interaction strength
    interp_init(T_init, T_final, N_nu);
    reaction_init();
    // initial condition
    time = expansion_time(T_init);
    y = 0.;
    y[n_index] = 1/(exp(Q/T_init) + 1);// neutron
    y[p_index] = 1/(1 + exp(-Q/T_init));// proton
}

void BBN::diff_eq(double t, const Vec_DP& y, Vec_DP& f)
// input:
//   t = time since T=T0 / sec
//   y = number density of elements / that of nucleons
// output:
//   f = dy/dt (right hand side of differential equation)
{
    int i,j;
    double T, p_n, n_p, N, x[4], X[N_element], a;
    Vec_DP r1(N_reaction), r2(N_reaction);
    const int *id;

    T = temperature(t);
    p_n = proton_to_neutron(T);
    n_p = neutron_to_proton(T);
    reaction_rate(r1, r2, T);
    // number density of nucleons / cm^-3
    N = n0*pow(neutrino_temperature(T), 3);
    // number density of elements / cm^-3
    for(i=0; i<N_element; i++) X[i] = N*y[i];
    a = y[p_index]*p_n - y[n_index]*n_p;
    f = 0.;
    f[n_index] += a;
    f[p_index] -= a;
    for(i=0; i<N_reaction; i++) {
        id = index[i];
        for(j=0; j<4; j++) x[j] = (id[j]<0 ? 1 : X[id[j]]);
        a = (r1[i]*x[0]*x[1] - r2[i]*x[2]*x[3])/N;
        // halve rate for identical particles
        if(id[0] == id[1]) a *= 0.5;
        for(j=0; j<4; j++)
            if(id[j]>=0) f[id[j]] -= (j<2 ? a : -a);
    }
}

void BBN::jac(double t, const Vec_DP& y, Vec_DP& fx, Mat_DP& fy)
// jacobian of diff_eq, passed to stiff equation solver
// input: same as diff_eq
// output: fx = df/dt, fy = df/dy
{
    static double EPS(1e-8);// for finite difference
    int i,j,k;
    double T, p_n, n_p, N, x[4], X[N_element], a, h(EPS*t);
    Vec_DP r1(N_reaction), r2(N_reaction), f(N_element);
    const int *id;

    // numerical derivative
    diff_eq(t-h,y,f);
    diff_eq(t+h,y,fx);
    for(i=0; i<N_element; i++) fx[i] = (fx[i]-f[i])/h/2;
    // analytic jacobian
    T = temperature(t);
    p_n = proton_to_neutron(T);
    n_p = neutron_to_proton(T);
    reaction_rate(r1, r2, T);
    // number density of nucleons / cm^-3
    N = n0*pow(neutrino_temperature(T), 3);
    // number density of elements / cm^-3
    for(i=0; i<N_element; i++) X[i] = N*y[i];
    fy = 0.;
    fy[n_index][n_index] -= n_p;
    fy[n_index][p_index] += p_n;
    fy[p_index][n_index] += n_p;
    fy[p_index][p_index] -= p_n;
    for(i=0; i<N_reaction; i++) {
        id = index[i];
        for(j=0; j<4; j++) x[j] = (id[j]<0 ? 1 : X[id[j]]);
        // halve rate for identical particles
        a = (id[0] == id[1] ? 0.5 : 1);
        for(j=0; j<4; j++)
            for(k=0; k<4; k++)
                if(id[j]>=0 && id[k]>=0)
                    fy[id[j]][id[k]] -= (j<2 ? a : -a)*(k<2 ? r1[i] : -r2[i])*x[k^1];
    }
}

void BBN::set_temperature(double T) {
    double t(expansion_time(T));
    odeint(y, diff_eq, jac, time, t, 1e-6);
    time = t;
}