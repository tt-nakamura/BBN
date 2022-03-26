#ifndef __BBN_h__
#define __BBN_h__

#include "nr.h"
#include "constants.h"

void expansion_init(double, double=3);
void expansion(double&, double&, double);
double expansion_rate(double, double);
void weak_rate(double&, double&, double, double, double=tau_n);

void spline(const Vec_DP &x, const Vec_DP &y, double yp1, double ypn, Vec_DP &y2);
double splint(const Vec_DP &xa, const Vec_DP &ya, const Vec_DP &y2a, double x);

struct particle {
    double mass;
    int spin;// statistical weight
    int A;// mass numbers
    std::string name;
    particle(double, int, const char*);
    inline bool operator==(const particle& p) { return name==p.name; }
};

struct BBN {
    static int N_element;// number of elements
    static int N_reaction;// number of nuclear reactions
    static int n_index;// index of neutron
    static int p_index;// index of proton
    static Mat_INT index;// table of particles
    static particle element[];// list of synthesized elements
    static particle reaction[][4];// list of nuclear reactions
    static void reaction_init();
    static void reaction_rate(Vec_DP&, double);
    static void reaction_rate(Vec_DP&, Vec_DP&, double);
    // interpolation variables
    static Vec_DP x0,x1,y0,y1,y2,y3,y4,dy0,dy1,dy2,dy3,dy4;
    static void interp_init(double, double, double);

    static Vec_DP y;// abundance of elements (dependent variable)
    static double time;// time since T=T0 / sec
    static double n0;// number density of nucleons at T_nu=1MeV
    static double weak;// weak interaction strength (1 if tau=tau_n)
    static void init(double, double, double, double=3, double=tau_n);
    static void diff_eq(double, const Vec_DP&, Vec_DP&);
    static void jac(double, const Vec_DP&, Vec_DP&, Mat_DP&);
    static void set_temperature(double);

    inline static double temperature(double t)
    // given time t / sec, return temperature / MeV
    { return splint(x0,y0,dy0,t); }
    inline static double expansion_time(double T)
    // given temperature T / MeV, return time since T=T_init / sec
    { return splint(x1,y1,dy1,T); }
    inline static double neutrino_temperature(double T)
    // given temperature T / MeV, return neutrino temperature / MeV
    { return splint(x1,y2,dy2,T)*T; }
    inline static double proton_to_neutron(double T)
    // given temperature T / MeV, return weak interaction rate p_n / sec^-1
    { return splint(x1,y3,dy3,T)*weak; }
    inline static double neutron_to_proton(double T)
    // given temperature T / MeV, return weak interaction rate n_p / sec^-1
    { return splint(x1,y4,dy4,T)*weak; }
    inline static double mass_fraction(int i)
    // given index i, return mass fraction of element i
    { return element[i].A * y[i]; }
};

void gaulag(Vec_DP &x, Vec_DP &w, double alf);

void odeint(Vec_DP& y, void f(double, const Vec_DP& y, Vec_DP& f),
            double a, double b, double eps);

void odeint(Vec_DP& y, void f(double, const Vec_DP& y, Vec_DP& f),
            void jac(double, const Vec_DP&, Vec_DP&, Mat_DP&),
            double a, double b, double eps);

#endif // _BBN_h__
