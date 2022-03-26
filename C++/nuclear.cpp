// configure nuclear reactions

#include<cmath>
#include "BBN.h"

particle::particle(double m, int s, const char *n)
: mass(m), spin(s), A(int(round(m/amu))), name(n) {;}

particle photon(0, 2, "photon");
particle electron(0.510998911, 2, "electron");
particle proton(938.272013, 2, "proton");
particle neutron(939.565346, 2, "neutron");
particle deutron(1875.612793, 3, "deutron");
particle tritium(2808.920906, 2, "tritium");
particle helium3(2808.391383, 2, "helium3");
particle helium4(3727.379109, 1, "helium4");
particle lithium7(6533.833166, 4, "lithium7");
particle beryllium7(6534.184060, 4, "beryllium7");

particle BBN::element[] = {
    neutron,
    proton,
    deutron,
    tritium,
    helium3,
    helium4,
    lithium7,
    beryllium7
};

// reaction[][0:2] = destroyed particles
// reaction[][2:4] = created particles
particle BBN::reaction[][4] = {
    { neutron, proton,  deutron, photon },
    { deutron, proton,  helium3, photon },
    { deutron, deutron, helium3, neutron },
    { deutron, deutron, tritium, proton },
    { helium3, neutron, tritium, proton },
    { tritium, deutron, helium4, neutron },
    { helium3, deutron, helium4, proton },
    { helium3, helium4, beryllium7, photon },
    { helium4, tritium, lithium7, photon },
    { beryllium7, neutron, lithium7, proton },
    { lithium7, proton, helium4, helium4 }
};

static int N = sizeof(BBN::element)/sizeof(particle);
static int M = sizeof(BBN::reaction)/sizeof(particle)/4;

void BBN::reaction_rate(Vec_DP& r, double T)
// input: T = temperature / MeV
// output: r = <(cross section)(relative velocity)> / cm^3/sec
//             averaged over Maxwellian distribution of v
// reverence: M. S. Smith, L. H. Kawano and R. A. Malaney
//   The Astrophysical Journal Supplement 85 (1993) 219
{
    static double T9K(1.e9*kB);// 10^9 Kelvin / MeV
    double t9(T/T9K), t912(sqrt(t9)), t932(t9*t912);
    double t913(pow(t9, 1./3)), t923(pow(t913, 2)), t943(pow(t923,2)), t953(t9*t923);
    double t9f(t9/(1.0+0.1071*t9)), t9f13(pow(t9f, 1./3)), t9f56(pow(t9f, 5./6));
    double t9e(t9/(1.0+0.1378*t9)), t9e13(pow(t9e, 1./3)), t9e56(pow(t9e, 5./6));
    double t9a(t9/(1.0+13.076*t9)), t9a32(pow(t9a, 1.5));
    double t9d(t9/(1.0+0.759*t9)), t9d13(pow(t9d, 1./3)), t9d56(pow(t9d, 5./6));
    // n + p -> d
    r[0] = 4.742e+4*(1.-.8504*t912+.4895*t9-.09623*t932+8.471e-3*t9*t9-2.80e-4*t9*t932);
    // p + d -> 3He + gamma
    r[1] = 2.65e+3/t923*exp(-3.720/t913)
    *(1.+.112*t913+1.99*t923+1.56*t9+.162*t943+.324*t953);
    // d + d -> n + 3He
    r[2] = 3.95e+8/t923*exp(-4.259/t913)
    *(1.+.098*t913+.765*t923+.525*t9+9.61e-3*t943+.0167*t953);
    // d + d -> p + t
    r[3] = 4.17e+8/t923*exp(-4.258/t913)
    *(1.+.098*t913+.518*t923+.355*t9-.010*t943-.018*t953);
    // n + 3He -> p + t
    r[4] = 7.21e+8*(1.-.508*t912+.228*t9);
    // d + t -> n + 4He
    r[5] = 1.063e+11/t923*exp(-4.559/t913-pow(t9/.0754,2))
    *(1.+.092*t913-.375*t923-.242*t9+33.82*t943+55.42*t953)
    + 8.047e+8/t923*exp(-0.4857/t9);
    // 3He + d -> 4He + p
    r[6] = 5.021e+10/t923*exp(-7.144/t913-pow(t9/.270,2))
    *(1.+.058*t913+.603*t923+.245*t9+6.97*t943+7.19*t953)
    + 5.212e+8/t912*exp(-1.762/t9);
    // 3He + 4He -> 7Be + gamma
    r[7] = 4.817e+6/t923*exp(-14.964/t913)
    *(1.+.0325*t913-1.04e-3*t923-2.37e-4*t9-8.11e-5*t943-4.69e-5*t953)
    + 5.938e+6*t9f56/t932*exp(-12.859/t9f13);
    // 4He + t -> Li7 + gamma
    r[8] = 3.032e+5/t923*exp(-8.090/t913)
    *(1.+.0516*t913+.0229*t923+8.28e-3*t9-3.28e-4*t943-3.01e-4*t953)
    + 5.109e+5*t9e56/t932*exp(-8.068/t9e13);
    // 7Be + n -> 7Li + p
    r[9] = 2.675e+9*(1.-.560*t912+.179*t9-.0283*t932 + 2.214e-3*t9*t9-6.851e-5*t9*t932)
    + 9.391e+8*t9a32/t932 + 4.467e+7/t932*exp(-0.07486/t9);
    // Li7 + p -> 4He + 4He
    r[10] = 1.096e+9/t923*exp(-8.472/t913) - 4.830e+8*t9d56/t932*exp(-8.472/t9d13)
    + 1.06e+10/t932*exp(-30.442/t9) + 1.56e+5/t923*exp((-8.472/t913)-pow(t9/1.696,2))
    *(1.+.049*t913-2.498*t923+.860*t9+3.518*t943+3.08*t953)
    + 1.55e+6/t932*exp(-4.478/t9);
    for(int i=0; i<M; i++) r[i] /= NA;
}

int BBN::N_element(N);// number of elements synthesized
int BBN::N_reaction(M);// number of reactions
int BBN::n_index(-1);// index of neutron
int BBN::p_index(-1);// index of proton
Mat_INT BBN::index(-1,M,4);// table of particles
Vec_DP BBN::y(N);// abundance of elements (dependent variable)

static Vec_DP bind(M);// binding energy / MeV
static Vec_DP balance(M);// balancing factor

void BBN::reaction_init()
{
    static double hc3(pow(hbar*c,3));
    int i,j,k;
    double m[4],g[4];
    
    if(n_index>=0) return;// execute only once
    // search proton and neutron in elements
    for(i=0; i<N; i++) {
        if(element[i] == neutron) n_index = i;
        else if(element[i] == proton) p_index = i;
        if(p_index>=0 && n_index>=0) break;
    }
    if(i==N) nrerror("proton or neutron absent");
    
    // precompute binding energy and balancing factor
    for(i=0; i<M; i++) {
        for(j=0; j<4; j++) {
            m[j] = reaction[i][j].mass;
            g[j] = reaction[i][j].spin;
            for(k=0; k<N; k++)// make table of particles
                if(reaction[i][j] == element[k])
                { index[i][j] = k; break; }
        }
        bind[i] = m[0] + m[1] - m[2] - m[3];
        if(index[i][3] < 0)// in case of photon creation
            balance[i] = g[0]*g[1]/g[2]*pow(m[0]*m[1]/m[2]/2/PI, 1.5)/hc3;
        else if(index[i][2] < 0)// in case of photon creation
            balance[i] = g[0]*g[1]/g[3]*pow(m[0]*m[1]/m[3]/2/PI, 1.5)/hc3;
        else
            balance[i] = g[0]*g[1]/g[2]/g[3]*pow(m[0]*m[1]/m[2]/m[3], 1.5);
    }
}

void BBN::reaction_rate(Vec_DP& r1, Vec_DP& r2, double T)
// input: T = temperature / MeV
// output: r1 = forward reaction rate / cm^3/sec
//         r2 = backward reaction rate / cm^3/sec
// in case of photon creation, r2 is in unit of sec^-1
{
    double T32(pow(T, 1.5));
    reaction_rate(r1,T);
    for(int i=0; i<M; i++) {
        r2[i] = r1[i]*balance[i]*exp(-bind[i]/T);
        if(index[i][3] < 0 || index[i][2] < 0)
            r2[i] *= T32;// in case of photon creation
    }
}