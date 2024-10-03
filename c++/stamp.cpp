#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <string>

std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(-1.0, 1.0);
std::uniform_int_distribution<int> int_distribution(0, 1);

/** INVERT NORMAL CDF *************************************************************************************************/


float my_logf (float);

double erf_inv_approx(double a)
{
    float p, r, t;
    t = fmaf (a, 0.0f - a, 1.0f);
    t = my_logf (t);
    if (fabsf(t) > 6.125f) { // maximum ulp error = 2.35793
        p =              3.03697567e-10f; //  0x1.4deb44p-32 
        p = fmaf (p, t,  2.93243101e-8f); //  0x1.f7c9aep-26 
        p = fmaf (p, t,  1.22150334e-6f); //  0x1.47e512p-20 
        p = fmaf (p, t,  2.84108955e-5f); //  0x1.dca7dep-16 
        p = fmaf (p, t,  3.93552968e-4f); //  0x1.9cab92p-12 
        p = fmaf (p, t,  3.02698812e-3f); //  0x1.8cc0dep-9 
        p = fmaf (p, t,  4.83185798e-3f); //  0x1.3ca920p-8 
        p = fmaf (p, t, -2.64646143e-1f); // -0x1.0eff66p-2 
        p = fmaf (p, t,  8.40016484e-1f); //  0x1.ae16a4p-1 
    } else { // maximum ulp error = 2.35002
        p =              5.43877832e-9f;  //  0x1.75c000p-28 
        p = fmaf (p, t,  1.43285448e-7f); //  0x1.33b402p-23 
        p = fmaf (p, t,  1.22774793e-6f); //  0x1.499232p-20 
        p = fmaf (p, t,  1.12963626e-7f); //  0x1.e52cd2p-24 
        p = fmaf (p, t, -5.61530760e-5f); // -0x1.d70bd0p-15 
        p = fmaf (p, t, -1.47697632e-4f); // -0x1.35be90p-13 
        p = fmaf (p, t,  2.31468678e-3f); //  0x1.2f6400p-9 
        p = fmaf (p, t,  1.15392581e-2f); //  0x1.7a1e50p-7 
        p = fmaf (p, t, -2.32015476e-1f); // -0x1.db2aeep-3 
        p = fmaf (p, t,  8.86226892e-1f); //  0x1.c5bf88p-1 
    }
    r = a * p;
    return r;
}

float my_logf (float a)
{
    float i, m, r, s, t;
    int e;

    m = frexpf (a, &e);
    if (m < 0.666666667f) { // 0x1.555556p-1
        m = m + m;
        e = e - 1;
    }
    i = (float)e;
    /* m in [2/3, 4/3] */
    m = m - 1.0f;
    s = m * m;
    /* Compute log1p(m) for m in [-1/3, 1/3] */
    r =             -0.130310059f;  // -0x1.0ae000p-3
    t =              0.140869141f;  //  0x1.208000p-3
    r = fmaf (r, s, -0.121484190f); // -0x1.f19968p-4
    t = fmaf (t, s,  0.139814854f); //  0x1.1e5740p-3
    r = fmaf (r, s, -0.166846052f); // -0x1.55b362p-3
    t = fmaf (t, s,  0.200120345f); //  0x1.99d8b2p-3
    r = fmaf (r, s, -0.249996200f); // -0x1.fffe02p-3
    r = fmaf (t, m, r);
    r = fmaf (r, m,  0.333331972f); //  0x1.5554fap-2
    r = fmaf (r, m, -0.500000000f); // -0x1.000000p-1
    r = fmaf (r, s, m);
    r = fmaf (i,  0.693147182f, r); //  0x1.62e430p-1 // log(2)
    if (!((a > 0.0f) && (a <= 3.40282346e+38f))) { // 0x1.fffffep+127
        r = a + a;  // silence NaNs if necessary
        if (a  < 0.0f) r = ( 0.0f / 0.0f); //  NaN
        if (a == 0.0f) r = (-1.0f / 0.0f); // -Inf
    }
    return r;
}




int main(){
    int const N = 200;
    int const M = 2;
    double tau = 30.0;
    double g = 1.5;
    double eta = 0.1;
    double alpha = 0.5;
    double dt = 1.0;
    int T = 1000;
    int const n_trials = 10000;
    int tottime = T * n_trials;
    double b = dt / tau;
    int const maxdelay=1;
    double beta=0.05;
    double th=1e-3;

    double J[N*N];
    double deltaJ[N*N];
    double x[N*(1+maxdelay)];  // X is a matrix (N,5)      X[:,0] --> X(t)     X[:,1] --> X(t-1)    X[:,2] --> X(t-2)   etx.
    double x_ave[N];
    double r[N*(1+maxdelay)];
     	
    double E[N*N];
    double B[N*M];

	

	

    double z=0;
 
    int seed=13;   /// NOT ALL SEED WORK EQUALLY WELL
   
      
    //generator.seed(seed);

	
	// RAMDOMLY INITIALIZE J AND B
	
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
        	double uni = 0.5*(1+distribution(generator));
        	double norm= std::sqrt(2) * erf_inv_approx(2*uni-1);
        	
            J[N*i+j] = norm * g / std::sqrt(N);
            
        }
    }

	
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            B[M*i+j] = distribution(generator);
            std::cout << B[M*i+j] << " ";
        }
    }
    
    // entries of read-out neuron are zero
    for (int j = 0; j < M; j++) {
        B[M*(N-1)+j] = 0;
    }
 

    double rate = 1.0 / 1000.0;
    int n_fixed = 4;


    
    int reading_neuron = N-1;

    bool corr = true;

    //generator.seed(0);

	double rprev[N];
    

 


	
 
    
    
    std::ofstream fB;
    fB.open("B_input.txt");



    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            fB << B[M*i+j] << "\t";
        }
        fB<<std::endl;
    }
    fB.close();

    


    return 0;
}