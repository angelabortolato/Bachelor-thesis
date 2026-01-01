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



/********************************************************************************************************************/


/***** GET TRIALS ***************************************************************************************************************/



struct trial{

  double u[2000];
  int ty;
  double y;
};



trial get_trial(int stim1, int stim2) {
    
    trial this_trial;
    
    int T=1000;
      
    double u[2000];
    
    //int stim1 = int_distribution(generator);
    //int stim2 = int_distribution(generator);
    
    int ty = 2*stim1 + stim2;

    for (int i = 0; i < 2000; i++) {
		 u[i] = 0.0;
	}

    for (int i = 0; i < 200; i++) {
        u[T*stim1+i] = 1.0;
    }
    for (int i = 400; i < 600; i++) {
        u[T*stim2+i] = 1.0;
    }
    
    double y = -1.0 + 2.0 * (double)((stim1 + stim2)%2);

    //printf("%d %d %d %f\n",stim1,stim2,ty,y);

 	   
	std::copy(std::begin(u), std::end(u), std::begin(this_trial.u)); // Copy local array to struct array
	
	this_trial.ty=ty;
	this_trial.y=y;

    return this_trial;
    
}


/********************************************************************************************************************/

/* all matrices are implemented as  unidimensional vectors

Let A of dimension (N,M)
Implemented as vector of size N*M
Entry A[i][j] -->   A[M*i+j]
*/

int main() {

    int N = 200;
    int M = 2;
    double tau = 30.0;
    double g = 1.5;
    double eta = 0.1;
    double alpha = 0.5;
    double dt = 1.0;
    int T = 1000;
    int n_trials = 10000;
    int tottime = T * n_trials;
    double b = dt / tau;
    int maxdelay=1;
    double beta=0.05;
    double th=1e-3;

    double J[N*N];
    double deltaJ[N*N];
    double x[N*(1+maxdelay)];  // X is a matrix (N,5)      X[:,0] --> X(t)     X[:,1] --> X(t-1)    X[:,2] --> X(t-2)   etx.
    double x_ave[N];
    double r[N*(1+maxdelay)];
     	
    double E[N*N];
    double B[N*M];

    double R_ave[4];
    for (int i = 0; i < 4; i++) {
	R_ave[i]=0;
    }
	

	
    double err[n_trials];	

    double z=0;
 
    int seed=13;   /// NOT ALL SEED WORK EQUALLY WELL
   
      
    generator.seed(seed);

	
	// RAMDOMLY INITIALIZE J AND B
	
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
        	double uni = 0.5*(1+distribution(generator));
        	double norm= std::sqrt(2) * erf_inv_approx(2 * uni - 1);
        	
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

    generator.seed(0);

	double deltaE[N*N];
    double deltax[N];
	double rprev[N];
    
	trial thistrial;

	// INITIALIZE X AND R
    
    for (int i = n_fixed; i < N; i++) {
                x[i] = -0.1 + 0.2 * distribution(generator);
                r[i] = 0.0;
                
                for(int k=1;k<=maxdelay;k++) {
	                x[k*N+i] = 0.0;
	                r[k*N+i] = 0.0;
	            }
                
     }

	 // first 4 neurons are 1
     for (int i = 0; i < n_fixed; i++) {
                for(int k=0;k<=maxdelay;k++) {
	                x[k*N+i] = 1.0;
	            }
    
     }
 
	for (int i = 0; i < N; i++) {
			x_ave[i] = x[i];
	}


	// INITIALIZE ELIGIBILITY TRACE 

    for (int i = 0; i < N; i++) {
       for (int j = 0; j < N; j++) {
           E[N*i+j] = 0.0;
           }
    }

	
    std::string SUFFIX = "_1g" + std::to_string(g) + "_th" + std::to_string(th) + "_eta" + std::to_string(eta) + "_alpha" + std::to_string(alpha) + "_rate" + std::to_string(rate) + "_beta" + std::to_string(beta) + "_seed" + std::to_string(seed);
 
	std::ofstream fcurrerr;
    fcurrerr.open("currerr.txt");  


    std::ofstream ferr;
    ferr.open("err" + SUFFIX + ".txt");  
    
	
    std::ofstream fJ;
    fJ.open("J" + SUFFIX + ".txt");
    
    
    std::ofstream fJ0;
    fJ0.open("J0" + SUFFIX + ".txt");

    for (int i = 0; i < N; i++) {
    	for (int j = 0; j < N; j++) {
        	fJ0 << J[N*i+j] << "\t";
        }
     	fJ0 << std::endl;
    }
    fJ0.close();



	/// ITERATIONS s

    for (int ntime = 0; ntime < tottime; ntime++) {
        if (ntime % (1 * T) == 0) {
            std::cout << ntime/T << std::endl;
            for (int i = 0; i < 4; i++) {
                std::cout << R_ave[i] << " ";
            }
            std::cout << std::endl;
        }
        
		        

        int trialtime = ntime % T;   // time within trial
        int curr_trial = ntime / T;   // index of current trial

		// get trial 

        if (trialtime==0) {
        	int stim2=curr_trial%2;
        	int stim1=(curr_trial%4)/2;
        	thistrial=get_trial(stim1, stim2);
			//printf("%d %f  %f  %f %f\n", thistrial.ty, thistrial.u[0],thistrial.u[400],thistrial.u[1000],thistrial.u[1400]);
        }
	
       
       	// used to update x_ave
       /*
        int timerunning = 5;
        int short_term;
        if (trialtime < timerunning) {
            short_term = trialtime+1;
        } 
        else {
            short_term = timerunning;
        }
        */

		// readout        
        if (trialtime >= 800 && trialtime < 1000) {
            z += r[N-1];
    		//printf("%f\n",z); 
        }
                
        
        // get r
        
        for (int i = 0; i < N; i++) {
            r[i] = std::tanh(x[i]);
        }
        
        
        // Update Equation
        		
        double sumJ[N];
        for (int i = 0; i < N; i++) {
        	sumJ[i]=0;
            for (int j = 0; j < N; j++) {
                sumJ[i] += J[N*i+j] * r[j];
            }
        
        }
        double sumB[N];
        
        for (int i = 0; i < N; i++) {
        	sumB[i]=0;
            for (int j = 0; j < M; j++) {
            	double uu= thistrial.u[T*j+trialtime];
                sumB[i] += B[M*i+j] * uu;
            }
        }
        
        // copy X[:,0] into X[:,1]

        for (int i = 0; i < N; i++) {
        	for(int k=maxdelay;k>0;k--) {
	            x[k*N+i] = x[(k-1)*N+i];
			    r[k*N+i] = r[(k-1)*N+i];

			}
        }
 
        
        for (int i = n_fixed; i < N; i++) {
            x[i] = x[N+i] / (1.0 + b) + b / (1.0 + b) * (sumJ[i] + sumB[i]);
            //x[i] = x[N+i] * (1.0 - b) + b  * (sumJ[i] + sumB[i]);
        }
 
 
 		// add noise       
        for (int i = n_fixed; i < N; i++) {
        
        	double n1=0.5*(1+distribution(generator));
        	
        	if(n1<rate) {
       			double noise1=0.5*(distribution(generator));
       			x[i] += noise1;
       			//count++;
        	}
        }


		/*

        for (int i = 0; i < N; i++) {
            x_ave[i] = 0.0;
            for (int j = 0; j < short_term; j++) {
                x_ave[i] += x[j*N+i];
            }
            x_ave[i] /= short_term;
        }
		*/

		
        for (int i = 0; i < N; i++) {
            deltax[i] = x[i] - x_ave[i];
        }
        
    
        
        for (int i = 0; i < N; i++) {
			x_ave[i] =  beta * x_ave[i] + (1.0 - beta) * x[i];
		}


		// update eligibility trace        
        
        for (int i = 0; i < N; i++) {
            rprev[i] = r[N+i];
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                deltaE[N*i+j] = std::pow(deltax[i] * rprev[j], 3);
            }
        }
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                E[N*i+j] += deltaE[N*i+j];
            }
        }


        
        
        if (trialtime == T - 1) {
                
			
            int curr_ty = thistrial.ty;
            z /= 200.0;
                       
            err[curr_trial] = std::abs(thistrial.y - z);
            double R = -err[curr_trial];
            
            printf("type %d Y %f Z %f R %f\n", curr_ty, thistrial.y, z, R);   
            fcurrerr <<err[curr_trial] <<std::endl;         
            
            
            if (corr) {
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        //deltaJ[N*i+j] = eta * E[N*i+j]*R_ave[curr_ty];  // as in the paper
                       
                        deltaJ[N*i+j] = eta * (-R + R_ave[curr_ty]) * E[N*i+j]*R_ave[curr_ty]; // modulate eta

						// weigth clipping
                        if (deltaJ[N*i+j] > th) {
                            deltaJ[N*i+j] =th;
                        }
                        if (deltaJ[N*i+j] < -th) {

                            deltaJ[N*i+j] =-th;
                        }
                        J[N*i+j] += deltaJ[N*i+j];
                    }
                }
                
                R_ave[curr_ty] = alpha * R_ave[curr_ty] + (1.0 - alpha) * R;
                
               // re-initialize E and deltaJ
                
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        E[N*i+j] = 0.0;
                    }
                }
                for (int i = 0; i < N; i++) {
                    for (int j = 0; j < N; j++) {
                        deltaJ[N*i+j] = 0.0;
                    }
                }
            }
 
            
            // RE-INITIALIZE X, R 
            
            for (int i = n_fixed; i < N; i++) {
                x[i] = -0.1 + 0.2 * distribution(generator);
                r[i] = 0.0;
                
                for(int k=1;k<=maxdelay;k++) {
	                x[k*N+i] = 0.0;
	                r[k*N+i] = 0.0;
	            }
                
     		}

		   	for (int i = 0; i < n_fixed; i++) {
                for(int k=1;k<=maxdelay;k++) {
	                x[k*N+i] = 1.0;
	            }
    		}

            
            for (int i = 0; i < N; i++) {
				x_ave[i] = x[i];
			}
            
            z = 0.0;
            
        }
    	
	}  /// end iterations loop



    for (double e : err) {
        ferr << e << std::endl;
    }
    ferr.close();
	
    for (int i = 0; i < N; i++) {
    	for (int j = 0; j < N; j++) {
        	fJ << J[N*i+j] << "\t";
        }
     	fJ << std::endl;
    }
    fJ.close();
    
	

    return 0;
}



