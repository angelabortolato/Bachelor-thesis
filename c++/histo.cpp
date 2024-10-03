#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <string>

int main() {
    std::ifstream infJstep[21];
    std::ofstream fhistJstep[21];
    for (int i = 0; i < 21; i++)
    {
        infJstep[i].open("Jstep"+std::to_string(i)+".txt");
        fhistJstep[i].open("histJstep"+std::to_string(i)+".txt");
    }




    int step=0;
    double width;
    int const n_bin=20;

    double dat[40000];
    double max[21]={0.};
    double min[21]={0.};

    for (int step = 0; step < 21; step++)
    {
        for (int i = 0; i < 40000; i++)
        {
            infJstep[step]>>dat[i];
            if (dat[i]>max[step])
            {
                max[step]=dat[i];
            }
            if (dat[i]<min[step])
            {
                min[step]=dat[i];
            }

        }

            double mean=0;
            double s2=0, s=0;
        for (int i = 0; i < 40000; i++)
        {
            mean+=dat[i];
        }
            mean/=40000;

        for (int i = 0; i < 40000; i++)
        {
            s2+=pow(dat[i]-mean, 2);
        }
            s2/=(40000-1);
            s=sqrt(s2);
        
            std::cout<<step <<"\t" <<mean <<"\t" <<s <<std::endl;

        width=(max[step]-min[step])/n_bin;
        int freq[n_bin]={0};

        for (int i = 0; i < 40000; i++)
        {
           //fhistJstep[step]<< dat[i] <<std::endl;
           int f=floor(dat[i]/width)-floor(min[step]/width);
           //std::cout <<f <<std::endl;
           freq[f]++;
        }
        
        for (int i = 0; i < n_bin; i++)
        {
            fhistJstep[step] <<(i+floor(min[step]/width))*width <<"\t" <<freq[i] <<std::endl;
        }
        
    }
    

    return 0;
}