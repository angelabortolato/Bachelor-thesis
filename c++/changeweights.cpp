#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
#include <string>

int main() {
    std::ifstream infJ0;
    infJ0.open("plainJstep0.txt");

    std::ifstream infJfin;
    infJfin.open("plainJfin.txt");

    std::ofstream fchweights;
    fchweights.open("index_changeweights.txt");

    double dat0[40000];
    double datfin[40000];

    int ch_weights[10];
    double diffweights[10]={0.};

    for (int i = 0; i < 40000; i++)
        {
            infJ0>>dat0[i];
            infJfin>>datfin[i];
            //std::cout << dat0[i] <<"\t" <<datfin[i] <<"\t" <<std::abs(dat0[i]-datfin[i]) <<std::endl;
        }

    for (int i = 0; i < 40000; i++)
    {
      //  std::cout<<abs(dat0[i]-datfin[i]);
        for (int w = 0; w < 10; w++)
            {
                if (std::abs(dat0[i]-datfin[i])>diffweights[w])
                {
                    diffweights[w]=std::abs(dat0[i]-datfin[i]);
                    ch_weights[w]=i;
                    w=10;                           //una volta sostituito uno dei pesi collezionati, non continuo a valutare sostituzioni
                }
            
            }

    }
    

    std::cout <<"weights that change the most: " ;
    for (int w = 0; w < 10; w++)
    {
        std::cout <<floor(ch_weights[w]/200) <<"," <<ch_weights[w]-floor(ch_weights[w]/200)*200 <<"\t" <<diffweights[w] <<std::endl;
        fchweights<<ch_weights[w] <<"\t";
    }
    
    std::cout <<std::endl;
    

    return 0;
}