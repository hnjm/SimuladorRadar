#include <iostream>

#include "include/Antena.h"
#include "include/Objetos.h"
#include "include/Receptor.h"

/*Main function of the SimRadar simulator.   */

/*Created by Arturo Collado Rosell 16/10/2018. First version. */

double Tfun = 5.0 ; // Time of radar work [seconds]
double Fc = 5e9; // carry frequency [Hz]
double c = 3e8; //speed of light [m/s]
double lambda = c/Fc ; // wavelength [meters]
double sigma_sct = 1.0; // scattering coeficient
double fs = 50e-6; // sampling frequency [Hz]
double PRF = 1000; //[HZ]

int main()
{
    std::vector<double> reflector1(3,0);
    reflector1[0] = 10;
    reflector1[1] = 10;
    reflector1[2] = 10;

    Antena A1 = Antena(1.0,12,13,5,6);
    Objetos O1 = Objetos(reflector1, sigma_sct, std::vector<double>(3,0));

    std::vector<double> aux( O1.get_pos());
    std::cout<<aux[0] <<"  "<<aux[1]<<"  "<<aux[2]<<std::endl;

    //Creo un receptor
    Receptor R1 = Receptor(1000,500,4.2,2.3, O1);

    R1.CreaM_datos();


    return 0;
}
