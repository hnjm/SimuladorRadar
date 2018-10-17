#include <iostream>

#include "include/Antena.h"
#include "include/Objetos.h"
#include "include/Receptor.h"



/*Main function of the SimRadar simulator.   */

/*Created by Arturo Collado Rosell 16/10/2018. First version. */

double Tfun = 1.0 ; // Time of radar work [seconds]
double Fc = 5e9; // carry frequency [Hz]
double c = 3e8; //speed of light [m/s]
double lambda = c/Fc ; // wavelength [meters]
double sigma_sct = 1.0; // scattering coeficient
double fs = 50e6; // sampling frequency [Hz]
double PRF = 1000; //[HZ]
double te = 10e-6; // [s]tiempo entre que se lanza el pulso y se abre la ventana de recepcion
double Tventana = 1e-5 ; // [s] tiempo que esta abierta la ventana de recepción
double tau = 2.5e-6; // [s] duracion del pulso transmitido
int main()
{
    int numObj = 2;
    std::vector<std::vector<double>> reflector1(numObj,std::vector<double>(3));

    reflector1[0][0] = 2500;
    reflector1[0][1] = 0;
    reflector1[0][2] = 0;

    reflector1[1][0] = 2000;
    reflector1[1][1] = 0;
    reflector1[1][2] = 0;

    std::vector<double> sigma_relf(numObj);
    sigma_relf[0] = sigma_sct;
    sigma_relf[1] = sigma_sct;

    Antena A1 = Antena(0.3,90,90,2,0.75);
    Objetos O1 = Objetos(reflector1, sigma_relf, std::vector<std::vector<double>>(numObj,std::vector<double>(3,0)));

    std::vector<double> aux( O1.get_pos(0));
    std::cout<<aux[0] <<"  "<<aux[1]<<"  "<<aux[2]<<std::endl;

    //Creo un receptor
    int N = int(Tfun*PRF),M =int(Tventana*fs) ;
    Receptor R1 = Receptor(N,M,fs,te, tau, Fc, O1);

    R1.CreaM_datos();


    return 0;
}
