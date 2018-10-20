#include <iostream>

#include "include/Antena.h"
#include "include/Objetos.h"
#include "include/Receptor.h"

//const double PI = 3.141592653589793;

/*Main function of the SimRadar simulator.   */

/*Created by Arturo Collado Rosell 16/10/2018. First version. */

double Tfun = 6.0 ; // Time of radar work [seconds]
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
    std::vector<std::vector<double>> velocidadR(numObj,std::vector<double>(3));
    reflector1[0][0] = 0;
    reflector1[0][1] = 1800;
    reflector1[0][2] = 0;

    reflector1[1][0] = 0;
    reflector1[1][1] = 2500;
    reflector1[1][2] = 0;


    velocidadR[0][0] = 0; velocidadR[0][1] = 5; velocidadR[0][2] = 0;
    velocidadR[1][0] = 0 ; velocidadR[1][1] = 10; velocidadR[1][2] = 0;


    std::vector<double> sigma_relf(numObj);
    sigma_relf[0] = sigma_sct;
    sigma_relf[1] = sigma_sct;

    Antena A1 = Antena(2.0*PI/5.0,PI/2.0 - 10*PI/180,PI/2.0,1,0.75);

    Objetos O1 = Objetos(reflector1, sigma_relf, velocidadR);

    std::vector<double> aux( O1.get_pos(0));
    std::vector<double> aux2( O1.get_pos(1));
    std::cout<<aux[0] <<"  "<<aux[1]<<"  "<<aux[2]<<std::endl;
    std::cout<<aux2[0] <<"  "<<aux2[1]<<"  "<<aux2[2]<<std::endl;

    double d1 = O1.get_dist(0),d2 = O1.get_dist(1);

    std::cout<<"Lad distancia de los reflectores son: "<<d1<<" y "<<d2<<std::endl;

    //Creo un receptor
    int N = int(Tfun*PRF),M =int(Tventana*fs) ;
    Receptor R1 = Receptor(N,M,fs,te, tau, Fc, PRF, c, O1, A1);

    R1.CreaM_datos();


    return 0;
}
