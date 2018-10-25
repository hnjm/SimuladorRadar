#include <iostream>

#include "include/Antena.h"
#include "include/Objetos.h"
#include "include/Receptor.h"

//const double PI = 3.141592653589793;

/*Main function of the SimRadar simulator.   */

/*Created by Arturo Collado Rosell 16/10/2018. First version. */
/*
static double Tfun = 5.0 ; // Time of radar work [seconds]
static double Fc = 5e9; // carry frequency [Hz]
static double c = 3e8; //speed of light [m/s]
static double lambda = c/Fc ; // wavelength [meters]
static double sigma_sct = 1.0; // scattering coeficient
static double fs = 50e6; // sampling frequency [Hz]
static double PRF = 1000; //[HZ]
static double te = 10e-6; // [s]tiempo entre que se lanza el pulso y se abre la ventana de recepcion
static double Tventana = 1e-5 ; // [s] tiempo que esta abierta la ventana de recepción
static double tau = 2.5e-6; // [s] duracion del pulso transmitido
*/


 double Tfun  ; // Time of radar work [seconds]
 double Fc ; // carry frequency [Hz]
 double c ; //speed of light [m/s]
 double lambda  ; // wavelength [meters]
 double sigma_sct = 1.0; // scattering coeficient
 double fs ; // sampling frequency [Hz]
 double PRF = 1000; //[HZ]
 double te ; // [s]tiempo entre que se lanza el pulso y se abre la ventana de recepcion
 double Tventana ; // [s] tiempo que esta abierta la ventana de recepción
 double tau ; // [s] duracion del pulso transmitido


//Función que lee los parámetros desde el archivo de texto.
void lee_parametros(std::string filename)
{
    std::string line;
    std::fstream archivo;
    std::vector<double> parametros;
    try{
        archivo.open("parametros.txt", std::fstream::in);
        if(!archivo.is_open())
            throw 4;

    }catch(...){std::cout<<"Hubo un problema abriendo el archivo";}

    int i=0;
    while (getline (archivo,line) )
    {
     //std::cout<<line<<std::endl;
     i++;
    parametros.push_back(std::stod (line)) ;
    }


    Tfun = parametros[0] ; // Time of radar work [seconds]
    Fc = parametros[1]; // carry frequency [Hz]
    c = parametros[2]; //speed of light [m/s]
    lambda = c/Fc ; // wavelength [meters]
    fs = parametros[3]; // sampling frequency [Hz]
    PRF = parametros[4]; //[HZ]
    te = parametros[5]; // [s]tiempo entre que se lanza el pulso y se abre la ventana de recepcion
    Tventana = parametros[6] ; // [s] tiempo que esta abierta la ventana de recepción
    tau = parametros[7]; // duracion del pulso transmitido



    archivo.close();
}



int main()
{
    lee_parametros("parametros.txt");

    std::cout<<"Tfun =  "<<Tfun<<std::endl;
    std::cout<<"Fc =  "<<Fc<<std::endl;
    std::cout<<"c =  "<<c<<std::endl;
    std::cout<<"lambda =  "<<lambda<<std::endl;
    std::cout<<"fs =  "<<fs<<std::endl;
    std::cout<<"PRF =  "<<PRF<<std::endl;
    std::cout<<"te =  "<<te<<std::endl;
    std::cout<<"Tventana =  "<<Tventana<<std::endl;
    std::cout<<"tau =  "<<tau<<std::endl;





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

    std::cout<<"La distancia de los reflectores son: "<<d1<<" y "<<d2<<std::endl;

    //Creo un receptor
    int N = int(Tfun*PRF),M =int(Tventana*fs) ;

    Receptor R1 = Receptor(N,M,fs,te, tau, Fc, PRF, c, O1, A1);



    R1.CreaM_datos();


    return 0;
}
