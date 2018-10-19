#ifndef RECEPTOR_H
#define RECEPTOR_H

#include <fstream>
#include <complex>
#include <vector>
#include <iostream>
#include <cmath>
#include "Objetos.h"
#include "Antena.h"

const double PI = 3.141592653589793;

class Receptor
{
    public:
        Receptor(int , int , double , double, double, double, double , Objetos Objc, Antena Antc);
        virtual ~Receptor();

        void CreaM_datos() ;
        std::vector<std::vector<std::complex<double>>> Simula() ;

        double Potencia_receptor(double t, double rk , double c);
        void imprimeVentana(std::vector<std::complex<double>>);
        double Pot_recibida(int n);
        double patron_sinc(std::vector<double> angulosAntena, std::vector<double> angulosReflector );


    protected:

    private:
        double fs, te, tau, Fc, PRF;
        int N, M;
        Objetos Obj;
        Antena Ant;


};

#endif // RECEPTOR_H
