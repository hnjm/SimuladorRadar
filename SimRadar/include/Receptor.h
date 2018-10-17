#ifndef RECEPTOR_H
#define RECEPTOR_H

#include <fstream>
#include <complex>
#include <vector>
#include <iostream>
#include <cmath>
#include "Objetos.h"

const double PI = 3.141592653589793;

class Receptor
{
    public:
        Receptor(int , int , double , double, double, double , Objetos Objc);
        virtual ~Receptor();

        void CreaM_datos() ;
        std::vector<std::vector<std::complex<double>>> Simula() ;

        double Potencia_receptor(double t, double rk , double c);


    protected:

    private:
        double fs, te, tau, Fc;
        int N, M;
        Objetos Obj;


};

#endif // RECEPTOR_H
