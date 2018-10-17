#ifndef RECEPTOR_H
#define RECEPTOR_H

#include <fstream>
#include <complex>
#include <vector>
#include <iostream>
#include "Objetos.h"

class Receptor
{
    public:
        Receptor(int , int , double , double , Objetos Objc);
        virtual ~Receptor();

        void CreaM_datos()const ;
        std::vector<std::vector<std::complex<double>>> Simula() const;


    protected:

    private:
        double fs, te;
        int N, M;
        Objetos Obj;


};

#endif // RECEPTOR_H
