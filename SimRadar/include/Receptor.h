#ifndef RECEPTOR_H
#define RECEPTOR_H

#include <fstream>

#include "Objetos.h"

class Receptor
{
    public:
        Receptor(int , int , double , double , Objetos Objc);
        virtual ~Receptor();

        void CreaM_datos()const ;

    protected:

    private:
        double fs, te;
        int N, M;
        Objetos Obj;


};

#endif // RECEPTOR_H
