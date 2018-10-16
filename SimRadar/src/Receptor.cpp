#include "Receptor.h"


Receptor::Receptor(int Nc, int Mc, double fsc, double tec, Objetos Objc): N(Nc), M(Mc), fs(fsc), te(tec), Obj(Objc)
{

}

Receptor::~Receptor()
{
    //dtor
}


void Receptor::CreaM_datos() const
{
    std::fstream archivo;
    archivo.open("data.txt", std::fstream::out);

    // body here





    archivo.close();
}
