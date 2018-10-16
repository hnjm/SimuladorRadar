#include <iostream>

#include "Antena.h"
#include "Objetos.h"

/*Main function of the SimRadar simulator.   */

/*Created by Arturo Collado Rosell 16/10/2018. First version. */



int main()
{
    std::vector<int> reflector1(3,0);
    reflector1[0] = 10;
    reflector1[1] = 10;
    reflector1[2] = 10;

    Antena A1 = Antena(1.0,12,13,5,6);
    Objetos O1 = Objetos(reflector1);

    std::vector<int> aux( O1.get());
    std::cout<<aux[0] <<"  "<<aux[1]<<"  "<<aux[2]<<std::endl;
    return 0;
}
