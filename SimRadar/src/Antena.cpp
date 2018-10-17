#include "../include/Antena.h"
#include <iostream>

//Constructor
Antena::Antena(double Wc, double titac, double phic, double Dc, double Gc): W(Wc),tita(titac),phi(phic),D(Dc),G(Gc)
{
    std::cout<<"Se construye la antena"<<std::endl;
}

//Destructor
Antena::~Antena()
{
    std::cout<<"Se destruye la antena"<<std::endl;
}
