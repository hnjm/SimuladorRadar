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

void Antena::Antena_gira(double Ts)
{
    this->tita = this->tita + this->W*Ts;
    //phi se mantiene constante.
}

void Antena::set_angles(int titas, int phis)
{
    this->tita = titas;
    this->phi = phis;
}

std::vector<double> Antena::get_angles()
{
    std::vector<double> angles(2);
    angles[0] = this->tita;
    angles[1] = this->phi;

    return angles;
}
