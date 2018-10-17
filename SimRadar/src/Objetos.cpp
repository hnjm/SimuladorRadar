#include "../include/Objetos.h"
#include <iostream>

//Constructors
Objetos::Objetos()
{
    //ctor
}

Objetos::Objetos(std::vector<std::vector<double>> reflectorc, std::vector<double> srcc, std::vector<std::vector<double>> velc):reflector(reflectorc),src(srcc), vel(velc)
{
    std::cout<<"Los objetos se construyeron "<<std::endl;
}

//Destructor
Objetos::~Objetos()
{
    std::cout<<"Los objetos se destruyeron"<<std::endl;
}


// Other methods

std::vector<double> Objetos::get_pos(int n)
{
    return this->reflector[n];
}

void Objetos::set_pos(std::vector<double> vectorI, int n)
{
    this->reflector[n] = vectorI;
}

double Objetos::get_dist(int n)
{
    return std::sqrt(reflector[n][0]*reflector[n][0] + reflector[n][1]*reflector[n][1] + reflector[n][2]*reflector[n][2]);
}


int Objetos::numbReflec()
{
    return reflector.size();
}




