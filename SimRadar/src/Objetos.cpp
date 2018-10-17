#include "../include/Objetos.h"
#include <iostream>

//Constructors
Objetos::Objetos()
{
    //ctor
}

Objetos::Objetos(std::vector<double>reflectorc, double srcc, std::vector<double> velc):reflector(reflectorc),src(srcc), vel(velc)
{
    std::cout<<"El objeto se construyo"<<std::endl;
}

//Destructor
Objetos::~Objetos()
{
    std::cout<<"El objeto se construyo"<<std::endl;
}


// Other methods

std::vector<double> Objetos::get_pos()
{
    return this->reflector;
}

void Objetos::set_pos(std::vector<double> vectorI)
{
    this->reflector = vectorI;
}
