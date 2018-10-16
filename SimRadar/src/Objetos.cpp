#include "Objetos.h"
#include <iostream>

//Constructors
Objetos::Objetos()
{
    //ctor
}

Objetos::Objetos(std::vector<int>reflectorc):reflector(reflectorc)
{
    std::cout<<"El objeto se construyo"<<std::endl;
}

//Destructor
Objetos::~Objetos()
{
    std::cout<<"El objeto se construyo"<<std::endl;
}


// Other methods

std::vector<int> Objetos::get()
{
    return this->reflector;
}

void Objetos::set(std::vector<int> vectorI)
{
    this->reflector = vectorI;
}
