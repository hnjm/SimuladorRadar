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


void Objetos::actualiza_pos(int n, double Ts)
{
    std::vector<double> auxpos(this->reflector[n]);
    auxpos[0] = auxpos[0] + Ts*this->vel[n][0];
    auxpos[1] = auxpos[1] + Ts*this->vel[n][1];
    auxpos[2] = auxpos[2] + Ts*this->vel[n][2];

    set_pos(auxpos,n);

}

std::vector<double> Objetos::get_angles(int n)
{
    std::vector<double> angles_n(2);
    angles_n[0] = std::atan(this->reflector[n][1]/this->reflector[n][0]);
    double aux = std::sqrt(this->reflector[n][0]*this->reflector[n][0] + this->reflector[n][1]*this->reflector[n][1]);
    angles_n[1] = std::atan(aux/this->reflector[n][2]);

    return angles_n;



}



