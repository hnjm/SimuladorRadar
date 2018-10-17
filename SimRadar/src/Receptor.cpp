#include "../include/Receptor.h"


Receptor::Receptor(int Nc, int Mc, double fsc, double tec, Objetos Objc): N(Nc), M(Mc), fs(fsc), te(tec), Obj(Objc)
{

}

Receptor::~Receptor()
{
    //dtor
}

std::vector<std::vector<std::complex<double>>> Receptor::Simula() const
{
    std::vector<std::vector<std::complex<double>>> matrix (M, std::vector<std::complex<double>> (N) );

    std::cout<<"Las dimensiones de la matriz de datos son: "<< matrix.size()<<" x "<< matrix[0].size()<<std::endl;

    for(int i=0; i<matrix[0].size(); ++i)
    {
        for(int j=0 ; j< matrix.size() ; ++j)
        {
            matrix[j][i] = std::complex<double>(j+i, 1 );

        }

    }




    return matrix;
}




void Receptor::CreaM_datos() const
{
    std::fstream archivo_real;
    std::fstream archivo_imag;

    archivo_real.open("data_real.txt", std::fstream::out);
    archivo_imag.open("data_imag.txt", std::fstream::out);

    std::vector<std::vector<std::complex<double>>> matrix(Simula());

    // body here
    for(int i=0; i<matrix[0].size(); ++i)
    {
        for(int j=0 ; j< matrix.size() ; ++j)
        {
          archivo_real<< std::real(matrix[j][i])<<" " ;
          archivo_imag<< std::imag(matrix[j][i])<<" " ;

        }
        archivo_real<<std::endl;
        archivo_imag<<std::endl;

    }


    archivo_real.close();
    archivo_imag.close();
}
