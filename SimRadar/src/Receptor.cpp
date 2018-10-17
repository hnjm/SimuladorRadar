#include "../include/Receptor.h"


Receptor::Receptor(int Nc, int Mc, double fsc, double tec, double tauc, double Fcc, Objetos Objc): N(Nc), M(Mc), fs(fsc), te(tec), tau(tauc), Fc(Fcc), Obj(Objc)
{

}

Receptor::~Receptor()
{
    //dtor
}


double Receptor::Potencia_receptor(double t, double rk, double c)
{
    if(t-2.0*rk/c >0 && t-2.0*rk/c< this->tau )
        return 1.0;
    return 0;
}




std::vector<std::vector<std::complex<double>>> Receptor::Simula()
{
    std::vector<std::vector<std::complex<double>>> matrix (N, std::vector<std::complex<double>> (M) );

    std::cout<<"Las dimensiones de la matriz de datos son: "<< matrix.size()<<" x "<< matrix[0].size()<<std::endl;

   /* for(int i=0; i<matrix[0].size(); ++i)
    {
        for(int j=0 ; j< matrix.size() ; ++j)
        {
            matrix[j][i] = std::complex<double>(j+i, 1 );

        }

    }*/

    // tiempo de adquisicion relativo a cada ventana de adquisicion
    std::vector<double> t(M);
    for(int tt=0; tt<M ; ++tt)
    {
        t[tt] = te + tt/fs ;

    }



    int K = Obj.numbReflec(); // Por el momento un unico reflector
    std::vector<std::vector<std::complex<double>>> P (K, std::vector<std::complex<double>> (M) );
    std::vector<std::complex<double>> Pdef (M,0);
    for(int n=0; n<N; ++n)
    {


        for(int k=0 ; k<K ; ++k)
        {
            double rk = Obj.get_dist(k);

            for(int i=0 ; i< M; ++i)
            {
                double auxPot = Potencia_receptor(t[i],rk, 3e8);
                P[k][i] = std::complex<double>( auxPot*std::cos(2* PI * Fc* (t[i] - 2.0*rk/3e8)) ,auxPot* std::sin(2* PI * Fc* (t[i] - 2.0*rk/3e8)));

            }
        }


        //sumo todas las contribuciones de las potencias

        for(int i=0 ; i<M; ++i)
        {


            for(int k=0 ; k<K ; ++k)
            {
              Pdef[i] = +P[k][i] ;

            }
        }

        //Armo el vector complejo de la ventana con la fase del retardo de la onda y las fases correspondientes
        for(int i = 0 ; i< M ; ++i)
        {
            matrix[n][i] = Pdef[i];

        }




    }



    return matrix;
}




void Receptor::CreaM_datos()
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
