#include "../include/Receptor.h"


Receptor::Receptor(int Nc, int Mc, double fsc, double tec, double tauc, double Fcc, double PRFc, double cc, Objetos Objc, Antena Antc): N(Nc), M(Mc), fs(fsc), te(tec), tau(tauc), Fc(Fcc), PRF(PRFc), c(cc), Obj(Objc), Ant(Antc)
{

}

Receptor::~Receptor()
{
    //dtor
}


double Receptor::Potencia_receptor(double t, double rk, double c)
{
    //if(t-2.0*rk/c >0 && t-2.0*rk/c< this->tau )
    if(t>=2.0*rk/c && t<=2.0*rk/c + this->tau )
        return 1.0;
    return 0;
}

std::vector<std::vector<std::complex<double>>> Receptor::Simula()
{
    std::vector<std::vector<std::complex<double>>> matrix (N, std::vector<std::complex<double>> (M) );


    std::cout<<"Las dimensiones de la matriz de datos son: "<< matrix.size()<<" x "<< matrix[0].size()<<std::endl;

    // tiempo de adquisicion relativo a cada ventana de adquisicion
    std::vector<double> t(M);
    for(int tt=0; tt<M ; ++tt)
    {
        t[tt] = te + tt/fs ;

    }
    int K = Obj.numbReflec(); // Por el momento un unico reflector
    std::cout<<"El numero de reflectores es : "<<K<<std::endl;

    std::vector<std::vector<std::complex<double>>> P (K, std::vector<std::complex<double>> (M) );
    //std::vector<std::complex<double>> Pdef (M,0);

    double auxPot;
    double Pot_r;
    for(int n=0; n<N; ++n)
    {
        Ant.Antena_gira(1.0/this->PRF);
       /* std::vector<double> angulo_anten(Ant.get_angles());
        std::cout<<angulo_anten[0]<<"  "<< angulo_anten[1] <<std::endl;*/
        for(int k=0 ; k<K ; ++k)
        {
            Obj.actualiza_pos(k, 1.0/PRF);

            double rk = Obj.get_dist(k);


            for(int i=0 ; i< M; ++i)
            {
                auxPot = Potencia_receptor(t[i],rk, this->c);
                if(auxPot != 0)
                {
                     Pot_r = Pot_recibida(k);
                     //std::cout<<Pot_r<<std::endl;
                }
                else
                {
                     Pot_r = 1.0;
                }
                auxPot *= Pot_r;

                P[k][i] = auxPot*signal_recibida_pulso_rectangular( rk,  this->Fc,  0,  this->c);
                //P[k][i] = std::complex<double>( auxPot*std::cos(2* PI * Fc* (t[i] - 2.0*rk/3e8)) ,auxPot* std::sin(2* PI * Fc* (t[i] - 2.0*rk/3e8)));


            }
        }

        //sumo todas las contribuciones de las potencias
        std::vector<std::complex<double>> Pdef (M,0);
        for(int i=0 ; i<M; ++i)
        {


            for(int k=0 ; k<K ; ++k)
            {
              Pdef[i] = Pdef[i] + P[k][i] ;

            }
        }

        //esto es solo una forma de mirar que se esta almacenando en cada ventana de recepcion.
        /*if(n==1)
        {
            imprimeVentana(P[0]);
            //imprimeVentana(Pdef);
        }*/

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

void Receptor::imprimeVentana(std::vector<std::complex<double>> ventanaRecepcion)
{
    std::fstream archivo_ventana_real, archivo_ventana_imag;
    archivo_ventana_real.open("ventana_real.txt", std::fstream::out);
    archivo_ventana_imag.open("ventana_imag.txt", std::fstream::out);

    for(int i=0; i<ventanaRecepcion.size() ; ++i)
    {
        archivo_ventana_real<< std::real(ventanaRecepcion[i])<<" ";
        archivo_ventana_imag<< std::imag(ventanaRecepcion[i])<<" ";

    }

    archivo_ventana_real.close();
    archivo_ventana_imag.close();
}

double Receptor::Pot_recibida(int n)
{
    //esta parte busca la potencia en funcion del angulo usando un patron de antena tipo sinc.
    std::vector<double> angles_reflector_n(Obj.get_angles(n));
    std::vector<double> angles_antena(Ant.get_angles());
    double R_reflector = Obj.get_dist(n);

    double P_radiacion = patron_sinc(angles_antena,angles_reflector_n);
    // por ahora solo pondre la dependecia con la distancia y con el patron de radiacion
    double Potencia_recibida = P_radiacion*P_radiacion/std::pow(R_reflector,4);

    //std::cout<<angles_reflector_n[0]<<"  "<<angles_reflector_n[1] <<std::endl;
    //std::cout<<Potencia_recibida<<std::endl;
    return Potencia_recibida;


}

double Receptor::patron_sinc(std::vector<double>angulosAntena, std::vector<double>angulosReflector)
{
    double lambda = this->c/this->Fc;
    double D = Ant.get_diamter();

    double Ant_tita = angulosAntena[0];
    double Ant_phi = angulosAntena[1];
    double reflec_tita = angulosReflector[0]; if(reflec_tita < 0) reflec_tita = 2.0*PI+ reflec_tita;
    double reflec_phi = angulosReflector[1];
    double diff_tita;
    double aux;
    double E_phi,E_tita;

    if(std::abs(reflec_tita-Ant_tita)>PI)
    {
        aux = std::abs(reflec_tita-Ant_tita)-PI;
        diff_tita = PI-aux;

    }
    else
    {
        diff_tita = std::abs(reflec_tita-Ant_tita);
    }





    if(diff_tita >=PI/2.0) return 0;

    if(diff_tita == 0 && reflec_phi-Ant_phi!=0 )
    {
         E_tita = 1.0 ;
         E_phi = std::sin(PI * (D/lambda)*std::sin(reflec_phi-Ant_phi) )/ (PI * (D/lambda)*std::sin(reflec_phi-Ant_phi)) ;
    }
    else if(diff_tita != 0 && reflec_phi-Ant_phi==0)
    {
         E_tita = std::sin(PI * (D/lambda)*std::sin(diff_tita) )/ (PI * (D/lambda)*std::sin(diff_tita)) ;
         E_phi = 1.0 ;


    }
    else if(diff_tita == 0 && reflec_phi-Ant_phi==0)
    {
         E_tita = 1.0 ;
         E_phi = 1.0 ;
    }
    else
    {
         E_tita = std::sin(PI * (D/lambda)*std::sin(diff_tita) )/ (PI * (D/lambda)*std::sin(diff_tita)) ;
         E_phi = std::sin(PI * (D/lambda)*std::sin(reflec_phi-Ant_phi) )/ (PI * (D/lambda)*std::sin(reflec_phi-Ant_phi)) ;
    }
    //std::cout<<reflec_tita-Ant_tita<<std::endl;
    //std::cout<<reflec_phi-Ant_phi<<std::endl;
    //std::cout<<E_phi<<std::endl;
    return E_tita*E_phi;


}


std::complex<double> Receptor::signal_recibida_pulso_rectangular(double rk, double Fc, double fase_0, double c)
{
    //señal en banda base
   return std::complex<double>( std::cos(2.0* PI * Fc* ( - 2.0*rk/c) + fase_0) , std::sin(2.0* PI * Fc* ( - 2.0*rk/c)+ fase_0));
}

