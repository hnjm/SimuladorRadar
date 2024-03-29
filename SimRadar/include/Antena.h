#ifndef ANTENA_H
#define ANTENA_H

#include <vector>
const double PI = 3.141592653589793;
class Antena
{
    public:
        Antena(double Wc, double titac, double phic, double Dc, double Gc);
        virtual ~Antena();

        void Antena_gira(double Ts);
        void set_angles(int titas , int phis);
        std::vector<double> get_angles();
        double get_diamter();

    protected:

    private:
        double W;
        double tita,phi;
        double D;
        double G;
};

#endif // ANTENA_H
