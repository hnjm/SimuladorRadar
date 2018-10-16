#ifndef ANTENA_H
#define ANTENA_H


class Antena
{
    public:
        Antena(double Wc, double titac, double phic, double Dc, double Gc);
        virtual ~Antena();

    protected:

    private:
        double W;
        double tita,phi;
        double D;
        double G;
};

#endif // ANTENA_H
