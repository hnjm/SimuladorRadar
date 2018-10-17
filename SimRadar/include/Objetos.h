#ifndef OBJETOS_H
#define OBJETOS_H
#include <vector>
#include <cmath>

class Objetos
{
    public:
        Objetos();
        Objetos(std::vector<std::vector<double>> reflectorc, std::vector<double> src, std::vector<std::vector<double>> vel);
        virtual ~Objetos();

        //get and set
        std::vector<double> get_pos(int n);
        void set_pos(std::vector<double>, int n);

        double get_dist(int n);
        int numbReflec();

    protected:

    private:
        std::vector<std::vector<double>> reflector;
        std::vector<double> src;
        std::vector<std::vector<double>> vel;
};

#endif // OBJETOS_H
