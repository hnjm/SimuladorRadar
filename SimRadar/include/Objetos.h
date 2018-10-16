#ifndef OBJETOS_H
#define OBJETOS_H
#include <vector>

class Objetos
{
    public:
        Objetos();
        Objetos(std::vector<double> reflectorc, double src, std::vector<double> vel);
        virtual ~Objetos();

        //get and set
        std::vector<double> get_pos();
        void set_pos(std::vector<double>);

    protected:

    private:
        std::vector<double> reflector;
        double src;
        std::vector<double> vel;
};

#endif // OBJETOS_H
