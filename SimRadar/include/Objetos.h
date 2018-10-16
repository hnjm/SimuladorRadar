#ifndef OBJETOS_H
#define OBJETOS_H
#include <vector>

class Objetos
{
    public:
        Objetos();
        Objetos(std::vector<int> reflectorc);
        virtual ~Objetos();

        //get and set
        std::vector<int> get();
        void set(std::vector<int>);

    protected:

    private:
        std::vector<int> reflector;
};

#endif // OBJETOS_H
