#include "higher-order.hpp"
#include "pairwise-cover.hpp"

int main(int argc, char** argv) {
    typedef int REAL;
    PairwiseCover<REAL, 4> pc;
    pc.AddVars(4);
    pc.AddUnaryTerm(0, 5);
    pc.AddUnaryTerm(1, -3);
    pc.AddUnaryTerm(2, 1);
    pc.AddUnaryTerm(3, 15);
    std::vector<int> vars = {0, 1, 2, 3};
    std::vector<REAL> energyTable 
        = {0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 0, 0,
           0, 0, 1, 1};
    pc.AddClique(vars, energyTable);


    QPBO<REAL> qpbo(0, 0);
    pc.ToQuadratic(qpbo);
    qpbo.Save("qpbo.dat");
    return 0;
}
