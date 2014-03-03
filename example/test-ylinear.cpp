#include "higher-order.hpp"
#include "y-linear-energy.hpp"

int main(int argc, char** argv) {
    typedef int REAL;
    {
        YLinearEnergy<REAL, 4> pc;
        pc.AddVars(4);
        std::vector<int> vars = {0, 1, 2, 3};
        std::vector<REAL> energyTable 
            = {1, 2, 3, 4,
               5, 6, 7, 8,
               9, 10, 11, 12,
               13, 14, 15, 16};
        pc.AddClique(vars, energyTable);


        QPBO<REAL> qpbo(0, 0);
        pc.ToQuadratic(qpbo);
        qpbo.Save("qpbo.dat");
    }
    {
        PairwiseCoverGrid<REAL, 4> pc;
        pc.AddVars(4);
        pc.SetWidth(2);
        pc.SetHeight(2);
        std::vector<int> vars = {0, 1, 2, 3};
        std::vector<REAL> energyTable 
            = {1, 2, 3, 4,
               5, 6, 7, 8,
               9, 10, 11, 12,
               13, 14, 15, 16};
        pc.AddClique(vars, energyTable);


        QPBO<REAL> qpbo(0, 0);
        pc.ToQuadratic(qpbo);
        qpbo.Save("qpbo-grid.dat");

    }
    return 0;
}
