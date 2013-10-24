#include "higher-order.hpp"
#include "pairwise-cover.hpp"

int main(int argc, char** argv) {
    typedef int REAL;
    {
        PairwiseCover<REAL, 4> pc;
        pc.AddVars(4);
        std::vector<int> vars = {0, 2, 1, 3};
        std::vector<REAL> energyTable 
            = {2414, 6437, 3801, 5912,
               5227, 7887, 4934, 7496,
               1840, 6194, 2689, 5803,
               5521, 7523, 5347, 7149};
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
            = {2414, 6437, 5227, 7887,
               3801, 5912, 4934, 7496,
               1840, 6194, 5521, 7523,
               2689, 5803, 5347, 7149};
        pc.AddClique(vars, energyTable);


        QPBO<REAL> qpbo(0, 0);
        pc.ToQuadratic(qpbo);
        qpbo.Save("qpbo-grid.dat");

    }
    return 0;
}
