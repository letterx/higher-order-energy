#ifndef _Y_LINEAR_ENERGY_HPP_
#define _Y_LINEAR_ENERGY_HPP_

/*
 * y-linear.hpp
 *
 * Copyright 2012 Alexander Fix
 * See LICENSE.txt for license information
 *
 * A representation for higher order pseudo-boolean energy functions
 *
 * Immediately transforms cliques added into their quadratic equivalent
 */

#include <vector>
#include <list>
#include <boost/foreach.hpp>

template <typename R, int D>
class YLinearEnergy {
    public:
        typedef int VarId;
        typedef VarId NodeId;

        // Constructs empty YLinearEnergy with no variables or terms
        YLinearEnergy();

        // Adds variables to the YLinearEnergy. Variables must be added 
        // before any terms referencing them can be added
        VarId AddVar();
        VarId AddVars(int n);
        VarId NumVars() const { return _varCounter; }
        NodeId AddNode(int n = 1) { return AddVars(n); }

        // Adds a clique defined by an table of energies, from 0 to 1 << k
        // where k is the size of the clique.
        void AddClique(const std::vector<VarId>& vars, const std::vector<R>& energy);

        void AddUnaryTerm(VarId v, R coeff);
        void AddUnaryTerm(VarId v, R E0, R E1) { AddUnaryTerm(v, E1 - E0); } 

        // Calculates attractive partitions
        void CalculateAttractivePartitions();
 
        // Utility - TODO(inrwherrmann): move this somewhere else
        int Count(uint32_t x, int n);
   
        // Reduces the YLinearEnergy to quadratic form
        // NOTE: THIS IS A DESTRUCTIVE OPERATION, so do not rely on the 
        // YLinearEnergy being in any useful state after this operation.
        //
        // This is a templated function, so it will work with any class 
        // implementing the necessary interface. See quadratic-rep.hpp for
        // the minimal requirements.
        template <typename QR>
        void ToQuadratic(QR &qr);

        template <typename QR>
        inline void AddReducedClique(const std::vector<VarId>& vars,
                const std::vector<R>& energyTable,
                QR& qr);
    private:
        struct Clique {
            std::vector<VarId> vars;
            std::vector<R> energyTable; // Keep all the 2^degree energy values of the clique. 
            // Note: energy is indexed so that energy[0] is all zero energy,
            // energy[1] is energy with x_0 = 1, energy[2] is with x_1 = 1, 
            // energy[3] is energy with x_0 = x_1 = 1, etc.
        };

        R _constantTerm;

        VarId _varCounter;

        typedef std::vector<Clique> CliqueVec;
        std::vector<R> _unaryTerms;
        CliqueVec _cliques;

        typedef std::vector<uint32_t> Partition;
        std::vector<Partition> _partitions;
};

template <typename R, int D>
inline YLinearEnergy<R, D>::YLinearEnergy()
    : _constantTerm(0), _varCounter(0), _cliques()
{
    CalculateAttractivePartition(D);
}

/**
* Calculate attractive partitions
*/
template <typename R, int D>
inline void YLinearEnergy<R, D>::CalculateAttractivePartition(int n) {
    typedef std::map<uint32_t, int> SubsetMap;
    std::vector<SubsetMap> subsetMaps(n+1);

   for (int i = 1; i <= n; ++i) {
        uint32_t max_subset = 1 << n;
        int subsets_size_i = 0;
        for (uint32_t subset = 0; subset < max_subset; ++subset) {
            int bits_set = 0;
            for (int b = 0; b < n; ++b) {
                if (subset & (1 << b)) {
                    bits_set++;
                }
            }
            if (bits_set == i) {
                subsetMaps[i][subset] = subsets_size_i;
                subsets_size_i++;
            }
        }
        subsetCounts[i] = subsets_size_i;
    }

    for (int i = 1; i < n; ++i) {
        auto& subsets = subsetMaps[i];
        auto& subsetsNext = subsetMaps[i+1];

        try {

        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);
        model.set(GRB_IntAttr_ModelSense, 1.0);

        std::vector<GRBVar> vars;
        // TODO(irwinherrmann): adds more variables than needed?
        // They're probably discarded in the minimization?
        for (const auto& p : subsets) {
            vars.push_back(model.addVar(0.0, 1.0, 1.0, GRB_BINARY));
        }
        model.update();
        for (const auto& p : subsetsNext) {
            GRBLinExpr constr{};
            auto subset = p.first;
            // Ensures constraint (3) is met
            for (int b = 0; b < n; ++b) {
                if (subset & (1 << b)) {
                    uint32_t neighbor = subset ^ (1 << b);
                    constr += vars[subsets[neighbor]];
                }
            }
            model.addConstr(constr, GRB_GREATER_EQUAL, 1.0);
        }
        model.update();
        model.optimize();
        partitionCounts[i] = model.get(GRB_DoubleAttr_ObjVal);

        // TODO(irwinherrmann): finish this part
        //_partitions[i].push_back(A_{i})

        } catch (GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch (...) {
            cout << "Exception during optimization" << endl;
            throw;
        }

    }
}


/**
* Term manipulation
*/

template <typename R, int D>
inline typename YLinearEnergy<R, D>::VarId 
YLinearEnergy<R, D>::AddVar() {
    _unaryTerms.push_back(0);
    return _varCounter++;
}

template <typename R, int D>
inline typename YLinearEnergy<R, D>::VarId 
YLinearEnergy<R, D>::AddVars(int n) {
    VarId firstVar = _varCounter;
    for (int i = 0; i < n; ++i)
        this->AddVar();
    return firstVar;


/**
* Clique manipulations
*/

template <typename R, int D>
inline void YLinearEnergy<R, D>::AddClique(const std::vector<VarId>& vars,
        const std::vector<R>& energyTable) {
    const unsigned int size = vars.size();
    const unsigned int numAssignments = 1 << size;
    assert(energyTable.size() == numAssignments);
    _cliques.push_back(Clique(vars, energyTable));
}



template <typename R, int D>
template <typename QR>
inline int YLinearEnergy<R, D>::AddReducedClique(const std::vector<VarId>& vars,
        const std::vector<R>& energyTable,
        QR& qr) {
    int n = vars.size();
    int numAssign = 1 << n;
    // Note: energyTable is indexed by the boolean assignment. Index always reads from left to right.
    // Ie. For vars = [0, 1], then energyTable[1] refers to var 1 being set and var 0 not being set.
    if (n == 1) {
        qr.AddUnaryTerm(vars[0], energyTable[0], energyTable[1]);
        return 0;
    } else if (n == 2) {
        qr.AddPairwiseTerm(vars[0], vars[1], energyTable[0], energyTable[1], energyTable[2], energyTable[3]);
        return 0;
    }

    // Start the y-linear reduction
    int coeff[n];

    // Indexed by k.
    int beta[n];
    // f(x)
    std::vector<R> sigma_x[n];
    std::vector<R> h_x[n]
    // f(a)
    R M_a[n][numAssign];    
    // f(a, x) linear in i, slope of each index
    std::vector<std::vector<R>> l_ai; // indexed by a, i

    // Initialize
    for (int i = 0; i < numAssign; i++) {
        sigma[0][i] = 0;
    }

    // TODO(irwinherrmann): abstract things
    for (int k = 1; k < vars.size(); k++) {

        // for all a in L_k
        for (uint32_t a : _partitions[k]) {
            // M_a calculation
            M_a[k][a] = 1;
            for (int i = 0; i < n; i++) {
                if ((a >> i) & 1 == 0) {
                    int neighbor = a | (1 << i);
                    M_a[k][a] += sigma_x[k-1][neighbor] - energyTable[neighbor];
                }
            }

            // l_ax calculation
            for (int i = 0; i < n; i++) {
                if ((a >> i) & 1) {
                    l_ai[a][i] += -M_a[k][a];
                } else {
                    int neighbor = a | (1 << i);
                    l_ai[a][i] -= (sigma_x[k-1][neighbor] - energyTable[neighbor]); 
                }
                //
            }
        }

        // h_x calculation
        for (uint32_t x = 0; x < numAssign; x++) {
            h_x[k][x] = sigma_x[k-1][x];
            for (uint32_t a : _partitions[k]) {
                R l_ax = 0;
                // TODO(irwinherrmann): abstract out
                for (int i = 0; i < n; i++) {
                    if (x & (1 << i)) {
                        l_ax += l_ai[a][i];
                    }
               }
               if (l_ax < 0) {
                    h_x[k][x] += l_ax;
               }
            }
        }
    
        // B_k calculation
        beta[k] = INT_MIN;
        for (uint32_t x = 0; x < numAssign; x++) {
            if (Count(x, n) == k + 2) {
                int value = energyTable[x] - h_x[k][x];
                if (value > max) {
                    beta[k] = value;
                }
            }
        }
        
        // sigma_x calculation
        for (uint32_t x = 0; x < numAssign; x++) {
            sigma[k][x] = h_x[k][x];
            if (Count(x, n) == k+2) {
                sigma[k][x] += beta[k];
            }
        }
    }

    // Add Terms 

    return 1;

}

inline int YLinearEnergy<R, D>::Count(uint32_t x, int n) {
    int count = 0;
    for (int i = 0; i < n; i++) {
        if ((x >> i) & 1) {
            count++;
        } 
    }
    return count;
}

template <typename R, int D>
template <typename QR>
inline void YLinearEnergy<R, D>::ToQuadratic(QR &qr) {
    for (const auto& c : _cliques)
        AddReducedClique(c.vars, c.energyTable, qr);
}


#endif
