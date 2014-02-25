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
#include "gurobi_c++.h"

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

        void AddUnaryTerm(VarId v, R coeff) { _unaryTerms[v] += coeff; }
        void AddUnaryTerm(VarId v, R E0, R E1) { AddUnaryTerm(v, E1 - E0); } 

        // Calculates attractive partitions
        void CalculateAttractivePartition(int n);
 
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
        inline int AddReducedClique(const std::vector<VarId>& vars,
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
        typedef std::map<uint32_t, std::vector<uint32_t>> DeltaMap;
        std::vector<DeltaMap> _deltas;
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
    std::vector<std::vector<uint32_t>> subsets;
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
                subsets[i].push_back(subset);
            }
        }
    }

    for (int i = 1; i < n; ++i) {
        auto& subsetsThis = subsetMaps[i];
        auto& subsetsNext = subsetMaps[i+1];

        try {

        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);
        model.set(GRB_IntAttr_ModelSense, 1.0);

        std::vector<GRBVar> vars;
        // TODO(irwinherrmann): adds more variables than needed?
        // They're probably discarded in the minimization?
        for (const auto& p  : subsetsThis) {
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
                    constr += vars[subsetsThis[neighbor]];
                }
            }
            model.addConstr(constr, GRB_GREATER_EQUAL, 1.0);
        }
        model.update();
        model.optimize();

        // Add to partition
        for (int j = 0; j < vars.size(); j++) {
            if (vars[j].get(GRB_DoubleAttr_X) == 1) {
                _partitions[i].push_back(j);
            }
        }
        // Calculate delta
        // Assumes it doesn't matter which A_{i+1} goes to which A_i as long as assigned correctly 
        // TODO(irwinherrmann): abstract
        for (uint32_t a_next : subsets[i+1]) {
            for (uint32_t a : _partitions[i]) {
                int diff = 0;
                for (int i = 0; i < n; i++) {
                    diff += ((a_next << i) & 1) ^ ((a << i) & 1);
                    if (diff > 1) {
                        break; // no need to keep looking
                    }
                }
                if (diff == 1) {
                    _deltas[i][a].push_back(a_next);
                    break; // been assigned don't need to find another one
                }
            } 
        }


        } catch (GRBException e) {
            std::cout << "Error code = " << e.getErrorCode() << std::endl;
            std::cout << e.getMessage() << std::endl;
        } catch (...) {
            std::cout << "Exception during optimization" << std::endl;
            throw;
        }

    }
}

/**
* Clique manipulations
*/

template <typename R, int D>
inline void YLinearEnergy<R, D>::AddClique(const std::vector<VarId>& vars,
        const std::vector<R>& energyTable) {
    const unsigned int size = vars.size();
    const unsigned int numAssignments = 1 << size;
    assert(energyTable.size() == numAssignments);
    _varCounter += vars.size();
    _cliques.push_back(Clique{vars, energyTable});
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
    std::vector<std::vector<R>> sigma_x(n);
    std::vector<std::vector<R>> h_x(n);
    // f(a)
    R M_a[n][numAssign];    
    // f(a, x) linear in i, slope of each index
    std::vector<std::vector<R>> l_ai; // indexed by a, i

    // Initialize
    for (int i = 0; i < numAssign; i++) {
        sigma_x[0][i] = 0;
    }

    // TODO(irwinherrmann): abstract things
    for (int k = 1; k < vars.size(); k++) {

        // for all a in L_k
        for (uint32_t a : _partitions[k]) {
            // M_a calculation
            M_a[k][a] = 1;
            std::vector<uint32_t> delta_a = _deltas[k][a];
            for (uint32_t ai : delta_a) {
                M_a[k][a] += sigma_x[k-1][ai] - energyTable[ai];
            }

            // l_ax calculation
            for (int i = 0; i < n; i++) {
                uint32_t ai = a | (1 << i);
                std::vector<uint32_t> delta_a = _deltas[k][a];
                if ((a >> i) & 1) {
                    l_ai[a][i] += -M_a[k][a];
                } else if (std::find(delta_a.begin(), delta_a.end(), ai) != delta_a.end()) {
                    l_ai[a][i] -= (sigma_x[k-1][ai] - energyTable[ai]); 
                } else {
                    l_ai[a][i] += M_a[k][a];
                }
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
                if (value > beta[k]) {
                    beta[k] = value;
                }
            }
        }
        
        // sigma_x calculation
        for (uint32_t x = 0; x < numAssign; x++) {
            sigma_x[k][x] = h_x[k][x];
            if (Count(x, n) == k+2) {
                sigma_x[k][x] += beta[k];
            }
        }
    }

    // Quadralization of symmetric term (sums of beta_k values)
    // g(0) = g(1) = 0 in this case
    int C = 0;
    for (int k = 0; k < n; k++) {
        C += beta[k]; 
    }

    int C_i[n];
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            int c_ik = 0;
            if (i == k) {
                c_ik = 3;
            } else if (i == k + 1 || i == k - 1) {
                c_ik = 0;
            } else {
                c_ik = 1;
            }
            C_i[i] += beta[k] * c_ik;
        }
    } 

    for (int j = 0; j < n; j++) {
        for (int i = 0; i < j; i++) {
            qr.AddPairwiseTerm(vars[i], vars[j], 0, 0, 0, C);
        }
    }

    // TODO(irwinherrmann): check! esp pairwise term
    for (int i = 0; i < n - 1; i++) {
        VarId wi = ++_varCounter;
        qr.AddUnaryTerm(wi, 0, C_i[i]*i);
        for (int j = 0; j < n; j++) {
            qr.AddPairwiseTerm(wi, vars[j], 0, 0, 0, C_i[i]);
        } 
    }

    // Quadralization of min {l_a(x), 0} term. n+1 is new variable
    VarId y = ++_varCounter; // new variable

    for (int i = 0; i < n; i++) {
        R xi_coeff = 0;
        // all a\in A
        for (std::vector<uint32_t> partitions : _partitions) {
            for (uint32_t a : partitions) {
                xi_coeff += l_ai[a][i];
            }
        }
        qr.AddPairwiseTerm(y, vars[i], 0, 0, 0, xi_coeff);
    }
    // no unary term due to construction of l_ai

    return 1;

}

template <typename R, int D>
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
