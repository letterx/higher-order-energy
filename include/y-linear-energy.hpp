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
};

template <typename R, int D>
inline YLinearEnergy<R, D>::YLinearEnergy()
    : _constantTerm(0), _varCounter(0), _cliques()
{ }

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
}

template <typename R, int D>
inline void YLinearEnergy<R, D>::AddUnaryTerm(VarId var, R coeff) {
    _unaryTerms[var] += coeff;
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
    _cliques.push_back(Clique{vars, energyTable});
}

template <typename R, int D>
template <typename QR>
inline void YLinearEnergy<R, D>::AddReducedClique(const std::vector<VarId>& vars,
        const std::vector<R>& energyTable,
        QR& qr) {
    // TODO(irwinherrmann): finish this

    // find attractive partions
    
    /* Commenting out to make compile
    std::vector<double> sigma;
    sigma.push_back(energyTable[

    

    // iteratively calculate terms

    // add terms to class

    */
}

template <typename R, int D>
template <typename QR>
inline void YLinearEnergy<R, D>::ToQuadratic(QR &qr) {
    for (const auto& c : _cliques)
        AddReducedClique(c.vars, c.energyTable, qr);
}


#endif
