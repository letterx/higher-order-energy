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

        // Adds a monomial to the YLinearEnergy. degree must be <= D
        // vars is an array of length at least degree, with the indices of 
        // the corresponding variables or literals in the monomial
        void AddTerm(R coeff, int degree, const VarId vars[]);

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

        void Clear();

    private:
        struct Term
        {
            R coeff;
            int degree;
            VarId vars[D];

            Term(R _coeff, int _degree, const VarId _vars[])
                : coeff(_coeff), degree(_degree)
            {
                for (int i = 0; i < degree; ++i)
                    vars[i] = _vars[i];
            }

            bool operator<(const Term& t) const;
            bool operator==(const Term& t) const;

            std::string ToString() const;

        };
        static int Compare(int d1, 
                const VarId vars1[], 
                int d2, 
                const VarId vars2[]);

        struct VarRecord {
            VarRecord(VarId id) 
                : _id(id), _positiveTerms(0), _higherOrderTerms(0), 
                _quadraticTerms(0), _sumDegrees(0), _terms(), _coeff(0) { }
            VarId _id;
            int _positiveTerms;
            int _higherOrderTerms;
            int _quadraticTerms;
            int _sumDegrees;
            std::list<Term> _terms;
            R _coeff;

            void PrintTerms() const;
        };

        R _constantTerm;

        void RemoveTerm(Term* tp);
        void _ReportMultilinearStats();

        size_t NumTerms() const {
            size_t numTerms = 0;
            BOOST_FOREACH(const VarRecord& vr, _varRecords)
                numTerms += vr._terms.size();
            return numTerms;
        }

        VarId _varCounter;

        typedef std::vector<VarRecord> varRecordVec_t;
        typedef std::vector<varRecordVec_t> cliqueRecord_t;
        cliqueRecord_t _cliqueRecords;

        typedef std::map<uint32_t, int> SubsetMap;
        std::vector<SubsetMap> subsetMaps;
};

template <typename R, int D>
inline YLinearEnergy<R, D>::YLinearEnergy()
    : _constantTerm(0), _varCounter(0), _varRecords()
{ }

/**
* Term manipulation
*/

template <typename R, int D>
inline typename YLinearEnergy<R, D>::VarId 
YLinearEnergy<R, D>::AddVar() {
    VarRecord vr(_varCounter);
    _varRecords.push_back(vr);
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
inline void 
YLinearEnergy<R, D>::AddTerm(int cliqueIndex, R coeff, int d, const VarId vars[]) {
    // TODO(irwinherrmann): change from class/global to handle each clique

    if 

    if(coeff == 0) {
        return;
    } else if (d == 0) {
        _constantTerm += coeff;
        return;
    } else if (d == 1) {
        _varRecords[vars[0]]._coeff += coeff;
        return;
    } else {
        VarRecord& smallestVarRec = _varRecords[vars[0]];
        typename std::list<Term>::iterator it = smallestVarRec._terms.begin();
        int compareVars = 1;
        while (it != smallestVarRec._terms.end()) {
            compareVars = Compare(d, vars, it->degree, it->vars); 
            if (compareVars == 0) {
                break;
            } else if (compareVars < 0) {
                break;
            } else {
                ++it;
            }
        }
        if (compareVars == 0) {
            it->coeff += coeff;
        } else {
            if (d > 2) {
                smallestVarRec._higherOrderTerms++;
                smallestVarRec._sumDegrees += d;
            } else {
                smallestVarRec._quadraticTerms++;
            }
            smallestVarRec._terms.insert(it, Term(coeff, d, vars));
        }
        if (coeff > 0)
            smallestVarRec._positiveTerms++;
        return;
    }
}

template <typename R, int D>
inline void YLinearEnergy<R, D>::AddUnaryTerm(VarId var, R coeff) {
    _varRecords[var]._coeff += coeff;
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

    cliqueRecordsVec_t::iterator end = _cliqueRecords.end();

    R coeffs[numAssignments];
    for (unsigned int subset = 1; subset < numAssignments; ++subset) {
        coeffs[subset] = 0;
    }
    // For each boolean assignment, add that to the corresponding monomials
    // with the correct parity
    for (unsigned int assignment = 0;
            assignment < numAssignments;
            ++assignment)
    {
        const R energy = energyTable[assignment];
        for (unsigned int subset = 1;
                subset < numAssignments;
                ++subset)
        {
            if (assignment & ~subset) {
                continue;
            } else {
                int parity = 0;
                for (unsigned int b = 0; b < size; ++b) {
                    parity ^=
                        (((assignment ^ subset) & (1 << b)) != 0);
                }
                coeffs[subset] += parity ? -energy : energy;
            }
        }
    }
    VarId subsetVars[D];
    for (unsigned int subset = 1; subset < numAssignments; ++subset) {
        int degree = 0;
        for (unsigned int b = 0; b < size; ++b) {
            if (subset & (1 << b)) {
                subsetVars[degree++] = vars[size-1-b];
            }
        }
        std::sort(subsetVars, subsetVars+degree);
        // TODO(irwinherrmann): make local
        AddTerm(end, coeffs[subset], degree, subsetVars);
    }

    AddReducedClique(coeff, energyTable);
}

template <typename R, int D>
template <typename QR>
inline void YLinearEnergy<R, D>::AddReducedClique(const std::vector<VarId>& vars,
        const std::vector<R>& energyTable) {
    // TODO(irwinherrmann): finish this

    // find attractive partions
    
    std::vector<double> sigma;
    sigma.push_back(energyTable[

    

    // iteratively calculate terms

    // add terms to class

}



template <typename R, int D>
inline int YLinearEnergy<R, D>::Compare(int d1, const VarId vars1[], int d2, const VarId vars2[]) {
    if (d1 < d2)
        return -1;
    if (d1 > d2)
        return 1;
    for (int index = 0; index < d1; ++index) {
        if (vars1[index] != vars2[index])
            return (vars1[index] < vars2[index]) ? -1 : 1;
    }
    return 0;
}

template <typename R, int D>
inline void YLinearEnergy<R, D>::Clear() {
    _varRecords.clear();
};

#endif
