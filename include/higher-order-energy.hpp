#ifndef _HIGHER_ORDER_ENERGY_HPP_
#define _HIGHER_ORDER_ENERGY_HPP_

/*
 * higher-order-energy.hpp
 *
 * Copyright 2012 Alexander Fix
 * See LICENSE.txt for license information
 *
 * A representation for higher order pseudo-boolean energy functions
 *
 * Built up by calls to AddVar and AddTerm
 * Main operation is the reduction to quadratic form, ToQuadratic()
 * ToQuadratic() impements the reduction described in 
 *  A Graph Cut Algorithm for Higher Order Markov Random Fields
 *  Alexander Fix, Artinan Gruber, Endre Boros, Ramin Zabih, ICCV 2011
 *
 * Template parameters:
 *  Energy   The result type of the energy function
 *  D        The maximum degree of any monomial 
 *
 * Example usage:
 *
 * To represent the function
 *      f(x_0, ..., x_3) = 7x_0 + 9x_1 - 8x_2 + 13x_3 + 27x_1x_2x_3 
 *                       - 31x_1x_3x_4 + 18x_1x_2x_3x_4
 * we can do the following:
 *
 * HigherOrderEnergy<int, 4> f;
 * f.AddVars(4); 
 * f.AddUnaryTerm(0, 7);
 * f.AddUnaryTerm(1, 9);
 * f.AddUnaryTerm(2, -8);
 * f.AddUnaryTerm(3, 13);
 *
 * VarId term1[] = {1, 2, 3};
 * VarId term2[] = {1, 3, 4};
 * VarId term3[] = {1, 2, 3, 4};
 * f.AddTerm(27, 3, term1);
 * f.AddTerm(-31, 3, term2);
 * f.AddTerm(18, 4, term3);
 *
 *
 * Then, we can get an equivalent quadratic form (by adding additional 
 * variables and performing transformations) by calling ToQuadratic. We can
 * solve this using QPBO, or any other quadratic optimizer that has the same
 * interface.
 *
 * QPBO<int> qr;
 * f.ToQuadratic(qr);
 * qr.Solve();
 *
 *
 */

#include <vector>
#include <set>
#include <unordered_map>
#include <list>
#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>

template <typename R, int D>
class HigherOrderEnergy {
    public:
        typedef int VarId;
        typedef VarId NodeId;
        typedef std::set<VarId> VarIdSet_t;

        // Constructs empty HigherOrderEnergy with no variables or terms
        HigherOrderEnergy();

        // Adds variables to the HigherOrderEnergy. Variables must be added 
        // before any terms referencing them can be added
        VarId AddVar();
        VarId AddVars(int n);
        VarId NumVars() const { return _varCounter; }
        NodeId AddNode(int n = 1) { return AddVars(n); }

        // Adds a monomial to the HigherOrderEnergy. degree must be <= D
        // vars is an array of length at least degree, with the indices of 
        // the corresponding variables or literals in the monomial
        void AddTerm(R coeff, int degree, const VarId vars[]);

        // Adds a clique defined by an table of energies, from 0 to 1 << k
        // where k is the size of the clique.
        void AddClique(const std::vector<VarId>& vars, const std::vector<R>& energy);

        void AddUnaryTerm(VarId v, R coeff);
        void AddUnaryTerm(VarId v, R E0, R E1) { AddUnaryTerm(v, E1 - E0); } 
    
        // Reduces the HigherOrderEnergy to quadratic form
        // NOTE: THIS IS A DESTRUCTIVE OPERATION, so do not rely on the 
        // HigherOrderEnergy being in any useful state after this operation.
        //
        // This is a templated function, so it will work with any class 
        // implementing the necessary interface. See quadratic-rep.hpp for
        // the minimal requirements.
        template <typename QR>
        void ToQuadratic(QR &qr);

        void SetWidth(int n);
        void SetHeight(int n);

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

        struct CoverRecord {
            CoverRecord(VarId id, VarIdSet_t set)
                : _id(id), _bkValue(0), _cover(set), _size(set.size()) { }

            VarId _id;
            double _bkValue;
            VarIdSet_t _cover;
            int _size;
        };

        struct VectorHash {
            std::size_t operator()(VarIdSet_t const& c) const {
                return boost::hash_range(c.begin(), c.end());
            }
        };

        R _constantTerm;

        void RemoveTerm(Term* tp);
        template <typename QR>
        void _SetupPairwiseCover(QR& qr);
        void _CalculateBkValues();
        template <typename QR>
        void _EliminateTerms(QR& qr);
        void _ReportMultilinearStats();
        void AddBkValue(VarIdSet_t set, int bkValue);
        CoverRecord GetCoverRecord(VarIdSet_t set);

        size_t NumTerms() const {
            size_t numTerms = 0;
            BOOST_FOREACH(const VarRecord& vr, _varRecords)
                numTerms += vr._terms.size();
            return numTerms;
        }

        VarId _varCounter;
        int _width;
        int _height;

        typedef std::vector<VarRecord> VarRecordVec_t;
        VarRecordVec_t _varRecords;
        std::vector<std::list<Term> > _termsByDegree;
        typedef std::unordered_map<VarIdSet_t, CoverRecord, VectorHash> Cover_t;
        Cover_t _coverRecords;
};

template <typename R, int D>
inline HigherOrderEnergy<R, D>::HigherOrderEnergy()
    : _constantTerm(0), _varCounter(0), _varRecords(), _termsByDegree(), _coverRecords() { 
    _termsByDegree.resize(D+1);
}

template <typename R, int D>
inline typename HigherOrderEnergy<R, D>::VarId 
HigherOrderEnergy<R, D>::AddVar() {
    VarRecord vr(_varCounter);
    _varRecords.push_back(vr);
    return _varCounter++;
}

template <typename R, int D>
inline typename HigherOrderEnergy<R, D>::VarId 
HigherOrderEnergy<R, D>::AddVars(int n) {
    VarId firstVar = _varCounter;
    for (int i = 0; i < n; ++i)
        this->AddVar();
    return firstVar;
}

template <typename R, int D>
inline void
HigherOrderEnergy<R, D>::SetWidth(int n) {
    _width = n;
}

template <typename R, int D>
inline void
HigherOrderEnergy<R, D>::SetHeight(int n) {
    _height = n;
}

template <typename R, int D>
inline void 
HigherOrderEnergy<R, D>::AddTerm(R coeff, int d, const VarId vars[]) {
    if(coeff == 0) {
        return;
    } else if (d == 0) {
        _constantTerm += coeff;
        return;
    } else if (d == 1) {
        _varRecords[vars[0]]._coeff += coeff;
    } else {
        VarRecord& smallestVarRec = _varRecords[vars[0]];
        typename std::list<Term>::iterator it = smallestVarRec._terms.begin();
        int compareVars = 1;
        // TODO(irwinherrmann) - abstract this out
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
    }

    _termsByDegree[d].push_back(Term(coeff, d, vars));
    return;
}

template <typename R, int D>
inline void HigherOrderEnergy<R, D>::AddClique(const std::vector<VarId>& vars,
        const std::vector<R>& energyTable) {
    const unsigned int size = vars.size();
    const unsigned int numAssignments = 1 << size;
    assert(energyTable.size() == numAssignments);

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
                subsetVars[degree++] = vars[b];
            }
        }
        std::sort(subsetVars, subsetVars+degree);
        AddTerm(coeffs[subset], degree, subsetVars);
    }
}

template <typename R, int D>
inline void HigherOrderEnergy<R, D>::AddUnaryTerm(VarId var, R coeff) {
    _varRecords[var]._coeff += coeff;
}


template<typename R, int D>
template<typename QR>
void HigherOrderEnergy<R, D>::_SetupPairwiseCover(QR& qr) {
    // TODO(irwinherrmann) more abstracted in the future. also consider 3x3 case.
    // 2x1 vertical block
    for (int varIndex = 0; varIndex < _width * (_height - 1); ++varIndex) {
        VarId newVar = AddVar();
        VarIdSet_t set;
        VarRecord& vr = _varRecords[varIndex];
        set.insert(vr._id);
        set.insert(vr._id + _width);
        CoverRecord cr(newVar, set);
        _coverRecords.insert(std::make_pair(set, cr));
    }

    // 1x1 unary block
    for (int varIndex = 0; varIndex < _width * _height; ++varIndex) {
        VarIdSet_t set;
        VarRecord& vr = _varRecords[varIndex];
        set.insert(vr._id);
        CoverRecord cr(vr._id, set);
        _coverRecords.insert(std::make_pair(set, cr));
    }
}

template <typename R, int D>
void HigherOrderEnergy<R, D>::_CalculateBkValues() {
    // TODO(irwinherrmann) more abstracted
    int degree = D;
    while (degree > 0) {
        typename std::list<Term> terms = _termsByDegree[degree];
        typename std::list<Term>::iterator termIt = terms.begin();
        while (termIt != terms.end()) {
            Term& t = *termIt;
            VarIdSet_t A;
            VarIdSet_t B;
            VarIdSet_t H;
            A.insert(t.vars[0]);
            H.insert(t.vars[0]);
            for (int i = 1; i < t.degree; ++i) {
                if ((t.vars[0] - t.vars[i]) % _width == 0) {
                    A.insert(t.vars[i]);
                } else {
                    B.insert(t.vars[i]);
                }
                H.insert(t.vars[i]);
            }
 
            double bK = 2 * std::abs(t.coeff);

            typename Cover_t::iterator coverIt = _coverRecords.find(H);
            if (coverIt != _coverRecords.end()) {
                bK += coverIt->second._bkValue;
            }

            AddBkValue(A, bK);
            AddBkValue(B, bK);

            ++termIt;
        }
        --degree;
    }
}

template <typename R, int D>
void HigherOrderEnergy<R, D>::AddBkValue(VarIdSet_t set, int bkValue) {
    if (set.empty()) {
        return;
    }

    typename Cover_t::iterator coverIt = _coverRecords.find(set);
    if (coverIt != _coverRecords.end()) {
        coverIt->second._bkValue += bkValue;
    } else {
        exit(1); // Error
    }
}

template <typename R, int D>
inline typename HigherOrderEnergy<R, D>::CoverRecord HigherOrderEnergy<R, D>::GetCoverRecord(VarIdSet_t set) {
    typename Cover_t::iterator coverIt = _coverRecords.find(set);
    if (coverIt != _coverRecords.end()) {
        return coverIt->second;
    } else {
        exit(1); // Error
    }
}

// specifically written for 2x2 clique size
template <typename R, int D>
template <typename QR>
void HigherOrderEnergy<R, D>::_EliminateTerms(QR& qr) {
    
    int expectedVars = _varCounter;
    int expectedEdges = 0;
    BOOST_FOREACH(const VarRecord& vr, _varRecords) {
        expectedEdges += vr._quadraticTerms;
    }

    qr.SetMaxEdgeNum(expectedEdges);
    qr.AddNode(expectedVars);

    size_t numPixels = _width * _height;
    // all H in curlyH
    for (size_t varIndex = 0; varIndex < numPixels; ++varIndex) {
        VarRecord& vr = _varRecords[varIndex];
        typename std::list<Term>::iterator termIt = vr._terms.begin();
        while (termIt != vr._terms.end()) {
            Term& t = *termIt;
            typename std::list<Term>::iterator currIt = termIt;
            ++termIt;

            if (t.degree <= 2) {
                qr.AddPairwiseTerm(t.vars[0], t.vars[1], 0, 0, 0, t.coeff);
            } else {
                VarIdSet_t s;
                for (int i = 0; i < t.degree; ++i) {
                    s.insert(t.vars[i]);
                }

                VarId newVars[2];
                VarIdSet_t A;
                VarIdSet_t B;
                A.insert(t.vars[0]);
                for (int i = 1; i < t.degree; ++i) {
                    if ((t.vars[0] - t.vars[i]) % _width == 0) {
                        A.insert(t.vars[i]);
                    } else {
                        B.insert(t.vars[i]);
                    }
                }
                
                newVars[0] = GetCoverRecord(A)._id;
                newVars[1] = GetCoverRecord(B)._id;

                qr.AddPairwiseTerm(newVars[0], newVars[1], 0, 0 , 0, t.coeff);
            }
            vr._terms.erase(currIt);
        }
    }

    typename Cover_t::iterator coverIt = _coverRecords.begin();
    while (coverIt != _coverRecords.end()) {
        CoverRecord cover = coverIt -> second;
        ++coverIt;

        if (cover._size < 2) {
            continue;
        }

        int newCoeff = cover._bkValue * (cover._size - .5);
        _varRecords[cover._id]._coeff = newCoeff;

        typename VarIdSet_t::iterator coverElemsIt = cover._cover.begin();
        while (coverElemsIt != cover._cover.end()) {
            VarId var = *coverElemsIt;
            qr.AddPairwiseTerm(var, cover._id, 0, 0, 0, -cover._bkValue);
            ++coverElemsIt;
        }

        // only for 2x2 case
        qr.AddPairwiseTerm(*cover._cover.begin(), *(++cover._cover.begin()), 0, 0, 0, .5*cover._bkValue);
    }

    BOOST_FOREACH(VarRecord& vr, _varRecords) {
        qr.AddUnaryTerm(vr._id, 0, vr._coeff);
    }
}

template <typename R, int D>
template <typename QR>
inline void HigherOrderEnergy<R, D>::ToQuadratic(QR& qr) {
    _SetupPairwiseCover(qr);
    _CalculateBkValues();
    _EliminateTerms(qr);
}

template <typename R, int D>
inline int HigherOrderEnergy<R, D>::Compare(int d1, const VarId vars1[], int d2, const VarId vars2[]) {
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


/**
 * Takes in a vector of node locations and returns a pattern which given a location uniquely identifies
 * the pattern of nodes. The pattern is defined as follows. Start at a given location and iterate
 * through each node in the figure. If the node is included in the list, use a 1. If the node is not
 * included in the list, use a 0. If the remaining nodes on a given row are 0, then remove the sequence
 * of 0s and write a 2. Eg. [0, width+1] -> 121. [0, width+2, width+4] -> 120101
 * Function assumes that the vector is sorted.
 *
template <typename R, int D>
inline string HigherOrderEnergy<R, D>::ToPattern(std::vector<int> nodes) {
    std::string pattern = "";
    int prev = ToLocation(nodes);
    for (int& i : nodes) {
        int dist = i - prev;
        if (i / width > prev / width) {
            pattern += "2";
            dist = i % width;
        }
        for (int j = 0; j < dist; ++j)
            pattern += "0";
        pattern += "1";
        prev = i;
    }
    return atoi(pattern.c_str());
}

template <typename R, int D>
inline int HigherOrderEnergy<R, D>::ToLocation(std::vector<int> nodes) {
    int minRow = nodes[0] / width;
    int minCol = nodes[0] % width;
    for (int& i : nodes) {
        int row = i / width;
        int col = i % width;
        if (minRow < row) {
            minRow = row;
        }
        if (minCol < col) {
            minCol = col;
        } 
    }

    return minRow * width + minCol;
}
*/

template <typename R, int D>
inline void HigherOrderEnergy<R, D>::Clear() {
    _varRecords.clear();
};

#endif
