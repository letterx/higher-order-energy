#ifndef _HIGHER_ORDER_ENERGY_GRAPH_HPP_
#define _HIGHER_ORDER_ENERGY_GRAPH_HPP_

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
 * HigherOrderEnergyGraph<int, 4> f;
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

#include <list>
#include <memory>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/heap/binomial_heap.hpp>
#include <boost/intrusive/set.hpp>

template <typename R, int D>
class HigherOrderEnergyGraph {
    public:
        typedef int VarId;
        typedef VarId NodeId;

        // Constructs empty HigherOrderEnergyGraph with no variables or terms
        HigherOrderEnergyGraph();

        // Adds variables to the HigherOrderEnergyGraph. Variables must be added 
        // before any terms referencing them can be added
        VarId AddVar();
        VarId AddVars(int n);
        VarId NumVars() const { return _varCounter; }
        NodeId AddNode(int n = 1) { return AddVars(n); }

        // Adds a monomial to the HigherOrderEnergyGraph. degree must be <= D
        // vars is an array of length at least degree, with the indices of 
        // the corresponding variables or literals in the monomial
        void AddTerm(R coeff, int degree, const VarId vars[]);

        // Adds a clique defined by an table of energies, from 0 to 1 << k
        // where k is the size of the clique.
        void AddClique(const std::vector<VarId>& vars, const std::vector<R>& energy);

        void AddUnaryTerm(VarId v, R coeff);
        void AddUnaryTerm(VarId v, R E0, R E1) { AddUnaryTerm(v, E1 - E0); } 
    
        // Reduces the HigherOrderEnergyGraph to quadratic form
        // NOTE: THIS IS A DESTRUCTIVE OPERATION, so do not rely on the 
        // HigherOrderEnergyGraph being in any useful state after this operation.
        //
        // This is a templated function, so it will work with any class 
        // implementing the necessary interface. See quadratic-rep.hpp for
        // the minimal requirements.
        template <typename QR>
        void ToQuadratic(QR &qr);

        void Clear();

    private:
        typedef boost::intrusive::set_base_hook<
            boost::intrusive::link_mode<boost::intrusive::normal_link>>
            TermSetHook;
        struct Term : public TermSetHook
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

            std::string ToString() const;

        };
        typedef std::unique_ptr<Term> TermPtr;
        struct TermComp {
            bool operator()(const Term& t1, const Term& t2) const {
                if (t1.degree != t2.degree) {
                    return t1.degree < t2.degree;
                } else {
                    for (int i = 0; i < t1.degree; ++i) {
                        if (t1.vars[i] != t2.vars[i])
                            return t1.vars[i] < t2.vars[i];
                    }
                    return false;
                }
            }
        };

        typedef boost::intrusive::set<
            Term,
            boost::intrusive::base_hook<TermSetHook>,
            boost::intrusive::compare<TermComp>> TermSet;

        typedef std::pair<VarId, int> VarPriority;
        struct PriorityComp {
            bool operator()(const VarPriority& p1, const VarPriority& p2) const {
                return p1.second < p2.second;
            }
        };
        typedef boost::heap::binomial_heap<
                VarPriority, 
                boost::heap::compare<PriorityComp>> 
            PositiveTermsPQ;
        typedef typename PositiveTermsPQ::handle_type PQHandle;

        struct VarRecord {
            VarRecord(VarId id) 
                : _id(id)
                , _higherOrderTerms(0)
                , _sumDegrees(0)
                , _terms()
                , _coeff(0) 
            { }
            VarId _id;
            int _higherOrderTerms;
            int _sumDegrees;
            std::vector<int> _terms; // Indices of positive terms
            R _coeff;
            PQHandle _pqHandle;

            void PrintTerms() const;
        };
        typedef std::vector<VarRecord> varRecordVec_t;


        void RemoveTerm(Term* tp);
        void AddTerm(Term& t);
        void _EliminatePositiveTerms();
        template <typename QR>
        void _ReduceNegativeTerms(QR& qr);
        void _ReportMultilinearStats();

        size_t NumTerms() const {
            return _terms.size();
        }

        R _constantTerm;
        std::vector<TermPtr> _terms;
        VarId _varCounter;
        varRecordVec_t _varRecords;
        bool _usePriority = true;
        PositiveTermsPQ _positiveTermsPQ;
        TermSet _termSet;
};

template <typename R, int D>
inline HigherOrderEnergyGraph<R, D>::HigherOrderEnergyGraph()
    : _constantTerm(0)
    , _terms()
    , _varCounter(0)
    , _varRecords()
    , _positiveTermsPQ()
    , _termSet()
{ }


template <typename R, int D>
inline typename HigherOrderEnergyGraph<R, D>::VarId 
HigherOrderEnergyGraph<R, D>::AddVar() {
    VarRecord vr(_varCounter);
    _varRecords.push_back(vr);
    return _varCounter++;
}

template <typename R, int D>
inline typename HigherOrderEnergyGraph<R, D>::VarId 
HigherOrderEnergyGraph<R, D>::AddVars(int n) {
    VarId firstVar = _varCounter;
    for (int i = 0; i < n; ++i)
        this->AddVar();
    return firstVar;
}

template <typename R, int D>
inline void 
HigherOrderEnergyGraph<R, D>::AddTerm(R coeff, int d, const VarId vars[]) {
    if(coeff == 0) {
        return;
    } else if (d == 0) {
        _constantTerm += coeff;
        return;
    } else if (d == 1) {
        _varRecords[vars[0]]._coeff += coeff;
        return;
    } else {
        _terms.emplace_back(new Term{coeff, d, vars});
        auto p = _termSet.insert(*_terms.back());
        if (!p.second) { // Then term already exists
            Term& t = *p.first;
            t.coeff += coeff;
            _terms.pop_back();
        } else {
            int idx = _terms.size()-1;
            for (int i = 0; i < d; ++i) {
                auto& vr = _varRecords[vars[i]];
                vr._terms.push_back(idx);
                if (d > 2) {
                    vr._higherOrderTerms++;
                    vr._sumDegrees += d;
                }
            }
        }
        return;
    }
}

template <typename R, int D>
inline void 
HigherOrderEnergyGraph<R, D>::AddTerm(Term& t) {
    if (t.coeff == 0) {
        return;
    } else if (t.degree == 0) {
        _constantTerm += t.coeff;
        t.coeff = 0;
        return;
    } else if (t.degree == 1) {
        _varRecords[t.vars[0]]._coeff += t.coeff;
        t.coeff = 0;
        return;
    } else {
        auto p = _termSet.insert(t);
        if (!p.second) { // Then term already exists
            Term& existing = *p.first;
            assert(&existing != &t);
            assert(t.degree == existing.degree);
            for (int i = 0; i < t.degree; ++i) {
                assert(t.vars[i] == existing.vars[i]);
            }
            existing.coeff += t.coeff;
            t.coeff = 0;
        }
        return;
    }
}

template <typename R, int D>
inline void HigherOrderEnergyGraph<R, D>::AddClique(const std::vector<VarId>& vars,
        const std::vector<R>& energyTable) {
    const unsigned int size = vars.size();
    const unsigned int numAssignments = 1 << size;
    assert(energyTable.size() == numAssignments);

    std::vector<R> coeffs(numAssignments);
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
        AddTerm(coeffs[subset], degree, subsetVars);
    }
}

template <typename R, int D>
inline void HigherOrderEnergyGraph<R, D>::AddUnaryTerm(VarId var, R coeff) {
    _varRecords[var]._coeff += coeff;
}

template <typename R, int D>
void HigherOrderEnergyGraph<R, D>::_EliminatePositiveTerms() {
    size_t numVars = _varRecords.size();
    for (size_t varIndex = 0; varIndex < numVars; ++varIndex) {
        R positiveSum = 0;
        VarId newPosVar = AddVar();

        VarRecord& vr = _varRecords.at(varIndex);
        assert(vr._id == varIndex);

        std::sort(vr._terms.begin(), vr._terms.end(), 
                [&](int i1, int i2) { return TermComp()(*_terms[i1], *_terms[i2]); }
                );

        for (auto termIdx : vr._terms) {
            assert(termIdx < _terms.size());
            Term& t = *_terms.at(termIdx);
            if (_termSet.find(t) == _termSet.end()) {
                assert(t.coeff == 0);
                continue;
            }
            assert(t.degree > 1);
            assert(t.degree <= D);
            assert(_termSet.find(t) != _termSet.end()); // May be too conservative
            //std::cout << "\t" << t.ToString() << std::endl;
            if (t.coeff <= 0)
                continue;

            positiveSum += t.coeff;
            _termSet.erase(_termSet.iterator_to(t));
            auto newEnd = std::remove(t.vars, t.vars + t.degree, varIndex);
            assert(newEnd+1 == t.vars+t.degree);
            t.degree--;
            auto coeff = t.coeff;
            AddTerm(t);
            t.vars[t.degree] = newPosVar;
            AddTerm(-coeff, t.degree+1, t.vars);
        }
        VarId quadratic[2];
        quadratic[0] = vr._id;
        quadratic[1] = newPosVar;
        AddTerm(positiveSum, 2, quadratic);
    }
}

template <typename R, int D>
template <typename QR>
void HigherOrderEnergyGraph<R, D>::_ReduceNegativeTerms(QR& qr) {
    // Estimate expected size of quadratic problem. Only nodes/edges are
    // created below, so we can count them ahead of time
    int expectedEdges = 0;
    for (const auto& t : _termSet) {
        expectedEdges += (t.degree > 2) ? t.degree : 1;
    }
    //std::cout << "\tExpected Vars: " << expectedVars << "\tExpected Edges: " << expectedEdges << std::endl;

    qr.SetMaxEdgeNum(expectedEdges);
    qr.AddNode(_varCounter);

    // Term-by-term reduction from Friedman & Drineas
    for (const auto& t : _termSet) {
        if (t.degree == 2) {
            qr.AddPairwiseTerm(t.vars[0], t.vars[1], 0, 0, 0, t.coeff);
        } else {
            assert(t.degree > 2);
            assert(t.coeff <= 0);
            if (t.coeff == 0)
                continue;
            auto w = qr.AddNode();
            for (int i = 0; i < t.degree; ++i) {
                qr.AddPairwiseTerm(t.vars[i], w, 0, 0, 0, t.coeff);
            }
            qr.AddUnaryTerm(w, 0, t.coeff*(1-t.degree));
        }
    }
    BOOST_FOREACH(VarRecord& vr, _varRecords) {
        qr.AddUnaryTerm(vr._id, 0, vr._coeff);
    }
}

template <typename R, int D>
template <typename QR>
inline void HigherOrderEnergyGraph<R, D>::ToQuadratic(QR& qr) {
    _EliminatePositiveTerms();
    _ReduceNegativeTerms(qr);
}

template <typename R, int D>
inline void HigherOrderEnergyGraph<R, D>::Clear() {
    _varRecords.clear();
};

#endif
