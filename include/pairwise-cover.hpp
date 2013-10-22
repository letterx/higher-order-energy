#ifndef _PAIRWISE_COVER_HPP_
#define _PAIRWISE_COVER_HPP_

#include <cstdint>
#include <vector>
#include <list>
#include <boost/foreach.hpp>
#include <cassert>

template <typename R, int D>
class PairwiseCover {
    public:
        typedef int VarId;
        typedef VarId NodeId;

        // Constructs empty PairwiseCover with no variables or terms
        PairwiseCover();

        // Adds variables to the PairwiseCover. Variables must be added 
        // before any terms referencing them can be added
        VarId AddVar();
        VarId AddVars(int n);
        VarId NumVars() const { return _varCounter; }
        NodeId AddNode(int n = 1) { return AddVars(n); }

        // Adds a clique defined by an table of energies, from 0 to 1 << k
        // where k is the size of the clique.
        void AddClique(const std::vector<VarId>& vars, const std::vector<R>& energy);

        void AddConstantTerm(R coeff) { _constantTerm += coeff; }
        void AddUnaryTerm(VarId v, R coeff);
        void AddUnaryTerm(VarId v, R E0, R E1) { AddUnaryTerm(v, E1 - E0); } 
    
        // Reduces the PairwiseCover to quadratic form
        // NOTE: THIS IS A DESTRUCTIVE OPERATION, so do not rely on the 
        // PairwiseCover being in any useful state after this operation.
        //
        // This is a templated function, so it will work with any class 
        // implementing the necessary interface. See quadratic-rep.hpp for
        // the minimal requirements.
        template <typename QR>
        void ToQuadratic(QR &qr);

    private:
        typedef uint32_t Assgn;
        struct Clique {
            int size;
            VarId vars[D];
            std::vector<R> subset_coeffs;
            int a_size;
            std::vector<R> a_coeffs;
            std::vector<R> b_coeffs;
        };

        void ComputeBeta();

        R _constantTerm;
        VarId _varCounter;
        std::vector<R> _unaryTerms;
        std::vector<Clique> _cliques;
};

template <typename R, int D>
inline PairwiseCover<R, D>::PairwiseCover()
    : _constantTerm(0), _varCounter(0), _unaryTerms(), _cliques()
{ }


template <typename R, int D>
inline typename PairwiseCover<R, D>::VarId 
PairwiseCover<R, D>::AddVar() {
    _unaryTerms.push_back(0);
    return _varCounter++;
}

template <typename R, int D>
inline typename PairwiseCover<R, D>::VarId 
PairwiseCover<R, D>::AddVars(int n) {
    VarId firstVar = _varCounter;
    for (int i = 0; i < n; ++i)
        this->AddVar();
    return firstVar;
}

/*
template <typename R, int D>
template <typename QR>
R PairwiseCover<R, D>::ComputeEnergy(QR& qr) {
    R energy = _constantTerm;
    for (NodeId i = 0; i < _varCounter; ++i) {
        if (qr.GetLabel(i) == 1)
            energy += _unaryTerms[i];
    }
    for (const Clique& c : _cliques) {
        Assgn max_subset = (1 << c.size);
        for (Assgn subset = 0; subset < max_subset; ++subset) {
            bool all_set = true;
            for (int i = 0; i < c.size; ++i) {
                if ((subset & (1 << i)) && (qr.GetLabel(c.vars[i]) != 1))
                    all_set = false;
            }
            if (all_set)
                energy += c.subset_coeffs[subset];
        }
    }
    return energy;
}
*/

template <typename R, int D>
inline void PairwiseCover<R, D>::AddClique(const std::vector<VarId>& vars,
        const std::vector<R>& energyTable) {
    const unsigned int size = vars.size();
    assert(size > 1);
    const unsigned int numAssignments = 1 << size;
    assert(energyTable.size() == numAssignments);

    _cliques.push_back(Clique());
    Clique& c = _cliques.back();
    c.size = size;
    for (size_t i = 0; i < size; ++i)
        c.vars[i] = vars[c.size - i - 1];
    c.subset_coeffs = std::vector<R>(numAssignments, 0);
    c.a_size = (size+1)/2;
    c.a_coeffs = std::vector<R>(1 << c.a_size, 0);
    c.b_coeffs = std::vector<R>(1 << (size - c.a_size), 0);

    // For each boolean assignment, add that to the corresponding monomials
    // with the correct parity
    for (unsigned int assignment = 0; 
            assignment < numAssignments; 
            ++assignment) 
    {
        const R energy = energyTable[assignment];
        for (unsigned int subset = 0; 
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
                c.subset_coeffs[subset] += parity ? -energy : energy;
            }
        }
    }
}

template <typename R, int D>
inline void PairwiseCover<R, D>::AddUnaryTerm(VarId var, R coeff) {
    _unaryTerms[var] += coeff;
}

template <typename R, int D>
inline void PairwiseCover<R, D>::ComputeBeta() {
    for (Clique& c: _cliques) {
        Assgn everything = (1 << c.size) - 1;
        Assgn a_vars = (1 << c.a_size) - 1;
        Assgn b_vars = everything & ~a_vars;
        for (Assgn term = everything; term > 0; --term) {
            Assgn term_a = term & a_vars;
            Assgn term_b = term & b_vars;
            assert(term_a <= a_vars);
            assert(term_b <= b_vars);
            assert((term_a | term_b) == term);
            if (term_a == 0) {
                // Then term is subset of B, hence is in pairwise cover
                // Split it into a singleton and rest of term, singleton
                // is least set bit
                assert(term_b == term);
                // term_a = don't care, it's a singleton, so beta doesn't matter
                term_b = (term - 1) & term;
                if (term_b == 0) continue; // Then term was just a singleton anyways
                c.b_coeffs[term_b >> c.a_size] += 2*abs(c.subset_coeffs[term]) 
                    + c.b_coeffs[term >> c.a_size];
            } else if (term_b == 0) {
                // Then term is subset of A, hence is in pairwise cover
                // Split it into a singleton and rest of term, singleton
                // is least set bit
                assert(term_a == term);
                term_a = (term - 1) & term;
                // term_b = don't care, it's a singelton, so beta doesn't matter
                if (term_a == 0) continue; // Then term was just a singleton anyways
                c.a_coeffs[term_a] += 2*abs(c.subset_coeffs[term]) 
                    + c.a_coeffs[term];
            } else {
                c.a_coeffs[term_a] += 2*abs(c.subset_coeffs[term]);
                c.b_coeffs[term_b >> c.a_size] += 2*abs(c.subset_coeffs[term]);
            }
        }
    }
}

template <typename R, int D>
template <typename QR>
inline void PairwiseCover<R, D>::ToQuadratic(QR& qr) {
    typedef typename QR::NodeId NodeId;
    ComputeBeta();
    qr.AddNode(_varCounter);
    // REMEMBER: SCALE ALL COSTS BY 2 WHEN ADDING TO QR
    for (NodeId i = 0; i < _varCounter; ++i)
        qr.AddUnaryTerm(i, 0, 2*_unaryTerms[i]);
    for (const Clique& c: _cliques) {
        // TODO(afix): note, there's some waste here, since we
        // don't actually need nodes for the singeltons / empty
        // subsets. Shouldn't really affect things too much
        NodeId a_nodes = qr.AddNode(c.a_coeffs.size());
        NodeId b_nodes = qr.AddNode(c.b_coeffs.size());
        Assgn everything = (1 << c.size) - 1;
        Assgn a_vars = (1 << c.a_size) - 1;
        Assgn b_vars = everything & ~a_vars;
        for (Assgn term = 1; term <= everything; ++term) {
            Assgn term_a = term & a_vars;
            Assgn term_b = term & b_vars;
            if (term_a == 0) {
                // Then term is subset of B, hence is in pairwise cover
                // Split it into a singleton and rest of term, singleton
                // is least set bit
                assert(term_b == term);
                term_b = (term - 1) & term;
                if (term_b == 0) { // Then term was singelton
                    int bit_set = __builtin_ctz(term); // Find which one it is
                    assert(term == Assgn(1 << bit_set));
                    qr.AddUnaryTerm(c.vars[bit_set], 0, 2*c.subset_coeffs[term]);
                } else {
                    R beta_H = c.b_coeffs[term >> c.a_size];
                    R alpha_H = c.subset_coeffs[term];
                    term_a = term & ~term_b;
                    assert(term_a != 0);
                    assert(term_b < b_vars && term_a < b_vars);
                    assert((term_a | term_b) == term);
                    NodeId a_node = b_nodes + (term_a >> c.a_size);
                    NodeId b_node = b_nodes + (term_b >> c.a_size);
                    if (__builtin_popcount(term_a) == 1) // It's a singelton, identify with original var
                        a_node = c.vars[__builtin_ctz(term_a)];
                    if (__builtin_popcount(term_b) == 1) 
                        b_node = c.vars[__builtin_ctz(term_b)];
                    qr.AddPairwiseTerm(a_node, b_node, 0, 0, 0, 2*alpha_H + beta_H);

                    // Also, add unary and connections to original vars
                    NodeId var = b_nodes + (term >> c.a_size);
                    qr.AddUnaryTerm(var, 0, beta_H*(2*__builtin_popcount(term) - 1));
                    for (int i = 0; i < c.size; ++i) {
                        if (term & (1 << i))
                            qr.AddPairwiseTerm(var, c.vars[i], 0, 0, 0, -2*beta_H);
                    }
                }
            } else if (term_b == 0) {
                // Then term is subset of A, hence is in pairwise cover
                // Split it into a singleton and rest of term, singleton
                // is least set bit
                assert(term_a == term);
                term_a = (term - 1) & term;
                if (term_a == 0) { // Then term was singelton
                    int bit_set = __builtin_ctz(term); // Find which one it is
                    assert(term == Assgn(1 << bit_set));
                    qr.AddUnaryTerm(c.vars[bit_set], 0, 2*c.subset_coeffs[term]);
                } else {
                    R beta_H = c.a_coeffs[term];
                    R alpha_H = c.subset_coeffs[term];
                    term_b = term & ~term_a;
                    assert(term_a != 0);
                    assert(term_b < a_vars && term_a < a_vars);
                    assert((term_a | term_b) == term);
                    NodeId a_node = a_nodes + term_a;
                    NodeId b_node = a_nodes + term_b;
                    if (__builtin_popcount(term_a) == 1) // It's a singelton, identify with original var
                        a_node = c.vars[__builtin_ctz(term_a)];
                    if (__builtin_popcount(term_b) == 1) 
                        b_node = c.vars[__builtin_ctz(term_b)];
                    qr.AddPairwiseTerm(a_node, b_node, 0, 0, 0, 2*alpha_H + beta_H);

                    // Also, add unary and connections to original vars
                    NodeId var = a_nodes + term;
                    qr.AddUnaryTerm(var, 0, beta_H*(2*__builtin_popcount(term) - 1));
                    for (int i = 0; i < c.size; ++i) {
                        if (term & (1 << i))
                            qr.AddPairwiseTerm(var, c.vars[i], 0, 0, 0, -2*beta_H);
                    }
                }
            } else {
                R alpha_H = c.subset_coeffs[term];
                NodeId a_node = a_nodes + term_a;
                if (__builtin_popcount(term_a) == 1) // It's a singelton, identify with original var
                    a_node = c.vars[__builtin_ctz(term_a)];
                NodeId b_node = b_nodes + (term_b >> c.a_size);
                if (__builtin_popcount(term_b) == 1) 
                    b_node = c.vars[__builtin_ctz(term_b)];
                qr.AddPairwiseTerm(a_node, b_node, 0, 0, 0, 2*alpha_H);
            }
        }
    }

    //std::cout << "Constant term: " << _constantTerm << "\n";
    //std::cout << "Before Energy: " << ComputeEnergy(qr) << "\n";
    //std::cout << "Before Energy: " << ComputeEnergy2(qr) << "\n";
    //std::cout << "QPBO Energy:   " << qr.ComputeTwiceEnergy() << "\n";
    qr.Solve();

    /*
    int one_labels = 0;
    int aux_labels = 0;
    NodeId node_base = _varCounter;
    for (const Clique& c : _cliques) {
        NodeId a_nodes = node_base;
        node_base += c.a_coeffs.size();
        NodeId b_nodes = node_base;
        node_base += c.b_coeffs.size();
        Assgn everything = (1 << c.size) - 1;
        Assgn a_vars = (1 << c.a_size) - 1;
        Assgn b_vars = everything & ~a_vars;
        for (Assgn a = 0; a <= a_vars; ++a) {
            if (__builtin_popcount(a) > 1) {
                int label = qr.GetLabel(a_nodes + a);
                if (label >= 0) aux_labels++;
                if (label == 1) {
                    one_labels++;
                    for (int i = 0; i < c.size; ++i) {
                        if (a & (1 << i)) {
                            assert(qr.GetLabel(c.vars[i]) != 0);
                            qr.SetLabel(c.vars[i], 1);
                            assert(qr.GetLabel(c.vars[i]) == 1);
                        }
                    }
                } else {
                    bool sublabeled = true;
                    for (int i = 0; i < c.size; ++i) {
                        if ((a & (1 << i)) && (qr.GetLabel(a_nodes+a) < 1))
                            sublabeled = false;
                    }
                    assert(!sublabeled);
                }
            }
        }
        for (Assgn b = 0; b <= (b_vars >> c.a_size); ++b) {
            Assgn shifted_b = b << c.a_size;
            if (__builtin_popcount(shifted_b) > 1) {
                int label = qr.GetLabel(b_nodes + b);
                if (label >= 0) aux_labels++;
                if (label == 1) {
                    one_labels++;
                    for (int i = 0; i < c.size; ++i) {
                        if (shifted_b & (1 << i)) {
                            assert(qr.GetLabel(c.vars[i]) != 0);
                            qr.SetLabel(c.vars[i], 1);
                            assert(qr.GetLabel(c.vars[i]) == 1);
                        }
                    }
                } else {
                    bool sublabeled = true;
                    for (int i = 0; i < c.size; ++i) {
                        if ((shifted_b & (1 << i)) && (qr.GetLabel(b_nodes+b) < 1))
                            sublabeled = false;
                    }
                    assert(!sublabeled);
                }
            }
        }
    }
    */
    //std::cout << "Auxilliary labels: " << aux_labels << "\n";
    //std::cout << "Auxilliary swaps: " << one_labels << "\n";
    qr.ComputeWeakPersistencies();
    //std::cout << "After Energy:  " << ComputeEnergy(qr) << "\n";
    //std::cout << "After Energy:  " << ComputeEnergy2(qr) << "\n";
    //std::cout << "QPBO Energy:   " << qr.ComputeTwiceEnergy(1) << "\n";
}


#endif
