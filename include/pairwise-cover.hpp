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

        struct Clique {
            int size;
            VarId vars[D];
            R subset_coeffs[1<<D];
            int a_size;
            R a_coeffs[1 << ((D+1)/2)];
            R b_coeffs[1 << (D/2)];
        };

        VarId NumNodes() const { return _varCounter; }
        size_t NumCliques() const { return _cliques.size(); }

    private:
        typedef uint32_t Assgn;
        R _constantTerm;
        VarId _varCounter;
        std::vector<R> _unaryTerms;
        std::vector<Clique> _cliques;
};

template <typename Clique>
inline void ComputeBeta(Clique& c);

template <typename R, int D, typename QR>
void CliqueToQuadratic(const PairwiseCover<R, D>& pc, const typename PairwiseCover<R, D>::Clique& c, QR& qr);

template <typename R, typename QR>
void CliqueToQuadratic(const PairwiseCover<R, 4>& pc, const typename PairwiseCover<R, 4>::Clique& c, QR& qr);

template <typename R, int D>
int NumQRNodes(const PairwiseCover<R, D>& pc) {
    return pc.NumNodes() + 2*(1 << (D/2))*pc.NumCliques();
}

template <typename R>
int NumQRNodes(const PairwiseCover<R, 4>& pc) {
    return pc.NumNodes() + 2*pc.NumCliques();
}

template <typename R, int D>
int NumQREdges(const PairwiseCover<R, D>& pc) {
    return pc.NumCliques();
}

template <typename R>
int NumQREdges(const PairwiseCover<R, 4>& pc) {
    return 15*pc.NumCliques();
}

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
    Assgn numAssignments = 1 << size;
    assert(energyTable.size() == numAssignments);

    _cliques.push_back(Clique());
    Clique& c = _cliques.back();
    c.size = size;
    for (size_t i = 0; i < size; ++i)
        c.vars[i] = vars[c.size - i - 1];
    for (size_t i = 0; i < numAssignments; ++i)
        c.subset_coeffs[i] = 0;
    c.a_size = (size+1)/2;
    for (int i = 0; i < (1 << c.a_size); ++i)
        c.a_coeffs[i] = 0;
    for (int i = 0; i < (1 << (c.size - c.a_size)); ++i)
        c.b_coeffs[i] = 0;

    // For each boolean assignment, add that to the corresponding monomials
    // with the correct parity
    for (Assgn assignment = 0; assignment < numAssignments; ++assignment) {
        const R energy = energyTable[assignment];
        Assgn assignment_bar = (numAssignments-1)^assignment;
        for (Assgn subset = assignment_bar; subset > 0; subset = (subset-1)&assignment_bar) {
            int parity = __builtin_parity(subset);
            c.subset_coeffs[subset|assignment] += parity ? -energy : energy;
        }
        c.subset_coeffs[assignment] += energy;
    }
    ComputeBeta(c);
}

template <typename R, int D>
inline void PairwiseCover<R, D>::AddUnaryTerm(VarId var, R coeff) {
    _unaryTerms[var] += coeff;
}

template <typename Clique>
inline void ComputeBeta(Clique& c) {
    typedef uint32_t Assgn;
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

template <typename R, int D, typename QR>
void CliqueToQuadratic(const PairwiseCover<R, D>& pc, const typename PairwiseCover<R, D>::Clique& c, QR& qr) {
    // Exact nodes and cliques worked out by hand for size 4 clique

    typedef uint32_t Assgn;
    auto a_nodes = qr.AddNode(1 << c.a_size);
    auto b_nodes = qr.AddNode(1 << (c.size - c.a_size));
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
                auto alpha_H = c.subset_coeffs[term];
                term_a = term & ~term_b;
                assert(term_a != 0);
                assert(term_b < b_vars && term_a < b_vars);
                assert((term_a | term_b) == term);
                auto a_node = b_nodes + (term_a >> c.a_size);
                auto b_node = b_nodes + (term_b >> c.a_size);
                if (__builtin_popcount(term_a) == 1) // It's a singelton, identify with original var
                    a_node = c.vars[__builtin_ctz(term_a)];
                if (__builtin_popcount(term_b) == 1) 
                    b_node = c.vars[__builtin_ctz(term_b)];
                qr.AddPairwiseTerm(a_node, b_node, 0, 0, 0, 2*alpha_H + beta_H);

                // Also, add unary and connections to original vars
                auto var = b_nodes + (term >> c.a_size);
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
                auto beta_H = c.a_coeffs[term];
                auto alpha_H = c.subset_coeffs[term];
                term_b = term & ~term_a;
                assert(term_a != 0);
                assert(term_b < a_vars && term_a < a_vars);
                assert((term_a | term_b) == term);
                auto a_node = a_nodes + term_a;
                auto b_node = a_nodes + term_b;
                if (__builtin_popcount(term_a) == 1) // It's a singelton, identify with original var
                    a_node = c.vars[__builtin_ctz(term_a)];
                if (__builtin_popcount(term_b) == 1) 
                    b_node = c.vars[__builtin_ctz(term_b)];
                qr.AddPairwiseTerm(a_node, b_node, 0, 0, 0, 2*alpha_H + beta_H);

                // Also, add unary and connections to original vars
                auto var = a_nodes + term;
                qr.AddUnaryTerm(var, 0, beta_H*(2*__builtin_popcount(term) - 1));
                for (int i = 0; i < c.size; ++i) {
                    if (term & (1 << i))
                        qr.AddPairwiseTerm(var, c.vars[i], 0, 0, 0, -2*beta_H);
                }
            }
        } else {
            auto alpha_H = c.subset_coeffs[term];
            auto a_node = a_nodes + term_a;
            if (__builtin_popcount(term_a) == 1) // It's a singelton, identify with original var
                a_node = c.vars[__builtin_ctz(term_a)];
            auto b_node = b_nodes + (term_b >> c.a_size);
            if (__builtin_popcount(term_b) == 1) 
                b_node = c.vars[__builtin_ctz(term_b)];
            qr.AddPairwiseTerm(a_node, b_node, 0, 0, 0, 2*alpha_H);
        }
    }
}

template <typename R, typename QR>
void CliqueToQuadratic(const PairwiseCover<R, 4>& pc, const typename PairwiseCover<R, 4>::Clique& c, QR& qr) {
    assert(c.size == 4);
    typedef uint32_t Assgn;
    auto a_node = qr.AddNode(1);
    auto b_node = qr.AddNode(1);
    Assgn a_var = 3;
    Assgn b_var = 12;
    R beta_a = c.a_coeffs[3];
    R beta_b = c.b_coeffs[3];
    for (int i = 0; i < c.size; ++i)
        qr.AddUnaryTerm(c.vars[i], 0, 2*c.subset_coeffs[1 << i]);

    // Edges involving a_node
    qr.AddPairwiseTerm(c.vars[0], c.vars[1], 0, 0, 0, 2*c.subset_coeffs[a_var] + beta_a);
    qr.AddUnaryTerm(a_node, 0, beta_a*3);
    qr.AddPairwiseTerm(a_node, 0, 0, 0, 0, -2*beta_a);
    qr.AddPairwiseTerm(a_node, 1, 0, 0, 0, -2*beta_a);

    // Edges involving b_node
    qr.AddPairwiseTerm(c.vars[2], c.vars[3], 0, 0, 0, 2*c.subset_coeffs[b_var] + beta_b);
    qr.AddUnaryTerm(b_node, 0, beta_b*3);
    qr.AddPairwiseTerm(b_node, 2, 0, 0, 0, -2*beta_b);
    qr.AddPairwiseTerm(b_node, 3, 0, 0, 0, -2*beta_b);
    
    // Rest of the edges
    qr.AddPairwiseTerm(c.vars[0], c.vars[2], 0, 0, 0, 2*c.subset_coeffs[5]);
    qr.AddPairwiseTerm(c.vars[0], c.vars[3], 0, 0, 0, 2*c.subset_coeffs[9]);
    qr.AddPairwiseTerm(c.vars[1], c.vars[2], 0, 0, 0, 2*c.subset_coeffs[6]);
    qr.AddPairwiseTerm(c.vars[1], c.vars[3], 0, 0, 0, 2*c.subset_coeffs[10]);
    qr.AddPairwiseTerm(a_node, c.vars[2], 0, 0, 0, 2*c.subset_coeffs[7]);
    qr.AddPairwiseTerm(a_node, c.vars[3], 0, 0, 0, 2*c.subset_coeffs[11]);
    qr.AddPairwiseTerm(b_node, c.vars[0], 0, 0, 0, 2*c.subset_coeffs[13]);
    qr.AddPairwiseTerm(b_node, c.vars[1], 0, 0, 0, 2*c.subset_coeffs[14]);
    qr.AddPairwiseTerm(a_node, b_node, 0, 0, 0, 2*c.subset_coeffs[15]);
}

template <typename R, int D>
template <typename QR>
inline void PairwiseCover<R, D>::ToQuadratic(QR& qr) {
    typedef typename QR::NodeId NodeId;

    qr.AddNode(_varCounter);
    qr.SetMaxEdgeNum(NumQREdges(*this));
    // REMEMBER: SCALE ALL COSTS BY 2 WHEN ADDING TO QR
    for (NodeId i = 0; i < _varCounter; ++i)
        qr.AddUnaryTerm(i, 0, 2*_unaryTerms[i]);
    for (const Clique& c: _cliques) {
        CliqueToQuadratic(*this, c, qr);
    }

    //std::cout << "Constant term: " << _constantTerm << "\n";
    //std::cout << "Before Energy: " << ComputeEnergy(qr) << "\n";
    //std::cout << "Before Energy: " << ComputeEnergy2(qr) << "\n";
    //std::cout << "QPBO Energy:   " << qr.ComputeTwiceEnergy() << "\n";
    qr.Solve();
    qr.ComputeWeakPersistencies();
    //qr.Improve();

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
    //std::cout << "After Energy:  " << ComputeEnergy(qr) << "\n";
    //std::cout << "After Energy:  " << ComputeEnergy2(qr) << "\n";
    //std::cout << "QPBO Energy:   " << qr.ComputeTwiceEnergy(1) << "\n";
}


#endif
