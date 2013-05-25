#ifndef _GENERIC_HIGHER_ORDER_HPP_
#define _GENERIC_HIGHER_ORDER_HPP_

#include "higher-order-energy.hpp"
#include "HOCR.h"

enum class OptType {
    Fix,
    HOCR
};

template <typename Optimizer>
void AddVars(Optimizer& opt, size_t numVars) {
    opt.AddVars(numVars);
}

template <typename REAL, int D>
void AddVars(PBF<REAL, D>& opt, size_t numVars) {
    // noop
}

template <typename Optimizer, typename Energy>
void AddUnaryTerm(Optimizer& opt, int v, Energy coeff) {
    opt.AddUnaryTerm(v, coeff);
}

template <typename REAL, int D>
void AddUnaryTerm(PBF<REAL, D>& opt, int v, REAL coeff) {
    opt.AddUnaryTerm(v, 0, coeff);
}

template <typename Optimizer, typename Energy>
void AddTerm(Optimizer& opt, Energy coeff, int d, const int vars[]) {
    opt.AddTerm(coeff, d, vars);
}

template <typename REAL, int D>
void AddTerm(PBF<REAL, D>& opt, REAL coeff, int d, const int vars[]) {
    std::vector<REAL> coeffs(1 << d, 0);
    coeffs[(1 << d) - 1] = coeff;
    opt.AddHigherTerm(d, const_cast<int*>(vars), coeffs.data());
}

#endif
