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
void AddClique(Optimizer& opt, int d, const Energy *coeffs, const int *vars) {
    std::vector<int> vec_vars(vars, vars+d);
    std::vector<Energy> vec_coeffs(coeffs, coeffs+(1 << d));
    opt.AddClique(vec_vars, vec_coeffs);
}

template <typename REAL, int D>
void AddClique(PBF<REAL, D>& opt, int d, const REAL *coeffs, const int *vars) {
    opt.AddHigherTerm(d, const_cast<int*>(vars), const_cast<REAL*>(coeffs));
}

#endif
