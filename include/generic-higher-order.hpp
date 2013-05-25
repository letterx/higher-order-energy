#ifndef _GENERIC_HIGHER_ORDER_HPP_
#define _GENERIC_HIGHER_ORDER_HPP_

#include "higher-order-energy.hpp"

template <typename Optimizer>
void AddVars(Optimizer& opt, size_t numVars) {
    opt.AddVars(numVars);
}

template <typename Optimizer, typename Energy>
void AddUnaryTerm(Optimizer& opt, int v, Energy coeff) {
    opt.AddUnaryTerm(v, coeff);
}

template <typename Optimizer, typename Energy>
void AddTerm(Optimizer& opt, Energy coeff, int d, const int vars[]) {
    opt.AddTerm(coeff, d, vars);
}

#endif
