#ifndef _FOE_CLIQUES_HPP_
#define _FOE_CLIQUES_HPP_
/*
 * foe-cliques.hpp
 *
 * Copyright 2012 Alexander Fix
 * See LICENSE.txt for license information
 *
 * Example of deriving from CliqueEnergy (defined in clique.hpp) to implement 
 * the Field of Experts energy used for the denoising algorithm
 */

#include "clique.hpp"
#include <math.h>
#include <iostream>
#include <cstdint>

typedef int REAL;
const double DoubleToREAL = 10000;

/*
 * The Field of Experts energy, defined for a 2x2 patch of the image.
 * Note that the only thing we really need to override from the abstract 
 * base class is operator(), which actually calculates the FoE energy of a 
 * 2x2 patch.
 */
class FoEEnergy : public CliqueEnergy<REAL, unsigned char, 4> {
    public:
        FoEEnergy(int size, int nbd[])
            : CliqueEnergy<REAL, unsigned char, 4>(size, nbd) { }
        virtual REAL operator()(const unsigned char buf[]) const;
        virtual void AddGradient(double grad[], const unsigned char image[]) const;
};

/*
 * The unary energy is defined for a single pixel. It penalizes the squared
 * distance from the original observed value. Note that we've added a new data
 * member _orig, which keeps track of the originally observed value. 
 */
extern double FoEUnarySigma;
template <int D>
class FoEUnaryEnergy : public CliqueEnergy<REAL, unsigned char, D> {
    public:
        FoEUnaryEnergy(int *index, unsigned char originalImagePixel) 
            : CliqueEnergy<REAL, unsigned char, D>(1, index),
              _orig(originalImagePixel) { }

        virtual REAL operator()(const unsigned char buf[]) const {
            double dist = (double)_orig - (double)buf[0];
            double e = dist*dist / (FoEUnarySigma*FoEUnarySigma * 2);
            return DoubleToREAL * e;
        }
        virtual void AddGradient(double grad[], const unsigned char image[]) const {
            grad[this->_neighborhood[0]] += ((double)image[this->_neighborhood[0]] - (double)_orig) / (FoEUnarySigma * FoEUnarySigma);
        }

    private:
        unsigned char _orig;
};


class FoE2x3Energy : public CliqueEnergy<REAL, unsigned char, 6> {
    public:
        FoE2x3Energy(int size, int nbd[])
            : CliqueEnergy<REAL, unsigned char, 6>(size, nbd) { }
        virtual REAL operator()(const unsigned char buf[]) const;
        virtual void AddGradient(double grad[], const unsigned char image[]) const;
};

class FoE3x3Energy : public CliqueEnergy<REAL, unsigned char, 9> {
    public:
        FoE3x3Energy(int size, int nbd[])
            : CliqueEnergy<REAL, unsigned char, 9>(size, nbd) { }
        virtual REAL operator()(const unsigned char buf[]) const;
        virtual void AddGradient(double grad[], const unsigned char image[]) const;

        static void InitFoE3x3();
};
#endif
