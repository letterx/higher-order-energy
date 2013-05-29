/*
 * foe-cliques.hpp
 * 
 * Copyright 2012 Alexander Fix
 * See LICENSE.txt for license information
 *
 * Implementation of the Field of Experts cliques
 */
#include "foe-cliques.hpp"

double alpha[3] = {0.586612685392731, 1.157638405566669, 0.846059486257292};
double expert[3][4] = {
	{-0.0582774013402734, 0.0339010363051084, -0.0501593018104054, 0.0745568557931712},
	{0.0492112815304123, -0.0307820846538285, -0.123247230948424, 0.104812330861557},
	{0.0562633568728865, 0.0152832583489560, -0.0576215592718086, -0.0139673758425540}
};


REAL FoEEnergy::operator()(const unsigned char buf[]) const {
    double energy = 0.0;
    if (_size != 4) {
        throw "Wrong size for FoeEnergy";
    }
    for (int i = 0; i < 3; ++i) {
        double dot = 0.0;
        for (int j = 0; j < 4; ++j) {
            dot += expert[i][j] * buf[j];
        }
        energy += alpha[i] * log(1 + 0.5 * dot * dot);
    }
    return energy*DoubleToREAL;
}

double FoEUnaryEnergy::sigma = 20.0;

REAL FoEUnaryEnergy::operator()(const unsigned char buf[]) const {
        double dist = (double)_orig - (double)buf[0];
        double e = dist*dist / (sigma*sigma * 2);
        return DoubleToREAL * e;
}

double alpha9[8] = {1.201143e-01, 7.520515e-02, 9.078330e-02, 1.280545e-01, 6.276734e-02, 1.201840e-01, 1.092460e-01, 1.217102e-01};
double expert9[8][8] = {
    {1.359507e+00, 8.525110e-01, 2.002059e-01, -5.937472e-01, -2.325712e+01, 7.502941e-01, 1.233328e-01, -4.874199e+00},
    {4.772056e+00, 1.361361e+00, 2.460155e+01, -6.766019e-02, -5.826954e-01, -3.316523e-01, -1.474663e-01, 7.793239e-01},
    {-6.690177e-01, 6.026398e-01, -5.105003e+00, -4.109973e+00, 9.656851e-01, -4.186919e+00, -2.120063e+01, 4.745476e-01},
    {1.943910e+00, -1.227122e+01, 1.042264e+01, 7.818959e-01, -4.940315e+00, -1.424995e+00, -1.004415e+01, -1.650005e+00},
    {-1.050778e+00, -2.146662e+00, -1.256774e+00, -1.434304e+00, 1.455968e+00, -2.401311e+01, 8.438753e+00, 1.829643e-01},
    {2.957102e+01, 7.639969e+00, -2.264681e+00, -2.922909e+00, -1.970511e+00, 9.113470e-01, 1.624417e-01, -5.807883e-01},
    {-3.108021e+00, -3.945750e+00, 3.026307e-01, 2.343125e-01, 5.325069e-01, 2.187780e-01, -2.584377e-01, 2.931325e+01},
    {-1.331707e+00, 1.828494e+01, -3.705114e-01, 2.913379e+01, -7.288777e-01, -5.627995e+00, -5.043575e+00, 7.468907e-01}
};
double basis9[8][9] = {
    {1.398454e-02, 2.065362e-04, -1.400086e-02, 1.616972e-02, 5.170091e-05, -1.625187e-02, 1.396836e-02, -1.474055e-04, -1.401394e-02},
    {1.492367e-02, 1.757808e-02, 1.487793e-02, -1.380950e-04, 2.190938e-05, 2.063829e-04, -1.498414e-02, -1.763298e-02, -1.486115e-02},
    {-2.998559e-02, 2.270054e-02, -2.722440e-03, -9.887762e-03, 3.939620e-02, -9.649374e-03, -2.777473e-03, 2.256777e-02, -2.994718e-02},
    {2.332423e-02, 8.093303e-03, -3.993632e-02, -1.693912e-03, 1.974521e-02, -1.491368e-03, -3.982216e-02, 8.283091e-03, 2.336171e-02},
    {-1.408454e-02, -2.905109e-02, -8.821682e-03, 3.730378e-02, 2.846007e-02, 3.776222e-02, -8.717803e-03, -2.833641e-02, -1.470643e-02},
    {-3.994405e-02, 6.108758e-02, -3.279410e-02, 5.985998e-03, 3.985209e-04, -5.208773e-03, 3.253787e-02, -6.191047e-02, 3.983567e-02},
    {3.358739e-02, 6.142535e-03, -4.096458e-02, -6.460519e-02, -6.896593e-04, 6.434845e-02, 4.118725e-02, -5.500375e-03, -3.352551e-02},
    {-3.923639e-02, 6.483810e-02, -3.836991e-02, 6.635343e-02, -1.063839e-01, 6.540974e-02, -3.862160e-02, 6.462672e-02, -3.872724e-02}
};
double **filter9;

void FoE3x3Energy::InitFoE3x3() {
    filter9 = new double *[8];
    for (int i = 0; i < 8; ++i)
        filter9[i] = new double[9];
    for (int i = 0; i < 8; ++i) {
        for (int j = 0; j < 9; ++j) {
            double mult = 0;
            for (int k = 0; k < 8; ++k) {
                mult += basis9[k][j] * expert9[k][i];
            }
            filter9[i][j] = mult;
        }
    }
    /*
    std::cout << "filter9[8][9] = {" << std::endl;
    for (int i = 0; i < 8; ++i) {
        std::cout << "\t{";
        for (int j = 0; j < 9; ++j) {
            std::cout << filter9[i][j] << ", ";
        }
        std::cout << "}," << std::endl;
    }
    std::cout << "};" << std::endl;
    */
}

REAL FoE3x3Energy::operator()(const unsigned char buf[]) const {
    double energy = 0.0;
    if (_size != 9) {
        throw "Wrong size for FoE3Energy";
    }
    for (int i = 0; i < 8; ++i) {
        double dot = 0.0;
        for (int j = 0; j < 9; ++j) {
            dot += filter9[i][j]*buf[j];
        }
        energy += alpha9[i] * log(1 + 0.5 * dot * dot);
    }
    return energy*DoubleToREAL;
}

double FoE3x3UnaryEnergy::sigma = 20.0;

REAL FoE3x3UnaryEnergy::operator()(const unsigned char buf[]) const {
        double dist = (double)_orig - (double)buf[0];
        double e = dist*dist / (sigma*sigma * 2);
        return DoubleToREAL * e;
}
