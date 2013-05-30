/*
 * higher-order-example.cpp
 *
 * Copyright 2012 Alexander Fix
 * See LICENSE.txt for license information
 *
 * Provides an example implementation of the Field of Experts denoising
 * algorithm, using a blur-and-random fusion move algorithm
 */

#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/program_options.hpp>
#include "foe-cliques.hpp"
#include "image.hpp"
#include "higher-order.hpp"

double sigma = 20.0;
int kernelRadius;
std::vector<double> gaussianKernel;
REAL threshold = 100.0 * DoubleToREAL;
int thresholdIters = 20;

// Set up the RNGs
static boost::mt19937 rng;
static boost::uniform_int<> uniform255(0, 255);
static boost::uniform_int<> uniform3sigma(-1.5*sigma, 1.5*sigma);
static boost::variate_generator<boost::mt19937&, boost::uniform_int<> > noise255(rng, uniform255);
static boost::variate_generator<boost::mt19937&, boost::uniform_int<> > noise3sigma(rng, uniform3sigma);


Image_uc GetProposedImage(const Image_uc& im, unsigned int iteration, Image_uc& blur);
Image_uc GetProposedImageGrad(const Image_uc& im, unsigned int iteration, CliqueSystem<REAL, unsigned char, 4> &cs, double eta);
CliqueSystem<REAL, unsigned char, 4> SetupCliques(const Image_uc& im);
void InitGaussKernel(double sigma, int& radius, std::vector<double>& kernel);
Image_uc ApplyGaussBlur(const Image_uc& im, int radius, const std::vector<double>& kernel);
double getPSNR(const Image_uc& d, const Image_uc& o);

int main(int argc, char **argv) {
    namespace po = boost::program_options;
    // Parse command arguments
    std::string basename;
    std::string infilename;
    std::string outfilename;
    std::vector<std::string> param_methods;
    std::vector<OptType> methods;
    bool lockstep;
    bool computePSNR;
    std::string original_name;
    bool grad_descent;
    double eta;


    int iterations;

    po::options_description desc("Example options");
    desc.add_options()
        ("help", "Display this help message")
        ("image", po::value<std::string>(&basename)->required(), "Base name of image file (without extension)")
        ("method,m", po::value<std::vector<std::string>>(&param_methods), "[fix|hocr] -> Method to use for higher-order reduction")
        ("iters,i", po::value<int>(&iterations)->default_value(300), "Number of iterations to run")
        ("lockstep", po::value<bool>(&lockstep)->default_value(false), "Run in lockstep mode: all methods have same energy problems to solve, determined by first method specified")
        ("original,o", po::value<std::string>(&original_name)->default_value(""), "Filename for original image -- for computing PSNR")
        ("grad,g", po::value<bool>(&grad_descent)->default_value(false), "Flag for using gradient descent proposals")
        ("eta", po::value<double>(&eta)->default_value(60.0), "Scale for gradient descent proposals")
    ;
    po::positional_options_description popts;
    popts.add("image", 1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
            options(desc).positional(popts).run(), vm);

    try {
        po::notify(vm);
        if (param_methods.empty()) {
            methods.push_back(OptType::Fix);
        } else {
            for (const std::string& m : param_methods) {
                if (m == std::string("hocr")) {
                    methods.push_back(OptType::HOCR);
                } else if (m == std::string("grd")) {
                    methods.push_back(OptType::GRD);
                } else if (m == std::string("grd-heur")) {
                    methods.push_back(OptType::GRD_Heur);
                } else if (m == std::string("fix")) {
                    methods.push_back(OptType::Fix);
                } else {
                    std::cout << "Unrecognized method type: " << m << "\n";
                    exit(-1);
                }
            }
        }
    } catch (std::exception& e) {
        std::cout << "Parsing error: " << e.what() << "\n";
        std::cout << "Usage: higher-order-example [options] basename\n";
        std::cout << desc;
        exit(-1);
    }
    infilename = basename + ".pgm";


    // Initialize the gaussian kernel for blurring the image
    InitGaussKernel(2.0, kernelRadius, gaussianKernel);

    Image_uc in = ImageFromFile(infilename.c_str());
    Image_uc current, blur;

    Image_uc original;
    if (original_name != std::string("")) {
        original = ImageFromFile(original_name.c_str());
        computePSNR = true;
    }


    // Set up the clique system, which defines the energy to be minimized
    // by fusion move
    CliqueSystem<REAL, unsigned char, 4> cliques = SetupCliques(in);

    // energies keeps track of last [thresholdIters] energy values to know
    // when we reach convergence
    REAL energies[thresholdIters];

    std::map<OptType, std::vector<FusionStats>> allStats;

    std::vector<OptType> full_run_methods; //Methods to run the full optimization on
    if (lockstep)
        full_run_methods.push_back(methods[0]);
    else
        full_run_methods = methods;
        
    for (OptType ot : full_run_methods) {
        current.Copy(in);
        rng.seed(0);
        std::chrono::system_clock::time_point startTime = std::chrono::system_clock::now();
        for (int i = 0; i < iterations; ++i) {
            std::cout << "Iteration " << i+1 << "\t\t";
            FusionStats stats;
            stats.iter = i;

            REAL energy  = cliques.Energy(current.Data()); 
            stats.initialEnergy = energy;
            // check if we've reached convergence
            if (i > thresholdIters 
                    && energies[i%thresholdIters] - energy < threshold) {
                break;
            }
            // Do some statistic gathering
            energies[i%thresholdIters] = energy;
            std::cout << "Current Energy: " << (double)energy / DoubleToREAL << std::endl;

            // Real work here: get proposed image, then fuse it with current image
            Image_uc proposed;
            if (grad_descent)
                proposed = GetProposedImageGrad(current, i, cliques, eta);
            else
                proposed = GetProposedImage(current, i, blur);

            if (lockstep) {
                std::vector<Image_uc> outputs;
                for (OptType lockstep_ot : methods) {
                    std::cout << "\t" << ToString(lockstep_ot) << "...\t";
                    std::cout.flush();
                    Image_uc tmp(current.Height(), current.Width());
                    FusionMove(stats, current.Height()*current.Width(), current.Data(), proposed.Data(), tmp.Data(), cliques, lockstep_ot);
                    REAL e = cliques.Energy(tmp.Data()); 
                    std::cout << e << "\n";

                    stats.finalEnergy = e;
                    if (computePSNR)
                        stats.psnr = getPSNR(tmp, original);
                    std::chrono::duration<double> cumulativeTime = (std::chrono::system_clock::now() - startTime);
                    stats.cumulativeTime = cumulativeTime.count();
                    allStats[lockstep_ot].push_back(stats);
                    outputs.push_back(tmp);
                }
                current.Copy(outputs[0]);
            } else {
                FusionMove(stats, current.Height()*current.Width(), current.Data(), proposed.Data(), current.Data(), cliques, ot);
                stats.finalEnergy = cliques.Energy(current.Data());
                if (computePSNR)
                    stats.psnr = getPSNR(current, original);
                std::chrono::duration<double> cumulativeTime = (std::chrono::system_clock::now() - startTime);
                stats.cumulativeTime = cumulativeTime.count();
                allStats[ot].push_back(stats);
            }
        }
        std::string optMethod = ToString(ot);
        outfilename = basename + "-" + optMethod + ".pgm";
        ImageToFile(current, outfilename.c_str());
        REAL energy  = cliques.Energy(current.Data());
        std::cout << "Final Energy: " << energy << std::endl;
    }

    for (OptType ot : methods) {
        std::string optMethod = ToString(ot);
        std::string statsName = basename + "-" + optMethod + ".stats";
        std::ofstream statsFile(statsName);
        for (const FusionStats& s : allStats[ot]) {
            statsFile << s.iter << " ";
            statsFile << s.numVars << " ";
            statsFile << s.additionalVars << " ";
            statsFile << s.labeled << " ";
            statsFile << s.swaps << " ";
            statsFile << s.time << " ";
            statsFile << s.cumulativeTime << " ";
            statsFile << double(s.initialEnergy) / DoubleToREAL  << " ";
            statsFile << double(s.finalEnergy) / DoubleToREAL << " ";
            statsFile << s.psnr << " ";
            statsFile << "\n";
        }
    }


    return 0;
}

Image_uc GetProposedImage(const Image_uc& im, unsigned int iteration, Image_uc& blur) {
    Image_uc proposed(im.Height(), im.Width());
    if (iteration % 2 == 0) {
        // On even iterations, proposal is a gaussian-blurred version of the 
        // current image, plus a small amount of gaussian noise
        blur = ApplyGaussBlur(im, kernelRadius, gaussianKernel);
        for (int i = 0; i < im.Height(); ++i) {
            for (int j = 0; j < im.Width(); ++j) {
                int p = (int)blur(i, j) + (int)noise3sigma();
                if (p > 255) p = 255;
                if (p < 0) p = 0;
                proposed(i, j) = (unsigned char)p;
            }
        }
    } else {
        // On odd iterations, proposal is a uniform random image
        for (int i = 0; i < im.Height(); ++i) {
            for (int j = 0; j < im.Width(); ++j) {
                proposed(i, j) = (unsigned char)(noise255());
            }
        }
    }
    return proposed;
}

Image_uc GetProposedImageGrad(const Image_uc& im, unsigned int iteration, CliqueSystem<REAL, unsigned char, 4> &cs, double eta) {
    int size = im.Height() * im.Width();
    boost::shared_array<double> grad(new double[size]);
    for (int i = 0; i < size; ++i)
        grad[i] = 0.0;
    for(const auto& cp : cs.GetCliques())
        cp->AddGradient(grad.get(), im.Data());
    Image_uc proposed(im.Height(), im.Width());
    for (int i = 0; i < size; ++i) {
        double value = ((double)im.At(i) - grad[i] * eta/(iteration + 1));
        if (value < 0)
            value = 0;
        if (value > 255)
            value = 255;
        proposed.At(i) = (unsigned char) value;
        //proposed.At(i) = (unsigned char)(double)im.At(i);
    }    
    return proposed;
}

CliqueSystem<REAL, unsigned char, 4> SetupCliques(const Image_uc& im) {
    CliqueSystem<REAL, unsigned char, 4> cs;
    int height = im.Height();
    int width = im.Width();
    // For each 2x2 patch, add in a Field of Experts clique
    for (int i = 0; i < height - 1; ++i) {
        for (int j = 0; j < width - 1; ++j) {
            int buf[4];
            int bufIdx = 0;
            buf[bufIdx++] = i*width + j;
            buf[bufIdx++] = (i+1)*width + j;
            buf[bufIdx++] = i*width + j+1;
            buf[bufIdx++] = (i+1)*width + j+1;
            cs.AddClique(CliqueSystem<REAL, unsigned char, 4>::CliquePointer(new FoEEnergy(4, buf)));
        }
    }
    // Add the unary terms
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j) {
            int buf[1];
            buf[0] = i*width + j;
            cs.AddClique(CliqueSystem<REAL, unsigned char, 4>::CliquePointer(new FoEUnaryEnergy(buf, im(i, j))));
        }
    }
    return cs;
}

// Calculate the gaussian kernel, given a std-dev sigma
void InitGaussKernel(double sigma, int& radius, std::vector<double>& kernel) {
    radius = ceil(3*sigma); 
    kernel.reserve(2*radius + 1);
    double pi = 4.0 * atan(1.0);
    double oneOverSqrt2PiSigmaSquared = 1.0 / (sqrt(2.0 * pi) * sigma);
    double oneOverTwoSigmaSquared = 1.0 / (2.0* sigma * sigma);
    for (int i = 0; i <= radius; ++i) {
        double value = oneOverSqrt2PiSigmaSquared 
            * exp(-(i*i)*oneOverTwoSigmaSquared);
        kernel[radius+i] = value;
        kernel[radius-i] = value;
    }
    double sum = 0.0;
    for (int i = 0; i < 2*radius + 1; ++i) {
        sum += kernel[i];
    }
    for (int i = 0; i < 2*radius + 1; ++i) {
        kernel[i] = kernel[i] / sum;
    }
}

// Blur an image, given a kernel and its size
Image_uc ApplyGaussBlur(const Image_uc& im, int radius, const std::vector<double>& kernel) {
    Image_uc vertical(im.Height(), im.Width());
    for (int i = 0; i < im.Height(); ++i) {
        for (int j = 0; j < im.Width(); ++j) {
            double acc = 0.0;
            for (int k = 0; k < 2*radius+1; ++k) {
                acc += kernel[k] * im(i + k - radius, j);
            }
            vertical(i, j) = (unsigned char)acc;
        }
    }
    Image_uc horizontal(im.Height(), im.Width());
    for (int i = 0; i < im.Height(); ++i) {
        for (int j = 0; j < im.Width(); ++j) {
            double acc = 0.0;
            for (int k = 0; k < 2*radius+1; ++k) {
                acc += kernel[k] * vertical(i, j + k - radius);
            }
            horizontal(i, j) = (unsigned char)acc;
        }
    }
    return horizontal;
}

double getPSNR(const Image_uc& d, const Image_uc& o)
{
    int N = o.Height() * o.Width();
	int s = 0;
	for (int i = 0; i < N; i++) {
        int diff = d.At(i) - o.At(i);
        s += diff * diff;
    }
	return 20 * log10(255.0 / sqrt((double)s / N));
}

