/*
 * random-graph.cpp
 *
 * Copyright 2012 Alexander Fix
 * See LICENSE.txt for license information
 *
 * Generates a random binary hypergraph instance, and solves it using
 * one of several optimization methods
 */

#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <random>
#include <boost/program_options.hpp>
#include "generic-higher-order.hpp"
#include "QPBO.h"

/* 
 * RandomInstance
 *
 * Setup a random graph with n variables and m hyper-edges, of size k
 * Seed the random generator with the same seed to get consistent results
 */
template <typename Opt>
void SetupRandomInstance(Opt& opt, int seed, int n, int m, int k, int unary_range, int term_range);

/*
 * RandomGrid
 *
 * Setup a grid-shaped random graph, with given width and height, where
 * the cliques are rectangles of specified width and height, tiled across grid
 */
template <typename Opt>
void SetupRandomGrid(Opt& opt, int seed, int width, int height, int clique_width, int clique_height, int unary_range, int term_range);

struct Params {
    std::string method;
    int size = 0;
    int num_cliques = 0;
    int clique_width = 0;
    int clique_height = 0;
    int unary_range = 0;
    int term_range = 0;
    int seed = 0;
    bool grid = false;
};

template <typename Opt>
void RandomSolve(Opt& opt, const Params& params);

int main(int argc, char **argv) {
    namespace po = boost::program_options;
    
    Params params;

    po::options_description desc("random-graph options");
    desc.add_options()
        ("help", "Display this help message")
        ("method,m", po::value<std::string>(&params.method)->required(), "Method to use for higher-order reduction")
        ("grid,g", "Set if graph should be a grid")
        ("size,n", po::value<int>(&params.size)->required(), "Number of variables (nongrid) or length of grid (grid)")
        ("cliques", po::value<int>(&params.num_cliques)->default_value(-1), "Number of cliques (nongrid only)")
        ("cw", po::value<int>(&params.clique_width)->default_value(2), "Width of cliques")
        ("ch", po::value<int>(&params.clique_height)->default_value(2), "Height of cliques")
        ("unary,u", po::value<int>(&params.unary_range)->default_value(1000), "Strength of unary terms")
        ("term,t", po::value<int>(&params.term_range)->default_value(1000), "Strength of higher-order terms")
        ("seed,s", po::value<int>(&params.seed)->default_value(0), "Seed for RNG")
    ;

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv).
                options(desc).run(), vm);
        
        if (vm.count("help")) {
            std::cout << desc;
            exit(0);
        }
        if (vm.count("grid"))
            params.grid = true;

        po::notify(vm);

        if (params.num_cliques <= 0)
            params.num_cliques = params.size;
    } catch (std::exception& e) {
        std::cout << "Parsing error: " << e.what() << "\n";
        std::cout << "Usage: random-graph -m <method> -n <size> [options]\n";
        std::cout << desc;
        exit(-1);
    }
    std::cout<<"HERE"<<std::endl;
    if (params.method == std::string("fix")) {
        switch (params.clique_width*params.clique_height) {
            // Why did I have to make clique size compile time dependent??!?!
#define CASE(k) case k: { HigherOrderEnergy<int, k> hoe; RandomSolve(hoe, params); } break
            CASE(1);
            CASE(2);
            CASE(3);
            CASE(4);
            CASE(5);
            CASE(6);
            CASE(7);
            CASE(8);
            CASE(9);
#undef CASE
            default: 
            std::cout << "Clique size bigger than we coded for! Add more cases...\n";
            exit(-1);
            break;
        }
    } else if (params.method == std::string("pc")) {
        switch (params.clique_width*params.clique_height) {
            // Why did I have to make clique size compile time dependent??!?!
#define CASE(k) case k: { PairwiseCover<int, k> hoe; RandomSolve(hoe, params); } break
            CASE(1);
            CASE(2);
            CASE(3);
            CASE(4);
            CASE(5);
            CASE(6);
            CASE(7);
            CASE(8);
            CASE(9);
#undef CASE
            default: 
            std::cout << "Clique size bigger than we coded for! Add more cases...\n";
            exit(-1);
            break;
        }
    } else if (params.method == std::string("pc-grid")) {
        if (!params.grid) {
            std::cout << "Cannot run pc-grid on non-grid! What did you expect?\n";
            exit(-1);
        }
        switch (params.clique_width*params.clique_height) {
            // Why did I have to make clique size compile time dependent??!?!
#define CASE(k) case k: { PairwiseCoverGrid<int, k> hoe; RandomSolve(hoe, params); } break
            CASE(1);
            CASE(2);
            CASE(3);
            CASE(4);
            CASE(5);
            CASE(6);
            CASE(7);
            CASE(8);
            CASE(9);
#undef CASE
            default: 
            std::cout << "Clique size bigger than we coded for! Add more cases...\n";
            exit(-1);
            break;
        }
    } else if (params.method == std::string("hocr")) {
        switch (params.clique_width*params.clique_height) {
            // Why did I have to make clique size compile time dependent??!?!
#define CASE(k) case k: { PBF<int, k> hoe; RandomSolve(hoe, params); } break
            CASE(1);
            CASE(2);
            CASE(3);
            CASE(4);
            CASE(5);
            CASE(6);
            CASE(7);
            CASE(8);
            CASE(9);
#undef CASE
            default: 
            std::cout << "Clique size bigger than we coded for! Add more cases...\n";
            exit(-1);
            break;
        }
    } else {
        std::cout << "Unrecognized option!\n";
        exit(-1);
    }

    return 0;
}

class DummyOpt {
    public:
        struct DummyClique {
            std::vector<int> vars;
            std::vector<int> coeffs;
            DummyClique(const std::vector<int> vars_, const std::vector<int> coeffs_)
                : vars(vars_),
                coeffs(coeffs_) { }
        };

        DummyOpt()
            : m_nvars(0) { }
        void AddVars(int n) { m_nvars = n; m_unary.resize(m_nvars, 0); }
        void AddUnaryTerm(int i, int coeff) { m_unary[i] += coeff; }
        void AddClique(const std::vector<int>& vars, const std::vector<int>& coeffs) { m_cliques.push_back(DummyClique(vars, coeffs)); }
        int Energy(const std::vector<int>& labels) const {
            int e = 0;
            for (int i = 0; i < m_nvars; ++i) {
                if (labels[i] == 1)
                    e += m_unary[i];
            }
            for (const DummyClique& c : m_cliques) {
                uint32_t assgn = 0;
                for (int i = 0; i < c.vars.size(); ++i)
                    if (labels[c.vars[i]] == 1)
                        assgn |= (1 << i);
                e += c.coeffs[assgn];
            }
            return e;
        }

        int m_nvars;
        std::vector<int> m_unary;
        std::vector<DummyClique> m_cliques;
};

template <typename Opt>
void RandomSolve(Opt& opt, const Params& params) {
    int size;
    std::chrono::system_clock::time_point startTime = std::chrono::system_clock::now();
    if (params.grid) {
        size = params.size * params.size;
        SetupRandomGrid(opt, params.seed, params.size, params.size, 
                params.clique_width, params.clique_height, 
                params.unary_range, params.term_range);    
    } else {
        size = params.size;
        SetupRandomInstance(opt, params.seed, params.size, 
                params.num_cliques, params.clique_width,
                params.unary_range, params.term_range);
    }
    QPBO<int> qr(0,0);
    ToQuadratic(opt, qr);
    qr.Solve();
    qr.ComputeWeakPersistencies();

    std::chrono::duration<double> time = std::chrono::system_clock::now() - startTime;

    std::vector<int> labels(size);
    int labeled = 0;
    int ones = 0;
    for (int i = 0; i < size; ++i) {
        labels[i] = qr.GetLabel(i);
        if (labels[i] >= 0) labeled++;
        if (labels[i] == 1) ones++;
    }

    int num_edges = 0;
    QPBO<int>::EdgeId e;
    for (e = qr.GetNextEdgeId(-1); e >= 0; e = qr.GetNextEdgeId(e))
        num_edges++;

    // Set up identical energy with dummy, to compute result
    DummyOpt dummy;
    if (params.grid) 
        SetupRandomGrid(dummy, params.seed, params.size, params.size, 
                params.clique_width, params.clique_height, 
                params.unary_range, params.term_range);    
    else
        SetupRandomInstance(dummy, params.seed, params.size, 
                params.num_cliques, params.clique_width,
                params.unary_range, params.term_range);
    int energy = dummy.Energy(labels);
    std::cout << "Energy:           " << energy << "\n";
    std::cout << "Original vars:    " << size << "\n";
    std::cout << "Labeled:          " << labeled << "\n";
    std::cout << "Ones:             " << ones << "\n";
    std::cout << "Time:             " << time.count() << "\n";
    std::cout << "Additional vars:  " << qr.GetNodeNum() - size << "\n";
    std::cout << "Total edges:      " << num_edges << "\n";
}

template <typename Opt>
void SetupRandomInstance(Opt& opt, int seed, int n, int m, int k, int unary_range, int term_range) {
    std::mt19937 generator;
    generator.seed(seed);

    std::uniform_int_distribution<int> unary_dist(-unary_range, unary_range);
    std::uniform_int_distribution<int> term_dist(-term_range, term_range);
    std::uniform_int_distribution<int> var_dist(0, n-1);
    auto unary_rand = std::bind(unary_dist, generator);
    auto term_rand = std::bind(term_dist, generator);
    auto var_rand = std::bind(var_dist, generator);

    AddVars(opt, n);
    for (int i = 0; i < m; ++i) {
        std::set<int> var_set;
        int added = 0;
        while (added < k) {
            int new_var = var_rand();
            bool inserted = false;
            std::tie(std::ignore, inserted) = var_set.insert(new_var);
            if (inserted) added++;
        }
        std::vector<int> vars(var_set.begin(), var_set.end());
        std::vector<int> energy_table(1 << k);
        for (int j = 0; j < (1 << k); ++j)
            energy_table[j] = term_rand();
        AddClique(opt, k, energy_table.data(), vars.data());
    }
    for (int i = 0; i < n; ++i) {
        AddUnaryTerm(opt, i, unary_rand());
    }

}

template <typename Opt>
void SetupRandomGrid(Opt& opt, int seed, int width, int height, int clique_width, int clique_height, int unary_range, int term_range) {
    std::mt19937 generator;
    generator.seed(seed);

    std::uniform_int_distribution<int> unary_dist(-unary_range, unary_range);
    std::uniform_int_distribution<int> term_dist(-term_range, term_range);
    auto unary_rand = std::bind(unary_dist, generator);
    auto term_rand = std::bind(term_dist, generator);

    AddVars(opt, width*height);
    for (int i = 0; i < width; ++i) {
        for (int j = 0; j < height; ++j) {
            AddUnaryTerm(opt, i*height+j, unary_rand());
        }
    }

    int clique_size = clique_width*clique_height;
    for (int i = 0; i < width-clique_width+1; ++i) {
        for (int j = 0; j < height-clique_height+1; ++j) {
            int vars[clique_size];
            int coeffs[1 << clique_size];
            int *vp = vars;
            for (int k = 0; k < clique_width; ++k) {
                for (int l = 0; l < clique_height; ++l) {
                    *(vp++) = (i+k)*height+j+l;
                }
            }
            for (int assgn = 0; assgn < (1 << clique_size); ++assgn)
                coeffs[assgn] = term_rand();
            AddClique(opt, clique_size, coeffs, vars);
        }
    }
}
