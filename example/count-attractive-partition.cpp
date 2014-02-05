#include "gurobi_c++.h"
#include <stdint.h>
#include <iostream>
#include <vector>
#include <map>
using namespace std;

int main (int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Usage: count-attractive-partition n\n";
        exit(-1);
    }
    int n = atoi(argv[1]);
    std::cout << "Counting partitions for size " << n << "\n";

    typedef std::map<uint32_t, int> SubsetMap;
    std::vector<SubsetMap> subsetMaps(n+1);
    std::vector<int> subsetCounts(n+1, 0);
    std::vector<double> partitionCounts(n+1, 0);

    for (int i = 1; i <= n; ++i) {
        uint32_t max_subset = 1 << n;
        int subsets_size_i = 0;
        for (uint32_t subset = 0; subset < max_subset; ++subset) {
            int bits_set = 0;
            for (int b = 0; b < n; ++b) {
                if (subset & (1 << b)) {
                    bits_set++;
                }
            }
            if (bits_set == i) {
                subsetMaps[i][subset] = subsets_size_i;
                subsets_size_i++;
            }
        }
        subsetCounts[i] = subsets_size_i;
    }

    for (int i = 1; i < n; ++i) {
        auto& subsets = subsetMaps[i];
        auto& subsetsNext = subsetMaps[i+1];

        try {

        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);
        model.set(GRB_IntAttr_ModelSense, 1.0);

        std::vector<GRBVar> vars;
        for (const auto& p : subsets)
            vars.push_back(model.addVar(0.0, 1.0, 1.0, GRB_BINARY));
        model.update();
        for (const auto& p : subsetsNext) {
            GRBLinExpr constr{};
            auto subset = p.first;
            for (int b = 0; b < n; ++b) {
                if (subset & (1 << b)) {
                    uint32_t neighbor = subset ^ (1 << b);
                    constr += vars[subsets[neighbor]];
                }
            }
            model.addConstr(constr, GRB_GREATER_EQUAL, 1.0);
        }
        model.update();
        model.optimize();
        partitionCounts[i] = model.get(GRB_DoubleAttr_ObjVal);

        } catch (GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch (...) {
            cout << "Exception during optimization" << endl;
            throw;
        }

    }

    int totalSubsets = 0;
    int totalPartitions = 0;
    for (int i = 0; i <= n; ++i) {
        std::cout << "Subsets size " << i << ": " << subsetCounts[i] << "\tPartition size: " << partitionCounts[i] << "\n";
        totalSubsets += subsetCounts[i];
        totalPartitions += partitionCounts[i];
    }
    std::cout << "\nTotal subsets: " << totalSubsets << "\tTotal partition size: " << totalPartitions << "\n";

    return 0;
}

