#ifndef _FUSION_STATS_HPP_
#define _FUSION_STATS_HPP_


struct FusionStats {
    size_t iter;
    size_t numVars;
    size_t additionalVars;
    size_t labeled;
    size_t swaps;
    double time;
    double cumulativeTime;
    double initialEnergy;
    double finalEnergy;

    FusionStats()
        : iter(0), numVars(0), additionalVars(0), labeled(0), swaps(0),
        time(0), cumulativeTime(0), initialEnergy(0), finalEnergy(0)
    { }
};


#endif
