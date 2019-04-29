#ifndef HYDRO_H
#define HYDRO_H

#include "riverprofile.h"
#include "tinyxml2/tinyxml2.h"
using namespace std;

class hydro
{
public:
    double preissTheta;                        // Theta constant for Preissmann scheme
    double hydUpw;
    unsigned int regimeCounter;                // Regime cross section geometry; march upstream

    vector<double> Qw_Ct;                      // Current discharge, [0] main channel, and [1..] tribs
    vector< vector < TS_Object > > Qw;         // 2D Vector; 1st is sources along grid; 2nd is entries over time.
                                               //   --> Q[Coord][TimeStep]
    vector<double> Fr2;                        // Froude #, squared
    vector<double> QwCumul;
    vector<double> bedSlope;                   // Bedslope

    hydro(RiverProfile *r, XMLElement *params_root);                                   // Constructor

    void backWater(RiverProfile *r);           // Principal Hydro routine: calculate water surface profile

    void initHydro(unsigned int nodes, XMLElement *params_root);

    void setQuasiSteadyNodalFlows(RiverProfile *r);

    void xsCritDepth(unsigned int n, RiverProfile *r, double Q);         // Critical depth at a cross-section for a given Qw

    int energyConserve(unsigned int node, RiverProfile *r);                 // Energy conservation between two nodes

    int quasiNormal(unsigned int node, RiverProfile *r);                    // Quasi-normal approximation of water-surface profile

    void fullyDynamic(RiverProfile *r);                            // Preissmann Scheme approximation of water-surface profile

    vector<double> matsol(int N, vector<vector<double> > EQN);      // Matrix solver

    void regimeModel(unsigned int n, RiverProfile *r);                           // Compute Millar-Eaton equilibrium channel width

    void channelState(unsigned int n, RiverProfile *r);

    void findStable(unsigned int n, RiverProfile *r);

    void setRegimeWidth( RiverProfile *r );

    void findQ(unsigned int n, RiverProfile *r);

    double interp1();
};

#endif // HYDRO_H
