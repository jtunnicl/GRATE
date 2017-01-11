#ifndef RIVERPROFILE_H
#define RIVERPROFILE_H

#include <vector>
#include <QVector>
#include <QDateTime>
#include <cmath>
using namespace std;

double gammln2(double xx);

class NodeGSDObject
{

    // Object intended to hold grain size info at each node

public:

    NodeGSDObject();

    vector <double> abrasion;                        // abrasion value for each lithology type (3)
    vector < double > psi;                            // psi (base 2) grain size categories
    vector < vector < double > > pct;                         // Grain-size  (ngsz x nlith)
    float dsg;                                 // Geometric mean grain size
    float d84;
    float d90;
    float stdv;                                // Standard deviation in GSD
    float sand_pct;                            // Percentage of sand (< 2 mm) in GSD

    void norm_frac();

    void dg_and_std();                         // Calculate D50, sand%, geometric (log2)
};

class NodeCHObject
{

    // Object intended to hold info on channel configuration - multiple channels may exist within a reach

public:

    NodeCHObject();

    double flowProp;                           // Proportion of total flow going into this channel
    double depth;                              // Flow depth (m from channel bottom)
    double width;                              // Channel width (m) at bottom of trapezoid for each node
    double b2b;                                // Bank-to-bank width (top of in-channel flow section)
    double flowArea;                           // Flow area within the channel
    double flowPerim;                          // Perimeter, within the channel

    double ustar;                              // Shear velocity
    int ovBank;                                // Flow has gone overbank
    double Tbed;                               // Shear stress acting on the channel bed (Pa)
    double Tbank;                              // Shear stress acting on the channel banks (Pa)
    double Qb_cap;                             // Transport capacity (m3/s)
    double comp_D;                             // The largest grain that the flow can move
    double K;                                  // Estimated division between key stones and bed material load
    double bankHeight;                         // Characteristic bank height above channel bottom (m)
    double Hmax;                               // Bank strength as a vertical upper bank section (m)
    double mu;                                 // Bank strength, relative to bed (afer Millar, 2005)
    double theta;                              // Bank sideslope angle (degrees)

    void chGeom();                             // Calculate x-sec area for a given depth

    void chCentr();                            // Elevation of xsec centre of mass
};

class NodeXSObject
{
    // This object holds info on the cross-section as a whole;
    // channels and grain size are encompassed within this parent object.

public:

    NodeXSObject();

    int node;
    int numChannels;                           // Number of channels
    double wsl;                                // Water surface level (m above sea level)
    vector<NodeCHObject> CHList;               // Vector containing channel characteristics
    double fpWidth;                            // Floodplain width (m)
    double chSinu;                             // Sinuosity (>1, channel length/valley length)
    double topW;                               // Total width of water surface, across all channels
    double velocity;                           // Mean velocity (m/s) at each node
    double xsDepth;                            // Total (maximum) flow depth, including overbank, in reach
    int mainChannel;                           // Channel with deepest flow
    double xsFlowArea[3];                      // [0] Channel [1] Floodplain [2] Total area
    double xsFlowPerim[3];                     // [0] Channel [1] Floodplain [2] Total perimeter
    double hydRadius;                          // Hydraulic radius
    double critdepth;                          // Critical depth
    double centr;                              // Vertical centroid of flow
    double rough;                              // Grain roughness height
    double omega;                              // Reciprocal of Dingman's Omega (~prop u*/U), Eqn. 6.17
    double k_mean;                             // Conveyance coefficient
    double eci;                                // Energy coefficient related to channel form drag

    void xsGeom();                             // Calculate x-sec area for a given depth

    void xsCentr();                            // Elevation of xsec centre of mass

    void xsECI(NodeGSDObject F);               // Energy coefficient

};

class TS_Object
{
                                               // Generic Time Series Object for either water or sediment inputs
public:

    TS_Object();

    QDateTime date_time;                       // Date and time of inputs (see QDateTime doc)
    double Q;                                  // m3/s (Qw) for water, m3/s (Qs) for sediment
    int Coord;                                 // Stream-wise coordinate of input (m)
    int GRP;                                   // For sediment: # of sed group
};

class RiverProfile
{

public:

    RiverProfile();                            // Constructor
    // Profile Elements

    int nnodes;                                // No. of points in the computational grid
    int npts;                                  // No. of points in the long-profile supplied (later interpolated to nnodes, if necessary)
    QDateTime cTime;                           // Current model time
    QDateTime startTime;
    QDateTime endTime;
    int counter;
    int yearCounter;
    int dt;                                    // Delta t in seconds
    float dx;                                  // Delta x - distance between cross-sections
    vector<double> xx;                         // Chainage (m) at each node (ordered, increasing)
    vector<double> eta;                        // Elevation (m) at each node (high z to low)

    // Sedimentary Elements

    int ngsz;                                  // 18
    int nlith;                                 // 3
    int ngrp;                                  // 31
    int nlayer;                                // 20; No. of sublayers
    double poro;                               // 0.6
    double default_la;                         // Active layer default thickness at each node (0.5)
    double layer;                              // Storage layer default thickness (5 m)

    vector< vector < NodeGSDObject > > storedf;// Array of subsurface GSD elements [nnodes][#store layers]
    vector<NodeGSDObject> grp;                 // 'Library' of grain size distributions
    vector<NodeGSDObject> F;                   // Array of surface GSD elements [nnodes]
    vector<double> la;                         // Thickness of the active layer (~2 D90)
    vector<int> algrp;                         // Active layer group #
    vector<int> ntop;                          // Top storage layer number (indicates remaining layers beneath current one, '0' means bedrock);
    vector<int> stgrp;                         // Storage layer group #
    vector<float> toplayer;                    // Thickness of top storage layer
    vector<float> bedrock;                     // Elevation of bedrock at each node (m)
    vector<float> rand_nums;                   // 10 random nums for Monte-Carlo run. Uses rand1().

    vector<NodeXSObject> RiverXS;              // Array of river cross-section objects

    // Randomizers
    vector<float> tweakArray;
    float qsTweak;                             // Augment the rate of tributary Qs, Qw inputs
    float qwTweak;
    float substrDial;
    float feedQw;
    float feedQs;
    float HmaxTweak;
    float randAbr;

    vector<float> N;                           // Transition matrix for coarsening or fining mixtures

    void initData();
    
    vector<float> hydroGraph();

    void readData();

    const char *getNextParam(ifstream &openFile, const char *nextParam);

    void getLongProfile(std::ifstream &openFile);

    void getGrainSizeLibrary(std::ifstream &openFile);

    void getLibraryLith(std::ifstream &openFile);

    void getStratigraphy(std::ifstream &openFile);

};


#endif // RIVERPROFILE_H
