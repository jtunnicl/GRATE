#ifndef SED_H
#define SED_H

#include <vector>
#include "riverprofile.h"
#include "hydro.h"

using namespace std;

class sed
{
public:

    vector< vector < TS_Object > > Qs_series;  // 2D Vector; 1st is sources along grid; 2nd is entries over time.
    vector < TS_Object >Qs_bc;                 // Current discharge boundary conditions, [0] main channel, and [1..] tribs

    NodeGSDObject fpp;                         // Bedload temp item
    vector <NodeGSDObject> p;                  // Array of bedload GSD elements [nnodes]
    vector <NodeGSDObject> df;

    vector <double> Qs;                        // Sediment discharge (m3/s)
    vector <double> deta;                      // Delta bed elevation change
    vector <double> dLa_over_dt;

    sed(RiverProfile *r);

    void initSedSeries(unsigned int nodes);                      // Set inputs

    void setNodalSedInputs(RiverProfile *r);

    NodeGSDObject multiplyGSD(NodeGSDObject &M, NodeGSDObject &N, double weight, RiverProfile *r);

    void computeTransport(RiverProfile *r);

    void exner(RiverProfile *r);
};


#endif // SED_H
