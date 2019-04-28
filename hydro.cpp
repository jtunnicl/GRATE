/*******************
 *
 *
 *  GRATE 9
 *
 *  Hydraulic parameters and flow routing algorithms
 *
 *
 *
*********************/

#include "hydro.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
#include "tinyxml2/tinyxml2.h"
#include "tinyxml2_wrapper.h"
using namespace std;

#define PI 3.14159265
#define G 9.80665
#define RHO 1000  // water density
#define Gs 1.65   // submerged specific gravity

hydro::hydro(RiverProfile *r, XMLElement *params_root)
{
    preissTheta = 0.7;
    hydUpw = 0.3;
    regimeCounter = (r->nnodes-2);

    initHydro(r->nnodes, params_root);
}

void hydro::initHydro(unsigned int nodes, XMLElement *params_root)
{
    double currentCoord = 0.;
    GrateTime NewDate;
    vector< TS_Object > tmp;
    TS_Object NewEntry;

    // get hydro_series element from XML file
    XMLElement *hydro_series = params_root->FirstChildElement("hydro_series");
    if (hydro_series == NULL) {
        throw std::string("Error getting hydro_series element from XML file");
    }

    // loop over all "STEP" elements in the XML file
    for (XMLElement* e = hydro_series->FirstChildElement("STEP"); e != NULL; e = e->NextSiblingElement("STEP")) {
        int year = getIntValue(e, "year");
        int month = getIntValue(e, "month");
        int day = getIntValue(e, "day");
        int hour = getIntValue(e, "hour");
        int minute = getIntValue(e, "minute");
        int second = getIntValue(e, "second");
        NewDate.setDate(year, month, day);
        NewDate.setTime(hour, minute, second);

        NewEntry.date_time = NewDate;
        NewEntry.Q = getDoubleValue(e, "Qw");
        NewEntry.Coord = getIntValue(e, "loc");
        NewEntry.GRP = 1;  // this could be added to the XML file if desired...

        if (NewEntry.Coord > currentCoord) {            // Have we moved to a new source coordinate?
            Qw.push_back( tmp );
            tmp.clear();
            currentCoord = NewEntry.Coord;
            tmp.push_back(NewEntry);  // Start new tmp
        }
        else {
            tmp.push_back(NewEntry);
        }
    }

    Qw.push_back( tmp );                  // Final tmp loaded into Qw array
    Fr2.resize(nodes);
    QwCumul.resize(nodes);
    bedSlope.resize(nodes);
}

void hydro::backWater(RiverProfile *r)
{

    double g = 9.81;
    double FrN2 = 0.8 * 0.8;                     // Threshold for critical flow - 0.8 (squared)
    int iret = 0;
    unsigned int n = 0;
    bool bQuasiNormal = 0;
    unsigned int lastNode = r->nnodes-1;

    setQuasiSteadyNodalFlows(r);

    // Divide QwCumul by Number of Channels !!
//    for ( n = 0; n < r->nnodes; n++)
//    {
//        QwCumul[n] = QwCumul[n] / r->RiverXS[n].noChannels;
//    }


    // update stats for last (down-stream) node:
    r->RiverXS[lastNode].xsArea();               // update area
    r->RiverXS[lastNode].velocity = QwCumul[lastNode] / r->RiverXS[lastNode].flow_area[2];
    r->RiverXS[lastNode].xsPerim();              // update perim
    r->RiverXS[lastNode].xsCentr();              // update centr
    r->RiverXS[lastNode].xsECI(r->F[lastNode]);                // update eci
    Fr2[lastNode] = r->RiverXS[lastNode].eci * r->RiverXS[lastNode].velocity
                  * r->RiverXS[lastNode].velocity / ( g * r->RiverXS[lastNode].depth );


    // update bed slope array
    for ( n = r->nnodes-2; n > 0 ; n-- )
        bedSlope[n] = (hydUpw * ( r->eta[n-1] - r->eta[n] ) / r->dx
               + (1 - hydUpw) * ( r->eta[n] - r->eta[n+1] ) / r->dx)
                   / r->RiverXS[n].chSinu; // note inclusion of sinuosity

    bedSlope[0] = ( r->eta[0] - r->eta[1] ) / r->dx;             // Slope at upstream/downstream nodes
    bedSlope[r->nnodes-1] = ( r->eta[r->nnodes-2] - r->eta[r->nnodes-1] ) / r->dx;

    // Boundary nodes: fixed or computed (default)
    quasiNormal(0, r);
    //quasiNormal(lastNode, r);

    r->RiverXS[lastNode].depth = 0.3 * pow( QwCumul[lastNode],0.3 );
    r->RiverXS[lastNode].wsl = r->eta[lastNode] + r->RiverXS[lastNode].depth;

    for (unsigned int n = r->nnodes-2; n > 0 ; n--)
    {
        xsCritDepth( n, r, QwCumul[n] );       // Calculate critical depth

        // Initial guess at depth
        r->RiverXS[n].depth = 0.3 * pow( QwCumul[n], 0.3 );

        // Flat profile if the bed steps up
        if ( bedSlope[n] < 0 )
            r->RiverXS[n].depth = r->RiverXS[n+1].depth - bedSlope[n] * r->dx;

        r->RiverXS[n].xsArea();                // update area
        r->RiverXS[n].velocity = QwCumul[n] / r->RiverXS[n].flow_area[2];
        r->RiverXS[n].xsPerim();               // update perim
        r->RiverXS[n].xsCentr();               // update centr
        r->RiverXS[n].xsECI(r->F[n]);          // update eci
        Fr2[n] = r->RiverXS[n].eci * r->RiverXS[n].velocity *
                r->RiverXS[n].velocity / ( g * r->RiverXS[n].depth );

        if ( ( Fr2[n] < FrN2 ) || ( bedSlope[n] <= 0 ) || ( n == 0 ) ) // Not super-crit; use energy eqn
            iret = energyConserve(n, r);

        else
        {                                                              // else - recalculate using quasi-normal assumption
            if (bQuasiNormal == 0)
            {
                iret = quasiNormal(n+1, r);                            // Recalculate i+1'th node   <JMW 20080303>
                // if ((iret > 0) || (n < (r->nnodes-3)))
                //     r->RiverXS[n+1].depth = r->RiverXS[n+2].depth;
            }

            iret = quasiNormal(n, r);

            if (iret > 0)
                r->RiverXS[n].depth = r->RiverXS[n+1].depth;

            bQuasiNormal = 1;
        };

        if ( ( iret > 0 ) || ( r->RiverXS[n].depth < r->RiverXS[n].critdepth ) )
            r->RiverXS[n].depth =  r->RiverXS[n].critdepth;

        if ( ( r->RiverXS[n].depth > 0 ) && ( bedSlope[n] > 0 ) )
            r->RiverXS[n].ustar = sqrt( 9.81 * r->RiverXS[n].depth * bedSlope[n] );
        else
            r->RiverXS[n].ustar = 1e-3;

        r->RiverXS[n].wsl = r->eta[n] + r->RiverXS[n].depth;       // Update water surface level at n
    };
}

void hydro::setQuasiSteadyNodalFlows(RiverProfile *r){

    unsigned int j = 0;
    unsigned int i = 0;

    if (Qw[0][0].date_time.secsTo(r->cTime) < 1)                   // Start of run?
        for (i = 0; i < Qw.size(); i++)                            // Qw.size is the # of tribs/sources
            Qw_Ct.push_back( Qw[i][0].Q );                         // Qw_Ct is effectively initialized, here
    else
    {
        j = 0;
        while( Qw[0][j].date_time.secsTo(r->cTime) > 0 )
            j++;
        for (i = 0; i < Qw.size(); i++)
        {
            Qw_Ct[i] = ( Qw[i][j-1].Q + ( Qw[i][j-1].date_time.secsTo(r->cTime) ) *
                       ( Qw[i][j].Q - Qw[i][j-1].Q ) /
                       ( Qw[i][j-1].date_time.secsTo(Qw[i][j].date_time) ));

            Qw_Ct[i] *= r->tweakArray[r->yearCounter];             // Flood = 0.8 to 1.8 mean flow
        }
        Qw_Ct[0] *= r->feedQw;                                     // Feed randomizer
    }

    // Accumulate QwCumul Array.

    QwCumul[0] = Qw_Ct[0];

    i = 0;
    for (j = 1; j < QwCumul.size(); j++)
    {
        QwCumul[j] = QwCumul[j-1];
        if ( i < (Qw.size()-1) && ( r->xx[j] > Qw[i + 1][0].Coord ) )
        {
            i++;
            QwCumul[j] += Qw_Ct[i];
        }
    }
}

void hydro::xsCritDepth(unsigned int n, RiverProfile *r, double Q){

    // Compute critical depth, given a flow

    NodeXSObject& xs = r->RiverXS[n];
    int it = 0;
    int itmax = 50;                        // Max iterations
    double toler = 0.0005;                 // Convergence criteria
    double ymin = 0.15;
    double ymax = xs.bankHeight + 1.5;           //
    double orig_depth = xs.depth;
    double dy = 0;
    double y1 = 0;
    double y2 = 0;                         // dummy variables

    // Make sure ymax is subcritical (ff <= 0); keep increasing ymax until this is so.

    double ff = 1;

    while (ff > 0)
    {
        if (it > 0)
            ymax *= 1.25;

        xs.depth = ymax;
        xs.xsArea();                          // Update statistics at node
        xs.xsPerim();
        xs.xsCentr();                         // re-calc topW
        xs.xsECI(r->F[n]);
        xs.velocity = Q / xs.flow_area[2];
        ff = r->RiverXS[n].eci * r->RiverXS[n].velocity / ( 9.81 * xs.depth ) - 1.0;

        it++;
        if (it > itmax)
        {
            cout << "Unable to initialise max depth for critical depth calculation at xc = \n";
            exit(1);
        }
    }

    y1 = ( ymin + ymax ) / 2.;                   //Initial trial depth

    // Solve by bisection

    while (it < itmax)
    {
        xs.depth = y1;
        xs.xsArea();
        xs.xsPerim();
        xs.xsCentr();
        xs.xsECI(r->F[n]);
        xs.velocity = Q / xs.flow_area[2];
        ff = r->RiverXS[n].eci * r->RiverXS[n].velocity / ( 9.81 * r->RiverXS[n].depth ) - 1.0;

        if (ff < 0)
            ymax = y1;
        else
            ymin = y1;

        y2 = ( ymin + ymax ) / 2.;
        dy = y2 - y1;

        if ( abs( dy / y2 ) < toler )
                break;

        y1 = y2;
        it++;

        if (it > itmax)
        {
            cout << "Critical depth did not converge \n";
            exit(1);
        }
    }
    //Success ... return critical depth
    xs.critdepth = y2;

    // Job done, return appropriate depth to XS
    xs.depth =orig_depth;

    xs.xsArea();
    xs.xsPerim();
    xs.xsCentr();
    xs.xsECI(r->F[n]);
    xs.velocity = Q / xs.flow_area[2];
}

int hydro::energyConserve(unsigned int n, RiverProfile *r)
    {

    // Energy conservation between the two nodes - using bisection algorithm
    // developed by JMW --> v.3.3 energy_conserve2()

    double ff;                  // Objective function
    double h1, h2, hu2;         // Straddle points in bisection
    double Sf, Sf2, Sfx;        // Friction slopes
    double Vhd;                 // Velocity, velocity head at downstream node
    double Vhu;                 // Velocity, velocity head at upstream node
    int flag, iter, itermax;    // Return flag, iteration counter, counter max
    double error;               // error during iterations
    double qm, km;              // Mean Qw and K (conveyance) between nodes

    NodeXSObject& XSu = r->RiverXS[n];           // "upstream" cross-secction  [n]--> Objective
    NodeXSObject& XSd = r->RiverXS[n+1];         // "downstream" cross-section [n+1]--> already computed

    flag = 0;                   // Function flag to be returned
    itermax = 300;

    XSu.xsArea();               // update area
    XSu.velocity = QwCumul[n] / XSu.flow_area[2];
    XSu.xsPerim();              // update perim
    XSu.xsECI(r->F[n]);         // update eci

    //Sf2 = bedSlope[n+1];       // slope gradient between n and n+1
    Vhd = XSd.eci * XSd.velocity * XSd.velocity / (2 * 9.81);
                                 // Velocity head, downtream
    //Sf2 = pow( ( QwCumul[n+1] / XSd.k_mean ), 2);   // Dingman 9B2.4
    Sf2 = ( QwCumul[n+1] ) * ( QwCumul[n+1] ) / ( XSd.k_mean * XSd.k_mean);
                                  // Friction slope, downstream
    h1 = r->RiverXS[n].critdepth;          // Bisection: lower straddle point
    h2 = max(10 * r->RiverXS[n].critdepth, (XSd.depth - bedSlope[n+1] * r->dx) * 2 );  // upper straddle point

    ff = -1;
    while (ff <= 0)
    {
        XSu.depth = h2;             // Update section data based on new depth
        XSu.xsArea();
        XSu.velocity = QwCumul[n] / XSu.flow_area[2];
        XSu.xsPerim();
        XSu.xsECI(r->F[n]);

        Sf = QwCumul[n] * QwCumul[n] / ( XSu.k_mean * XSu.k_mean );
        Vhu = XSu.eci * XSu.velocity * XSu.velocity / (2 * 9.81);

        ff = (XSu.depth + Vhu) - (XSd.depth + Vhd) + ( (bedSlope[n+1] + bedSlope[n]) / 2. - Sf) * r->dx;

        h2 = 2 * XSu.depth;
    }

    h2 = XSu.depth;
    XSu.depth = (h1 + h2) / 2;

    error = 1;
    iter = 0;
    while (error > 5e-4)
    {
        XSu.xsArea();    // Update section data based on new depth
        XSu.velocity = QwCumul[n] / XSu.flow_area[2];
        XSu.xsPerim();
        if ( XSu.depth > 0 )
            XSu.xsECI(r->F[n]);

        //Sf1 = Sf2;        // Initial approximations
        Sf = Sf2;

        Vhu = XSu.eci * XSu.velocity * XSu.velocity / ( 2. * 9.81 );

        if (iter > 1)
        {
            qm = ( QwCumul[n] + QwCumul[n+1] ) / 2.;
            km = ( XSu.k_mean + XSd.k_mean) / 2.;
            Sfx = qm / km;
            Sf = Sfx * Sfx;
        }

        ff = (XSu.depth + Vhu) - (XSd.depth + Vhd) + (bedSlope[n] - Sf) * r->dx;

        if (ff > 0)
            h2 = XSu.depth;
        else
            h1 = XSu.depth;

        if (h2 > XSu.critdepth)
        {
            hu2 = (h1 + h2) / 2.;
            error = abs(hu2 - XSu.depth) / XSu.depth;
            XSu.depth = hu2;
        }
        else
        {
            //  Limit XSu.depth to critical depth
            XSu.depth = XSu.critdepth;
            break;
        }

        iter ++;

        if (iter > itermax)
        {
            cout << "energy_conserve: std step backwater calculation failed to converge \n";
            exit(1);
            flag = 8;
        }

        if (XSu.depth < 0)
        {
            cout << "energy_conserve: negative depth results \n";
            exit(1);
            flag = 16;
        }
    }
    r->RiverXS[n] = XSu;
    return flag;
}

int hydro::quasiNormal(unsigned int n, RiverProfile *r){

    double ff, fp;                             // Objective function and derivative
    int iter, maxiter;
    double error;

    NodeXSObject& XS = r->RiverXS[n];           // Cross-section object to be calculated
    NodeGSDObject& f = r->F[n];
    f.norm_frac();
    f.dg_and_std();

    error = 1;
    iter = 0;
    maxiter = 900;

    while ( error > 0.0001 )
    {
        XS.xsArea();    // Update section data based on new depth
        XS.xsPerim();
        XS.xsCentr();
        if ( XS.depth > 0 )
            XS.xsECI(f);

        ff = QwCumul[n] / XS.topW - XS.depth *
                sqrt( 9.81 * abs( XS.depth ) * bedSlope[n] ) / XS.omega;

        fp = -2.5 * sqrt( 9.81 * abs( XS.depth ) * bedSlope[n] )
                * ( 1.5 * log(11.0 * abs ( XS.depth ) / XS.rough) + 1.0 );

        error = -ff / fp;
        XS.depth += error / 2.;
        error = abs(error / XS.depth);
        ++iter;

        if (iter> maxiter)
        {
            cout << "Iteration Count exceeded in routine quasiNormal at " << n << "\n";
            return 8;
        }
    }

    XS.xsArea();    // Update section data based on new depth
    XS.xsPerim();
    XS.xsCentr();
    XS.xsECI(f);

    return 0;
}

void hydro::fullyDynamic(RiverProfile *r){

    unsigned int i, j, idx, K, iflag, iter, fadj, NNODES;
    double ARI, ARIP1, KI, KIP1, BI, BIP1, ECI, ECIP1, HRI2;
    double AM, DTX2, SF1, SF2, SUMM, TOL, QM, THETA, Ybc;
    double DAY1, DAY2, DCDY1, DCDY2, DSDQ1, DSDQ2, DSDY1, DSDY2;
    double PERI, PERIP1, HRIP12, DPDY1, DPDY2;
    double TERM1, TERM2, TERM3;
    double FR2T, FD_FR_MIN, FD_FR_MAX, FR2_TRIG1, FR2_TRIG2;            // Trigger levels for transition to critical flow

    vector<double> Q, Y, C1, C2, C2A, DF, tmp;
    vector<vector<double> > EQN;

    setQuasiSteadyNodalFlows(r);

    Q.resize(r->nnodes);
    Y.resize(r->nnodes);
    C1.resize(r->nnodes);                       // Matrix elements, used below
    C2.resize(r->nnodes);
    C2A.resize(r->nnodes);                  

    for (i = 0; i < 5; i++)
        tmp.push_back(0.0);

    for (i = 0; i < r->nnodes; i++){
        EQN.push_back(tmp);                    // 5 x nnodes matrix
        EQN.push_back(tmp);                    // this vector is 2 * NNODES in size
        DF.push_back(0.0);                     // Solution matrix gets sent to 'matsol'
        DF.push_back(0.0);
        Q[i] = QwCumul[i];
        Y[i] = r->eta[i] + r->RiverXS[i].depth; // W.S. Elevation
    }

    iflag = 0;
    THETA = preissTheta;                         // Preissmann Weighting Coefficient
    TOL = 0.001;                                 // Tolerance for interactions
    NNODES = r->nnodes;

    quasiNormal(0, r);
    quasiNormal(NNODES-1, r);
    Ybc = r->eta[NNODES-1] + 2.2; // r->RiverXS[NNODES-1].depth;   // d/s boundary condition

    FD_FR_MIN = 0.8;
    FD_FR_MAX = 0.9;

    // COMPUTE TRANSIENT CONDITIONS C
    iter = 0;
    i = 0;
    idx = 0;

    //Variables for transitioning to critical flow by neglecting the
    //spatial derivative of area part of the convective momentum term.
    //FR2_TRIG1 is the trigger level, in terms of Froude No squared.
    //The adjustment transitions limearly between no adjustment at FR2_TRIG1
    //and full adjustment at FR2_TRIG2 where the trigger levels are
    //in terms of Froude No. squared.

    FR2_TRIG1 = FD_FR_MIN * FD_FR_MIN;
    FR2_TRIG2 = FD_FR_MAX * FD_FR_MAX;


    // GENERATE SYSTEM OF EQUATIONS C
    while ( i < NNODES - 1 ) {

        if (i==0){
            r->RiverXS[i].xsArea();
            r->RiverXS[i].velocity = Q[i] / r->RiverXS[i].flow_area[2];
            r->RiverXS[i].xsPerim();
            r->RiverXS[i].xsCentr();
            r->RiverXS[i].xsECI(r->F[i]);
        }

        r->RiverXS[i+1].xsArea();                         // update area, cross-section params for d/s node
        r->RiverXS[i+1].velocity = Q[i+1] / r->RiverXS[i+1].flow_area[2];
        r->RiverXS[i+1].xsPerim();
        r->RiverXS[i+1].xsCentr();
        r->RiverXS[i+1].xsECI(r->F[i+1]);

        ARI = r->RiverXS[i].flow_area[2];                 // Statement flow area
        ARIP1 = r->RiverXS[i+1].flow_area[2];
        AM = ( ARI + ARIP1 ) / 2.;

        KI = r->RiverXS[i].k_mean;
        KIP1 = r->RiverXS[i+1].k_mean;

        ECI = r->RiverXS[i].eci;
        ECIP1 = r->RiverXS[i+1].eci;

        DTX2 = 2 * r->dt / r->dx;

        if (i == 0) Q[i] = QwCumul[0];           //  Upstream boundary condition

        FR2T = ECI * r->RiverXS[i].velocity * r->RiverXS[i].velocity * r->RiverXS[i].topW / ( G * ARI );

        if (FR2T >= FR2_TRIG2) fadj = 0;
            else
            if (FR2T <= FR2_TRIG1) fadj = 1.0;
                else
                fadj = ( FR2_TRIG2-FR2T ) / ( FR2_TRIG2-FR2_TRIG1 );

        //if (Qw[i].Coord = i)..  Conditional clause for trib inputs (not yet implemented!)

        C1[i] = DTX2 * ( 1 - THETA ) * ( Q[i+1] - Q[i] ) - ARI - ARIP1;

        SF1 = abs(Q[i]) * Q[i] / (KI * KI);
        SF2 = abs(Q[i+1]) * Q[i+1] / (KIP1 * KIP1);

        TERM1 = r->dt * (1 - THETA) * G * (ARIP1 * SF2 + ARI * SF1);
        TERM2 = -( Q[i] + Q[i+1] );
        TERM3 = DTX2 * ( 1 - THETA ) *
                ( ECIP1 * Q[i+1] * r->RiverXS[i+1].velocity - ECI * Q[i] * r->RiverXS[i].velocity +
                  G * fadj * ( r->RiverXS[i+1].centr - r->RiverXS[i].centr ) );

        C2[i] = TERM1 + TERM2 + TERM3;
        C2A[i] = -TERM2 * ( 1 - THETA );

        i++;
    }

    SUMM = TOL + 10;

    // 'Line 100'

    while ( SUMM > TOL ){
        i = 0;
        while( i < 2 * NNODES ){                 // Zero out EQN array
            j = 0;
            while( j < 5 ){
                EQN[i][j] = 0.0;
                j++;
            }
            i++;
        }

        // BOUNDARY EQUATIONS C

        EQN[0][1] = 1.0;
        EQN[0][4] = -( Q[0] - QwCumul[0] );
        EQN[2*NNODES-1][2] = 1.0;
        EQN[2*NNODES-1][4] = -( Y[NNODES-1] - Ybc );

        // INTERIOR NODES
        DTX2 = 2 * r->dt / r->dx;

        i = 0;
        while ( i < NNODES - 1 ) {

            // CROSS SECTION UPDATE
            if (i==0){
                r->RiverXS[i].xsArea();
                r->RiverXS[i].velocity = Q[i] / r->RiverXS[i].flow_area[2];
                r->RiverXS[i].xsPerim();
                r->RiverXS[i].xsCentr();
                r->RiverXS[i].xsECI(r->F[i]);
            }

            r->RiverXS[i+1].xsArea();                     // update area, cross-section params for d/s node
            r->RiverXS[i+1].velocity = Q[i+1] / r->RiverXS[i+1].flow_area[2];
            r->RiverXS[i+1].xsPerim();
            r->RiverXS[i+1].xsCentr();
            r->RiverXS[i+1].xsECI(r->F[i+1]);

            ARI = r->RiverXS[i].flow_area[2];             // Statement flow area
            ARIP1 = r->RiverXS[i+1].flow_area[2];
            AM = ( ARI + ARIP1 ) / 2.;

            PERI = r->RiverXS[i].flow_perim[2];           // Wetted Perimeter at node I
            PERIP1 = r->RiverXS[i+1].flow_perim[2];       // Wetted Perimeter at node I+1

            HRI2 = pow( r->RiverXS[i].hydRadius, 0.667 );     // Hyd Radius^0.667 at node I
            HRIP12 = pow( r->RiverXS[i+1].hydRadius, 0.667 ); // Hyd Radius^0.667 at node I+1

            KI = r->RiverXS[i].k_mean;                    // Conveyance at node I
            KIP1 = r->RiverXS[i+1].k_mean;                // Conveyance at node I+1

            ECI = r->RiverXS[i].eci;                      // Energy Coefficient at node I
            ECIP1 = r->RiverXS[i+1].eci;                  // Energy Coefficient at node I+1

            BI = r->RiverXS[i].topW;                      // Top Width at node I
            BIP1 = r->RiverXS[i+1].topW;                  // Top Width at node I+1

            DPDY1 = r->RiverXS[i].centr;                  // Derivative of Wetted Perimeter wrt y at node I
            DPDY2 = r->RiverXS[i+1].centr;                // Derivative of Wetted Perimeter wrt y at node I+1

            SF1 = abs( Q[i] ) * Q[i] / ( KI * KI );
            SF2 = abs( Q[i+1]) * Q[i+1] /( KIP1 * KIP1 );

            FR2T = ECI * r->RiverXS[i].velocity * r->RiverXS[i].velocity * r->RiverXS[i].topW / ( G * ARI );

            if (FR2T >= FR2_TRIG2) fadj = 0;
                else
                if (FR2T <= FR2_TRIG1) fadj = 1.0;
                    else
                    fadj = ( FR2_TRIG2-FR2T ) / ( FR2_TRIG2-FR2_TRIG1 ); //Linearly adjust fadj factor between the trigger levels

            //if (Qw[i].Coord = i)..  Conditional clause for trib inputs (not yet implemented!)

            K = 2 * i + 1;                       // EQN array Index [1,3,5..]

            if ( FR2T < 0.9 ){

                EQN[K][4]= -( ARI + ARIP1 + DTX2 * THETA * ( Q[i+1] - Q[i] ) + C1[i] );

                //Following term modified by JMW to avoid including the
                //spatial derivative of area in the gAy term.
                //TERM1 = DTX2 * THETA *
                //      ( ( ECIP1 * Q[i+1] * Q[i+1] ) / ARIP1 + G * AM * Y[i+1]
                //       -( ECI   * Q[i] *   Q[i]   ) / ARI   - G * AM * Y[i]);
                TERM1 = DTX2  *THETA * ( ECIP1 * Q[i+1] * r->RiverXS[i+1].velocity - ECI * Q[i] * r->RiverXS[i].velocity +
                        G * fadj * ( r->RiverXS[i+1].centr - r->RiverXS[i].centr ) );

                TERM2 = THETA * r->dt * G * ( SF2 * ARIP1 + SF1 * ARI );

                EQN[K+1][4] = -( Q[i] + Q[i+1] + TERM1 + TERM2 + C2[i] );

                DAY1 = BI;
                DAY2 = BIP1;
                EQN[K][0] = DAY1;
                EQN[K][1] = -DTX2 * THETA;
                EQN[K][2] = DAY2;
                EQN[K][3] = DTX2 * THETA;

                //DCDY1 = DCENDY(I,D);
                //DCDY2 = DCENDY(I+1,D1);
                DCDY1 = ARI;
                DCDY2 = ARIP1;

                DSDQ1 = 2 * Q[i]   / (KI   * KI );
                DSDQ2 = 2 * Q[i+1] / (KIP1 * KIP1 );

                TERM1 = DPDY1 * ARI - DAY1 * PERI;
                TERM2 = sqrt(HRI2) * ARI * ARI;
                DSDY1 = Q[i] * abs(Q[i]) / ( KI * KI ) * ( 1.333 * TERM1 / TERM2 - 2 * DAY1 / ARI / sqrt(HRI2) );

                TERM1 = DPDY2 * ARIP1 - DAY2 * PERIP1;
                TERM2 = sqrt(HRIP12) * ARIP1 * ARIP1;
                DSDY2 = Q[i+1] * abs( Q[i+1] ) / ( KIP1 * KIP1 ) *
                            ( 1.333 * TERM1 / TERM2 - 2 * DAY2 / ARIP1 / sqrt(HRIP12) );

                TERM1 = DTX2 * THETA * ( ECI * abs(r->RiverXS[i].velocity) * r->RiverXS[i].velocity
                            * DAY1 - G * DCDY1 );
                TERM2 = G * r->dt * THETA * SF1 * DAY1;

                EQN[K+1][0] = TERM1 + TERM2 + G * r->dt * THETA * ARI * DSDY1;
                EQN[K+1][1] = 1.0 - DTX2 * THETA * 2 * ECI * Q[i] / ARI + G * r->dt * THETA * ARI * DSDQ1;

                TERM1 = -DTX2 * THETA * ( ECIP1 * r->RiverXS[i+1].velocity * abs( r->RiverXS[i+1].velocity )
                            * DAY2 - G * DCDY2 );
                TERM2 = G * r->dt * G * SF2 * DAY2;

                EQN[K+1][2] = TERM1 + TERM2 + THETA * r->dt * G * ARIP1 * DSDY2;
                EQN[K+1][3] = 1.0 + DTX2 * THETA * 2 * ECIP1 * Q[i+1] / ARIP1 + THETA * r->dt * G * ARIP1 * DSDQ2;
                }
            else
                {
                    //Special treatment, in case flow goes supercritical.  In this case
                    //force the flow at the node to be critical depth and keep Q constant
                QM = (Q[i+1]+Q[i])/2;

                xsCritDepth( i, r, Q[i] );       // Calculate critical depth

                EQN[K][1] = -1.0;
                EQN[K][3] = 1.0;
                EQN[K][4] = (Q[i]-Q[i+1]);  //Force continuity of Q
                EQN[K+1][0] = -1.0;
                EQN[K+1][4] = (Y[i] - ( r->eta[i] + r->RiverXS[i].critdepth ) );
                }

            i++;
            }

        // SOLVE SYSTEM OF EQUATIONS
        DF = matsol(NNODES, EQN);

        i = 0;
        SUMM = 0.0;
        while (i <= 2 * NNODES - 1 ){
            SUMM = abs(DF[i]) + SUMM;
            idx = floor(i/2);
            if((i % 2) == 0){                     // Even entries in DF
                Y[idx] = Y[idx] + DF[i];
                r->RiverXS[idx].depth = Y[idx] - r->eta[idx];
            }
            else
            {
                Q[idx] = Q[idx] + DF[i];
                QwCumul[idx] = Q[idx];
            }
            i++;
        }

        iter++;

        if (iter > 1500){
            cout << "Preiss1: Maximum number of iterations exceeded";
            break;
        }
    }
}

vector<double> hydro::matsol(int N, vector<vector<double> > A){

    int i, j, k, inode, M;
    double t1, t2, t3, t4, d;
    vector<double> X, C;

    for (i = 0; i < N * 2; i++){                 // Initialize C array
        C.push_back(0.0);
        X.push_back(0.0);
    }

//Perform first sweep
    C[0] = 0.0;
    C[1] = A[0][4];

    for (inode = 0; inode < N - 1; inode++){
        j = 2 * inode + 1;
        k = j + 1;
        t1 = A[j][0] + A[j][1] * C[k-2];
        t2 = A[j+1][0] + A[j+1][1] * C[k-2];
        t3 = A[j+1][4] - A[j+1][1] * C[k-1];
        t4 = A[j][4] - A[j][1] * C[k-1];
        d  = t1 * A[j+1][3] - t2 * A[j][3];

        if( abs(d) <= 1E-08 )
            cout << "SINGULAR MATRIX --> NO UNIQUE SOLUTION EXISTS";

        C[k] = ( -t1 * A[j+1][2] + t2 * A[j][2]) / d;
        C[k+1] = ( t1 * t3 - t2 * t4 ) / d;
    }

//Perform second sweep
    M = 2 * N - 2;
    X[M] = A[M+1][4];
    X[M+1] = C[M] * X[M] + C[M+1];

    for (inode = N - 1; inode > 0; inode--){
        j = 2 * inode - 1;
        k = j - 1;
        t4 = A[j][4] - A[j][1] * C[k+1];
        d = A[j][0] + A[j][1] * C[k];

        if( abs(d) <= 1E-08 )
        cout << "SINGULAR MATRIX --> NO UNIQUE SOLUTION EXISTS";

        X[k] = ( t4 - ( A[j][2] * X[k+2] + A[j][3] * X[k+3] ) ) / d;
        X[k+1] = C[k] * X[k] + C[k+1];
        }

    return X;

}

// New Routines:   *********************************************************************

void hydro::regimeModel( unsigned int n, RiverProfile *r )
{
    NodeXSObject& XS = r->RiverXS[n];
    NodeGSDObject& f = r->F[n];

    double Tol = 0.00001;
    double Q = QwCumul[n] / r->RiverXS[n].noChannels;
    double test_plus, test_minus = 0;
    double p, p1, p2, p_upper, p_lower = 0;
    double converg, gradient = 0;
    double gradient_1 = 0;
    double gradient_2 = 0;
    double old_width = XS.width;

    p = 3 * pow( Q, 0.5 );

#ifdef DEBUG_REGIME_MODEL
    // DEBUGGING - not for production
    std::ofstream plotf;
    plotf.open("plot_regimeModel.csv");
    double plotmax = p * 4;
    int plotnum = 1000;
    double plotstep = plotmax / static_cast<double>(plotnum);
    for (int ploti = 1; ploti <= plotnum; ploti++) {
        XS.width = ploti * plotstep;
        xsCritDepth( n, r, Q );    // Calculate critical depth
        findStable( n, r );        // Update section data [n.b. bed stress] based on new theta
        XS.xsWilcockTransport(f);  // Work out transport potential
        test_plus = XS.Qb_cap;
        plotf << XS.width << ", " << XS.Qb_cap << std::endl;
    }
    plotf.close();
    // END DEBUGGING
#endif

    XS.width = p * 1.001;
    xsCritDepth( n, r, Q );    // Calculate critical depth
    findStable( n, r );        // Update section data [n.b. bed stress] based on new theta
    XS.xsWilcockTransport(f);  // Work out transport potential
    test_plus = XS.Qb_cap;

    XS.width = p * 0.999;
    xsCritDepth( n, r, Q );       // Calculate critical depth
    findStable( n, r );
    XS.xsWilcockTransport(f);
    test_minus = XS.Qb_cap;

    gradient_1 = test_plus - test_minus;
    p1 = p;

    // Now move in the direction of the gradient
    if (gradient_1 > 0)
        p = p + 0.25 * p;
    else
        p = p - 0.25 * p;

    XS.width = p * 1.001;
    xsCritDepth( n, r, Q );       // Calculate critical depth
    findStable( n, r );
    XS.xsWilcockTransport(f);
    test_plus = XS.Qb_cap;

    XS.width = p * 0.999;
    xsCritDepth( n, r, Q );       // Calculate critical depth
    findStable( n, r );
    XS.xsWilcockTransport(f);
    test_minus = XS.Qb_cap;

    gradient_2 = test_plus - test_minus;
    p2 = p;

    while( gradient_1 / gradient_2 > 0 )
    {
        gradient_1 = gradient_2;
        p1 = p;

        if (gradient_2 > 0)
            p = p + 0.25 * p;
        else
            p = p - 0.25 * p;

        XS.width = p * 1.001;
        xsCritDepth( n, r, Q );       // Calculate critical depth
        findStable( n, r );
        XS.xsWilcockTransport(f);
        test_plus = XS.Qb_cap;

        XS.width = p * 0.999;
        xsCritDepth( n, r, Q );       // Calculate critical depth
        findStable( n, r );
        XS.xsWilcockTransport(f);
        test_minus = XS.Qb_cap;

        gradient_2 = test_plus - test_minus;
        p2 = p;
    }

    p_upper = max( p1, p2 );
    p_lower = min( p1, p2 );
    p = 0.5 * ( p_upper + p_lower );
    converg = ( p_upper - p_lower ) / p;

    while(converg > Tol)
    {
        XS.width = p * 1.001;
        xsCritDepth( n, r, Q );       // Calculate critical depth
        findStable( n, r );
        XS.xsWilcockTransport(f);
        test_plus = XS.Qb_cap;

        XS.width = p * 0.999;
        xsCritDepth( n, r, Q );       // Calculate critical depth
        findStable( n, r );
        XS.xsWilcockTransport(f);
        test_minus = XS.Qb_cap;

        gradient = test_plus - test_minus;

        if ( gradient > 0 )
            p_lower = p;
        else
            p_upper = p;
        p = 0.5 * ( p_upper + p_lower );
        converg = ( p_upper - p_lower ) / p;
    }

    // Update reach geometry, with newly optimsed variables

    // Allow no more than 2% width change per timestep
    // Keep at least 10m channel width

    XS.deltaW = 1 + (p - old_width) / old_width;

/*    if ( n < r->nnodes - 3 )           // Change is weighted by downstream regime results (deltaW) to prevent instabilities
        XS.width = max( 20., XS.width * (r->RiverXS[n+2].deltaW + r->RiverXS[n+1].deltaW + 2 * XS.deltaW) / 4 );
    else
        XS.width = max( 20., XS.width * XS.deltaW);
*/

//    if ( XS.depth < XS.Hmax )
//        XS.bankHeight = XS.Hmax;                      // If equilibrium depth is less than Hmax, XS is a rectangle
//    else
//        XS.bankHeight = XS.depth;                     // MAKE SURE BANKHEIGHT IS CORRECT, HERE

    // Final geometry calcs
    XS.width = old_width;
    energyConserve(n, r);    // Work out depth based on energy considerations

    XS.xsArea();
    XS.xsPerim();
    XS.xsECI(f);
    XS.xsStressTerms(f, bedSlope[n]);
    xsCritDepth( n, r, Q );       // Calculate critical depth


}

void hydro::findStable( unsigned int n, RiverProfile *r )
{
    // Find the stable channel shape for the specified Q and specified bank character
    // Iteratively vary cross-section until tau_bank = bank_crit

    NodeXSObject& XS = r->RiverXS[n];
    NodeGSDObject& f = r->F[n];

    // specify constants and set mu to an equivalent phi
    double phi = 40;
    // double mod_phi = atan( XS.mu * tan( phi * PI / 180) ) * 180 / PI;
    double Tol = 0.00001;
    double deltaX = 0.00001 * phi;
    double tau_star = 0.020;
    double D90 = pow( 2, f.d90 ) / 1000.;
    double converg, bank_crit;
    int iter = 0;
    double b_upper = phi - deltaX;   // set the upper and lower angle limits
    double b_lower = deltaX;

    XS.theta = 0.67 * phi;
    findQ(n, r);    //  Update depth, given new XS.theta. Routine should be called FindY..
    XS.xsStressTerms(f, bedSlope[n]);   // Update bed and bank stresses

    // calculate the bank stability index
    bank_crit = G * RHO * Gs * D90 * tau_star *
              pow( 1 - ( pow( sin ( XS.theta * PI / 180. ), 2) /
              pow( sin( phi * PI / 180. ), 2) ), 0.5 );
    converg = ( XS.Tbank - bank_crit ) / bank_crit;

    while( ( abs( converg ) > Tol ) && ( iter < 50 ) )
    {
        if( converg > 0 )
            b_upper = XS.theta;
        else
            b_lower = XS.theta;

        XS.theta = 0.5 * (b_upper + b_lower);
        findQ(n, r);
        XS.xsStressTerms(f, bedSlope[n]);

        bank_crit = G * RHO * Gs * D90 * tau_star *
                  pow( 1 - ( pow( sin ( XS.theta * PI / 180. ), 2) /
                  pow( sin( phi * PI / 180. ), 2) ), 0.5 );
        converg = ( XS.Tbank - bank_crit ) / bank_crit;

        iter++;
    }
}

void hydro::findQ( unsigned int n, RiverProfile *r ){

    // Find the appropriate depth to satisfy continuity with the imposed flow, QwCumul[n]
    NodeXSObject& XS = r->RiverXS[n];
    NodeGSDObject& f = r->F[n];
    double deltaX = 0.1 * pow ( QwCumul[n], 0.3);
    double tol = 0.005;
    double converg = 1;
    double testQ;
    int iter = 0;

    // Guess the initial channel depth, set up variables
    XS.depth = 0.3 * pow ( QwCumul[n], 0.3);

    // First test
    XS.xsArea();
    XS.xsPerim();
    XS.xsECI(f);
    XS.xsStressTerms(f, bedSlope[n]);        // Note empirical velocity estimate, not based on energy considerations

    testQ = XS.flow_area[2] * XS.velocity;
    converg = (testQ - QwCumul[n]) / QwCumul[n];

    while( ( abs(converg) > tol) && ( iter < 250 ) ){
        if ( converg > 0 )
            XS.depth -= converg * XS.depth;
        else
            XS.depth += abs(converg) * XS.depth;

        XS.xsArea();
        XS.xsPerim();
        XS.xsECI(f);
        XS.xsStressTerms(f, bedSlope[n]);

        testQ = XS.flow_area[2] * XS.velocity;
        converg = (testQ - QwCumul[n]) / QwCumul[n];

        iter++;
    }
}

void hydro::setRegimeWidth(RiverProfile *r)
{

    // Adjust channel regime one cross-section at a time, marching upstream

    double deltaArea, oldArea = 0.;            // Change in reach cross-section area
    double deltaEta = 0.;
    double reachDrop = 0.;                     // Drop in elevation over reach river length
    double oldBankHeight = 0.;
    double aspect = 0.;
//    double Tol = 50.;                          // Maximum allowed channel aspect (w/d)

    oldBankHeight = r->RiverXS[regimeCounter].bankHeight;
    oldArea = r->RiverXS[regimeCounter].flow_area[2];

    r->RiverXS[regimeCounter].noChannels = 1;  // Starts on assumption of just one channel
    regimeModel( regimeCounter, r );           // Call to the Regime Functions.
    aspect = r->RiverXS[regimeCounter].width / r->RiverXS[regimeCounter].depth;

    if ( (r->counter > 260 ) )                 // Update floodplain volume
    {
        deltaArea = oldArea - r->RiverXS[regimeCounter].flow_area[2];       // Change induced by floodplain erosion
        deltaEta = r->RiverXS[regimeCounter].bankHeight - oldBankHeight;    // Change induced by channel aggr/degr
        deltaEta += deltaArea / r->RiverXS[regimeCounter+1].fpWidth;        // Total lateral change is a product of the two

        reachDrop = bedSlope[regimeCounter] * r->dx * r->RiverXS[regimeCounter].chSinu;
        // new sinuosity
        r->RiverXS[regimeCounter].chSinu *= ( (reachDrop + deltaEta ) / reachDrop );
        if ( r->RiverXS[regimeCounter].chSinu < 1 )
            r->RiverXS[regimeCounter].chSinu = 1;
        if ( r->RiverXS[regimeCounter].chSinu > 2.6 )
            r->RiverXS[regimeCounter].chSinu = 2;
        // readjust downstream elevation to account for gains or losses during width adjustment
        // [ old area - new area ] - positive value if material removed, channel widening.

    }
    regimeCounter --;
    if (regimeCounter < 2)
        regimeCounter = (r->nnodes-2);

}
