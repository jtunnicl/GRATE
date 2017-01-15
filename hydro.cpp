/*******************
 *
 *
 *  GRATE 9
 *
 *  Hydraulic parameters, backwater and channel regime algorithms
 *
 *
 *
*********************/

#include "hydro.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <fstream>
using namespace std;

#define PI 3.14159265
#define G 9.80665
#define RHO 1000  // water density
#define Gs 1.65   // submerged specific gravity

hydro::hydro(RiverProfile *r)
{
    preissTheta = 0.7;
    hydUpw = 0.3;
    regimeCounter = (r->nnodes-2);

    initHydro(r->nnodes);
}

void hydro::initHydro(int nodes)
{
    ifstream inHydroFile;
    double currentCoord = 0.;
    int v2 = 0; // Yr, Mt, Day, Hr, Mn, Sec, GRP
    int v3 = 0;
    int v4 = 0;
    int v5 = 0;
    int v6 = 0;
    int v7 = 0;
    int v9 = 0;
    double v1 = 0.;
    double v8 = 0.;                    // Coord, Qw
    QDateTime NewDate;
    vector< TS_Object > tmp;
    TS_Object NewEntry;

    NewDate.setDate(QDate(2000,1,1));
    NewDate.setTime(QTime(0,0,0));
    NewEntry.date_time = NewDate;
    NewEntry.Q = 0;
    NewEntry.Coord = 0;
    NewEntry.GRP = 1;
    //tmp.push_back(NewEntry);

    inHydroFile.open("hydro_series.dat");
    if (!inHydroFile) cout << "hydro file couldn't open properly" << endl;

    while(!inHydroFile.eof() )
    {
        if(inHydroFile >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9)
        {
            NewDate.setDate(QDate(v2, v3, v4));
            NewDate.setTime(QTime(v5, v6, v7));
            NewEntry.date_time = NewDate;
            NewEntry.Q = v8;
            NewEntry.Coord = v1;
            NewEntry.GRP = v9;

            if (v1 > currentCoord)        // Have we moved to a new source coordinate?
            {
                Qw.push_back( tmp );
                tmp.clear();
                currentCoord = v1;
                tmp.push_back(NewEntry);  // Start new tmp
            }
            else
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
    bool bQuasiNormal = 0;
    int lastNode = r->nnodes-1;

    setQuasiSteadyNodalFlows(r);

    // update stats for last (down-stream) node:
    r->RiverXS[lastNode].xsGeom();               // update area
    r->RiverXS[lastNode].meanVeloc = QwCumul[lastNode] / r->RiverXS[lastNode].xsFlowArea[2];
    r->RiverXS[lastNode].xsCentr();              // update centr
    r->RiverXS[lastNode].xsECI(r->F[lastNode]);                // update eci
    Fr2[lastNode] = r->RiverXS[lastNode].eci * r->RiverXS[lastNode].meanVeloc
                  * r->RiverXS[lastNode].meanVeloc / ( g * r->RiverXS[lastNode].maxDepth );


    // update bed slope array
    for (int n = r->nnodes-2; n > 0 ; n--)
        bedSlope[n] = (hydUpw * ( r->eta[n-1] - r->eta[n] ) / r->dx
               + (1 - hydUpw) * ( r->eta[n] - r->eta[n+1] ) / r->dx)
                   / r->RiverXS[n].chSinu; // note inclusion of sinuosity

    // Set slope at upstream/downstream boundary nodes
    bedSlope[0] = ( r->eta[0] - r->eta[1] ) / r->dx;
    bedSlope[r->nnodes-1] = ( r->eta[r->nnodes-2] - r->eta[r->nnodes-1] ) / r->dx;

    // Boundary nodes: fixed or computed (default)
    quasiNormal(0, r);
    //quasiNormal(lastNode, r);

    r->RiverXS[lastNode].maxDepth = 0.3 * pow( QwCumul[lastNode],0.3 );
    r->RiverXS[lastNode].wsl = r->eta[lastNode] + r->RiverXS[lastNode].maxDepth;

    for (int n = r->nnodes-2; n > 0 ; n--)
    {
        xsCritDepth( n, r, QwCumul[n] );       // Calculate critical depth

        // Initial guess at depth
        r->RiverXS[n].maxDepth = 0.3 * pow( QwCumul[n],0.3 );

        // Flat profile if the bed steps up
        if ( bedSlope[n] < 0 )
            r->RiverXS[n].maxDepth = r->RiverXS[n+1].maxDepth - bedSlope[n] * r->dx;

        r->RiverXS[n].xsGeom();                // update area
        r->RiverXS[n].meanVeloc = QwCumul[n] / r->RiverXS[n].xsFlowArea[2];
        r->RiverXS[n].xsCentr();               // update centr
        r->RiverXS[n].xsECI(r->F[n]);          // update eci
        Fr2[n] = r->RiverXS[n].eci * r->RiverXS[n].meanVeloc *
                r->RiverXS[n].meanVeloc / ( g * r->RiverXS[n].maxDepth );

        if ( ( Fr2[n] < FrN2 ) || ( bedSlope[n] <= 0 ) || ( n == 0 ) ) // Not super-crit; use energy eqn
            iret = energyConserve(n, r);
        else
        {                                                              // else - recalculate using quasi-normal assumption
            if (bQuasiNormal == 0)
            {
                iret = quasiNormal(n+1, r);                            // Recalculate i+1'th node   <JMW 20080303>
                // if ((iret > 0) || (n < (r->nnodes-3)))
                //     r->RiverXS[n+1].maxDepth = r->RiverXS[n+2].maxDepth;
            }

            iret = quasiNormal(n, r);

            if (iret > 0)
                r->RiverXS[n].maxDepth = r->RiverXS[n+1].maxDepth;

            bQuasiNormal = 1;
        };

        if ( ( iret > 0 ) || ( r->RiverXS[n].maxDepth < r->RiverXS[n].critDepth ) )
            r->RiverXS[n].maxDepth =  r->RiverXS[n].critDepth;

        if ( ( r->RiverXS[n].maxDepth > 0 ) && ( bedSlope[n] > 0 ) )
            r->RiverXS[n].ustar = sqrt( 9.81 * r->RiverXS[n].maxDepth * bedSlope[n] );
        else
            r->RiverXS[n].ustar = 1e-3;

        r->RiverXS[n].wsl = r->eta[n] + r->RiverXS[n].maxDepth;       // Update water surface level at n
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
    for ( j = 1; j < QwCumul.size(); j++ )
    {
        QwCumul[j] = QwCumul[j-1];
        if ( i < (Qw.size()-1) && ( r->xx[j] > Qw[i + 1][0].Coord ) )
        {
            i++;
            QwCumul[j] += Qw_Ct[i];
        }
    }
}

void hydro::xsCritDepth(int n, RiverProfile *r, double Q){

    double iter,i,k;
    double y_star, y_c1, y_c2, y_c3;                 // Following Chaudhry's technique, Section 3-7, in 2nd Ed. 2008
    double Tol = 0.001;
    double converg = 1;
    double upper, lower, Cmax;
    vector<double> y_r_test = {1.0001, 1.0005, 1.001, 1.005, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.1, 1.2, 1.5, 2, 3};    // 16 elements
    vector<double> C_result = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};

    // Relative quantities
    // y_f = r-> RiverXS[n].xsBankHt
    double y_r = r->RiverXS[n].maxDepth / r->RiverXS[n].xsBankHt;
    double b_r = r->RiverXS[n].fpWidth / r->RiverXS[n].xsBedWidth;
    double b_f = r->RiverXS[n].fpWidth / r->RiverXS[n].xsBankHt;
    double n_r = 0.05;                         // Needs Keulegan equivalent
    double m = 1 / ( 1 + 2 * n_r *
            pow( r->RiverXS[n].xsFlowArea[1] / r->RiverXS[n].xsFlowArea[0], 1.6667 ) *
            pow( r->RiverXS[n].xsFlowPerim[0] / r->RiverXS[n].xsFlowPerim[1], 0.6667 ) );

    k = G * pow( r->RiverXS[n].xsBedWidth, 2 ) * pow ( r->RiverXS[n].xsBankHt, 3 ) / pow ( QwCumul[n], 2 );

    double C = 1 / ( y_r + 2 * b_r * ( y_r - 1 ) ) * ( pow ( m / y_r, 2 ) + pow ( ( 1 - m ) / ( y_r - 1 ), 2 ) * ( 1 / 2 * b_r ) )
            + 2 * m * ( 1 -  m ) / 3 * ( y_r + 2 * b_r * ( y_r - 1 ) ) * ( 5 / ( y_r * ( y_r - 1 )) - 2 / ( b_f + y_r - 1 ) ) *
            ( ( m / y_r ) - ( ( 1 - m ) / ( y_r - 1 ) ) * 1 / 2 * b_r );

    // Compute critical depth, given a flow

    if ( k < 1 )                           // Only one critical depth, and it is greater than the floodplain depth, y_f;
                                           // Critical depth cannot occur when the flow is only in the main channel
    {
        iter = 0;
        y_r = 1.1;                         // First guess
        upper = 5;                         // Bounds on likely y_r solution
        lower = 1.0001;

        while ( ( converg > Tol ) && ( iter < 50 ) )
        {
            if ( converg > 0 )
                upper = y_r;
            else
                lower = y_r;
            y_r = 0.5 * ( upper + lower );

            y_star = (2 * b_r) / (2 * b_r + 1) + 1 / C * ( 2 * b_r + 1 ) *
                    ( pow ( m / r->RiverXS[n].maxDepth, 2 ) + pow ( ( ( 1 - m ) / ( y_r -1) ), 2 ) *
                      ( 1 / 2 * b_r ) ) + ( 2 * m * ( 1 - m ) ) / 3 * C * ( 2 * b_r + 1 ) *
                    ( 5 / y_r * ( y_r - 1 ) - 2 / ( b_f + y_r - 1 ) ) * ( m / y_r -  ( 1 - m ) / (y_r - 1 ) * 1 / ( 2 * b_r ));

            converg = abs( y_r - y_star );
            iter++;
        }
        r->RiverXS[n].critDepth = y_r * r->RiverXS[n].xsBankHt;
    }

    else   // ( k >= 1 )
    {
        for ( i = 1; i < y_r_test.size(); i++ )                 //  Generate C vs y_r relation
            C_result[i] = 1 / ( y_r_test[i] + 2 * b_r * ( y_r_test[i] - 1 ) ) * ( pow ( m / y_r_test[i], 2 ) + pow ( ( 1 - m ) / ( y_r_test[i] - 1 ), 2 ) * ( 1 / 2 * b_r ) )
                            + 2 * m * ( 1 -  m ) / 3 * ( y_r_test[i] + 2 * b_r * ( y_r_test[i] - 1 ) ) * ( 5 / ( y_r_test[i] * ( y_r_test[i] - 1 )) - 2 / ( b_f + y_r_test[i] - 1 ) ) *
                            ( ( m / y_r_test[i] ) - ( ( 1 - m ) / ( y_r_test[i] - 1 ) ) * 1 / 2 * b_r );
        Cmax = *max_element( C_result.begin(), C_result.end() );     //  A better scheme to find the maximum would be desirable, here!

        if ( k > Cmax )
            r->RiverXS[n].critDepth = pow ( pow ( QwCumul[n], 2 ) / ( G * pow( r->RiverXS[n].xsBedWidth, 2 )), 0.334 );
        else
        {
            y_c1 = pow ( pow ( QwCumul[n], 2 ) / ( G * pow( r->RiverXS[n].xsBedWidth, 2 )), 0.334 );   //  Found Yc1

            y_r = 1.1;                         // First guess
            upper = 5;                         // Bounds on likely y_r solution
            lower = 1.0001;
            iter = 0;

            while ( ( converg > Tol ) && ( iter < 50 ) )
            {
                if ( converg > 0 )
                    upper = y_r;
                else
                    lower = y_r;
                y_r = 0.5 * ( upper + lower );

                y_star = (2 * b_r) / (2 * b_r + 1) + 1 / C * ( 2 * b_r + 1 ) *
                         ( pow ( m / r->RiverXS[n].maxDepth, 2 ) + pow ( ( ( 1 - m ) / ( y_r -1) ), 2 ) *
                         ( 1 / 2 * b_r ) ) + ( 2 * m * ( 1 - m ) ) / 3 * C * ( 2 * b_r + 1 ) *
                         ( 5 / y_r * ( y_r - 1 ) - 2 / ( b_f + y_r - 1 ) ) * ( m / y_r -  ( 1 - m ) / (y_r - 1 ) * 1 / ( 2 * b_r ));

                converg = abs( y_r - y_star );
                iter++;
            }

            y_c3 = y_r * r->RiverXS[n].xsBankHt;              // Found Yc3

            y_r = 1.001;
            C = 1.5;

            while ( C < k )
            {
                C = 1 / ( y_r + 2 * b_r * ( y_r - 1 ) ) * ( pow ( m / y_r, 2 ) + pow ( ( 1 - m ) / ( y_r - 1 ), 2 ) * ( 1 / 2 * b_r ) ) +
                        2 * m * ( 1 -  m ) / 3 * ( y_r + 2 * b_r * ( y_r - 1 ) ) * ( 5 / ( y_r * ( y_r - 1 )) - 2 / ( b_f + y_r - 1 ) ) *
                        ( ( m / y_r ) - ( ( 1 - m ) / ( y_r - 1 ) ) * 1 / 2 * b_r );
                y_r *= 1.01;               // Again, fairly unsophisticated marching solution.
            }

            y_c2 = y_r * r->RiverXS[n].xsBankHt;              // Founds Yc2

            if ( r->RiverXS[n].maxDepth > r->RiverXS[n].xsBankHt )
                r->RiverXS[n].critDepth = y_c3;
            else if ( r->RiverXS[n].maxDepth == r->RiverXS[n].xsBankHt )
                r->RiverXS[n].critDepth = y_c2;
            else
                r->RiverXS[n].critDepth = y_c1;
        }
    }


    int it = 0;
    int itmax = 50;                        // Max iterations
    double toler = 0.0005;                 // Convergence criteria
    double ymin = 0.15;
    double ymax = r->RiverXS[n].CHList[0].bankHeight + 1.5;           //
    double dy = 0;
    double y1 = 0;
    double y2 = 0;                         // dummy variables

    // Make sure ymax is subcritical (ff <= 0); keep increasing ymax until this is so.

    double ff = 1;
    while (ff > 0)
    {
        if (it > 0)
            ymax *= 1.5;

        r->RiverXS[n].maxDepth = ymax;
        r->RiverXS[n].xsGeom();                          // Update statistics at node
        r->RiverXS[n].hydRadius = r->RiverXS[n].xsFlowArea[2] / r->RiverXS[n].xsFlowPerim[2];
        r->RiverXS[n].xsCentr();                         // re-calc topW

        ff = Q / r->RiverXS[n].xsFlowArea[2] / sqrt( 9.81 * r->RiverXS[n].hydRadius ) - 1.0;
        it++;
        if (it > itmax)
        {
            // cout << "Unable to initialise max depth for critical depth calculation at xc = \n";
            exit(1);
        }
    }

    y1 = ( ymin + ymax ) / 2;                   //Initial trial depth

    // Solve by bisection

    while (it < itmax)
    {
        r->RiverXS[n].maxDepth = y1;
        r->RiverXS[n].xsGeom();
        r->RiverXS[n].hydRadius = r->RiverXS[n].xsFlowArea[2] / r->RiverXS[n].xsFlowPerim[2];
        r->RiverXS[n].xsCentr();

        ff = Q / r->RiverXS[n].xsFlowArea[2] / sqrt( 9.81 * r->RiverXS[n].hydRadius ) - 1.0;
        if (ff < 0)
            ymax = y1;
        else
            ymin = y1;

        y2 = ( ymin + ymax ) / 2;
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
    r->RiverXS[n].critDepth = y2;
}

int hydro::energyConserve(int n, RiverProfile *r)
    {

    // Energy conservation between the two nodes - using bisection algorithm
    // developed by JMW --> v.3.3 energy_conserve2()

    double ff;                  // Objective function
    double h1, h2, hu2;         // Straddle points in bisection
    double Sf, Sf2, Sfx;        // Friction slopes
    double Vhd;                 // Velocity head at downstream node
    double Vhu;                 // Velocity head at upstream node
    int flag, iter, itermax;    // Return flag, iteration counter, counter max
    double error;               // error during iterations
    double qm, km;              // Mean Qw and K (conveyance) between nodes

    NodeXSObject& XSu = r->RiverXS[n];           // "upstream" cross-secction  [n]--> Objective
    NodeXSObject& XSd = r->RiverXS[n+1];         // "downstream" cross-section [n+1]--> already computed

    flag = 0;                   // Function flag to be returned
    itermax = 300;

    XSu.xsGeom();               // update area
    r->RiverXS[n].meanVeloc = QwCumul[n] / XSu.xsFlowArea[2];

    XSu.xsECI(r->F[n]);       // update eci

    XSd.xsGeom();               // update area
    r->RiverXS[n+1].meanVeloc = QwCumul[n+1] / XSd.xsFlowArea[2];

    XSd.xsECI(r->F[n+1]);         // update eci

    //Sf2 = bedSlope[n+1];       // slope gradient between n and n+1
    Vhd = XSd.eci * r->RiverXS[n+1].meanVeloc * r->RiverXS[n+1].meanVeloc / (2 * 9.81);
                                 // Velocity head
    //Sf2 = pow( ( QwCumul[n+1] / XSd.k_mean ), 2);   // Dingman 9B2.4
    Sf2 = QwCumul[n+1] * QwCumul[n+1] / ( XSd.k_mean * XSd.k_mean);

    h1 = r->RiverXS[n].critDepth;          // Bisection: lower straddle point
    h2 = max(10 * r->RiverXS[n].critDepth, (XSd.maxDepth + bedSlope[n+1] * r->dx) * 2 );  // upper straddle point

    ff = -1;
    while (ff <= 0)
    {
        XSu.maxDepth = h2;             // Update section data based on new depth
        XSu.xsGeom();
        r->RiverXS[n].meanVeloc = QwCumul[n] / XSu.xsFlowArea[2];
        if ( XSu.maxDepth > 0 )
            XSu.xsECI(r->F[n+1]);
        Sf = QwCumul[n] / XSu.k_mean;
        Vhu = XSu.eci * r->RiverXS[n].meanVeloc * r->RiverXS[n].meanVeloc / (2 * 9.81);

        ff = (XSu.maxDepth + Vhu) - (XSd.maxDepth + Vhd) + ( (bedSlope[n+1] + bedSlope[n]) / 2 - Sf) * r->dx;

        h2 = 2 * XSu.maxDepth;
    }

    h2 = XSu.maxDepth;
    XSu.maxDepth = (h1 + h2) / 1.5;

    error = 1;
    iter = 0;
    while (error > 5e-4)
    {
        XSu.xsGeom();    // Update section data based on new depth
        r->RiverXS[n].meanVeloc = QwCumul[n] / XSu.xsFlowArea[2];
        if ( XSu.maxDepth > 0 )
            XSu.xsECI(r->F[n+1]);

        //Sf1 = Sf2;        // Initial approximations
        Sf = Sf2;

        Vhu = XSu.eci * r->RiverXS[n].meanVeloc * r->RiverXS[n].meanVeloc / 2 * 9.81;

        if (iter > 1)
        {
            qm = (QwCumul[n] + QwCumul[n+1] ) / 2;
            km = (XSu.k_mean + XSd.k_mean) / 2;
            Sfx = qm / km;
            Sf = Sfx * Sfx;
        }

        ff = (XSu.maxDepth + Vhu) - (XSd.maxDepth + Vhd) + (bedSlope[n] - Sf) * r->dx;

        if (ff > 0)
            h2 = XSu.maxDepth;
        else
            h1 = XSu.maxDepth;

        if (h2 > XSu.critDepth)
        {
            hu2 = (h1 + h2) / 2;
            error = abs(hu2 - XSu.maxDepth) / XSu.maxDepth;
            XSu.maxDepth = hu2;
        }
        else
        {
            //Limit XSu.maxDepth to critical depth
            XSu.maxDepth = XSu.critDepth;
            break;
        }

        iter ++;

        if (iter > itermax)
        {
            cout << "energy_conserve: std step backwater calculation failed to converge \n";
            exit(1);
            flag = 8;
        }

        if (XSu.maxDepth < 0)
        {
            cout << "energy_conserve: negative depth results \n";
            exit(1);
            flag = 16;
        }
    }
    r->RiverXS[n] = XSu;
    return flag;
}

int hydro::quasiNormal(int n, RiverProfile *r){

    double ff, fp;                             // Objective function and derivative
    int iter, maxiter;
    double error;

    NodeXSObject& XS = r->RiverXS[n];          // Reach cross-section object to be calculated
    NodeGSDObject& f = r->F[n];
    f.norm_frac();
    f.dg_and_std();

    error = 1;
    iter = 0;
    maxiter = 900;

    while ( error > 0.0001 )
    {
        XS.xsGeom();                           // Update section data based on new depth
        XS.xsCentr();
        if ( XS.hydRadius > 0 )
            XS.xsECI(f);

        ff = QwCumul[n] / XS.topW - XS.hydRadius * sqrt( 9.81 * abs( XS.hydRadius ) * bedSlope[n] ) * XS.omega;

        fp = -2.5 * sqrt( 9.81 * abs( XS.hydRadius ) * bedSlope[n] )
                * ( 1.5 * log(11.0 * abs ( XS.hydRadius ) / XS.rough) + 1.0 );

        error = -ff / fp;
        XS.maxDepth += error / 2;
        error = abs(error / XS.maxDepth);
        ++iter;

        if (iter> maxiter)
        {
            cout << "Iteration Count exceeded in routine quasiNormal at " << n << "\n";
            return 8;
            exit(1);
        }
    }

    XS.xsGeom();    // Update section data based on new depth
    XS.xsCentr();
    XS.xsECI(f);

    return 0;
}

void hydro::fullyDynamic(RiverProfile *r){

    int i, j, idx, iflag, iter, fadj, NNODES;
    double ARI, ARIP1, KI, KIP1, BI, BIP1, ECI, ECIP1, HRI2;
    double AM, DTX2, SF1, SF2, SUMM, TOL, QM, THETA, Ybc;
    double DAY1, DAY2, DCDY1, DCDY2, DSDQ1, DSDQ2, DSDY1, DSDY2;
    double PERI, PERIP1, HRIP12, DPDY1, DPDY2;
    double TERM1, TERM2, TERM3;
    double FR2T, FD_FR_MIN, FD_FR_MAX, FR2_TRIG1, FR2_TRIG2;            // Trigger levels for transition to critical flow
    float K;

    vector<double> Q, Y, C1, C2, C2A, DF, tmp;
    vector<vector<double>> EQN;

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
        Y[i] = r->eta[i] + r->RiverXS[i].maxDepth; // W.S. Elevation
    }

    iflag = 0;
    THETA = preissTheta;                         // Preissmann Weighting Coefficient
    TOL = 0.001;                                 // Tolerance for interactions
    NNODES = r->nnodes;

    quasiNormal(0, r);
    quasiNormal(NNODES-1, r);
    Ybc = r->eta[NNODES-1] + 2.2; // r->RiverXS[NNODES-1].maxDepth;   // d/s boundary condition

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
            r->RiverXS[i].xsGeom();
            r->RiverXS[i].meanVeloc = Q[i] / r->RiverXS[i].xsFlowArea[2];
            r->RiverXS[i].xsCentr();
            r->RiverXS[i].xsECI(r->F[i]);
        }

        r->RiverXS[i+1].xsGeom();                         // update area, cross-section params for d/s node
        r->RiverXS[i+1].meanVeloc = Q[i+1] / r->RiverXS[i+1].xsFlowArea[2];
        r->RiverXS[i+1].xsCentr();
        r->RiverXS[i+1].xsECI(r->F[i+1]);

        ARI = r->RiverXS[i].xsFlowArea[2];                 // Statement flow area
        ARIP1 = r->RiverXS[i+1].xsFlowArea[2];
        AM = ( ARI + ARIP1 ) / 2;

        KI = r->RiverXS[i].k_mean;
        KIP1 = r->RiverXS[i+1].k_mean;

        ECI = r->RiverXS[i].eci;
        ECIP1 = r->RiverXS[i+1].eci;

        DTX2 = 2 * r->dt / r->dx;

        if (i == 0) Q[i] = QwCumul[0];           //  Upstream boundary condition

        FR2T = ECI * r->RiverXS[i].meanVeloc * r->RiverXS[i].meanVeloc * r->RiverXS[i].topW / ( G * ARI );

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
                ( ECIP1 * Q[i+1] * r->RiverXS[i+1].meanVeloc - ECI * Q[i] * r->RiverXS[i].meanVeloc +
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
                r->RiverXS[i].xsGeom();
                r->RiverXS[i].meanVeloc = Q[i] / r->RiverXS[i].xsFlowArea[2];
                r->RiverXS[i].hydRadius = r->RiverXS[i].xsFlowArea[2] / r->RiverXS[i].xsFlowPerim[2];
                r->RiverXS[i].xsCentr();
                r->RiverXS[i].xsECI(r->F[i]);
            }

            r->RiverXS[i+1].xsGeom();                     // update area, cross-section params for d/s node
            r->RiverXS[i+1].meanVeloc = Q[i+1] / r->RiverXS[i+1].xsFlowArea[2];
            r->RiverXS[i+1].hydRadius = r->RiverXS[i+1].xsFlowArea[2] / r->RiverXS[i+1].xsFlowPerim[2];
            r->RiverXS[i+1].xsCentr();
            r->RiverXS[i+1].xsECI(r->F[i+1]);

            ARI = r->RiverXS[i].xsFlowArea[2];             // Statement flow area
            ARIP1 = r->RiverXS[i+1].xsFlowArea[2];
            AM = ( ARI + ARIP1 ) / 2;

            PERI = r->RiverXS[i].xsFlowPerim[2];           // Wetted Perimeter at node I
            PERIP1 = r->RiverXS[i+1].xsFlowPerim[2];       // Wetted Perimeter at node I+1

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

            FR2T = ECI * r->RiverXS[i].meanVeloc * r->RiverXS[i].meanVeloc * r->RiverXS[i].topW / ( G * ARI );

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
                TERM1 = DTX2  *THETA * ( ECIP1 * Q[i+1] * r->RiverXS[i+1].meanVeloc - ECI * Q[i] * r->RiverXS[i].meanVeloc +
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

                TERM1 = DTX2 * THETA * ( ECI * abs(r->RiverXS[i].meanVeloc) * r->RiverXS[i].meanVeloc
                            * DAY1 - G * DCDY1 );
                TERM2 = G * r->dt * THETA * SF1 * DAY1;

                EQN[K+1][0] = TERM1 + TERM2 + G * r->dt * THETA * ARI * DSDY1;
                EQN[K+1][1] = 1.0 - DTX2 * THETA * 2 * ECI * Q[i] / ARI + G * r->dt * THETA * ARI * DSDQ1;

                TERM1 = -DTX2 * THETA * ( ECIP1 * r->RiverXS[i+1].meanVeloc * abs( r->RiverXS[i+1].meanVeloc )
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
                EQN[K+1][4] = (Y[i] - ( r->eta[i] + r->RiverXS[i].critDepth ) );
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
                r->RiverXS[idx].maxDepth = Y[idx] - r->eta[idx];
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
            break;
            cout << "Preiss1: Maximum number of iterations exceeded";
        }
    }
}

vector<double> hydro::matsol(int N, vector<vector<double>> A){

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

void hydro::regimeModel(int n, RiverProfile *r )
{
    int ch_idx = 1;
    NodeCHObject& CH = r->RiverXS[n].CHList[ch_idx];
    double Tol = 0.00001;
    double Q = QwCumul[n] * CH.QProp;
    double test_plus, test_minus = 0;
    double p, p1, p2, p_upper, p_lower = 0;
    double converg, gradient = 0;
    double gradient_1 = 0;
    double gradient_2 = 0;

    p = 4 * pow( Q, 0.5 );
    CH.width = p * 1.001;
    energyConserve( n, r );                 // Update section data based on new theta
    findStable( n, ch_idx, r );
    test_plus = CH.Qb_cap;

    CH.width = p * 0.999;
    energyConserve( n, r );
    findStable( n, ch_idx, r );
    test_minus = CH.Qb_cap;
    gradient_1 = test_plus - test_minus;
    p1 = p;

    // Now move in the direction of the gradient
    if (gradient_1 > 0)
        p = p + 0.25 * p;
    else
        p = p - 0.25 * p;
    CH.width = p * 1.001;
    energyConserve( n, r );
    findStable( n, ch_idx, r );
    test_plus = CH.Qb_cap;

    CH.width = p * 0.999;
    energyConserve( n, r );
    findStable( n, ch_idx, r );
    test_minus = CH.Qb_cap;

    gradient_1 = test_plus - test_minus;
    p2 = p;


    while( gradient_1 / gradient_2 > 0 )
    {
        gradient_1 = gradient_2;
        p1 = p;
        if (gradient_1 > 0)
            p = p + 0.25 * p;
        else
            p = p - 0.25 * p;
        CH.width = p * 1.001;
        energyConserve( n, r );
        findStable( n, ch_idx, r );
        test_plus = CH.Qb_cap;

        CH.width = p * 0.999;
        energyConserve( n, r );
        findStable( n, ch_idx, r );
        test_minus = CH.Qb_cap;

        gradient_2 = test_plus - test_minus;
        p2 = p;
    }

    p_upper = max( p1, p2 );
    p_lower = min( p1, p2 );
    p = 0.5 * ( p_upper + p_lower );
    converg = ( p_upper - p_lower ) / p;

    while(converg > Tol)
    {
        CH.width = p * 1.001;
        CH.chGeom();
        findStable( n, ch_idx, r );
        test_plus = CH.Qb_cap;

        CH.width = p * 0.999;
        energyConserve( n, r );
        findStable( n, ch_idx, r );
        test_minus = CH.Qb_cap;

        gradient = test_plus - test_minus;
        if ( gradient > 0 )
            p_lower = p;
        else
            p_upper = p;
        p = 0.5 * ( p_upper + p_lower );
        converg = ( p_upper - p_lower ) / p;
    }

    r->RiverXS[n].CHList[ch_idx].width = p;
    r->RiverXS[n].CHList[ch_idx].bankHeight = r->RiverXS[n].CHList[ch_idx].Hmax +
            sin( r->RiverXS[n].CHList[ch_idx].theta * PI / 180 ) *
            ( (r->RiverXS[n].CHList[ch_idx].b2b - r->RiverXS[n].CHList[ch_idx].width) / 2 );
    energyConserve( n, r );
    findStable( n, ch_idx, r );
}

void hydro::channelState( int n, int ch_idx, RiverProfile *r )
{
    // Compute stresses, transport capacity
    NodeCHObject& CH = r->RiverXS[n].CHList[ch_idx];
    double X, arg;
    double SFbank = 0;
    double tau_star_ref, tau_ref, totstress, W_star;

    NodeGSDObject& f = r->F[n];
    //double D84 = pow( 2, f.d84 ) / 1000;
    double D50 = pow( 2, f.dsg ) / 1000;

    NodeGSDObject bankMaterial = r->storedf[n][r->ntop[n]];
    bankMaterial.dg_and_std();
    float theta_rad = CH.theta * PI / 180;

    // use Ferguson 2007 to calculate the stream velocity
    // Res = a1 * a2 * ( CH.hydRadius / D50 ) / pow( ( pow( a1, 2 ) + pow( a2, 2 ) *
    //    pow(( CH.hydRadius / D84 ),( 5 / 3 ))),( 1 / 2 ));
    // use the Keulegan Equation
    // Res = (1/0.4)*log(12.2*R/(D84))
    // CH.velocity = Res * pow((G * CH.hydRadius * bedSlope[n]),(1/2));

    // use the equations from Knight and others to partition stress

    arg =  -1.4026 * log10( CH.width / ( CH.flowPerim - CH.width ) + 1.5 ) + 0.247;
    SFbank = pow ( 10.0 , arg );    // partioning equation, M&Q93 Eqn.8, E&M04 Eqn.2
    totstress = G * RHO * CH.depth * bedSlope[n];
    CH.Tbed =  totstress * (1 - SFbank) *
            ( CH.b2b / (2 * CH.width) + 0.5 );           // bed_str = stress acting on the bed, M&Q93 Eqn.10, E&M04 Eqn.4
    CH.Tbank =  totstress * SFbank *
            ( CH.b2b + CH.width ) * sin( theta_rad ) / (4 * CH.depth );

    // estimate the largest grain that the flow can move
    CH.comp_D = CH.Tbed / (0.02 * G * RHO * Gs );

    // estimate the division between key stones and bed material load
    //   (this corresponds to the approximate limit of full mobility)
    CH.K = CH.Tbed / (0.04 * G * RHO * Gs);

    // use Wilcock and Crowe to estimate the sediment transport rate
    tau_star_ref = 0.021 + 0.015 * exp (-20 * f.sand_pct);
    tau_ref = tau_star_ref * G * RHO * Gs * D50;
    X = CH.Tbed / tau_ref;

    if (X < 1.35)
        W_star = 0.002 * pow( X, 7.5 );
    else
        W_star = 14 * pow( ( 1 - ( 0.894 / pow( X, 0.5 ) ) ), ( 4.5 ) );

    CH.Qb_cap = CH.width * ( W_star / ( 1.65 * G ) ) * pow( ( CH.Tbed / RHO ), ( 3 / 2 ) );
}

void hydro::findStable( int n, int ch_idx, RiverProfile *r )
{
    // find the stable channel shape for the specified Q and relative bank strength, mu
    NodeCHObject& CH = r->RiverXS[n].CHList[ch_idx];
    NodeGSDObject& f = r->F[n];

    // specify constants and set mu to an equivalent phi
    double phi = 40;
    // double mod_phi = atan( CH.mu * tan( phi * PI / 180) ) * 180 / PI;
    double Tol = 0.00001;
    double deltaX = 0.00001 * phi;
    double tau_star = 0.02;
    double D90 = pow( 2, f.d90 ) / 1000;
    double converg, bank_crit;
    int iter = 0;

    // set the upper and lower angle limits
    double b_upper = phi - deltaX;
    double b_lower = deltaX;
    CH.theta = 0.67 * phi;

    // calculate the bank stability index
    bank_crit = G * RHO * Gs * D90 * tau_star *
              pow( 1 - ( pow( sin ( CH.theta * PI / 180 ), 2) /
              pow( sin( phi * PI / 180 ), 2) ), 0.5 );

    energyConserve( n, r );                 // Update section data based on new theta
    channelState( n, ch_idx, r );                   // Update stresses, transport capacity
    converg = ( CH.Tbank - bank_crit ) / bank_crit;

    while( ( abs( converg ) > Tol ) && ( iter < 50 ) )
    {
        if( converg > 0 )
            b_upper = CH.theta;
        else
            b_lower = CH.theta;
        CH.theta = 0.5 * (b_upper + b_lower);

        energyConserve( n, r );            // Update section data based on new theta
        channelState( n, ch_idx, r );              // Update stresses, transport capacity
        converg = ( CH.Tbank - bank_crit ) / bank_crit;

        bank_crit = G * RHO * Gs * D90 * tau_star *
                  pow( 1 - ( pow( sin ( CH.theta * PI / 180 ), 2) /
                  pow( sin( phi * PI / 180 ), 2) ), 0.5 );

        converg = (CH.Tbank - bank_crit) / bank_crit;

        iter++;
    }
}

void hydro::setRegimeWidth(RiverProfile *r)
{

    // Adjust channel regime one cross-section at a time, marching upstream

    float deltaArea, oldArea = 0;            // Change in reach cross-section area
    float deltaEta = 0;
    float reachDrop = 0;                     // Drop in elevation over reach river length
    float oldBankHeight = 0;

    // oldBankHeight = r->RiverXS[regimeCounter].bankHeight;
    // oldArea = r->RiverXS[regimeCounter].xsFlowArea[2];

    regimeModel( regimeCounter, r );         // This is a call to the Regime Functions.

    if ( (r->counter > 260 ) )
    {
        deltaArea = oldArea - r->RiverXS[regimeCounter].xsFlowArea[2];       // Change induced by floodplain erosion
        deltaEta = r->RiverXS[regimeCounter].xsBankHt - oldBankHeight;    // Change induced by channel aggr/degr
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
