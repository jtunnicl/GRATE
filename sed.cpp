#include "sed.h"
#include <math.h>
#include<iostream>
#include<fstream>
using namespace std;

sed::sed(RiverProfile *r)
{
    initSedSeries(r->nnodes);
}

void sed::initSedSeries(int nodes)
{
    std::ifstream inSedFile;
    double currentCoord = 0.;
    int v2, v3, v4, v5, v6, v7, v9;        // Yr, Mt, Day, Hr, Mn, Sec, GRP
    double v1, v8;                         // Coord, Qs
    QDateTime NewDate;
    vector< TS_Object > tmp;
    TS_Object NewEntry;

    NewDate.setDate(QDate(2000,1,1));
    NewDate.setTime(QTime(0,0,0));

    inSedFile.open("sed_series.dat", ios::in);

    if (!inSedFile) cout << "file couldn't open sed file properly" << endl;
    
    while(! inSedFile.eof() )
    {
        if(inSedFile >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8 >> v9)
        {
            NewDate.setDate(QDate(v2, v3, v4));
            NewDate.setTime(QTime(v5, v6, v7));
            NewEntry.date_time = NewDate;
            NewEntry.Q = v8;
            NewEntry.Coord = v1;
            NewEntry.GRP = v9 - 1;

            if (v1 > currentCoord){            // Have we moved to a new source coordinate?
                Qs_series.push_back( tmp );
                tmp.clear();
                currentCoord = v1;
                tmp.push_back(NewEntry);  // Start new tmp
                }
            else
                tmp.push_back(NewEntry);
        }
    }

    Qs_series.push_back( tmp );                 // Final tmp loaded into Qs_series array
    Qs.resize(nodes);                                       // Bedload transport (m3/s) at each node
    deta.resize(nodes);                                     // Rate of vertical bed change (d-eta) with time (dt)
    dLa_over_dt.resize(nodes);
    p.resize(nodes);
    df.resize(nodes);
}

void sed::setNodalSedInputs(RiverProfile *r)
{
    int j = 0;
    unsigned int i = 0;

    if ( Qs_series[0][0].date_time.secsTo( r->cTime ) < 1 )        // Start of run?
        for ( i = 0; i < Qs_series.size(); i++ )                 // Qs.size is the # of tribs/sources
            Qs_bc.push_back( Qs_series[i][0] );
    else
    {
        j = 0;
        while( Qs_series[0][j].date_time.secsTo( r->cTime ) > 0 )
            j++;
        for ( i = 0; i < Qs_series.size(); i++ )
        {
            Qs_bc[i].Coord = Qs_series[i][j-1].Coord;
            Qs_bc[i].GRP = Qs_series[i][j-1].GRP;
            Qs_bc[i].date_time =  Qs_series[i][j-1].date_time;

            Qs_bc[i].Q = ( Qs_series[i][j-1].Q + ( Qs_series[i][j-1].date_time.secsTo(r->cTime) ) *
                   ( Qs_series[i][j].Q - Qs_series[i][j-1].Q ) /
                   ( Qs_series[i][j-1].date_time.secsTo(Qs_series[i][j].date_time) ));
            Qs_bc[i].Q *= r->qsTweak;
            Qs_bc[i].Q *= r->tweakArray[r->yearCounter];                           // Flood = 0.8 to 1.8 mean flow
        }
    Qs_bc[0].Q *= r->qwTweak;                                         // Feed randomizer
    }
}

NodeGSDObject sed::multiplyGSD(NodeGSDObject &M, NodeGSDObject &N, double weight, RiverProfile *r)
{

    // this routine is used to multiply two grain size distributions together,
    // with proportion 'weight' used as the weighting on 1st element, '1-weight' as the other.
    // Primary use is for chi constant (0.7). '0.5' otherwise.

    NodeGSDObject fi;                          // return GSD object

    for ( int j = 0; j < r->ngsz; j++ )
    {
        for ( int k = 0; k < r->nlith; k++ )
        {
            fi.pct[k][j] = weight * M.pct[k][j] + ( 1.0 - weight ) * N.pct[k][j];
        }
    }

    fi.norm_frac();

    return(fi);
}

void sed::computeTransport(RiverProfile *r)
{
    unsigned int bc;
    int i, j, k;
    int inode;
    NodeGSDObject qtemp;             // temporary, for storing grain size fractions
    float ngsz, nlith;
    float taussrg;                 // Wilcock - reference (median) shear
    float b;                       // b exponent for each size fraction
    float arg;                     // decision for G
    float phisgo;
    float dj;                      // grain size;
    float ds50;
    float specWt;                  // submerged specific weight of gravel
    float a0;
    float Wwc;                     // Wi* from Wilcock Crowe
    float FGSum;

    vector<float> ktot, ktotn;

    ngsz = r->F[0].psi.size() - 2;
    nlith = r->F[0].abrasion.size();
    ktot.resize(ngsz);
    ktotn.resize(ngsz);

    specWt = 0.65; //(2650 - 1000) / 1000 - 1.;
    setNodalSedInputs(r);                                      // Calculate inputs at each tributary

    for ( i = 0; i < r->nnodes; i++ )                          // iterate nodes
    {
        for ( j = 0; j < r->ngsz; j++ )                        // iterate grain size
            for ( k = 0; k < r->nlith; k++ )                   // iterate lithology
                fpp.pct[k][j] = r->F[i].pct[k][j];             // temp bl is extracted from the surface layer

        fpp.norm_frac();                                       // Normalize f fractions
        fpp.dg_and_std();

        if (r->eta[i] >= r->bedrock[i])
        {
            taussrg = 0.021 + 0.015 * exp( -20 * fpp.sand_pct );
            phisgo = ( ( r->RiverXS[i].ustar * r->RiverXS[i].ustar ) / specWt / 9.81 / (pow( 2, fpp.dsg ) / 1000)) / taussrg;
            FGSum = 1e-10;
            Wwc = 0.;

            for ( j = 0; j < ngsz; j++ )
            {
                ktot[j] = 0;
                a0 = ( 0.5 * ( fpp.psi[j] + fpp.psi[j+1] ) );
                ds50 = pow( 2.0, fpp.dsg ) / 1000;
                dj = pow( 2.0, a0 ) / 1000;
                b = 0.67 / (1 + exp( 1.5 - ( dj / ds50 ) ) );    // Wilcock eqn. (4)
                arg = phisgo * pow( ( dj / ds50 ), -b );

                if (arg < 1.35)
                    Wwc = 0.002 * pow( arg, 7.5 );                 // eqn.7a
                else
                    Wwc = 14 * pow( ( 1 - 0.894 / sqrt(arg) ), 4.5 );  //eqn. 7b

                for ( k = 0; k < nlith; k++ )
                {
                    fpp.pct[k][j] *= Wwc;
                    ktot[j] += fpp.pct[k][j];
                }

                FGSum += ktot[j];
            }

            // Normalize the bedload fractions
            fpp.norm_frac();

            if (FGSum > 0)
                Qs[i] = FGSum * pow( r->RiverXS[i].ustar, 3 ) / specWt / 9.81 * ( r->RiverXS[i].xsBedWidth);
            else
                Qs[i] = 0.0;
        }
        else
            Qs[i] = 0.0;

        for ( j = 0; j < r->ngsz; j++ )
            for ( k = 0; k < r->nlith; k++ )
                p[i].pct[k][j] = fpp.pct[k][j];                // revise the bedload grain-size fractions

        p[i].abrasion[0] = r->randAbr;
        p[i].abrasion[1] = r->randAbr;
        p[i].abrasion[2] = r->randAbr;
    }
    
    //Qs[r->nnodes - 1] = Qs[r->nnodes - 2];               // Equilibrium bottom node

    // Adjust bedload for sed inflow - revise bedload for main channel and tributary sediment inflows

    for ( bc = 0; bc < Qs_series.size(); bc++ )
    {
        inode = ceil( Qs_series[bc][0].Coord / r->dx);                         // Node where trib is entering

        for ( j = 0; j < r->ngsz; j++ )
        {
            for ( k = 0; k < r->nlith; k++ )
            {
                qtemp.pct[k][j] = p[inode].pct[k][j] * Qs[inode] + r->grp[Qs_bc[bc].GRP].pct[k][j] * Qs_bc[bc].Q;
            }
        }

        qtemp.norm_frac();
	
        if (bc > 0)
        {
            Qs[inode] += Qs_bc[bc].Q * 0.75;
            Qs[inode+1] += Qs_bc[bc].Q * 0.25;                    // Distribute trib material downstream
        }
        else
            Qs[inode] += Qs_bc[bc].Q;

        for ( j = 0; j < r->ngsz; j++ )
            for ( k = 0; k < r->nlith; k++ )
                p[inode].pct[k][j] = qtemp.pct[k][j];

    }     //  trib_input                   //  Add additional, storage load from major tributaries

    exner(r);
}

void sed::exner(RiverProfile *r)
{
    int i, j, k, m = 0;
    double upw = 0.75;                                             // Upwinding constant
    double chi = 0.7;                                              // weighting for interfacial exchange
    double dmy;
    vector<double> tmp(3);
    tmp[0] = 0;
    tmp[1] = 0;
    tmp[2] = 0;

    NodeGSDObject fi, Fprime;                                      // Temporary grain-size container

    for ( i = 1; i < (r->nnodes-1); i++ )                          // Calculate deta
    {
        deta[i] = r->dt * ( upw * ( ( Qs[i-1] - Qs[i] ) / ( r->xx[i] - r->xx[i-1] ) )
            + ( 1 - upw ) * ( ( Qs[i] - Qs[i+1] ) / ( r->xx[i+1] - r->xx[i] ) ) )
                      / (1.0 - r->poro) / r->RiverXS[i].fpWidth;
    }

    deta[0] = ( ( Qs_bc[0].Q - Qs[1] ) / ( r->xx[1] - r->xx[0] ) );

    for ( i = 1; i < r->nnodes-1; i++ )                            // Upstream boundary - if floating, i = 0
        r->eta[i] += deta[i];

    //r->eta[r->nnodes-1] += deta[r->nnodes-2];                    // Downstream boundary - uncomment if floating

    for ( i = 1; i < r->nnodes-1; i++ )      // Calculate grain size changes
    {
        if ( r->toplayer[i] <= 0.0 )
        {
            r->toplayer[i] += r->layer;
            r->ntop[i] = r->ntop[i] - 1;
        }

        if ( deta[i] >= 0.0 )                  // interface, aggradational case
            fi = multiplyGSD(p[i], r->F[i], chi, r);

        else                                   // interface, degradational case
        {
            for ( j = 0; j < r->ngsz; j++ )
                for ( k = 0; k < r->nlith; k++ )
                    fi.pct[k][j] = r->storedf[i][r->ntop[i]].pct[k][j];    // applied to all degrading nodes

            if ( -deta[i] > r->toplayer[i] )      // degrade more than one layer
            {
                for ( j = 0; j < r->ngsz; j++ )
                    for ( k = 0; k < r->nlith; k++ )
                        fi.pct[k][j] *= r->toplayer[i];                // applied to all degrading nodes
                fi.norm_frac();

                dmy = -deta[i] - r->toplayer[i] - r->layer;
                m = r->ntop[i] - 1;

                while (dmy > 0.0)
                {
                    if (m <= 0)
                        cout << "Erosion has reached the bottom of the lowest storage layer at node " <<  i;

                    for (  j = 0; j < r->ngsz; j++ )
                        for (  k = 0; k < r->nlith; k++ )
                            fi.pct[k][j] += r->layer * r->storedf[i][m].pct[k][j];   // applied to all degrading nodes
                    fi.norm_frac();

                    dmy = dmy - r->layer;
                    m = m - 1;
                }                              // end while loop

                if (m <= 0)
                {
                    cout << "Erosion has reached the bottom of the lowest storage layer at node " <<  i;
                    break;
                }

                for ( j = 0; j < r->ngsz; j++ )
                    for ( k = 0; k < r->nlith; k++ )
                        fi.pct[k][j] += ( r->layer + dmy ) * r->storedf[i][m].pct[k][j];
                fi.norm_frac();
            }                                  // end degrading more than 1 layer
        }                                      // end aggradational/degradational cases


    // Estimate F'

        for ( j = 0; j < r->ngsz; j++ )
            for ( k = 0; k < r->nlith; k++ )
                Fprime.pct[k][j] = r->F[i].pct[k][j] / sqrt( pow( 2, ( r->F[i].psi[j+1] - r->F[i].psi[j] ) / 2 ) );

        Fprime.norm_frac();
        Fprime.pct.push_back(tmp);            // add a 'j+1' category (zeros), to satisfy 'df' equation below.

        for ( j = 0; j < r->ngsz; j++ )
        {
            df[i].pct[0][j] = 0.0;
            df[i].pct[1][j] = 0.0;
            df[i].pct[2][j] = 0.0;
            for ( k = 0; k < r->nlith; k++ )
                df[i].pct[k][j] += -(r->dt / r->RiverXS[i].xsBedWidth ) *
                              ( ( upw * ( Qs[i] * p[i].pct[k][j] - Qs[i-1] * p[i-1].pct[k][j] ) / ( r->xx[i] - r->xx[i-1] )
                              + ( 1 - upw ) * ( Qs[i+1] * p[i+1].pct[k][j] - Qs[i] * p[i].pct[k][j] ) / (r->xx[i+1] - r->xx[i] ) )
                              - p[i].abrasion[k] * Qs[i] * ( p[i].pct[k][j] + Fprime.pct[k][j] )
                              + p[i].abrasion[k] * Qs[i] * ( 1 / ( 3 * log(2) ) ) * ( (p[i].pct[k][j] + Fprime.pct[k][j] )
                              / (r->F[i].psi[j+1] - r->F[i].psi[j]) - ( p[i].pct[k][j+1] + Fprime.pct[k][j+1] ) / ( r->F[i].psi[j+2] - r->F[i].psi[j+1] ) ) )
                              / ( 1.0 - r->poro ) - fi.pct[k][j] * deta[i] + ( fi.pct[k][j] - r->F[i].pct[k][j] ) * dLa_over_dt[i] * r->dt;
        }
    }  // end for loop

    // New loop - update bed grain size distribution
    for ( i = 2; i < r->nnodes-1; i++ )
    {
        for ( j = 0; j < r->ngsz; j++ )
            for ( k = 0; k < r->nlith; k++ )
                r->F[i].pct[k][j] += df[i].pct[k][j];
        r->F[i].norm_frac();
    }

    // Update storage layers
    for ( i = 2; i < r->nnodes-1; i++ )
    {
        if (deta[i] < 0.0)
        {
            dmy = -deta[i] - r->toplayer[i];
            while (dmy >= 0.0)
            {
                dmy -= r->layer;
                r->ntop[i]--;
                if (r->ntop[i] <= 0.0)            // Raise exception here; bedrock reached.
                {
                    cout << "Bedrock reached at node " <<  i;
                    break;
                }
            }
            r->toplayer[i] = -dmy;
        }                                      // end degradational case
        else                                   // begin aggradational case
        {
            if ((deta[i] + r->toplayer[i]) <= r->layer)
            {
                for ( j = 0; j < r->ngsz; j++ )
                {
                    for ( k = 0; k < r->nlith; k++ )
                    {
                        r->storedf[i][r->ntop[i]].pct[k][j] = deta[i] * (chi * p[i].pct[k][j] + ( 1.0 - chi ) *
                                                      r->F[i].pct[k][j]) + r->toplayer[i] * r->storedf[i][r->ntop[i]].pct[k][j];
                    }                          // aggraded material is a mixture of p and f.
                }

                r->storedf[i][r->ntop[i]].norm_frac();
                r->toplayer[i] += deta[i];
            }
            else
            {                                  //aggrade more than current layer
                for ( j = 0; j < r->ngsz; j++ )
                    for ( k = 0; k < r->nlith; k++ )
                        r->storedf[i][r->ntop[i]].pct[k][j] = ( r->layer - r->toplayer[i] ) * ( chi * p[i].pct[k][j] +
                                ( 1.0 - chi ) * r->F[i].pct[k][j] ) + r->toplayer[i] * r->storedf[i][r->ntop[i]].pct[k][j];
                                               // fill in additional stratigraphy w/ mixture of p and f.

                r->storedf[i][r->ntop[i]].norm_frac();

                dmy = deta[i] + r->toplayer[i] - r->layer;
                while (dmy > 0.0)
                    {
                        r->ntop[i]++;
                        if (r->ntop[i] > (r->nlayer - 2))      //raise Exception: 'not enough storage layers for aggradation.'
                        {
                            cout << "Note enough storage layers for aggradation at node " <<  i;
                            break;
                        }
                        for ( j = 0; j < r->ngsz; j++ )
                            for ( k = 0; k < r->nlith; k++ )
                                r->storedf[i][r->ntop[i]].pct[k][j] = chi * p[i].pct[k][j] + (1.0 - chi) * r->F[i].pct[k][j];
                        r->storedf[i][r->ntop[i]].norm_frac();
                        dmy -= r->layer;
                    }
                    r->toplayer[i] = dmy + r->layer;

            }        //aggrade 1 or more layers;
        }        //aggradational case;
    }        //update storage layers
}
