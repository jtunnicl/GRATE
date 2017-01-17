/*******************
 *
 *
 *  GRATE 9
 *
 *  Long profile parameters
 *
 *
 *
*********************/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "riverprofile.h"
using namespace std;

#define PI 3.14159265

double gammln2(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,    -86.50532032941677,
            24.01409824083091,    -1.231739572450155,
            0.1208650973866179e-2,-0.5395239384953e-5};
  int j;

  y = x = xx;
  tmp = x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser = 1.000000000190015;
  for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
  return -tmp + log( 2.5066282746310005 * ser / x);
}

NodeGSDObject::NodeGSDObject()
    {
         vector <double> tmp;                  // dummy array for grain sizes
         abrasion.push_back(0.0000060000);
         abrasion.push_back(0.0000060000);
         abrasion.push_back(0.0000060000);

         for (int j = 0; j < 15; j++)
         {
             tmp.push_back(0);
             psi.push_back(-3 + j );           // psi -3 to 11 .. should be 9, but ngsz+1,2 is required throughout
         }

         for (int k = 0; k < 3; k++)
             pct.push_back(tmp);

         dsg = 0.;
         stdv = 0.;
         sand_pct = 0.;
    }

void NodeGSDObject::norm_frac()
{
    float ngsz, nlith, cumtot;
    vector<float> ktot;
    ktot.resize(psi.size());

    ngsz = psi.size() - 2;
    nlith = abrasion.size();

    sand_pct = 0;

    // Normalize
    cumtot = 0.0;
    for ( int j = 0; j < ngsz; j++ )           // Sum mass fractions
    {
        ktot[j] = 0.0;
        for ( int k = 0; k < nlith; k++ )
            if (pct[k][j] > 0)                 // Solves problems with rounding
                ktot[j] += pct[k][j];
            else
                pct[k][j] = 0;
        cumtot += ktot[j];
    }

    for ( int j = 0; j < ngsz; j++ )
        for ( int k = 0; k < nlith; k++ )
        {
            if (pct[k][j] > 0)
                pct[k][j] /= cumtot;
            if (psi[j] <= 0) sand_pct += pct[k][j];   // sum sand fraction
        }
}

void NodeGSDObject::dg_and_std()
{
    float tdev;
    float ngsz, nlith;
    vector<float> ktot;
    ktot.resize(psi.size());

    ngsz = psi.size() - 2;
    nlith = abrasion.size();

    dsg = 0.0;
    d84 = 0.0;
    d90 = 0.0;

    for ( int j = 0; j < ngsz; j++ )
    {
        //if (psi[j] >= -1)   // ***!!!*** Dg50 is based on gsizes >= 2 mm (psi = -1 or less)
        //{
        ktot[j] = 0;
        for ( int k = 0; k < nlith; k++ )
            ktot[j] += pct[k][j];               // lithology values for each size fraction are summed.
        dsg += 0.50 * (psi[j] + psi[j+1]) * ktot[j];
        d84 += 0.84 * (psi[j] + psi[j+1]) * ktot[j];
        d90 += 0.90 * (psi[j] + psi[j+1]) * ktot[j];
        //}
    }

    stdv = 0.0;
    for ( int j = 0; j < ngsz; j++ )
    {
        ktot[j] = 0;
        for ( int k = 0; k < nlith; k++ )
            ktot[j] = ktot[j] + pct[k][j];

        tdev = 0.5 * (psi[j] + psi[j+1]) - dsg;
        stdv += 0.5 * tdev * tdev * ktot[j];
    }

    if (stdv > 0)
        stdv = sqrt(stdv);
}

NodeCHObject::NodeCHObject()
{
    QProp = 0.;                             // Proportion of total flow directed to this channel
    depth = 1.;                                // Given the proportion of flow in the channel, this is the computed depth - modified later in xsGeom()
    width = 0.;
    bankHeight = 3.;                           // Measured relative to channel bottom
    b2b = 0.0;
    flowArea = 0.0;
    flowPerim = 0.0;
    ovBank = 0;
    Tbed = 20.;
    Tbank = 20.;
    Qb_cap = 0.2;
    comp_D = 0.005;
    K = 0;

    Hmax = 0.5;
    mu = 1.5;                                 // Not used - perhaps in future versions
    theta = 30.;
}

void NodeCHObject::chGeom(double relDepth)
{   // Update channel cross-section area and perimeter
    // This is only for the channel itself. Overbank flows are handled in the NodeXSObject
    float theta_rad = theta * PI / 180;        // theta is always in degrees

    b2b = width + 2 * ( bankHeight - Hmax) / tan( theta_rad );
    if ( depth <= ( bankHeight - Hmax ) )      // w.s.l. is below sloping bottom edges
    {
        flowArea = width * depth + pow ( depth, 2 ) / tan( theta_rad );
        flowPerim = width + 2 * depth / sin ( theta_rad );
    }
    else
    {
        flowArea = b2b * depth - pow ( ( bankHeight - Hmax ), 2 ) / tan( theta_rad );
        flowPerim = width + 2 * ( bankHeight - Hmax ) / sin ( theta_rad ) + 2 * ( depth - (bankHeight - Hmax) );
    }

    if (depth > bankHeight)
        ovBank = 1;
    else
        ovBank = 0;
}

NodeXSObject::NodeXSObject()                   // Initialize object
{
     NodeCHObject tmp;
     int i;

     node = 0;
     numChannels = 1;
     fpWidth = 0.;
     maxBankHt = 1.;
     chSinu = 1.05;
     topW = 10.;
     xsBedWidth = 2;
     meanVeloc = 0.;
     ustar = 0.;
     maxDepth = 1.;
     ovBankFlag = 0;
     mainChannel = 0;
     hydRadius = 0;
     critDepth = 0.;
     centr = 0.;
     rough = 0.;
     omega = 0.;
     k_mean = 0.;
     eci = 0.;

     for (i = 0; i < 3; i++)
     {
         xsFlowArea.push_back(0);
         xsFlowPerim.push_back(0);
     }

     for (int i = 0; i < 10; i++)              // Max 10 dummy channel objects initiated, all with QProp '0'..
         CHList.push_back(tmp);
     CHList[0].QProp = 1;                      // ..except 1st element. Assume single-channel default.
}

void NodeXSObject::xsGeom()                    // Given maxDepth, update flow XS area at a given node, including overbank flows
{
    int i;
    double deltaWSL = maxDepth - maxBankHt;    // When deltaWSL = 0, flow depth is at bank top. +ve is overbank, -ve is in channel
    ovBankFlag = 0;
    xsBedWidth = 0.;
    topW = 0.;

    for (i = 0; i < 3; i++)                    // Clear out old data
    {
        xsFlowArea[i] = 0.;
        xsFlowPerim[i] = 0.;
    }

    for (i = 0; i < numChannels; i++)          // Add up in-channel flow, area, perimeter and bed width
    {
        CHList[i].depth = maxDepth;            // Send XS depth to channel object, returns area possibly including overbank flow
        CHList[i].chGeom(deltaWSL);
        xsFlowArea[0] += CHList[i].flowArea;
        xsFlowPerim[0] += CHList[i].flowPerim;
        xsBedWidth += CHList[i].width;
        topW += CHList[i].b2b;
        xsFlowArea[1] += (CHList[i].depth);
        if (CHList[i].ovBank == 1)             // If flows go overbank..
            ovBankFlag = 1;

        if (CHList[i].depth > maxDepth)
        {
            maxDepth = CHList[i].depth;         // identify deepest, 'main' channel
            mainChannel = i;
        }
    }

    if (ovBankFlag == 1)
    {
        maxDepth = maxDepth + (xsFlowArea[1] / fpWidth);  // maxDepth includes overbank flows
        xsFlowPerim[1] += fpWidth + 2 * (xsFlowArea[1] / fpWidth);
    }

    xsFlowArea[2] = xsFlowArea[0] + xsFlowArea[1];       // Sum total area and perim for the cross-section
    xsFlowPerim[2] = xsFlowPerim[0] + xsFlowPerim[1];
    hydRadius = xsFlowArea[2] / xsFlowPerim[2];
}

void NodeXSObject::xsCentr()
{
    // Compute centroid of flow
                                                                        // Elevation at point where floodplain meets valley wall
    /* The formula for a trapezoid with base (a), top width (b),
       and height (h) - arbitrary side slope length - is:

                      / 2a + b  \
           (h/3)  *  | --------  |
                      \  a + b  /                             */


    centr = 1;         // Need to figure out computation scheme, here.

    // NEED TO DIVIDE BY TOTAL AREA (e.g. n.xs_area()) AFTER GETTING THIS ARRAY (?)

}

void NodeXSObject::xsECI(NodeGSDObject F)
{

    F.norm_frac();
    F.dg_and_std();                                    // Update grain size statistics

    rough = 2 * pow( 2, F.dsg ) / 1000. * pow( F.stdv, 1.28 );              // roughness height, ks, 2*D90
    if (rough <= 0)         // indicates problems with previous F calcs
        rough = 0.01;
    omega = 2.5 * log( 11.0 * ( CHList[mainChannel].depth / rough ) );                          // Parker (1991), Dingman 6.25, p.224
    double K_ch = xsFlowArea[0] * omega * sqrt( 9.81 * CHList[mainChannel].depth );              // Dingman, (2009) 8B2.3C, p.300
    double K_fp = 0;
    k_mean = 0;
    double ovBank = CHList[mainChannel].depth - CHList[0].bankHeight;

    if (ovBank > 0)
    {
        K_fp = xsFlowArea[1] * rough * sqrt( 9.81 * ovBank * 0.5 );   // Depth halfway across the floodplain
        k_mean = K_ch + K_fp;
        eci = ( pow(K_ch,3) / pow(xsFlowArea[0],2) + pow(K_fp,3) / pow(xsFlowArea[1],2) ) /
                ( pow(k_mean,3) / pow(xsFlowArea[2],2) );                   // Dingman, (2009) 8B2.4, Chaudhry (2nd ed) 4-41
    }
    else
    {
        eci = 1;
        k_mean = K_ch;
    }
}

TS_Object::TS_Object()
{

    date_time.setDate(QDate(2000,1,1));
    date_time.setTime(QTime(12, 0, 0, 0));
    Q = 0.;
    Coord = 0.;
    GRP = 1;

}

RiverProfile::RiverProfile()
    {

    NodeXSObject tmp;                          // Initialize RiverXS Object
    nnodes = 0;
    npts = 0;
    dt = 0;
    dx = 0.0;
    ngsz = 0;
    nlith = 0;
    ngrp = 0;
    ngsz = 0;
    nlith = 0;
    ngrp = 0;
    nlayer = 0;
    poro = 0.0;
    default_la = 0.0;
    layer = 5.0;

    RiverXS.push_back(tmp);

    cTime.setDate(QDate(2002,12,5));           // Arbitrary date initialisation
    cTime.setTime(QTime(12,0,0,0));
    startTime.setDate(QDate(2002,12,5));
    startTime.setTime(QTime(12,0,0,0));
    endTime.setDate(QDate(2002,12,5));
    endTime.setTime(QTime(12,0,0,0));

    counter = 0;
    yearCounter = 0;                           // Model 'year' counter that resets every 5 days
                                               //          e.g. 5d * 24hr * 60 * 60 / ( dt )

    for (int i = 0; i < 5; i++)
        N.push_back(0);                        // Substrate shift matrix
    for (int i = 0; i < 10; i++)
	    rand_nums.push_back(0);
                                               // Random tweak variables are based on logarithmic (e) scaled values
                                               // Augment the rate of tributary Qs, Qw inputs
    tweakArray = hydroGraph();                 // A gamma-distribution that simulates hydrograph form
    qsTweak = 1;                               // rand_nums[1] * 1.5 + 0.5;        // qs between 0.5 and 2
    qwTweak = 1;                               // Hydrograph multiplier
    substrDial = 0;                            // rand_nums[3]  * 3.8 - 1.9;      // Positive (up to +2) makes finer mix, negative (down to -2) coarsens all grain groups
    feedQw =  1;                               // rand_nums[4]  * 0.5 + 0.75;     // between 0.75 and 1.25
    feedQs =  1;                               // rand_nums[5] + 0.5;                     // between 0.5 and 1.5
    HmaxTweak = 1;                             // See line ~750ff
    randAbr = 1e-5;                            // between 10^-4 and 10^-7

    // Set up substrate shift matrix

    if (substrDial > 0 && substrDial < 1)
    {
         N[2] = 1 - substrDial;
         N[3] = substrDial;
    }
    if (substrDial >= 1 && substrDial < 2)
    {
         N[3] = 1 - ( substrDial - 1 );
         N[4] = ( substrDial - 1 );
    }
    if (substrDial >= 2) N[4] = 1;
    if (substrDial == 0) N[2] = 1;
    if (substrDial < 0 && substrDial > -1)
    {
        substrDial = abs(substrDial);
        N[1] = abs(substrDial);
        N[2] = 1 - abs(substrDial);
    }
    if (substrDial <= -1 && substrDial > -2)
    {
        N[0] = abs( substrDial + 1 );
        N[1] = 1 - abs( substrDial + 1 );
    }
    if (substrDial <= -2) N[0] = 1;

    if ( ( N[0] + N[1] + N[2] + N[3] + N[4] ) > 1)
       cout << "Interpolation Array is over 1.0";

    initData();

    }

vector<float> RiverProfile::hydroGraph()
{
    /* This routine creates a hydrograph, based on the gamma distribution,
       meant to simulate the range of flow experienced over the course of
       1 year.
    */
  int i = 0;

  double max_flow = 1.6;  // Up to 1.6 * 50 = 80 m3/s
  double min_flow = 0.8;
  double elems = 900;
  double alpha = 4;
  double beta = 5;
  double delta = 0.0052;
  double factor = 0;
  vector<float> x, xx, fac;

  for (int i=0; i<elems; ++i)
  {
    x.push_back(0);
    xx.push_back(0);
    fac.push_back(0);
  }

  x[i] = -1.5;
  for (int i=1; i<elems; ++i)
  {
    x[i] = x[i-1] + delta;
    xx[i] = pow(10,x[i]);
    factor = alpha*log(beta)-gammln2(alpha);
    fac[i] = exp(-beta * xx[i] + ( alpha - 1. ) * log( xx[i] ) + factor);
    if (fac[i] < 0.001)
      fac[i] = 1e-3;
    fac[i] = (fac[i]) * (max_flow - min_flow) + min_flow;
    // std::cout << x[i] << ":\t" << fac[i] << endl;
  }
  fac[0] = fac[1];

  return fac;
}

void RiverProfile::initData()
    {

    std::ifstream inDatFile;

    NodeGSDObject tmp;
    vector < NodeGSDObject > tmp2;
    const char * f;
    int i = 0;
    char g[8];

    ofstream myfile ("example.txt");
    if (myfile.is_open())
    {
      myfile << "This is a line.\n";
      myfile << "This is another line.\n";
      myfile.close();
    }

    inDatFile.open("GRATEInputFile1.dat");

    if( !inDatFile )                           //test the file open.
    {
        cout<<"Error opening profile intput file.."<<endl;
        system("pause");
    }

    f = getNextParam(inDatFile, "NNODES");
    for (i = 0; i < 8; i++) g[i] = *(f++);
    nnodes = atoi(g);

    // Allocate vectors
    xx.resize(nnodes);
    eta.resize(nnodes);
    algrp.resize(nnodes);
    stgrp.resize(nnodes);
    bedrock.resize(nnodes);
    RiverXS.resize(nnodes);

    f = getNextParam(inDatFile, "LAYER");
    for (i = 0; i < 8; i++) g[i] = *(f++);
    layer = atof(g);
    toplayer.assign(nnodes, layer);            // Thickness of the top storage layer; starts at 5 and erodes down

    f = getNextParam(inDatFile, "LA");
    for (i = 0; i < 8; i++) g[i] = *(f++);
    default_la = atof(g);
    la.assign(nnodes, default_la);             // Default active layer thickness

    f = getNextParam(inDatFile, "NLAYER");
    for (i = 0; i < 8; i++) g[i] = *(f++);
    nlayer = atoi(g);
    ntop.assign(nnodes, nlayer-15);            // Indicates # of layers remaining, below current (couple of layers left for aggradation)

    f = getNextParam(inDatFile, "PORO");
    for (i = 0; i < 8; i++) g[i] = *(f++);
    poro = atof(g);                            // Default deposit porosity

    for (i = 0; i < nlayer; i++)               // Init storedf stratigraphy matrix
        tmp2.push_back(tmp);

    for (i = 0; i < nnodes; i++)
    {
        storedf.push_back(tmp2);
        F.push_back(tmp);
    }

    f = getNextParam(inDatFile, "NGSZ");
    for (i = 0; i < 8; i++) g[i] = *(f++);
    ngsz = atoi(g);

    f = getNextParam(inDatFile, "NLITH");
    for (i = 0; i < 8; i++) g[i] = *(f++);
    nlith = atoi(g);

    f = getNextParam(inDatFile, "NGRP");
    for (i = 0; i < 8; i++) g[i] = *(f++);
    ngrp = atoi(g);

    for (i = 0; i < ngrp; i++)
        grp.push_back(tmp);

    getGrainSizeLibrary(inDatFile);

    getLibraryLith(inDatFile);

    for (i = 0; i < nnodes; i++)
    {
        RiverXS[i].maxDepth = 1.5;
        RiverXS[i].wsl = eta[i] + RiverXS[i].maxDepth;
        RiverXS[i].node = i;
    }

    f = getNextParam(inDatFile, "NPTS");
    for (i = 0; i < 8; i++) g[i] = *(f++);
    npts = atoi(g);

    getLongProfile(inDatFile);

    getStratigraphy(inDatFile);

    dx = xx[1]-xx[0];                       // Assume uniform grid
    dt = 10;

}

const char *RiverProfile::getNextParam(std::ifstream &openFile, const char *nextParam)
{

    const char *token[40] = {};              // initialize to 0; 40 tokens max
    char buf[512];                           // Max 512 chars per line
    int Found = 0;
    int n = 0;

    while (Found == 0)
    {
        openFile.getline(buf, 512);
        token[0] = strtok( buf, " " );       // first token
        if (token[0] == NULL || strcmp(token[0], "!") == 0)       // zero if line is blank or "!"
            continue;
        else
            for (n = 1; n < 10; n++)
            {
                token[n] = strtok(0, " ");   // subsequent tokens
                if (!token[n]) break;        // no more tokens
            }
        if (strcmp(token[0], nextParam) == 0 )
            Found = 1;
    }

    return (token[2]);
}

void RiverProfile::getGrainSizeLibrary(ifstream &openFile)
{
    const char* token[40] = {};                // initialize to 0; 40 tokens max
    char buf[512];                             // Max 512 chars per line
    int Found = 0;
    int n, grpCount = 0;

    while (Found == 0)                         // Skip through comments
    {
        openFile.getline(buf, 512);
        token[0] = strtok( buf, " " );         // first token
        if (token[0] == NULL || strcmp(token[0], "!") == 0)       // zero if line is blank or "!"
            continue;
        else
        Found = 1;
    }

    for (int gsCount = 0; gsCount < ngsz; gsCount++)
    {
        for (n = 0; n < 50; n++)
        {
            token[n] = strtok(0, " ");
            if (!token[n]) break;              // no more tokens
        }

        for (grpCount = 0; grpCount < (ngrp); grpCount++)
            {
                grp[grpCount].pct[0][gsCount] = atof(token[grpCount]);
                grp[grpCount].pct[1][gsCount] = atof(token[grpCount]);
                grp[grpCount].pct[2][gsCount] = atof(token[grpCount]);
            }

        openFile.getline(buf, 512);            // proceed to next line
        token[0] = strtok( buf, " " );
    }
}

void RiverProfile::getLibraryLith(ifstream &openFile)
{
    NodeGSDObject qtemp;
    const char* token[40] = {};              // initialize to 0; 40 tokens max
    char buf[512];                           // Max 512 chars per line
    int Found = 0;
    int n, gsCount, grpCount, lithCount = 0;

    gsCount = 0;
    while (Found == 0)                       // Skip through comments
    {
        openFile.getline(buf, 512);
        token[0] = strtok( buf, " " );                    // first token
        if (token[0] == NULL || strcmp(token[0], "!") == 0)       // zero if line is blank or "!"
            continue;
        else
        {
            for (lithCount = 0; lithCount < nlith; lithCount++)
            {
                for (n = 1; n < 50; n++)           // tokenize string
                {
                    token[n] = strtok(0, " "); // subsequent tokens
                    if (!token[n]) break;            // no more tokens
                }

                for (grpCount = 0; grpCount < ngrp; grpCount++)
                    grp[grpCount].pct[lithCount][gsCount] *= (strtod(token[grpCount+1], NULL) / 100);    //

                openFile.getline(buf, 512);
                token[0] = strtok( buf, " " );         // start loop again
            }
        }

        gsCount++;
        if (gsCount >= ngsz)
            Found = 1;
    }

    // Take cumulative data and turn it into normalized fractions
    for (grpCount = 0; grpCount < (ngrp); grpCount++)
        for (int gsCount = ngsz; gsCount > 0; gsCount--)
        {
            grp[grpCount].pct[0][gsCount] -= grp[grpCount].pct[0][gsCount - 1];
            grp[grpCount].pct[1][gsCount] -= grp[grpCount].pct[1][gsCount - 1];
            grp[grpCount].pct[2][gsCount] -= grp[grpCount].pct[2][gsCount - 1];
        }
        
    // Carry out substrate shift; for randomization work

    for (grpCount = 0; grpCount < ngrp; grpCount++)
    {
        for (int gsCount = 0; gsCount < ngsz; gsCount++)
        {
            for (lithCount = 0; lithCount < nlith; lithCount++)                           // last term is a sand content addition
            {
                if ( gsCount == 0 )
                    qtemp.pct[lithCount][gsCount] = N[2] * grp[grpCount].pct[gsCount][lithCount] + N[3] * grp[grpCount].pct[lithCount][gsCount+1]
                        + N[4] * grp[grpCount].pct[lithCount][gsCount+2];
                else if ( gsCount == 1 )
                    qtemp.pct[lithCount][gsCount] = N[1]*grp[grpCount].pct[lithCount][gsCount-1]
                        + N[2]*grp[grpCount].pct[lithCount][gsCount] + N[3]*grp[grpCount].pct[lithCount][gsCount+1]
                        + N[4]*grp[grpCount].pct[lithCount][gsCount+2];
                else if ( gsCount == ngsz-2 )
                    qtemp.pct[lithCount][gsCount]= N[0]*grp[grpCount].pct[lithCount][gsCount-2]+ N[1]*grp[grpCount].pct[lithCount][gsCount-1]
                        + N[2]*grp[grpCount].pct[lithCount][gsCount] + N[3]*grp[grpCount].pct[lithCount][gsCount+1];
                else if ( gsCount == ngsz-1 )
                    qtemp.pct[lithCount][gsCount]= N[0]*grp[grpCount].pct[lithCount][gsCount-2]+ N[1]*grp[grpCount].pct[lithCount][gsCount-1]
                        + N[2]*grp[grpCount].pct[lithCount][gsCount];
                else
                    qtemp.pct[lithCount][gsCount]= N[0]*grp[grpCount].pct[lithCount][gsCount-2]+ N[1]*grp[grpCount].pct[lithCount][gsCount-1]
                        + N[2]*grp[grpCount].pct[lithCount][gsCount] + N[3]*grp[grpCount].pct[lithCount][gsCount+1]
                       + N[4]*grp[grpCount].pct[lithCount][gsCount+2];
            }
        }

        qtemp.norm_frac();
        for (int gsCount = 0; gsCount < ngsz; gsCount++)
            for (lithCount = 0; lithCount < nlith; lithCount++)
                grp[grpCount].pct[lithCount][gsCount] = qtemp.pct[lithCount][gsCount];
        grp[grpCount].dg_and_std();
    }
    
}

void RiverProfile::getLongProfile(ifstream &openFile)
{

    const char* token[40] = {};                // initialize to 0; 40 tokens max
    char buf[512];                             // Max 512 chars per line
    int m, n = 0;
    int Found = 0;

    while (Found == 0)
    {
        openFile.getline(buf, 512);
        token[0] = strtok( buf, " " );                    // first token
        if (token[0] == NULL || strcmp(token[0], "!") == 0)       // zero if line is blank or "!"
            continue;
        else
        {
            for (m = 0; m < npts; m++)                  // Loop through grid points
            {
                for (n = 1; n < 10; n++)
                {
                    token[n] = strtok(0, " ");
                        if (!token[n]) break;            // no more tokens
                }

                xx[m] = atof(token[0]);
                eta[m] = atof(token[1]);
                bedrock[m] = atof(token[2]);
                if (bedrock[m] > eta[m])
                        (bedrock[m] = eta[m]);       // bedrock must be at, or lower than, initial bed

                RiverXS[m].node = m;
                RiverXS[m].CHList[0].width = atof(token[3]);
                RiverXS[m].fpWidth = atof(token[4]);
/*                if ( HmaxTweak < 0.5 )
                    RiverXS[m].Hmax = atof(token[4]) + ( HmaxTweak * 2 - 0.5 );    // Add height in the range [-0.5 to +0.5]
                else
                    RiverXS[m].Hmax = (HmaxTweak - 0.5) * 3.5 + 0.75; */             // Uniform range from 0.75 to 2.5
                RiverXS[m].CHList[0].Hmax = atof(token[5]);
                RiverXS[m].CHList[0].bankHeight = RiverXS[m].CHList[0].Hmax + 1;  // initial guess
                RiverXS[m].CHList[0].theta = atof(token[6]);
                algrp[m] = atof(token[7]) - 1;
                stgrp[m] = atof(token[8]) - 1;      // Ignoring STGRP for now

                openFile.getline(buf, 512);         // proceed to next line
                token[0] = strtok( buf, " " );
            }
            Found = 1;
        }
    }
}

void RiverProfile::getStratigraphy(ifstream &openFile)
{
    const char* token[150] = {};              // initialize to 0; 40 tokens max
    char buf[512];                           // Max 512 chars per line
    int m, n, i, j, k, idx;
    int Found = 0;

    idx = 0;

    while (Found == 0)
    {
        openFile.getline(buf, 512);
        token[0] = strtok( buf, " " );                    // first token
        if (token[0] == NULL || strcmp(token[0], "!") == 0)       // zero if line is blank or "!"
            continue;
        else
            Found = 1;
    }

    for (m = 0; m < nlayer; m++)
    {

        for (n = 1; n < 150; n++)
        {
            token[n] = strtok(0, " ");
                if (!token[n]) break;            // no more tokens
        }

        for (i = 0; i < npts; i++)
        {
            idx = atof(token[i]);
            for (j = 0; j < ngsz; j++)
                for (k = 0; k < nlith; k++)
                    storedf[i][m].pct[k][j] = grp[idx-1].pct[k][j];  // 'idx-1' because of C++ indexing
            token[i] = NULL;
        }

        openFile.getline(buf, 512);         // proceed to next line
        token[0] = strtok( buf, " " );
    }

    for (i = 0; i < nnodes; i++)            // Populate initial active layer bed GSD
    {
        for (j = 0; j < ngsz; j++)
            for (k = 0; k < nlith; k++)
                F[i].pct[k][j] = grp[algrp[i]].pct[k][j];
        F[i].norm_frac();
        F[i].dg_and_std();
        F[i].abrasion[0] = randAbr;
        F[i].abrasion[1] = randAbr;
        F[i].abrasion[2] = randAbr;
    }

}

