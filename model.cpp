#include "model.h"
#include "riverprofile.h"
#include "hydro.h"
#include "sed.h"
#include "tinyxml2/tinyxml2.h"
#include <iostream>
#include <fstream>

using namespace tinyxml2;

Model::Model(XMLElement* params_root, string out1) :
    rn(nullptr), wl(nullptr), sd(nullptr)
{
    rn = new RiverProfile(params_root);  // Long profile, channel geometry
    wl = new hydro(rn, params_root);  // Channel hydraulic parameters
    sd = new sed(rn, params_root);

    // initialise
    rn->cTime = wl->Qw[0][0].date_time;
    rn->startTime = wl->Qw[0][0].date_time;
    rn->endTime = wl->Qw[0][wl->Qw[0].size() - 1].date_time;
    rn->writeInterval = 100;  // CDJS: set to something small to get output for checking results
    rn->outputFile = out1;
    writeResults(0);
}

Model::~Model() {
    delete rn;
    delete wl;
    delete sd;
}

void Model::iteration() {
    wl->backWater(rn);
    sd->computeTransport(rn);
    stepTime();
    rn->qwTweak = rn->tweakArray[rn->yearCounter];

    if ( ( rn->regimeFlag == 1 ) && (rn->counter % 4 == 0) && ( rn->qwTweak < 1 ) )
            wl->setRegimeWidth(rn);         // kick off regime restraints, once hydraulics are working

    if (rn->counter % rn->writeInterval == 0) {
        writeResults(rn->counter);
    }
}

void Model::stepTime(){
    rn->cTime.addSecs(rn->dt);
    rn->counter++;
    rn->yearCounter++;
    //if (rn->yearCounter > 899) {
    //    rn->yearCounter = 0;
    //}
}

void Model::writeResults(int count){

    int i = 0;

    if (count == 0)
    {
        ofstream outDatFile;
        outDatFile.open(rn->outputFile);

        outDatFile << "Output file for program Grate_NESI" << '\n' <<
        "there are twenty-four columns in the output.  they are:" << '\n' <<
        "column no. 1:  X coordinates in meters" << '\n' <<
        "column no. 2:  Bed elevation in meters"  << '\n' <<
        "column no. 3:  Flow depth in meters"  << '\n' <<
        "column no. 4:  Channel width (m)"  << '\n' <<
        "column no. 5:  Channel theta (deg)"  << '\n' <<
        "column no. 6:  Number of channels"  << '\n' <<
        "column no. 7:  Geometric mean grain size (mm) below the surface layer"  << '\n' <<
        "column no. 8:  Geometric mean grain size (mm) of the surface layer"  << '\n' <<
        "column no. 9:  Standard deviation at the same position."  << '\n' <<
        "column no. 10:  Sediment transport rate (m2/s)"  << '\n' <<
        "column no. 11:  Sand percentage (Fs)"  << '\n' <<
        "column no. 12-24: Surface grain size matrix (12 classes)"  << '\n' <<
        "qwTweak = " << rn->qwTweak << '\n' <<
        "qsTweak = " << rn->qsTweak << '\n' <<
        "substrDial = " << rn->substrDial << '\n' <<
        "feedQw = " << rn->feedQw << '\n' <<
        "feedQs = " << rn->feedQs << '\n' <<
        "HmaxTweak = " << rn->HmaxTweak << '\n' <<
        "randAbr = " << rn->randAbr << '\n' <<
        "" << '\n';
        
        outDatFile << '\n';
        outDatFile << "Count:  " << rn->counter << '\n';

        for ( i = 0; i < rn->nnodes; i++ )
        {
            outDatFile << rn->xx[i] << '\t' <<
            rn->eta[i] << '\t' <<
            rn->RiverXS[i].depth << '\t' <<
            rn->RiverXS[i].width << '\t' <<
            rn->RiverXS[i].theta << '\t' <<
            rn->RiverXS[i].noChannels << '\t' <<
            rn->storedf[i][rn->ntop[i]].dsg << '\t' <<
            rn->F[i].dsg << '\t' <<
            rn->F[i].stdv << '\t' <<
            sd->Qs[i] << '\t' <<
            rn->F[i].sand_pct << '\t' <<
            rn->F[i].pct[0][0] + rn->F[i].pct[1][0] + rn->F[i].pct[2][0] << '\t' <<
            rn->F[i].pct[0][1] + rn->F[i].pct[1][1] + rn->F[i].pct[2][1] << '\t' <<
            rn->F[i].pct[0][2] + rn->F[i].pct[1][2] + rn->F[i].pct[2][2] << '\t' <<
            rn->F[i].pct[0][3] + rn->F[i].pct[1][3] + rn->F[i].pct[2][3] << '\t' <<
            rn->F[i].pct[0][4] + rn->F[i].pct[1][4] + rn->F[i].pct[2][4] << '\t' <<
            rn->F[i].pct[0][5] + rn->F[i].pct[1][5] + rn->F[i].pct[2][5] << '\t' <<
            rn->F[i].pct[0][6] + rn->F[i].pct[1][6] + rn->F[i].pct[2][6] << '\t' <<
            rn->F[i].pct[0][7] + rn->F[i].pct[1][7] + rn->F[i].pct[2][7] << '\t' <<
            rn->F[i].pct[0][8] + rn->F[i].pct[1][8] + rn->F[i].pct[2][8] << '\t' <<
            rn->F[i].pct[0][9] + rn->F[i].pct[1][9] + rn->F[i].pct[2][9] << '\t' <<
            rn->F[i].pct[0][10] + rn->F[i].pct[1][10] + rn->F[i].pct[2][10] << '\t' <<
            rn->F[i].pct[0][11] + rn->F[i].pct[1][11] + rn->F[i].pct[2][11] << '\t' <<
            rn->F[i].pct[0][12] + rn->F[i].pct[1][12] + rn->F[i].pct[2][12] << '\n';
        }
        outDatFile << '\n';
    }
    else   // append records
    {
        ofstream outDatFile;
        outDatFile.open(rn->outputFile, ios::out | ios::app);

        outDatFile << "Count:  " << rn->counter << '\n';

		for ( i = 0; i < rn->nnodes; i++ )
		{
		    outDatFile << rn->xx[i] << '\t' <<
			rn->eta[i] << '\t' <<
			rn->RiverXS[i].depth << '\t' <<
			rn->RiverXS[i].width << '\t' <<
			rn->RiverXS[i].theta << '\t' <<
            rn->RiverXS[i].noChannels << '\t' <<
			rn->storedf[i][rn->ntop[i]].dsg << '\t' <<
			rn->F[i].dsg << '\t' <<
			rn->F[i].stdv << '\t' <<
			sd->Qs[i] << '\t' <<
			rn->F[i].sand_pct << '\t' <<
			rn->F[i].pct[0][0] + rn->F[i].pct[1][0] + rn->F[i].pct[2][0] << '\t' <<
			rn->F[i].pct[0][1] + rn->F[i].pct[1][1] + rn->F[i].pct[2][1] << '\t' <<
			rn->F[i].pct[0][2] + rn->F[i].pct[1][2] + rn->F[i].pct[2][2] << '\t' <<
			rn->F[i].pct[0][3] + rn->F[i].pct[1][3] + rn->F[i].pct[2][3] << '\t' <<
			rn->F[i].pct[0][4] + rn->F[i].pct[1][4] + rn->F[i].pct[2][4] << '\t' <<
			rn->F[i].pct[0][5] + rn->F[i].pct[1][5] + rn->F[i].pct[2][5] << '\t' <<
			rn->F[i].pct[0][6] + rn->F[i].pct[1][6] + rn->F[i].pct[2][6] << '\t' <<
			rn->F[i].pct[0][7] + rn->F[i].pct[1][7] + rn->F[i].pct[2][7] << '\t' <<
			rn->F[i].pct[0][8] + rn->F[i].pct[1][8] + rn->F[i].pct[2][8] << '\t' <<
			rn->F[i].pct[0][9] + rn->F[i].pct[1][9] + rn->F[i].pct[2][9] << '\t' <<
			rn->F[i].pct[0][10] + rn->F[i].pct[1][10] + rn->F[i].pct[2][10] << '\t' <<
			rn->F[i].pct[0][11] + rn->F[i].pct[1][11] + rn->F[i].pct[2][11] << '\t' <<
			rn->F[i].pct[0][12] + rn->F[i].pct[1][12] + rn->F[i].pct[2][12] << '\n';
		}
        outDatFile << '\n';
    }
}
