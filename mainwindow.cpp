/*******************
 *
 *
 *  GRATE 8
 *
 *  Graphical Interface
 *
 *
 *
*********************/


#include "mainwindow.h"
#include "ui_RwaveWin.h"
#include <iostream>
#include <fstream>

#define PI 3.14159265

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow),
    rn(new RiverProfile),                        // Long profile, channel geometry
    wl(new hydro(rn)),                               // Channel hydraulic parameters
    sd(new sed(rn))

{
    ui->setupUi(this);
    rn->cTime = wl->Qw[0][0].date_time;
    rn->startTime = wl->Qw[0][0].date_time;
    rn->endTime = wl->Qw[0][wl->Qw[0].size() - 1].date_time;
    setupChart();                                // Setup GUI graph
    setWindowTitle("Raparapaririki River");
    ui->textFileName->setText("Input_Rip1_equil_1938.dat");
    ui->VectorPlot->replot();
                                               // Use <Refresh Plot> button to start run
    connect(ui->refreshButton, SIGNAL(clicked()), this, SLOT(kernel()), Qt::QueuedConnection);
}

void MainWindow::setupChart(){

    int i, j, k, n = 0;
    float theta_rad, a, b;
    unsigned int bc = 0;
    QVector<double> x( rn->nnodes );
    QVector<double> eta( rn->nnodes );
    QVector<double> WSL( rn->nnodes );
    QVector<double> Froude( rn->nnodes );
    QVector<double> Bedload( rn->nnodes );
    QVector<double> Qw_Plot( rn->nnodes );
    QVector<double> LeftBankLower( rn->nnodes );
    QVector<double> RightBankLower( rn->nnodes );
    QVector<double> LeftBankTop( rn->nnodes );
    QVector<double> RightBankTop( rn->nnodes );
    QVector<double> XsPlotX( 11 );
    QVector<double> XsPlotY( 11 );
    QVector<double> wsXS_X( 2 );
    QVector<double> wsXS_Y( 2 );
    QVector<double> Bottom_X( 2 );
    QVector<double> Bottom_Y( 2 );
    QVector<double> Qw_TS( 900 );
    QVector<double> time( 900 );
    QVector<double> CursorX( 2 );
    QVector<double> CursorY( 2 );
    QVector<double> tmp( 8 );
    QVector<QVector<double> > GSD_Data, GSD_Cumul;

    rn->counter = 0;
    rn->yearCounter = 0;

    tmp.fill( 1, rn -> nnodes );
    GSD_Data.fill( tmp, 7);
    GSD_Cumul.fill(tmp, 13);

    for ( i = 0; i < rn -> nnodes; i++ )
    {
        x[i] = rn->RiverXS[i].node;
        eta[i] = rn->eta[i];
        WSL[i] = rn -> eta[i] + rn->RiverXS[i].depth * 8;    // x8 exaggeration for display
        Froude[i] = rn->ntop[i];                             // RiverXS[i].ustar;
        Bedload[i] = sd->Qs[i];
        Qw_Plot[i] = wl->QwCumul[i] / 100;
        for ( j = 0; j < rn->ngsz; j++ )                     // Make a cumulative dist
        {
            if ( j > 0 )
                GSD_Cumul[j][i] = GSD_Cumul[j-1][i];
            for ( k = 0; k < rn->nlith; k++  )
                GSD_Cumul[j][i] -= rn->F[i].pct[k][j];
            if (GSD_Cumul[j][i] < 0.0001)
                GSD_Cumul[j][i] = 0;
        }
    }

    // dt control
    ui->dt_disp->setValue(rn->dt);
    ui->deltaT->setValue(rn->dt);

    // Setup cross-section graph
    n = ui->spinNode->value();
    theta_rad = rn->RiverXS[n].theta * PI / 180.;
    a = rn->RiverXS[n].bankHeight - rn->RiverXS[n].Hmax;  // Vert & Horiz triangle segments at lower channel.
    b = a / tan(theta_rad);

    XsPlotX[2] = -1.5 * rn->RiverXS[n].fpSlope;
    XsPlotY[2] = 0;
    XsPlotX[3] = 0;
    XsPlotY[3] = -1.5;
    XsPlotX[4] = 0.001;
    XsPlotY[4] = -1.5 - rn->RiverXS[n].Hmax;
    XsPlotX[5] = b;
    XsPlotY[5] = -1.5 - rn->RiverXS[n].bankHeight;
    XsPlotX[6] = b + rn->RiverXS[n].width;
    XsPlotY[6] = XsPlotY[5];
    XsPlotX[7] = XsPlotX[6] + b;
    XsPlotY[7] = XsPlotY[4];
    XsPlotX[8] = XsPlotX[7] + 0.001;
    XsPlotY[8] = -1.5;
    XsPlotX[9] = XsPlotX[8] + 1.5;
    XsPlotY[9] = 0;
    XsPlotX[10] = XsPlotX[9] + 5 * rn->RiverXS[n].valleyWallSlp;
    XsPlotY[10] = 5;

    XsPlotX[1] = XsPlotX[9] - rn->RiverXS[n].fpWidth;
    XsPlotY[1] = 0;
    XsPlotX[0] = XsPlotX[1] - ( 5 * rn->RiverXS[n].valleyWallSlp );
    XsPlotY[0] = 5;

    wsXS_X[0] = 0;
    wsXS_X[1] = XsPlotX[8];
    wsXS_Y[0] = -1.5;
    wsXS_Y[1] = -1.5;

    Bottom_X[0] = XsPlotX[0];
    Bottom_X[1] = XsPlotX[10];
    Bottom_Y[0] = -30;
    Bottom_Y[1] = -30;

    ui->spinNode->setMaximum(rn->nnodes);

    for ( bc = 0; bc < 899; bc++ )
    {
        Qw_TS[bc] = rn->tweakArray[bc];
        time[bc] = bc;
    }

    for ( j = 0; j < 7 ; j++ )
        for ( i = 0; i < rn->nnodes; i++ )
                GSD_Data[j][i] = GSD_Cumul[j*2][i];

    CursorX[0] = rn->yearCounter;
    CursorX[1] = rn->yearCounter;
    CursorY[0] = 0;
    CursorY[1] = 1200;

    ui->VectorPlot->xAxis->setLabel("Distance Downstream");
    ui->VectorPlot->yAxis->setLabel("Elevation");
    ui->VectorPlot->xAxis->setRange(0, rn->nnodes + 1);
    ui->VectorPlot->yAxis->setRange(180, 450);

    ui->VectorPlot->addGraph();  // Orig channel bed elevation plot (static)
    ui->VectorPlot->graph(0)->setData( x, eta );
    ui->VectorPlot->graph(0)->setPen(QPen(Qt::gray));

    ui->VectorPlot->addGraph();  // Channel bed elevation plot
    ui->VectorPlot->graph(1)->setData( x, eta );
    ui->VectorPlot->graph(1)->setPen(QPen(Qt::black));
    ui->VectorPlot->graph(1)->setBrush(QColor(255, 161, 0, 50));

    ui->VectorPlot->addGraph();  // Water surface plot
    ui->VectorPlot->graph(2)->setData( x, WSL );
    ui->VectorPlot->graph(2)->setPen(QPen(Qt::blue));
    ui->VectorPlot->graph(2)->setBrush(QBrush(QColor(0, 0, 255, 20)));
    ui->VectorPlot->graph(2)->setChannelFillGraph(ui->VectorPlot->graph(1));
                                          //(rand() % 254, rand() % 254, rand() % 254)));

    ui->BedloadPlot->xAxis->setLabel("Distance Downstream");
    ui->BedloadPlot->yAxis->setLabel("Qs - Bedload Discharge");
    ui->BedloadPlot->xAxis->setRange(0,rn->nnodes+1);
    ui->BedloadPlot->yAxis->setRange(0, 40);

    ui->BedloadPlot->addGraph();  // Ustar plot
    ui->BedloadPlot->graph(0)->setData( x, Froude );
    ui->BedloadPlot->graph(0)->setPen(QPen(Qt::black));

    ui->BedloadPlot->addGraph();  // Bedload plot
    ui->BedloadPlot->graph(1)->setData( x, Bedload );
    ui->BedloadPlot->graph(1)->setPen(QPen(Qt::red));

    ui->BedloadPlot->addGraph();  // Qw plot
    ui->BedloadPlot->graph(2)->setData( x, Qw_Plot );
    ui->BedloadPlot->graph(2)->setPen(QPen(Qt::blue));

    // BankWidth
    ui->BankWidthPlot->xAxis->setLabel("Distance Downstream");
    ui->BankWidthPlot->yAxis->setLabel("Bank Width");
    ui->BankWidthPlot->xAxis->setRange(0,rn->nnodes + 1);
    ui->BankWidthPlot->yAxis->setRange(-100,200);

    ui->BankWidthPlot->addGraph();  // Left Lower Bank
    ui->BankWidthPlot->graph(0)->setData( x, LeftBankLower );
    ui->BankWidthPlot->graph(0)->setPen(QPen(Qt::black));

    ui->BankWidthPlot->addGraph();  // Right Lower Bank
    ui->BankWidthPlot->graph(1)->setData( x, RightBankLower );
    ui->BankWidthPlot->graph(1)->setPen(QPen(Qt::black));

    ui->BankWidthPlot->addGraph();  // Left Upper Bank
    ui->BankWidthPlot->graph(2)->setData( x, LeftBankTop );
    ui->BankWidthPlot->graph(2)->setPen(QPen(Qt::red));

    ui->BankWidthPlot->addGraph();  // Right Upper Bank
    ui->BankWidthPlot->graph(3)->setData( x, RightBankTop );
    ui->BankWidthPlot->graph(3)->setPen(QPen(Qt::blue));

    ui->BankWidthPlot->addGraph();  // Sinuosity
    ui->BankWidthPlot->graph(4)->setData( x, RightBankTop );
    ui->BankWidthPlot->graph(4)->setPen(QPen(Qt::green));

    ui->XSectPlot->xAxis->setLabel("Cross Section Width");
    ui->XSectPlot->yAxis->setLabel("Relative Elevation");
    ui->XSectPlot->xAxis->setRange(-(rn->RiverXS[n].fpWidth / 2), rn->RiverXS[n].fpWidth / 2 );
    ui->XSectPlot->yAxis->setRange(-10, 12);

    ui->XSectPlot->addGraph();  // Cross-sectional plot
    ui->XSectPlot->graph(0)->setData( XsPlotX, XsPlotY );
    ui->XSectPlot->graph(0)->setPen(QPen(Qt::black));

    ui->XSectPlot->addGraph();  // Water surface plot
    ui->XSectPlot->graph(1)->setData( wsXS_X, wsXS_Y );
    ui->XSectPlot->graph(1)->setPen(QPen(Qt::blue));
    ui->XSectPlot->graph(1)->setBrush(QBrush(QColor(0, 0, 255, 20)));
    ui->XSectPlot->graph(1)->setChannelFillGraph(ui->XSectPlot->graph(0));

    ui->XSectPlot->addGraph();  // Bottom of cross-section (fill)
    ui->XSectPlot->graph(2)->setData( Bottom_X, Bottom_Y );
    ui->XSectPlot->graph(2)->setPen(QPen(Qt::black));
    ui->XSectPlot->graph(2)->setBrush(QColor(255, 161, 0, 50));
    ui->XSectPlot->graph(2)->setChannelFillGraph(ui->XSectPlot->graph(0));

    ui->QwSeries->xAxis->setTickLabelType(QCPAxis::ltDateTime);
    ui->QwSeries->xAxis->setDateTimeFormat("h:m");
    ui->QwSeries->yAxis->setLabel("Discharge");
    ui->QwSeries->xAxis->setLabel("Time");
    ui->QwSeries->xAxis->setRange(0, 900);
    ui->QwSeries->yAxis->setRange(0, 3);

    ui->QwSeries->addGraph();
    ui->QwSeries->graph(0)->setData( time, Qw_TS );
    ui->QwSeries->addGraph();
    ui->QwSeries->graph(1)->setData(CursorX, CursorY);
    ui->QwSeries->graph(1)->setPen(QPen(Qt::red));

    ui->GSD_Dash->xAxis->setLabel("Distance Downstream");
    ui->GSD_Dash->yAxis->setLabel("Cum %");
    ui->GSD_Dash->xAxis->setRange(0, rn->nnodes + 1);
    ui->GSD_Dash->yAxis->setRange(0, 1);

    ui->GSD_Dash->addGraph();  // Grain size distributions
    ui->GSD_Dash->graph(0)->setData( x, GSD_Data[0] );
    ui->GSD_Dash->addGraph();
    ui->GSD_Dash->graph(1)->setData( x, GSD_Data[1] );
    ui->GSD_Dash->addGraph();
    ui->GSD_Dash->graph(2)->setData( x, GSD_Data[2] );
    ui->GSD_Dash->addGraph();
    ui->GSD_Dash->graph(3)->setData( x, GSD_Data[3] );
    ui->GSD_Dash->addGraph();
    ui->GSD_Dash->graph(4)->setData( x, GSD_Data[4] );
    ui->GSD_Dash->addGraph();
    ui->GSD_Dash->graph(5)->setData( x, GSD_Data[5] );
    ui->GSD_Dash->addGraph();
    ui->GSD_Dash->graph(6)->setData( x, GSD_Data[6] );
    ui->GSD_Dash->graph(0)->setPen(QPen(QColor(0,0,0)));  // 0 is the finest grain size category: black
    for ( i = 1; i < 7 ; i++ )
           ui->GSD_Dash->graph(i)->setPen(QPen(QColor((i * 30), 254 - (i * 30), 113)));

    ui->grateDateTime->setDateTime(rn->cTime);
    ui->reportStep->setValue(rn->counter);
    ui->reportYear->setValue(rn->yearCounter);

    ui->spinBankHt->setValue(rn->RiverXS[n].bankHeight);
    ui->spinTheta->setValue(rn->RiverXS[n].theta);
    ui->spinDepth->setValue(rn->RiverXS[n].depth);
    ui->spinNoChnl->setValue(rn->RiverXS[n].noChannels);
    ui->spinD50->setValue(pow(2, rn->F[n].dsg));
    ui->spinDcomp->setValue(pow(2, rn->RiverXS[n].comp_D));
    ui->spinHmax->setValue(rn->RiverXS[n].Hmax);

    ui->VectorPlot->replot();
    ui->BedloadPlot->replot();
    ui->BankWidthPlot->replot();

    writeResults(0);
}

void MainWindow::stepTime(){

    rn->cTime = rn->cTime.addSecs(rn->dt);
    rn->counter++;
    rn->yearCounter++;
    if ( rn->yearCounter > 899 )
        rn->yearCounter = 0;
}

void MainWindow::kernel(){

    ui->grateDateTime->setDateTime(rn->cTime);

    connect(&dataTimer, SIGNAL(timeout()), this, SLOT(modelUpdate()));
    dataTimer.start(0); // Interval 0 means to refresh as fast as possible    
}

void MainWindow::modelUpdate(){

    unsigned int i, j, k, n, inc = 0;
    float theta_rad, a, b, c;
    float topFp, ovBank, ovFp;
    QVector<double> WSL( rn->nnodes );
    QVector<double> x( rn->nnodes );
    QVector<double> eta( rn->nnodes );
    QVector<double> Froude( rn->nnodes );
    QVector<double> Bedload( rn->nnodes );
    QVector<double> Qw_Plot( rn->nnodes );
    QVector<double> LeftBankLower( rn->nnodes );
    QVector<double> RightBankLower( rn->nnodes );
    QVector<double> LeftBankTop( rn->nnodes );
    QVector<double> RightBankTop( rn->nnodes );
    QVector<double> Sinuosity( rn->nnodes );
    QVector<double> XsPlotX( 11 );
    QVector<double> XsPlotY( 11 );
    QVector<double> wsXS_X( 2 );
    QVector<double> wsXS_Y( 2 );
    QVector<double> CursorX( 2 );
    QVector<double> CursorY( 2 );
    QVector<double> PSI( 8 );
    QVector<double> tmp( 8 );
    QVector<QVector<double> > GSD_Data, GSD_Cumul;


    tmp.fill( 1, rn -> nnodes );
    GSD_Data.fill( tmp, 7);
    GSD_Cumul.fill(tmp, 13);

    int prog;   // Model run progress

    wl->backWater(rn);
    sd->computeTransport(rn);
    stepTime();
    rn->qwTweak = rn->tweakArray[rn->yearCounter];

    if ( (rn->counter % 4 == 0) && ( rn->qwTweak < 1 ) )
            wl->setRegimeWidth(rn);         // kick off regime restraints, once hydraulics are working

    // dt control
    rn->dt = ui->deltaT->value();
    ui->dt_disp->setValue(rn->dt);                    // Control dt with slider
    

    // Calculate water profile data points, bank widths
    for ( i = 0; i < rn -> nnodes; i++ )
    {
        x[i] = rn->RiverXS[i].node;
        eta[i] = rn->eta[i];
        WSL[i] = rn->eta[i] + rn->RiverXS[i].depth * 8;
        Froude[i] = rn->ntop[i];               // RiverXS[i].ustar;
        Bedload[i] = sd->Qs[i];

        LeftBankLower[i] = rn->RiverXS[i].width/2;
        RightBankLower[i] = -rn->RiverXS[i].width/2;
        theta_rad = rn->RiverXS[i].theta * PI / 180;
        LeftBankTop[i] = ( rn->RiverXS[i].width + (2 * ( rn->RiverXS[i].bankHeight - rn->RiverXS[i].Hmax) / tan( theta_rad ) ) ) / 2;
        RightBankTop[i] = rn->RiverXS[i].depth * 100;
        Sinuosity[i] = ( rn->RiverXS[i].chSinu - 1 ) * 100;

        Qw_Plot[i] = wl->QwCumul[i] / 100;
        for ( j = 0; j < rn->ngsz; j++ )            // Make a cumulative dist
        {
            if ( j > 0 )
                GSD_Cumul[j][i] = GSD_Cumul[j-1][i];
            for ( k = 0; k < rn->nlith; k++  )
                GSD_Cumul[j][i] -= rn->F[i].pct[k][j];
            if (GSD_Cumul[j][i] < 0.0001)
                GSD_Cumul[j][i] = 0;
        }
    }

    // Setup cross-section graph
    n = ui->spinNode->value();
    theta_rad = rn->RiverXS[n].theta * PI / 180;
    a = rn->RiverXS[n].bankHeight - rn->RiverXS[n].Hmax;  // Vert & Horiz triangle segments at lower channel.
    c = tan(theta_rad);
    b = a / c;              // aka dW, horizontal distance between bed and bank, under toe of channel edge
    inc = 0.5;            // small increment to help with ordering linework.

    XsPlotX[2] = -1.5 * rn->RiverXS[n].fpSlope;
    XsPlotY[2] = 0;
    XsPlotX[3] = 0;
    XsPlotY[3] = -1.5;
    XsPlotX[4] = 0.001;
    XsPlotY[4] = -1.5 - rn->RiverXS[n].Hmax;
    XsPlotX[5] = b;
    XsPlotY[5] = -1.5 - rn->RiverXS[n].bankHeight;
    XsPlotX[6] = b + rn->RiverXS[n].width;
    XsPlotY[6] = XsPlotY[5];
    XsPlotX[7] = XsPlotX[6] + b;
    XsPlotY[7] = XsPlotY[4];
    XsPlotX[8] = XsPlotX[7] + 0.001;
    XsPlotY[8] = -1.5;
    XsPlotX[9] = XsPlotX[8] + 1.5;
    XsPlotY[9] = 0;
    XsPlotX[10] = XsPlotX[9] + 5 * rn->RiverXS[n].valleyWallSlp;
    XsPlotY[10] = 5;

    XsPlotX[1] = XsPlotX[9] - rn->RiverXS[n].fpWidth;
    XsPlotY[1] = 0;
    XsPlotX[0] = XsPlotX[1] - ( 5 * rn->RiverXS[n].valleyWallSlp );
    XsPlotY[0] = 5;

    wsXS_Y[0] = rn->RiverXS[n].depth;


    topFp = rn->RiverXS[n].bankHeight + 1.5;
    if (rn->RiverXS[n].depth > topFp)
    {
        ovFp = rn->RiverXS[n].depth - topFp;
        ovBank = 1.5;
        wsXS_X[0] = XsPlotX[1] - ( ovFp * rn->RiverXS[n].valleyWallSlp );
        wsXS_X[1] = XsPlotX[9] + ( ovFp * rn->RiverXS[n].valleyWallSlp );
        wsXS_Y[0] = ovFp;               // Floodplain elev. is '0' datum
        wsXS_Y[1] = wsXS_Y[0];
    }
    else if (rn->RiverXS[n].depth > rn->RiverXS[n].bankHeight)
    {
        ovBank = rn->RiverXS[n].depth - rn->RiverXS[n].bankHeight;
        wsXS_X[0] = - ( ovBank * rn->RiverXS[n].fpSlope );
        wsXS_X[1] = XsPlotX[8] + ovBank;
        wsXS_Y[0] = -1.5 + ovBank;
        wsXS_Y[1] = wsXS_Y[0];
    }
    else if (rn->RiverXS[n].depth > a)   // 'a' is computed, above, as bottom of vertical banks
    {
        wsXS_X[0] = 0;
        wsXS_X[1] = XsPlotX[8];
        wsXS_Y[0] = -1.5 - (rn->RiverXS[n].bankHeight - rn->RiverXS[n].depth);
        wsXS_Y[1] = wsXS_Y[0];
    }
    else   // otherwise, within the lower trapezoid
    {
        wsXS_X[0] = c;
        wsXS_X[1] = XsPlotX[8] - c;
        wsXS_Y[0] = -1.5 - (rn->RiverXS[n].bankHeight - rn->RiverXS[n].depth);
        wsXS_Y[1] = wsXS_Y[0];
    }

    // Grain size fining plot
    for ( j = 0; j < 7 ; j++ )
        for ( i = 0; i < rn -> nnodes; i++ )
                GSD_Data[j][i] = GSD_Cumul[j*2][i];

    // t = ui->selectNode->value();

    CursorX[0] = rn->yearCounter;
    CursorX[1] = rn->yearCounter;
    CursorY[0] = 0;
    CursorY[1] = 3;

    ui->VectorPlot->graph(1)->clearData();
    ui->VectorPlot->graph(1)->setData(x, eta);
    ui->VectorPlot->graph(2)->clearData();
    ui->VectorPlot->graph(2)->setData(x, WSL);
    //ui->VectorPlot->graph(1)->setPen(QPen(QColor(rand() % 254, rand() % 254, rand() % 254)));
    ui->VectorPlot->xAxis->setLabel("Distance Downstream");
    ui->VectorPlot->yAxis->setLabel("Elevation");
    ui->VectorPlot->xAxis->setRange(0,rn->nnodes + 1);
    ui->VectorPlot->yAxis->setRange(180,450);
    //ui->VectorPlot->setRangeDrag(Qt::Horizontal|Qt::Vertical);
    //ui->VectorPlot->setRangeZoom(Qt::Horizontal|Qt::Vertical);

    ui->BedloadPlot->graph(0)->clearData();
    ui->BedloadPlot->graph(0)->setData( x, Froude );
    ui->BedloadPlot->graph(1)->clearData();
    ui->BedloadPlot->graph(1)->setData( x, Bedload );
    ui->BedloadPlot->graph(2)->clearData();
    ui->BedloadPlot->graph(2)->setData( x, Qw_Plot );

    ui->BedloadPlot->xAxis->setLabel("Distance Downstream");
    ui->BedloadPlot->yAxis->setLabel("Qs - Bedload Discharge");
    ui->BedloadPlot->xAxis->setRange(0, rn->nnodes + 1);
    ui->BedloadPlot->yAxis->setRange(0, 40);

    ui->BankWidthPlot->graph(0)->clearData();
    ui->BankWidthPlot->graph(0)->setData( x, LeftBankLower );
    ui->BankWidthPlot->graph(1)->clearData();
    ui->BankWidthPlot->graph(1)->setData( x, RightBankLower );
    ui->BankWidthPlot->graph(2)->clearData();
    ui->BankWidthPlot->graph(2)->setData( x, LeftBankTop );
    ui->BankWidthPlot->graph(3)->clearData();
    ui->BankWidthPlot->graph(3)->setData( x, RightBankTop );
    ui->BankWidthPlot->graph(4)->clearData();
    ui->BankWidthPlot->graph(4)->setData( x, Sinuosity );

    ui->BankWidthPlot->xAxis->setLabel("Distance Downstream");
    ui->BankWidthPlot->yAxis->setLabel("Bank Width");
    ui->BankWidthPlot->xAxis->setRange(0,rn->nnodes + 1);
    ui->BankWidthPlot->yAxis->setRange(-100,200);

    ui->XSectPlot->xAxis->setRange( -20, 70 );
    ui->XSectPlot->yAxis->setRange( -7, 7 );
    ui->XSectPlot->graph(0)->clearData();
    ui->XSectPlot->graph(0)->setData( XsPlotX, XsPlotY );
    ui->XSectPlot->graph(1)->clearData();
    ui->XSectPlot->graph(1)->setData( wsXS_X, wsXS_Y );

    ui->QwSeries->graph(1)->clearData();
    ui->QwSeries->graph(1)->setData(CursorX, CursorY);

    //ui->GSD_Dash->graph(0)->clearData();
    ui->GSD_Dash->graph(0)->setData( x, GSD_Data[0] );
    ui->GSD_Dash->graph(1)->setData( x, GSD_Data[1] );
    ui->GSD_Dash->graph(2)->setData( x, GSD_Data[2] );
    ui->GSD_Dash->graph(3)->setData( x, GSD_Data[3] );
    ui->GSD_Dash->graph(4)->setData( x, GSD_Data[4] );
    ui->GSD_Dash->graph(5)->setData( x, GSD_Data[5] );
    ui->GSD_Dash->graph(6)->setData( x, GSD_Data[6] );
    ui->GSD_Dash->graph(0)->setPen(QPen(QColor(0,0,0)));  // 0 is the finest grain size category: black
    for ( i = 1; i < 7 ; i++ )
           ui->GSD_Dash->graph(i)->setPen(QPen(QColor((i * 30), 254 - (i * 30), 113)));

    ui->GSD_Dash->xAxis->setRange(0, rn->nnodes + 1);
    ui->GSD_Dash->yAxis->setRange(0, 1);

    ui->VectorPlot->replot();
    ui->BankWidthPlot->replot();
    ui->XSectPlot->replot();
    ui->BedloadPlot->replot();
    ui->QwSeries->replot();
    ui->GSD_Dash->replot();
    ui->grateDateTime->setDateTime(rn->cTime);

    // Progress bar
    prog = (rn->startTime.secsTo(rn->cTime) * 100) / (rn->startTime.secsTo(rn->endTime));
    if (prog < 100)
        ui->runProgress->setValue(prog);

    ui->reportQw->setValue(wl->QwCumul[rn->nnodes-1]);
    ui->reportStep->setValue(rn->counter);
    ui->reportYear->setValue(rn->yearCounter);

    ui->spinBankHt->setValue(rn->RiverXS[n].bankHeight);
    ui->spinTheta->setValue(rn->RiverXS[n].theta);
    ui->spinDepth->setValue(rn->RiverXS[n].depth);
    ui->spinWidth->setValue(rn->RiverXS[n].width);
    ui->spinNoChnl->setValue(rn->RiverXS[n].noChannels);
    ui->spinD50->setValue(pow(2, rn->F[n].dsg));
    ui->spinDcomp->setValue(rn->RiverXS[n].comp_D);
    ui->spinD90->setValue(pow(2, rn->F[n].d90));
    ui->spinHmax->setValue(rn->RiverXS[n].Hmax);

    if (rn->counter % 100 == 0)
        writeResults(rn->counter);

    // temporarily limit to 50 steps for testing
    if (rn->counter == 800) {
        dataTimer.stop();
    }
}

void MainWindow::writeResults(int count){

    int i = 0;


    if (count == 0)
    {
        ofstream outDatFile;
        outDatFile.open("Run_Results.txt");

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
        outDatFile.open("Run_Results.txt", ios::out | ios::app);

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

MainWindow::~MainWindow()
{
    delete ui;
    delete rn;
}

