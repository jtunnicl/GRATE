/********************************************************************************
** Form generated from reading UI file 'RwaveWin.ui'
**
** Created by: Qt User Interface Compiler version 5.13.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_RWAVEWIN_H
#define UI_RWAVEWIN_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDateTimeEdit>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>
#include "qcustomplot.h"

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QAction *action_load_XML;
    QWidget *centralWidget;
    QCustomPlot *VectorPlot;
    QDateTimeEdit *grateDateTime;
    QPushButton *startButton;
    QTabWidget *tabWidget;
    QWidget *Bedload;
    QCustomPlot *BedloadPlot;
    QWidget *BankWidth;
    QCustomPlot *BankWidthPlot;
    QWidget *XSect;
    QCustomPlot *XSectPlot;
    QLabel *label_4;
    QDoubleSpinBox *spinWidth;
    QLabel *label_5;
    QDoubleSpinBox *spinBankHt;
    QSpinBox *spinNode;
    QLabel *label_6;
    QDoubleSpinBox *spinDepth;
    QLabel *label_7;
    QLabel *label_8;
    QLabel *label_9;
    QSpinBox *spinNoChnl;
    QDoubleSpinBox *spinDcomp;
    QDoubleSpinBox *spinD50;
    QLabel *label_10;
    QLabel *label_11;
    QDoubleSpinBox *spinD90;
    QDoubleSpinBox *spinTheta;
    QLabel *label_12;
    QDoubleSpinBox *spinHmax;
    QLabel *label_13;
    QProgressBar *runProgress;
    QDoubleSpinBox *reportQw;
    QLabel *label;
    QDoubleSpinBox *reportStep;
    QLabel *label_2;
    QCustomPlot *QwSeries;
    QCustomPlot *GSD_Dash;
    QGroupBox *groupBox;
    QSlider *deltaT;
    QDoubleSpinBox *dt_disp;
    QDoubleSpinBox *writeInt_disp;
    QLabel *label_16;
    QTextEdit *outputFileName;
    QLabel *label_17;
    QTextEdit *textFileName;
    QPushButton *RegimeButton;
    QPushButton *stopRun;
    QDoubleSpinBox *reportQs;
    QLabel *label_3;
    QLabel *label_14;
    QLabel *label_15;
    QLabel *loadingAdvice;
    QMenuBar *menuBar;
    QMenu *menuFile;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->setEnabled(true);
        MainWindow->resize(1102, 667);
        action_load_XML = new QAction(MainWindow);
        action_load_XML->setObjectName(QString::fromUtf8("action_load_XML"));
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        VectorPlot = new QCustomPlot(centralWidget);
        VectorPlot->setObjectName(QString::fromUtf8("VectorPlot"));
        VectorPlot->setGeometry(QRect(170, 50, 551, 321));
        grateDateTime = new QDateTimeEdit(centralWidget);
        grateDateTime->setObjectName(QString::fromUtf8("grateDateTime"));
        grateDateTime->setGeometry(QRect(170, 20, 161, 22));
        grateDateTime->setReadOnly(false);
        startButton = new QPushButton(centralWidget);
        startButton->setObjectName(QString::fromUtf8("startButton"));
        startButton->setEnabled(false);
        startButton->setGeometry(QRect(20, 50, 90, 22));
        tabWidget = new QTabWidget(centralWidget);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tabWidget->setGeometry(QRect(170, 380, 551, 231));
        Bedload = new QWidget();
        Bedload->setObjectName(QString::fromUtf8("Bedload"));
        BedloadPlot = new QCustomPlot(Bedload);
        BedloadPlot->setObjectName(QString::fromUtf8("BedloadPlot"));
        BedloadPlot->setGeometry(QRect(0, 0, 541, 191));
        tabWidget->addTab(Bedload, QString());
        BankWidth = new QWidget();
        BankWidth->setObjectName(QString::fromUtf8("BankWidth"));
        BankWidthPlot = new QCustomPlot(BankWidth);
        BankWidthPlot->setObjectName(QString::fromUtf8("BankWidthPlot"));
        BankWidthPlot->setGeometry(QRect(0, 0, 541, 191));
        tabWidget->addTab(BankWidth, QString());
        XSect = new QWidget();
        XSect->setObjectName(QString::fromUtf8("XSect"));
        XSectPlot = new QCustomPlot(XSect);
        XSectPlot->setObjectName(QString::fromUtf8("XSectPlot"));
        XSectPlot->setGeometry(QRect(120, 0, 421, 191));
        label_4 = new QLabel(XSect);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setGeometry(QRect(0, 120, 71, 16));
        spinWidth = new QDoubleSpinBox(XSect);
        spinWidth->setObjectName(QString::fromUtf8("spinWidth"));
        spinWidth->setGeometry(QRect(0, 54, 51, 22));
        spinWidth->setReadOnly(true);
        spinWidth->setButtonSymbols(QAbstractSpinBox::NoButtons);
        spinWidth->setMaximum(999.990000000000009);
        label_5 = new QLabel(XSect);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setGeometry(QRect(0, 40, 57, 14));
        spinBankHt = new QDoubleSpinBox(XSect);
        spinBankHt->setObjectName(QString::fromUtf8("spinBankHt"));
        spinBankHt->setGeometry(QRect(0, 135, 51, 23));
        spinBankHt->setReadOnly(true);
        spinBankHt->setButtonSymbols(QAbstractSpinBox::NoButtons);
        spinNode = new QSpinBox(XSect);
        spinNode->setObjectName(QString::fromUtf8("spinNode"));
        spinNode->setGeometry(QRect(0, 13, 51, 23));
        spinNode->setMaximum(130);
        spinNode->setValue(45);
        label_6 = new QLabel(XSect);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setGeometry(QRect(0, 0, 57, 14));
        spinDepth = new QDoubleSpinBox(XSect);
        spinDepth->setObjectName(QString::fromUtf8("spinDepth"));
        spinDepth->setGeometry(QRect(0, 95, 51, 23));
        spinDepth->setReadOnly(true);
        spinDepth->setButtonSymbols(QAbstractSpinBox::NoButtons);
        spinDepth->setMaximum(999.990000000000009);
        label_7 = new QLabel(XSect);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setGeometry(QRect(0, 80, 57, 14));
        label_8 = new QLabel(XSect);
        label_8->setObjectName(QString::fromUtf8("label_8"));
        label_8->setGeometry(QRect(60, 40, 57, 14));
        label_9 = new QLabel(XSect);
        label_9->setObjectName(QString::fromUtf8("label_9"));
        label_9->setGeometry(QRect(60, 80, 71, 16));
        spinNoChnl = new QSpinBox(XSect);
        spinNoChnl->setObjectName(QString::fromUtf8("spinNoChnl"));
        spinNoChnl->setGeometry(QRect(60, 14, 51, 22));
        spinNoChnl->setReadOnly(true);
        spinNoChnl->setButtonSymbols(QAbstractSpinBox::NoButtons);
        spinDcomp = new QDoubleSpinBox(XSect);
        spinDcomp->setObjectName(QString::fromUtf8("spinDcomp"));
        spinDcomp->setGeometry(QRect(60, 95, 51, 23));
        spinDcomp->setReadOnly(true);
        spinDcomp->setButtonSymbols(QAbstractSpinBox::NoButtons);
        spinDcomp->setMaximum(999.990000000000009);
        spinD50 = new QDoubleSpinBox(XSect);
        spinD50->setObjectName(QString::fromUtf8("spinD50"));
        spinD50->setGeometry(QRect(60, 55, 51, 23));
        spinD50->setReadOnly(true);
        spinD50->setButtonSymbols(QAbstractSpinBox::NoButtons);
        spinD50->setMaximum(999.990000000000009);
        label_10 = new QLabel(XSect);
        label_10->setObjectName(QString::fromUtf8("label_10"));
        label_10->setGeometry(QRect(60, 0, 57, 14));
        label_11 = new QLabel(XSect);
        label_11->setObjectName(QString::fromUtf8("label_11"));
        label_11->setGeometry(QRect(60, 120, 71, 16));
        spinD90 = new QDoubleSpinBox(XSect);
        spinD90->setObjectName(QString::fromUtf8("spinD90"));
        spinD90->setGeometry(QRect(60, 135, 51, 23));
        spinD90->setReadOnly(true);
        spinD90->setButtonSymbols(QAbstractSpinBox::NoButtons);
        spinD90->setMaximum(999.990000000000009);
        spinTheta = new QDoubleSpinBox(XSect);
        spinTheta->setObjectName(QString::fromUtf8("spinTheta"));
        spinTheta->setGeometry(QRect(60, 180, 51, 23));
        spinTheta->setReadOnly(true);
        spinTheta->setButtonSymbols(QAbstractSpinBox::NoButtons);
        spinTheta->setMaximum(999.990000000000009);
        label_12 = new QLabel(XSect);
        label_12->setObjectName(QString::fromUtf8("label_12"));
        label_12->setGeometry(QRect(60, 160, 71, 16));
        spinHmax = new QDoubleSpinBox(XSect);
        spinHmax->setObjectName(QString::fromUtf8("spinHmax"));
        spinHmax->setGeometry(QRect(0, 180, 51, 23));
        spinHmax->setReadOnly(true);
        spinHmax->setButtonSymbols(QAbstractSpinBox::NoButtons);
        spinHmax->setMaximum(999.990000000000009);
        label_13 = new QLabel(XSect);
        label_13->setObjectName(QString::fromUtf8("label_13"));
        label_13->setGeometry(QRect(0, 160, 71, 16));
        tabWidget->addTab(XSect, QString());
        runProgress = new QProgressBar(centralWidget);
        runProgress->setObjectName(QString::fromUtf8("runProgress"));
        runProgress->setGeometry(QRect(730, 20, 171, 23));
        runProgress->setValue(0);
        reportQw = new QDoubleSpinBox(centralWidget);
        reportQw->setObjectName(QString::fromUtf8("reportQw"));
        reportQw->setGeometry(QRect(30, 160, 81, 22));
        reportQw->setMaximum(9999.989999999999782);
        reportQw->setSingleStep(0.010000000000000);
        label = new QLabel(centralWidget);
        label->setObjectName(QString::fromUtf8("label"));
        label->setGeometry(QRect(30, 140, 91, 16));
        reportStep = new QDoubleSpinBox(centralWidget);
        reportStep->setObjectName(QString::fromUtf8("reportStep"));
        reportStep->setGeometry(QRect(30, 260, 81, 22));
        reportStep->setDecimals(0);
        reportStep->setMaximum(10000000.000000000000000);
        reportStep->setSingleStep(1.000000000000000);
        label_2 = new QLabel(centralWidget);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setGeometry(QRect(30, 240, 91, 16));
        QwSeries = new QCustomPlot(centralWidget);
        QwSeries->setObjectName(QString::fromUtf8("QwSeries"));
        QwSeries->setGeometry(QRect(730, 50, 351, 201));
        GSD_Dash = new QCustomPlot(centralWidget);
        GSD_Dash->setObjectName(QString::fromUtf8("GSD_Dash"));
        GSD_Dash->setGeometry(QRect(730, 260, 351, 201));
        groupBox = new QGroupBox(centralWidget);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        groupBox->setGeometry(QRect(730, 470, 351, 141));
        deltaT = new QSlider(groupBox);
        deltaT->setObjectName(QString::fromUtf8("deltaT"));
        deltaT->setGeometry(QRect(10, 20, 331, 16));
        deltaT->setMaximum(60);
        deltaT->setValue(12);
        deltaT->setSliderPosition(12);
        deltaT->setTracking(true);
        deltaT->setOrientation(Qt::Horizontal);
        deltaT->setTickPosition(QSlider::TicksBelow);
        deltaT->setTickInterval(10);
        dt_disp = new QDoubleSpinBox(groupBox);
        dt_disp->setObjectName(QString::fromUtf8("dt_disp"));
        dt_disp->setGeometry(QRect(10, 50, 62, 22));
        dt_disp->setReadOnly(true);
        writeInt_disp = new QDoubleSpinBox(groupBox);
        writeInt_disp->setObjectName(QString::fromUtf8("writeInt_disp"));
        writeInt_disp->setGeometry(QRect(10, 110, 62, 22));
        writeInt_disp->setReadOnly(false);
        writeInt_disp->setDecimals(0);
        writeInt_disp->setMaximum(10000.000000000000000);
        writeInt_disp->setValue(100.000000000000000);
        label_16 = new QLabel(groupBox);
        label_16->setObjectName(QString::fromUtf8("label_16"));
        label_16->setGeometry(QRect(10, 90, 91, 16));
        outputFileName = new QTextEdit(groupBox);
        outputFileName->setObjectName(QString::fromUtf8("outputFileName"));
        outputFileName->setGeometry(QRect(110, 70, 231, 61));
        outputFileName->setOverwriteMode(true);
        label_17 = new QLabel(groupBox);
        label_17->setObjectName(QString::fromUtf8("label_17"));
        label_17->setGeometry(QRect(110, 50, 221, 16));
        textFileName = new QTextEdit(centralWidget);
        textFileName->setObjectName(QString::fromUtf8("textFileName"));
        textFileName->setGeometry(QRect(340, 20, 381, 21));
        RegimeButton = new QPushButton(centralWidget);
        RegimeButton->setObjectName(QString::fromUtf8("RegimeButton"));
        RegimeButton->setEnabled(true);
        RegimeButton->setGeometry(QRect(20, 410, 91, 23));
        RegimeButton->setCheckable(true);
        RegimeButton->setChecked(false);
        stopRun = new QPushButton(centralWidget);
        stopRun->setObjectName(QString::fromUtf8("stopRun"));
        stopRun->setGeometry(QRect(20, 80, 90, 22));
        reportQs = new QDoubleSpinBox(centralWidget);
        reportQs->setObjectName(QString::fromUtf8("reportQs"));
        reportQs->setGeometry(QRect(30, 210, 81, 22));
        reportQs->setDecimals(5);
        reportQs->setMaximum(9999.989999999999782);
        reportQs->setSingleStep(0.010000000000000);
        label_3 = new QLabel(centralWidget);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setGeometry(QRect(30, 190, 91, 16));
        label_14 = new QLabel(centralWidget);
        label_14->setObjectName(QString::fromUtf8("label_14"));
        label_14->setGeometry(QRect(20, 440, 91, 16));
        label_15 = new QLabel(centralWidget);
        label_15->setObjectName(QString::fromUtf8("label_15"));
        label_15->setGeometry(QRect(20, 460, 111, 16));
        loadingAdvice = new QLabel(centralWidget);
        loadingAdvice->setObjectName(QString::fromUtf8("loadingAdvice"));
        loadingAdvice->setGeometry(QRect(20, 30, 111, 16));
        MainWindow->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(MainWindow);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1102, 21));
        menuFile = new QMenu(menuBar);
        menuFile->setObjectName(QString::fromUtf8("menuFile"));
        MainWindow->setMenuBar(menuBar);
        mainToolBar = new QToolBar(MainWindow);
        mainToolBar->setObjectName(QString::fromUtf8("mainToolBar"));
        MainWindow->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(MainWindow);
        statusBar->setObjectName(QString::fromUtf8("statusBar"));
        MainWindow->setStatusBar(statusBar);

        menuBar->addAction(menuFile->menuAction());
        menuFile->addAction(action_load_XML);

        retranslateUi(MainWindow);

        tabWidget->setCurrentIndex(1);


        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QCoreApplication::translate("MainWindow", "MainWindow", nullptr));
        action_load_XML->setText(QCoreApplication::translate("MainWindow", "Load XML...", nullptr));
        grateDateTime->setDisplayFormat(QCoreApplication::translate("MainWindow", "yyyy/MM/dd hh:mm AP", nullptr));
        startButton->setText(QCoreApplication::translate("MainWindow", "Start Run", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(Bedload), QCoreApplication::translate("MainWindow", "Bedload", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(BankWidth), QCoreApplication::translate("MainWindow", "Bank W", nullptr));
        label_4->setText(QCoreApplication::translate("MainWindow", "Bank Ht", nullptr));
        label_5->setText(QCoreApplication::translate("MainWindow", "Bed Width", nullptr));
        label_6->setText(QCoreApplication::translate("MainWindow", "Node", nullptr));
        label_7->setText(QCoreApplication::translate("MainWindow", "Depth", nullptr));
        label_8->setText(QCoreApplication::translate("MainWindow", "D50 (mm)", nullptr));
        label_9->setText(QCoreApplication::translate("MainWindow", "D Cmp (mm)", nullptr));
        label_10->setText(QCoreApplication::translate("MainWindow", "NoChannels", nullptr));
        label_11->setText(QCoreApplication::translate("MainWindow", "D90", nullptr));
        label_12->setText(QCoreApplication::translate("MainWindow", "Theta", nullptr));
        label_13->setText(QCoreApplication::translate("MainWindow", "Hmax", nullptr));
        tabWidget->setTabText(tabWidget->indexOf(XSect), QCoreApplication::translate("MainWindow", "Cross Section", nullptr));
        label->setText(QCoreApplication::translate("MainWindow", "Discharge m3/s", nullptr));
        label_2->setText(QCoreApplication::translate("MainWindow", "Time Step", nullptr));
        groupBox->setTitle(QCoreApplication::translate("MainWindow", "Delta t (time step in seconds)", nullptr));
        label_16->setText(QCoreApplication::translate("MainWindow", "Write Interval", nullptr));
        outputFileName->setPlaceholderText(QCoreApplication::translate("MainWindow", "RunResults.txt", nullptr));
        label_17->setText(QCoreApplication::translate("MainWindow", "Results File (.txt format)", nullptr));
#if QT_CONFIG(tooltip)
        RegimeButton->setToolTip(QString());
#endif // QT_CONFIG(tooltip)
        RegimeButton->setText(QCoreApplication::translate("MainWindow", "Regime Width", nullptr));
        stopRun->setText(QCoreApplication::translate("MainWindow", "Stop Run", nullptr));
        label_3->setText(QCoreApplication::translate("MainWindow", "Bedload m3/s", nullptr));
        label_14->setText(QCoreApplication::translate("MainWindow", "Auto-adjust bank", nullptr));
        label_15->setText(QCoreApplication::translate("MainWindow", "width (experimental)", nullptr));
        loadingAdvice->setText(QCoreApplication::translate("MainWindow", "Load file before start", nullptr));
        menuFile->setTitle(QCoreApplication::translate("MainWindow", "File", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_RWAVEWIN_H
