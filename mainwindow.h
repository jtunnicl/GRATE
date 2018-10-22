#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDateTime>
#include <QTimer>
#include "qcustomplot.h"
#include "riverprofile.h"
#include "hydro.h"
#include "sed.h"

#include <QtCore>
#include <QWidget>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void setupChart();
    void stepTime();
    void writeResults(int count);

public slots:
    void kernel();
    void modelUpdate();

signals:
    void drawChart();

private:
    Ui::MainWindow *ui;
    RiverProfile *rn;
    hydro *wl;
    sed *sd;
    QString demoName;
    QTimer dataTimer;
    QCPItemTracer *itemDemoPhaseTracer;
    int currentDemoIndex;
};


#endif // MAINWINDOW_H
