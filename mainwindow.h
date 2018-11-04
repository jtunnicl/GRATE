#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QDateTime>
#include <QTimer>
#include "qcustomplot.h"
#include "riverprofile.h"
#include "hydro.h"
#include "sed.h"
#include "model.h"

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

public slots:
    void kernel();
    void modelUpdate();

signals:
    void drawChart();

private:
    void showErrorMessage(const char *title, std::stringstream &msg_stream);
    Ui::MainWindow *ui;
    Model *model;
    QString demoName;
    QTimer dataTimer;
    QCPItemTracer *itemDemoPhaseTracer;
    int currentDemoIndex;
    bool initialised;
};


#endif // MAINWINDOW_H
