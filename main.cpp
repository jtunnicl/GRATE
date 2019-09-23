#include <QApplication>
#include <QWidget>
#include <QDir>
#include "mainwindow.h"
#include "riverprofile.h"

int main(int argc, char *argv[])
{
  //QApplication::setGraphicsSystem("raster");
  QApplication a(argc, argv);
  QCoreApplication::setAttribute(Qt::AA_DontUseNativeMenuBar);
  MainWindow w;
  w.show();

  return a.exec();
}
