#-------------------------------------------------
#
# Project created by QtCreator 2015-05-01T11:15:58
#
#-------------------------------------------------

QT       += core gui
QT       += printsupport
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = GrateRip
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    hydro.cpp \
    qcustomplot.cpp \
    riverprofile.cpp \
    sed.cpp

HEADERS  += mainwindow.h \
    hydro.h \
    qcustomplot.h \
    riverprofile.h \
    sed.h

FORMS    += \
    RwaveWin.ui