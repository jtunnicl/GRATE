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
    model.cpp \
    hydro.cpp \
    qcustomplot.cpp \
    riverprofile.cpp \
    sed.cpp \
    tinyxml2_wrapper.cpp \
    tinyxml2/tinyxml2.cpp

HEADERS  += mainwindow.h \
    model.h \
    hydro.h \
    qcustomplot.h \
    riverprofile.h \
    sed.h \
    tinyxml2_wrapper.h \
    tinyxml2/tinyxml2.h

FORMS    += \
    RwaveWin.ui
