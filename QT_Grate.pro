QT       += core gui
QT       += printsupport
QT       += widgets

CONFIG += c++11

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

TARGET = QT_Grate
TEMPLATE = app


SOURCES += \
    gratetime.cpp \
    hydro.cpp \
    main.cpp \
    mainwindow.cpp \
    model.cpp \
    qcustomplot.cpp \
    riverprofile.cpp \
    sed.cpp \
    tinyxml2/tinyxml2.cpp \
    tinyxml2_wrapper.cpp
    mainwindow.cpp
    gratetime.cpp
    hydro.cpp
    model.cpp
    qcustomplot.cpp
    riverprofile.cpp
    sed.cpp
    tinyxml2_wrapper.cpp

HEADERS += \
    gratetime.h \
    hydro.h \
    mainwindow.h \
    model.h \
    qcustomplot.h \
    riverprofile.h \
    sed.h \
    tinyxml2/tinyxml2.h \
    tinyxml2_wrapper.h \
    ui_RwaveWin.h
    gratetime.h
    hydro.h
    model.h
    qcustomplot.h
    riverprofile.h
    sed.h
    tinyxml2_wrapper.h

FORMS += \
    RwaveWin.ui \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
