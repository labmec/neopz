#-------------------------------------------------
#
# Project created by QtCreator 2013-05-06T17:18:39
#
#-------------------------------------------------

#Begin FAD
INCLUDEPATH += "/home/felps/FAD/"
INCLUDEPATH += "/home/felps/FAD/Fad"
INCLUDEPATH += "/home/felps/FAD/TinyFad"
INCLUDEPATH += "/home/felps/FAD/TinyFadET"
#End FAD

#Begin PZ
INCLUDEPATH += "/home/felps/neopz_lib/include"
LIBS += -L"/home/felps/neopz_lib/" -lpz
DEFINES += PZSOURCEDIR=\"/home/felps/neopz_code\" REFPATTERNDIR=\"/home/felps/neopz_code\"
DEFINES += USING_BOOST _AUTODIFF LOG4CXX USING_METIS BUILD_UNITTESTING BUILD_TUTORIAL REALdouble STATEdouble
#End PZ

#Begin LOG4CXX
INCLUDEPATH += "/usr/local/include/log4cxx"
LIBS += -L"/usr/local/include/" -llog4cxx
#End LOG4CXX

#Begin METIS
INCLUDEPATH += "/usr/local/include/"
LIBS += -L"/usr/local/include/" -lmetis
#End METIS

#Begin QWT
INCLUDEPATH += "/usr/local/include/"
LIBS += -L"/usr/local/include/" -lqwt
#End QWT


QT       += core gui

TARGET = Plasticity
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    plot.cpp \
    simulation.cpp \
    initialpointdock.cpp \
    canvaspicker.cpp \
    ../TPZPlasticityTest.cpp


HEADERS  += mainwindow.h \
    plot.h \
    simulation.h \
    initialpointdock.h \
    canvaspicker.h \
    common.h \
../TPZPlasticityTest.h

FORMS    += mainwindow.ui \
    initialpointdock.ui
