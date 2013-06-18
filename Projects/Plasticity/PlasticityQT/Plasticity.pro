#-------------------------------------------------
#
# Project created by QtCreator 2013-05-06T17:18:39
#
#-------------------------------------------------

#Begin FAD
INCLUDEPATH += "/dados/GOOGLE_PZ/externallibs/FAD/"
INCLUDEPATH += "/dados/GOOGLE_PZ/externallibs/FAD/Fad"
INCLUDEPATH += "/dados/GOOGLE_PZ/externallibs/FAD/TinyFad"
INCLUDEPATH += "/dados/GOOGLE_PZ/externallibs/FAD/TinyFadET"
#End FAD

#Begin PZ
#INCLUDEPATH += "C:/Users/Raul/Desktop/pzlib/include"
#LIBS += -L"C:/Users/Raul/Desktop/pzlib/lib" -lpz
INCLUDEPATH += "/dados/GOOGLE_PZ/PZLIB/include"
LIBS += -L"/dados/GOOGLE_PZ/PZLIB/" -lpz
DEFINES += PZSOURCEDIR=\"/dados/GOOGLE_PZ/neopz_build_teste\" REFPATTERNDIR=\"/dados/GOOGLE_PZ/neopz_build_teste/Refine/RefPatterns\"
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

QWT_ROOT = /usr/local/qwt-6.0.2-svn/
include( $${QWT_ROOT}/features/qwtconfig.pri )

#    QWT_CONFIG += QwtFramework


QMAKE_RPATHDIR *= $${QWT_ROOT}/lib

contains(QWT_CONFIG, QwtFramework) {

    LIBS      += -F$${QWT_ROOT}/lib/

}
else {

    LIBS      += -L$${QWT_ROOT}/lib
}

IPATH       = $${INCLUDEPATH}
qtAddLibrary(qwt)
INCLUDEPATH = $${IPATH}
INCLUDEPATH += $${QWT_ROOT}/lib/qwt.framework/Headers

#contains(QWT_CONFIG, QwtSvg) {

#    QT += svg
#}
#else {

    DEFINES += QWT_NO_SVG
#}


win32 {
    contains(QWT_CONFIG, QwtDll) {
        DEFINES    += QT_DLL QWT_DLL
    }
}

QT       += core gui

TARGET = Plasticity
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    plot.cpp \
    simulation.cpp \
    initialpointdock.cpp \
    canvaspicker.cpp

HEADERS  += mainwindow.h \
    plot.h \
    simulation.h \
    initialpointdock.h \
    canvaspicker.h \
    common.h

FORMS    += mainwindow.ui \
    initialpointdock.ui
