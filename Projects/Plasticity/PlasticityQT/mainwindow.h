#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <qwt_plot_zoomer.h>
#include "common.h"
#include <limits>

#include "simulation.h"

class QString;
class QModelIndex;
class QwtPicker;
class QwtEventPattern;
class QwtPlotCanvas;
class Plot;

// This struct contains all loaded lab curves, its name, size and X/Y values
class TXT {
public:
    TXT () : initialPnt (0), endPnt (0) {
    }

    //Copy constructor
    TXT (const TXT &copy) {
        name = QString (copy.name);
//        time = QVector<double> (copy.time);
        sigmaAxialTotal = QVector<double> (copy.sigmaAxialTotal);
//        sigmaAxialDesv = QVector<double> (copy.sigmaAxialDesv);
        sigmaConf = QVector<double> (copy.sigmaConf);
        defAxial = QVector<double> (copy.defAxial);
        defLateral = QVector<double> (copy.defLateral);
//        defVol = QVector<double> (copy.defVol);
        I1 = QVector<double> (copy.I1);
        SQRTJ2 = QVector<double> (copy.SQRTJ2);
//        sigmaVol = QVector<double> (copy.sigmaVol);
//        sigmaLateral = QVector<double> (copy.sigmaLateral);
        initialPnt = copy.initialPnt;
        endPnt = copy.endPnt;
        X = QVector<double> (copy.X);
        Y = QVector<double> (copy.Y);
        testType = testTypes (copy.testType);
    }

//    //Operator =
//    TXT & operator=(const TXT &copy) {
//        if(this == &copy)
//            return *this;
//    }

    //Destructor
    ~TXT () {

    }

    void cutFile (int initialPnt_n, int endPnt_n) {
        initialPnt = initialPnt_n;
        endPnt = endPnt_n;
    }

    int getInitialPoint () {
        return initialPnt;
    }
    int getEndPoint () {
        return endPnt;
    }
    void resetFile () {
        initialPnt = 0;
        endPnt = X.size() -1;
    }
    int size () {
        return (endPnt - initialPnt) + 1;
    }

    QVector<double> getX () {
        return X.mid(initialPnt, endPnt-initialPnt +1);
    }
    QVector<double> getY () {
        return Y.mid(initialPnt, endPnt-initialPnt +1);
    }
    QVector<double> getSigmaAxialTotal () {
        return sigmaAxialTotal.mid(initialPnt, endPnt-initialPnt +1);
    }
    QVector<double> getSigmaConf () {
        return sigmaConf.mid(initialPnt, endPnt-initialPnt +1);
    }
    QVector<double> getDefAxial () {
        return defAxial.mid(initialPnt, endPnt-initialPnt +1);
    }
    QVector<double> getDefLateral () {
        return defLateral.mid(initialPnt, endPnt-initialPnt +1);
    }

    void setTestType (testTypes typeTest) {
        testType = typeTest;
        X = defAxial;
        Y = sigmaAxialTotal;
    }

    testTypes getTestType () {
        return testType;
    }

    QString name;
//    QVector<double> time;
    QVector<double> sigmaAxialTotal;
    QVector<double> sigmaAxialDesv;
    QVector<double> sigmaConf;
    QVector<double> defAxial;
    QVector<double> defLateral;
    QVector<double> defVol;
    QVector<double> I1;
    QVector<double> SQRTJ2;
    QVector<double> sigmaVol;
private:
    int initialPnt;
    int endPnt;
    QVector<double> X;
    QVector<double> Y;
    testTypes testType;
};

//X	Sig Axial Total	Tempo	Sig Axial Desv	Tempo	Sig Conf	Tempo	Def Axial	Tempo	Def Lateral	Tempo	Def Volume
//Tempo	Sig Axial Desv	Tempo	Def Axial	Tempo	Def Lateral	Tempo	Def Volume


namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
private slots:
    void on_listWidget_clicked(const QModelIndex &index);
    void on_actionShowGraphList_triggered();
    void on_actionZoom_toggled(bool arg1);
    void on_actionOpenFile_triggered();
    void clickedOnPlot(Plot *plotTmp);
    void fullscreenOnPlot(Plot *plotTmp);
    void ShowListContextMenu(const QPoint& pos);
    void updateCurve(int indexCurves, int indexStartPoint, int indexEndPoint);
    void ChangePlotAxis(Plot *plot_ptr, QString PlotAxis);
    void cutCurve(int indexCurves, int indexStartPoint, int indexEndPoint);

    void on_actionExportPlot_triggered();

    void on_actionShow_Parameter_List_triggered();

    void on_pushButton_clicked();

private:
    Ui::MainWindow *ui;
    QHash <int, TXT > FilesList;

    Plot *currentPlot;
    void reloadCurvesList(Plot *plotTmp);

};

#endif // MAINWINDOW_H
