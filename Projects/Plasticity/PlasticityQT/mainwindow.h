#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <qwt_plot_zoomer.h>
#include "common.h"

class QString;
class QModelIndex;
class QwtPicker;
class QwtEventPattern;
class QwtPlotCanvas;
class Plot;

// This struct contains all loaded lab curves, its name, size and X/Y values
struct TXT {
            QString name; QVector<double> *X; QVector<double> *Y;
            QVector<double> *time; QVector<double> *sigmaAxialTotal; QVector<double> *sigmaAxialDesv; QVector<double> *sigmaConf;
            QVector<double> *defAxial; QVector<double> *defLateral; QVector<double> *defVol;
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

    void on_actionExportPlot_triggered();

private:
    Ui::MainWindow *ui;
    QHash <int, TXT > *FilesList;

    Plot *currentPlot;
    void reloadCurvesList(Plot *plotTmp);

};

#endif // MAINWINDOW_H
