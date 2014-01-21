#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <qwt_plot_zoomer.h>
#include "common.h"
#include <limits>
#include <qtreewidget.h>
#include <QTreeWidgetItem>

#include "TPBrDataControl.h"
#include "plot.h"


class QString;
class QModelIndex;
class QwtPicker;
class QwtEventPattern;
class QwtPlotCanvas;
//class Plot;

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
    void ShowListContextMenu(const QPoint& pos);
    void on_actionOpenFile_triggered();
    void clickedOnPlot(Plot *plotTmp);
    void fullscreenOnPlot(Plot *plotTmp);
    void on_treeWidget_itemClicked(QTreeWidgetItem *item, int column);
    void updateParameters(TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> obj);
    void setParameters(TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &obj);
    void setSymbAllPlots(int idx, int indexCurve, Plot::pointType ptnTyp);
    void generateEnvelopeAllPlots();
    void on_checkBox_toggled(bool checked);
    void on_actionZoom_toggled(bool on);

    void on_start_idx_slider_valueChanged(int value);

    void on_end_idx_slider_valueChanged(int value);

    void on_start_idx_value_valueChanged(double value);

    void on_end_idx_value_valueChanged(double value);

    void on_runSimBtn_clicked(bool checked);

    void on_comboBoxMed_currentIndexChanged(int index);

    void on_A_counter_valueChanged(double value);

    void on_B_counter_valueChanged(double value);

    void on_C_counter_valueChanged(double value);

    void on_comboBoxSim_currentIndexChanged(int index);

    void on_lockParamsBtn_toggled(bool checked);

    void on_elastic_trans_idx_value_valueChanged(double value);

    void on_elastic_trans_idx_slider_valueChanged(int value);

    void on_identifyBtn_clicked(bool checked);


private:
    Ui::MainWindow *ui;
    Plot *currentPlot;
    void reloadCurvesList(Plot *plotTmp);
};

#endif // MAINWINDOW_H
