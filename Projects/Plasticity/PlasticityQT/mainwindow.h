#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <qwt_plot_zoomer.h>
#include "common.h"
#include <limits>

#include "../PlasticityQT_CLI/TPBrDataControl.h"


class QString;
class QModelIndex;
class QwtPicker;
class QwtEventPattern;
class QwtPlotCanvas;
class Plot;

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

    void on_pushButton_clicked();

private:
    Ui::MainWindow *ui;

    Plot *currentPlot;

};

#endif // MAINWINDOW_H
