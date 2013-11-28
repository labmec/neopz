#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "initialpointdock.h"
#include "ui_initialpointdock.h"
#include "globalconfig.h"
#include <iostream>
#include <limits>
using namespace std;

#include <QMainWindow>
#include <QHash>
#include <QString>
#include <QListWidgetItem>
#include <QtGui>
#include <QSignalMapper>
#include <QMouseEvent>
#include <qwt_plot_grid.h>
#include <qwt_plot_curve.h>
#include <qwt_scale_engine.h>
#include <qwt_symbol.h>
#include <qwt_picker_machine.h>
#include <qwt_plot_panner.h>
#include <qwt_plot_renderer.h>
#include <qwt_plot_marker.h>
#include <qwt_counter.h>
#include <qwt_plot.h>
#include "plot.h"
#include "canvaspicker.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    setWindowTitle(tr("Plasticity"));

    //////////////////////////////// Interface stuff ////////////////////////////////

    // Test 1
    ui->Plot_1->setTitle("TEST 01");
    connect(ui->Plot_1->canvas_picker, SIGNAL(mouseLeftClicked(Plot *)),
            this, SLOT(clickedOnPlot(Plot *)));
    connect(ui->Plot_1->canvas_picker, SIGNAL(mouseDoubleClicked(Plot *)),
            this, SLOT(fullscreenOnPlot(Plot *)));
    connect(ui->Plot_1, SIGNAL(AxisChanged_signal (Plot*, QString)),
            this, SLOT(ChangePlotAxis (Plot*,QString)));
    // Initial Plot where graphs will be ploted
    currentPlot = ui->Plot_1;
    ui->Plot_1->setHighlighted(true);
    ui->Plot_1->canvas()->setCursor(Qt::CrossCursor);
    // Test 2
    ui->Plot_2->setTitle("TEST 02");
    connect(ui->Plot_2->canvas_picker, SIGNAL(mouseLeftClicked(Plot *)),
            this, SLOT(clickedOnPlot(Plot *)));
    connect(ui->Plot_2->canvas_picker, SIGNAL(mouseDoubleClicked(Plot *)),
            this, SLOT(fullscreenOnPlot(Plot *)));
    connect(ui->Plot_2, SIGNAL(AxisChanged_signal (Plot*, QString)),
            this, SLOT(ChangePlotAxis(Plot*,QString)));
    // Test 3
    ui->Plot_3->setTitle("TEST 03");
    connect(ui->Plot_3->canvas_picker, SIGNAL(mouseLeftClicked(Plot *)),
            this, SLOT(clickedOnPlot(Plot *)));
    connect(ui->Plot_3->canvas_picker, SIGNAL(mouseDoubleClicked(Plot *)),
            this, SLOT(fullscreenOnPlot(Plot *)));
    connect(ui->Plot_3, SIGNAL(AxisChanged_signal (Plot*, QString)),
            this, SLOT(ChangePlotAxis (Plot*,QString)));
    // Test 4
    ui->Plot_4->setTitle("TEST 04");
    connect(ui->Plot_4->canvas_picker, SIGNAL(mouseLeftClicked(Plot *)),
            this, SLOT(clickedOnPlot(Plot *)));
    connect(ui->Plot_4->canvas_picker, SIGNAL(mouseDoubleClicked(Plot *)),
            this, SLOT(fullscreenOnPlot(Plot *)));
    connect(ui->Plot_4, SIGNAL(AxisChanged_signal (Plot*, QString)),
            this, SLOT(ChangePlotAxis (Plot*,QString)));

//    on_actionZoom_toggled(false);

    // TXT List Context menu
    ui->treeWidget->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(ui->treeWidget, SIGNAL(customContextMenuRequested(const QPoint&)),
        this, SLOT(ShowListContextMenu(const QPoint&)));

    QSettings settings("myapp.ini", QSettings::IniFormat);
    settings.setValue("editor/wrapMargin", 68);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_actionOpenFile_triggered()
{
    QFileDialog Open_File(this);
    QStringList fileNames = Open_File.getOpenFileNames(this, "Select one or more files to open",
                                                 ".","Text files (*.txt);; All files (*)");

    foreach (QString fileName, fileNames) {

        if (fileName == "") {
            qDebug() << "Error: No filename specified.";
            return;
        }

        QFile curve_file;
        curve_file.setFileName (fileName);


        TPBrLaboratoryData newLabFile (fileName.toStdString());
        newLabFile.Set_start_idx(0);
        int med_idx = DADOS.InsertLaboratoryData(newLabFile);


        //Criando entrada na tabela de arquivos
        QFileInfo infofile = QFileInfo (fileName);
        QString fname = infofile.baseName()+"."+infofile.suffix();

        //Criando entrada na treeWidget
        QTreeWidgetItem *item;
        item = new QTreeWidgetItem();
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
        item->setCheckState(0, Qt::Unchecked);
        item->setText(0, fname);
        item->setData(0, 5, med_idx); //item "ID"
        ui->treeWidget->addTopLevelItem(item);
    }
}


void MainWindow::ShowListContextMenu(const QPoint& pos)
{
    QTreeWidgetItem *treeItem = ui->treeWidget->itemAt(pos);
    if (!treeItem) return;

    QPoint globalPos = ui->treeWidget->mapToGlobal(pos); // for most widgets
    // QPoint globalPos = myWidget->viewport()->mapToGlobal(pos); // for QAbstractScrollArea and derived classes you would use:

    QMenu myMenu;
    QAction *aSaveAs = myMenu.addAction("Save curve as...");
    QAction *aSelectPoint = myMenu.addAction("Select initial point");
    QAction *aResetPoint = myMenu.addAction("Reset initial point");
    QAction *aDeleteFile = myMenu.addAction("Unload file");
    QAction *aRunSimulation = myMenu.addAction("Run Simulation");

    if (treeItem->checkState(0) == 2) {
        aSelectPoint->setEnabled(1);
        aResetPoint->setEnabled(1);
    }
    else {
        aSelectPoint->setDisabled(1);
        aResetPoint->setDisabled(1);
    }

    QAction* selectedItem = myMenu.exec(globalPos);

//    if (aSaveAs == selectedItem) {

//        QFileDialog Save_File(this);

//        QString fileName = Save_File.getSaveFileName(this, "Save file as...",
//                                                     ".","Text files (*.txt);; All files (*)");

//        if (fileName == "") {
//            qDebug() << "Error: file cannot be created.";
//            return;
//        }

//        QFile curve_file;
//        curve_file.setFileName (fileName);

//        if (!curve_file.open(QIODevice::ReadWrite))
//        {
//            qDebug() << "Error: No filename specified.";
//            return;
//        }

//        QTextStream *out_curve = new QTextStream (&curve_file);

//        int indexCurve = ui->listWidget->item(ui->listWidget->indexAt(pos).row())->data(5).toInt();
//        for (int i=0; i<this->FilesList.value(indexCurve).X->size(); i++) {
//            *out_curve << this->FilesList.value(indexCurve).X->value(i) << " " << this->FilesList.value(indexCurve).Y->value(i) << endl;
//        }

//        //Atualizando entrada na tabela de arquivos / ListWidget
//        TXT arquivo_1 = this->FilesList.take(indexCurve);
//        QFileInfo infofile = QFileInfo (fileName);
//        arquivo_1.name = infofile.baseName()+"."+infofile.suffix();

//        this->FilesList.insert(indexCurve, arquivo_1);

//        listItem->setText(arquivo_1.name);

//        curve_file.close();
//    }

    if (aResetPoint == selectedItem) {

//        int indexCurve = ui->listWidget->item(ui->listWidget->indexAt(pos).row())->data(5).toInt();
//        //resetFile sets 'initial point' = first point of the file (0)
//        //and 'end point' = last point of the file (size - 1)
//        FilesList[indexCurve].resetFile();
//        int initialPnt = FilesList[indexCurve].getInitialPoint();
//        int endPnt = FilesList[indexCurve].getEndPoint();
//        updateCurve(indexCurve, initialPnt, endPnt);
    }

    if (aSelectPoint == selectedItem) {

//        int indexCurve = ui->listWidget->item(ui->listWidget->indexAt(pos).row())->data(5).toInt();

//        // VERIFICAR COMO NAO CRIAR 2 DIALOGOS PARA MESMA CURVA
//        initialpointdock *selectpointdock = new initialpointdock(this);

//        // connecting signals -> slots to show/hide symbols curve depending on dock visibility
//        connect(selectpointdock, SIGNAL(showSymbCurve(int)),
//                currentPlot, SLOT(showSymbCurve(int)));
//        connect(selectpointdock, SIGNAL(hideSymbCurve(int)),
//                currentPlot, SLOT(hideSymbCurve(int)));
//        // connecting signal/slot to perform curve redraw with new data
//        connect(selectpointdock, SIGNAL(cutCurve(int,int,int)),
//                this, SLOT(updateCurve(int,int,int)));
//        // connecting signals -> slots to sync information when points are changed
//        connect(selectpointdock, SIGNAL(SymbPointChanged(int,int,Plot::pointType)),
//                currentPlot->canvas_picker, SLOT(setSymbPoint(int,int,Plot::pointType)));
//        connect(currentPlot->canvas_picker, SIGNAL(SymbPointChanged(int,int,Plot::pointType)),
//                selectpointdock, SLOT(setSymbPoint(int,int,Plot::pointType)));

//        selectpointdock->setIndexCurve(indexCurve);
//        selectpointdock->setXdata(currentPlot->CurvesList[indexCurve].X);
//        selectpointdock->setYdata(currentPlot->CurvesList[indexCurve].Y);
//        int indexStartPoint = currentPlot->getSymbIndex(indexCurve, Plot::startPoint);
//        int indexEndPoint   = currentPlot->getSymbIndex(indexCurve, Plot::endPoint);
//        qDebug() << "aSelectPoint == selectedItem: [indexCurve] indexStartPoint / indexEndPoint" << "[" << indexCurve << "] " << indexStartPoint << " " << indexEndPoint ;
//        selectpointdock->setSlidersMinimum(indexStartPoint);
//        selectpointdock->setSlidersMaximum(indexEndPoint);
//        selectpointdock->setSymbPoint(indexStartPoint, indexCurve, Plot::startPoint);
//        selectpointdock->setSymbPoint(indexEndPoint, indexCurve, Plot::endPoint);
//        selectpointdock->setWindowTitle(listItem->text());
//        selectpointdock->setFloating(1);
//        selectpointdock->setGeometry(300,250,selectpointdock->width(),selectpointdock->height());
//        selectpointdock->show();
    }

    if (aDeleteFile == selectedItem) {

        int indexCurve = ui->treeWidget->indexAt(pos).row();
        QTreeWidgetItem *tmp_item =  ui->treeWidget->topLevelItem( indexCurve );

        indexCurve = -1;
        indexCurve = tmp_item->data(0, 5).toInt();
        qDebug() <<"INDEX CURVE!@%#$@!#$%$" <<indexCurve;
    }

    if (aRunSimulation == selectedItem) {

        int indexCurve = ui->treeWidget->indexAt(pos).row();
        QTreeWidgetItem *tmp_item =  ui->treeWidget->topLevelItem( indexCurve );

        indexCurve = -1;
        indexCurve = tmp_item->data(0, 5).toInt();
        qDebug() <<"aRunSimulation (" << indexCurve << ")";

        TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> sandlerObj;
        REAL poisson, E, A, B, C, R, D, W;
        E = 29269;
        poisson = 0.203;
        A = 616.67;
        B = 0.0036895;
        C = 111.48;
        D = 0.018768;
        R = 0.91969;
        W = 0.006605;
        sandlerObj.SetUp(poisson, E, A, B, C, R, D, W);
        DADOS.SetSandlerDimaggio(sandlerObj);

        DADOS.fMedicoes[indexCurve].Set_start_idx(100);
        int idx_sim = DADOS.fMedicoes[indexCurve].RunSimulation(sandlerObj);


        int idx_med = DADOS.fMedicoes[indexCurve].GlobalId();
        qDebug() << "Med idx: " << idx_med << " Sim idx: " << idx_sim;

        //Criando entrada na treeWidget
        QTreeWidgetItem *item;
        item = new QTreeWidgetItem();
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
        item->setCheckState(0, Qt::Unchecked);
        item->setText(0, "fname");
        item->setData(0, 5, idx_med); //item "ID"
        ui->treeWidget->addTopLevelItem(item);
    }
}


void MainWindow::on_pushButton_clicked()
{

}
