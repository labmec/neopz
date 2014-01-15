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

    on_actionZoom_toggled(false);

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

        int indexCurve = ui->treeWidget->indexAt(pos).row();
        QTreeWidgetItem *tmp_item =  ui->treeWidget->topLevelItem( indexCurve );
        indexCurve = tmp_item->data(0, 5).toInt();

        // VERIFICAR COMO NAO CRIAR 2 DIALOGOS PARA MESMA CURVA
        initialpointdock *selectpointdock = new initialpointdock(this);

        // connecting signals -> slots to show/hide symbols curve depending on dock visibility
        connect(selectpointdock, SIGNAL(showSymbCurve(int)),
                currentPlot, SLOT(showSymbCurve(int)));
        connect(selectpointdock, SIGNAL(hideSymbCurve(int)),
                currentPlot, SLOT(hideSymbCurve(int)));
        // connecting signal/slot to perform curve redraw with new data
        connect(selectpointdock, SIGNAL(cutCurve(int,int,int)),
                this, SLOT(updateCurve(int,int,int)));
        // connecting signals -> slots to sync information when points are changed
        connect(selectpointdock, SIGNAL(SymbPointChanged(int,int,Plot::pointType)),
                currentPlot->canvas_picker, SLOT(setSymbPoint(int,int,Plot::pointType)));
        connect(currentPlot->canvas_picker, SIGNAL(SymbPointChanged(int,int,Plot::pointType)),
                selectpointdock, SLOT(setSymbPoint(int,int,Plot::pointType)));

        selectpointdock->setGlobal_ID(indexCurve);

        selectpointdock->setWindowTitle(tmp_item->text(0));
        selectpointdock->setFloating(1);
        selectpointdock->setGeometry(300,250,selectpointdock->width(),selectpointdock->height());
        selectpointdock->show();
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
//        REAL poisson, E, A, B, C, R, D, W;
//        E = 29269;
//        poisson = 0.203;
//        A = 616.67;
//        B = 0.0036895;
//        C = 111.48;
//        D = 0.018768;
//        R = 0.91969;
//        W = 0.006605;
        setParameters(sandlerObj);
//        sandlerObj.SetUp(poisson, E, A, B, C, R, D, W);
        DADOS.SetSandlerDimaggio(sandlerObj);

        TPBrStrainStressDataBase *basedata = DADOS.getObj(indexCurve);
        TPBrLaboratoryData *labdata = dynamic_cast<TPBrLaboratoryData *>(basedata);
        if(!labdata) DebugStop();
        labdata->Set_start_idx(100);
        int idx_sim = labdata->RunSimulation(sandlerObj);
        int idx_med = labdata->GlobalId();
        if(idx_med != indexCurve) DebugStop();

        qDebug() << "Med idx: " << idx_med << " Sim idx: " << idx_sim;

        //Criando entrada na treeWidget
        QTreeWidgetItem *item;
        item = new QTreeWidgetItem(tmp_item);
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
        item->setCheckState(0, Qt::Unchecked);
        int count = DADOS.SizeLabData();
        qDebug() << "COUNTM = " <<count << "!!!!!!!!!!!";
        int counts = labdata->SizeSimData();
        qDebug() << "COUNTS = " <<counts << "!!!!!!!!!!!";
//        item->setText(0, QString ("Sim").append( QVariant(count).toString()) );
        item->setText(0, "Sim");
        item->setData(0, 5, idx_sim); //item "ID"
        ui->treeWidget->addTopLevelItem(item);
        ui->treeWidget->expandAll();
    }
}

void MainWindow::on_treeWidget_itemClicked(QTreeWidgetItem *item, int column)
{
    qDebug() << "CLICOU NO TREE WIDGET!!";
    int indexCurve = item->data(column, 5).toInt();
    qDebug() << "INDEX CURVE = "<<indexCurve;
    int checkStatus = item->checkState(column);
    qDebug() << "CHK STATUS = "<<checkStatus;

    if (DADOS.isMed(indexCurve) == 0) // eh uma simulacao
    {
//         TPBrSimulationData *simData = dynamic_cast<TPBrSimulationData*> (DADOS.getObj(indexCurve));
//         if(!simData)
//         {
//             DebugStop();
//         }
//         REAL poisson, E, A, B, C, R, D, W;
//         simData->fSandler.getParams(poisson, E, A, B, C, R, D, W);
//         qDebug () << "AAAAAAAHHHH: " << poisson << " " << E<< " " << A<< " " << B<< " " << C<< " " << R<< " " << D<< " " << W;
//         updateParameters(simData->getSandler());
    }
    //para medicao fazer oq?????????????????????

    if (checkStatus == 2) {
        currentPlot->createCurve(indexCurve, checkStatus);
    }
    else
    {
        currentPlot->deleteCurve(indexCurve);
    }
}

void MainWindow::clickedOnPlot(Plot *plotTmp)
{
    reloadCurvesList (plotTmp);

    //set plotTmp to highlighted
    ui->Plot_1->setHighlighted(false);
    ui->Plot_2->setHighlighted(false);
    ui->Plot_3->setHighlighted(false);
    ui->Plot_4->setHighlighted(false);
    plotTmp->setHighlighted(true);

    //set plotTmp cursor to cross
    ui->Plot_1->canvas()->setCursor(Qt::PointingHandCursor);
    ui->Plot_2->canvas()->setCursor(Qt::PointingHandCursor);
    ui->Plot_3->canvas()->setCursor(Qt::PointingHandCursor);
    ui->Plot_4->canvas()->setCursor(Qt::PointingHandCursor);
    plotTmp->canvas()->setCursor(Qt::CrossCursor);
}

void MainWindow::fullscreenOnPlot(Plot *plotTmp)
{
    if (ui->Plot_1->isHidden() || ui->Plot_2->isHidden() || ui->Plot_3->isHidden() || ui->Plot_4->isHidden()) {
        ui->Plot_1->show();
        ui->Plot_2->show();
        ui->Plot_3->show();
        ui->Plot_4->show();
    }
    else {
        ui->Plot_1->hide();
        ui->Plot_2->hide();
        ui->Plot_3->hide();
        ui->Plot_4->hide();

        plotTmp->show();
    }
}

void MainWindow::reloadCurvesList(Plot *plotTmp)
{
    currentPlot = plotTmp;

    // updating check-status
    for(int i = 0; i < ui->treeWidget->topLevelItemCount(); i++)
    {
       QTreeWidgetItem *item = ui->treeWidget->topLevelItem(i);
       int item_id = item->data(0,5).toInt();
       item->setCheckState(0, Qt::CheckState(plotTmp->CurvesList.value(item_id).chk_status));
       for( int j = 0; j < item->childCount(); j++ )
       {
           QTreeWidgetItem *item_child = item->child(j);
           int item_child_id = item_child->data(0,5).toInt();
           item_child->setCheckState(0, Qt::CheckState(plotTmp->CurvesList.value(item_child_id).chk_status));
        }
    }
}

void MainWindow::updateParameters(TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> obj)
{
//     REAL poisson, E, A, B, C, R, D, W;
//     obj.getParams(poisson, E, A, B, C, R, D, W);
//     ui->A_counter->setValue(A);
//     ui->B_counter->setValue(B);
//     ui->C_counter->setValue(C);
//     ui->D_counter->setValue(D);
//     ui->R_counter->setValue(R);
//     ui->W_counter->setValue(W);
//     ui->poisson_counter->setValue(poisson);
//     ui->young_counter->setValue(E);
}

void MainWindow::setParameters(TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &obj)
{
    REAL poisson, E, A, B, C, R, D, W;
    A = ui->A_counter->value();
    B = ui->B_counter->value();
    C = ui->C_counter->value();
    D = ui->D_counter->value();
    R = ui->R_counter->value();
    W = ui->W_counter->value();
    poisson = ui->poisson_counter->value();
    E = ui->young_counter->value();

    obj.SetUp(poisson, E, A, B, C, R, D, W);
}

void MainWindow::on_checkBox_toggled(bool checked)
{

    if (checked)
    {
        REAL A, B, C;
        A = ui->A_counter->value();
        B = ui->B_counter->value();
        C = ui->C_counter->value();

        //qDebug() <<"------------------------->"<< A << B << C;
	
	// a curva envelope deveria ser gerada cada vez que modifica os parametros A, B, C

        ui->Plot_1->Generate_Envelope(A, B, C);
        ui->Plot_2->Generate_Envelope(A, B, C);
        ui->Plot_3->Generate_Envelope(A, B, C);
        ui->Plot_4->Generate_Envelope(A, B, C);
    }

    else
    {
        ui->Plot_1->hide_envelope();
        ui->Plot_2->hide_envelope();
        ui->Plot_3->hide_envelope();
        ui->Plot_4->hide_envelope();
    }
}

void MainWindow::on_actionZoom_toggled(bool on)
{
    ui->Plot_1->panner->setEnabled(on);
    ui->Plot_1->zoomer->setEnabled(on);

    ui->Plot_2->panner->setEnabled(on);
    ui->Plot_2->zoomer->setEnabled(on);

    ui->Plot_3->panner->setEnabled(on);
    ui->Plot_3->zoomer->setEnabled(on);

    ui->Plot_4->panner->setEnabled(on);
    ui->Plot_4->zoomer->setEnabled(on);
}
