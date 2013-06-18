#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "initialpointdock.h"
#include "ui_initialpointdock.h"
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
#include "simulation.h"
#include "canvaspicker.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    setWindowTitle(tr("Plasticity"));

    // Creating objects and linking them to their pointers
    this->FilesList = new QHash <int, TXT > ;

    //////////////////////////////// Interface stuff ////////////////////////////////

    // Test 1
    ui->Plot_1->setTitle("TEST 01");
//    ui->Plot_1->setAxisTitle( QwtPlot::xBottom, "Deformation" );
//    ui->Plot_1->setAxisTitle( QwtPlot::yLeft, "Tension" );
    connect(ui->Plot_1->canvas_picker, SIGNAL(mouseLeftClicked(Plot *)),
            this, SLOT(clickedOnPlot(Plot *)));
    connect(ui->Plot_1->canvas_picker, SIGNAL(mouseDoubleClicked(Plot *)),
            this, SLOT(fullscreenOnPlot(Plot *)));
    // Initial Plot where graphs will be ploted
    currentPlot = ui->Plot_1;
    ui->Plot_1->setHighlighted(true);
    // Test 2
    ui->Plot_2->setTitle("TEST 02");
//    ui->Plot_2->setAxisTitle( QwtPlot::xBottom, "Deformation" );
//    ui->Plot_2->setAxisTitle( QwtPlot::yLeft, "Tension" );
    connect(ui->Plot_2->canvas_picker, SIGNAL(mouseLeftClicked(Plot *)),
            this, SLOT(clickedOnPlot(Plot *)));
    connect(ui->Plot_2->canvas_picker, SIGNAL(mouseDoubleClicked(Plot *)),
            this, SLOT(fullscreenOnPlot(Plot *)));
    // Test 3
    ui->Plot_3->setTitle("TEST 03");
//    ui->Plot_3->setAxisTitle( QwtPlot::xBottom, "Deformation" );
//    ui->Plot_3->setAxisTitle( QwtPlot::yLeft, "Tension" );
    connect(ui->Plot_3->canvas_picker, SIGNAL(mouseLeftClicked(Plot *)),
            this, SLOT(clickedOnPlot(Plot *)));
    connect(ui->Plot_3->canvas_picker, SIGNAL(mouseDoubleClicked(Plot *)),
            this, SLOT(fullscreenOnPlot(Plot *)));
    // Test 4
    ui->Plot_4->setTitle("TEST 04");
//    ui->Plot_4->setAxisTitle( QwtPlot::xBottom, "-" );
//    ui->Plot_4->setAxisTitle( QwtPlot::yLeft, "-" );
    connect(ui->Plot_4->canvas_picker, SIGNAL(mouseLeftClicked(Plot *)),
            this, SLOT(clickedOnPlot(Plot *)));
    connect(ui->Plot_4->canvas_picker, SIGNAL(mouseDoubleClicked(Plot *)),
            this, SLOT(fullscreenOnPlot(Plot *)));


    on_actionZoom_toggled(false);

    // TXT List Context menu
    ui->listWidget->setContextMenuPolicy(Qt::CustomContextMenu);
    connect(ui->listWidget, SIGNAL(customContextMenuRequested(const QPoint&)),
        this, SLOT(ShowListContextMenu(const QPoint&)));





    QSettings settings("myapp.ini", QSettings::IniFormat);
    settings.setValue("editor/wrapMargin", 68);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_listWidget_clicked(const QModelIndex &index)
{
    int pos = index.row();
    int checkStatus = ui->listWidget->item(pos)->checkState();

    if (checkStatus == 2) {
        currentPlot->createCurve(pos, this->FilesList->value(pos).X, this->FilesList->value(pos).Y, checkStatus);
    }
    else
    {
        currentPlot->deleteCurve(pos);
    }
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

void MainWindow::clickedOnPlot(Plot *plotTmp)
{
    reloadCurvesList (plotTmp);
    //set plotTmp to highlighted =)
    ui->Plot_1->setHighlighted(false);
    ui->Plot_2->setHighlighted(false);
    ui->Plot_3->setHighlighted(false);
    ui->Plot_4->setHighlighted(false);
    plotTmp->setHighlighted(true);
}

void MainWindow::reloadCurvesList(Plot *plotTmp)
{
    currentPlot = plotTmp;

    // repopulating listWidgetQwtPlot
    for (int i=0; i < FilesList->size(); i++) {
        QListWidgetItem *item;
        item = ui->listWidget->item(i);
        item->setCheckState(Qt::CheckState(plotTmp->CurvesList->value(i).chk_status));
    }
}

void MainWindow::on_actionShowGraphList_triggered()
{
    ui->dockWidget->show();
}

void MainWindow::updateCurve(int indexCurves, int indexStartPoint, int indexEndPoint) {
    qDebug() << "UpdateCurve MAIN";
    ui->Plot_1->updateCurve( indexCurves, indexStartPoint, indexEndPoint );
    ui->Plot_2->updateCurve( indexCurves, indexStartPoint, indexEndPoint );
    ui->Plot_3->updateCurve( indexCurves, indexStartPoint, indexEndPoint );
    ui->Plot_4->updateCurve( indexCurves, indexStartPoint, indexEndPoint );
}

void MainWindow::ShowListContextMenu(const QPoint& pos)
{
    QListWidgetItem *listItem = ui->listWidget->itemAt(pos);
    if (!listItem) return;

    QPoint globalPos = ui->listWidget->mapToGlobal(pos); // for most widgets
    // QPoint globalPos = myWidget->viewport()->mapToGlobal(pos); // for QAbstractScrollArea and derived classes you would use:

    QMenu myMenu;
    QAction *aSaveAs = myMenu.addAction("Save curve as...");
    QAction *aSelectPoint = myMenu.addAction("Select initial point");
    if (listItem->checkState() == 2)
        aSelectPoint->setEnabled(1);
    else
        aSelectPoint->setDisabled(1);

    QAction* selectedItem = myMenu.exec(globalPos);

    if (aSaveAs == selectedItem) {

        QFileDialog Save_File(this);

        QString fileName = Save_File.getSaveFileName(this, "Save file as...",
                                                     ".","Text files (*.txt);; All files (*)");

        if (fileName == "") {
            qDebug() << "Error: file cannot be created.";
            return;
        }

        QFile curve_file;
        curve_file.setFileName (fileName);

        if (!curve_file.open(QIODevice::ReadWrite))
        {
            qDebug() << "Error: No filename specified.";
            return;
        }

        QTextStream *out_curve = new QTextStream (&curve_file);

        int indexCurve = ui->listWidget->indexAt(pos).row();
        for (int i=0; i<this->FilesList->value(indexCurve).X->size(); i++) {
            *out_curve << this->FilesList->value(indexCurve).X->value(i) << " " << this->FilesList->value(indexCurve).Y->value(i) << endl;
        }



        //Atualizando entrada na tabela de arquivos / ListWidget
        TXT arquivo_1 = this->FilesList->take(indexCurve);
        QFileInfo infofile = QFileInfo (fileName);
        arquivo_1.name = infofile.baseName()+"."+infofile.suffix();

        this->FilesList->insert(indexCurve, arquivo_1);

        listItem->setText(arquivo_1.name);





        curve_file.close();
    }

    if (aSelectPoint == selectedItem) {
        int indexCurve = ui->listWidget->indexAt(pos).row();
        // VERIFICAR COMO NAO CRIAR 2 DIALOGOS PARA MESMA CURVA
        initialpointdock *selectpointdock = new initialpointdock(this);

        // connecting signals -> slots to show/hide symbols curve depending on dock visibility
        connect(selectpointdock, SIGNAL(showSymbCurve(int)),
                currentPlot, SLOT(showSymbCurve(int)));
        connect(selectpointdock, SIGNAL(hideSymbCurve(int)),
                currentPlot, SLOT(hideSymbCurve(int)));
        // connecting signal/slot to perform curve cut
        connect(selectpointdock, SIGNAL(cutCurve(int,int,int)),
                currentPlot, SLOT(cutCurve(int,int,int)));
        // connecting signal/slot to perform curve redraw with new data
        connect(selectpointdock, SIGNAL(cutCurve(int,int,int)),
                this, SLOT(updateCurve(int,int,int)));
        // connecting signals -> slots to sync information when points are changed
        connect(selectpointdock, SIGNAL(SymbPointChanged(int,int,Plot::pointType)),
                currentPlot->canvas_picker, SLOT(setSymbPoint(int,int,Plot::pointType)));
        connect(currentPlot->canvas_picker, SIGNAL(SymbPointChanged(int,int,Plot::pointType)),
                selectpointdock, SLOT(setSymbPoint(int,int,Plot::pointType)));

        selectpointdock->setIndexCurve(indexCurve);
        selectpointdock->setXdata(currentPlot->CurvesList->value(indexCurve).X);
        selectpointdock->setYdata(currentPlot->CurvesList->value(indexCurve).Y);
        int indexStartPoint = currentPlot->getSymbIndex(indexCurve, Plot::startPoint);
        int indexEndPoint   = currentPlot->getSymbIndex(indexCurve, Plot::endPoint);
        qDebug() << "aSelectPoint == selectedItem: [indexCurve] indexStartPoint / indexEndPoint" << "[" << indexCurve << "] " << indexStartPoint << " " << indexEndPoint ;
        selectpointdock->setSlidersMinimum(indexStartPoint);
        selectpointdock->setSlidersMaximum(indexEndPoint);
        selectpointdock->setSymbPoint(indexStartPoint, indexCurve, Plot::startPoint);
        selectpointdock->setSymbPoint(indexEndPoint, indexCurve, Plot::endPoint);
        selectpointdock->setWindowTitle(listItem->text());
        selectpointdock->setFloating(1);
        selectpointdock->show();

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

        if (!curve_file.open(QFile::ReadOnly | QFile::Text))
        {
            qDebug() << "Error: file not loaded.";
            return;
        }

//        QTextStream *in_curve = new QTextStream (&curve_file);
        double valueX, valueY;
        QVector<double> *Xvalues;
        QVector<double> *Yvalues;
        QVector<double> *time = new QVector<double>();
        QVector<double> *sigmaAxialTotal = new QVector<double>();
        QVector<double> *sigmaAxialDesv = new QVector<double>();
        QVector<double> *sigmaConf = new QVector<double>();
        QVector<double> *defAxial = new QVector<double>();
        QVector<double> *defLateral = new QVector<double>();
        QVector<double> *defVol = new QVector<double>();
        int countLines = 0;
        double Xsmallest=numeric_limits<double>::max(), Xbiggest=numeric_limits<double>::min(),
               Ysmallest=numeric_limits<double>::max(), Ybiggest=numeric_limits<double>::min();
        int typeTest = testTypes(All);

        QString line = curve_file.readLine(); // header line, ignoring it
        QStringList strings = line.split( QRegExp("[\\t]+") );

        if (strings.size() == 8) {
            qDebug() << "UCS";
            typeTest = testTypes(UCS); //0;
        }

        if (strings.size() == 12) {
            qDebug() << "Triaxial";
            typeTest = testTypes(Triaxial); //1;
        }

        while (!curve_file.atEnd()) {
            line = curve_file.readLine();
            line.chop(1); // to remove \n of the end

//            strings = line.split( QRegExp( "[ \\t]+") );
            strings = line.split( QRegExp("[\\t]+") );


            if (strings.size() < 2) continue;

            if (typeTest == testTypes(UCS)) {
//Tempo	SigAxialDesv	DefAxial	DefLateral	DefVolume
                time->insert( countLines, QVariant(strings.value(0)).toDouble() );
                sigmaAxialDesv->insert( countLines, QVariant(strings.value(1)).toDouble() );
                defAxial->insert( countLines, QVariant(strings.value(3)).toDouble() );
                defLateral->insert( countLines, QVariant(strings.value(5)).toDouble() );
                defVol->insert( countLines, QVariant(strings.value(7)).toDouble() );

                valueY = QVariant(strings.value(1)).toDouble(); //sigmaAxialDesv
                valueX = QVariant(strings.value(3)).toDouble(); //defAxial

                Xvalues = defAxial;
                Yvalues = sigmaAxialDesv;
            }
            else if (typeTest == testTypes(Triaxial))
            {
//Tempo	**SigAxialTotal**	SigAxialDesv	**SigConf**	DefAxial	DefLateral	DefVolume
                time->insert( countLines, QVariant(strings.value(0)).toDouble() );
                sigmaAxialTotal->insert( countLines, QVariant(strings.value(1)).toDouble() ); // Triaxial ONLY
                sigmaAxialDesv->insert( countLines, QVariant(strings.value(3)).toDouble() );
                sigmaConf->insert( countLines, QVariant(strings.value(5)).toDouble() );       // Triaxial ONLY
                defAxial->insert( countLines, QVariant(strings.value(7)).toDouble() );
                defLateral->insert( countLines, QVariant(strings.value(9)).toDouble() );
                defVol->insert( countLines, QVariant(strings.value(11)).toDouble() );

                valueY = QVariant(strings.value(1)).toDouble(); //sigmaAxialTotal
                valueX = QVariant(strings.value(7)).toDouble(); //defAxial

                Xvalues = defAxial;
                Yvalues = sigmaAxialDesv;
            }

            if (valueX < Xsmallest) Xsmallest = valueX;
            if (valueX > Xbiggest) Xbiggest = valueX;
            if (valueY < Ysmallest) Ysmallest = valueY;
            if (valueY > Ybiggest) Ybiggest = valueY;

//            qDebug() << "--->" << strings.size();
//            foreach (QString field, strings) {
//                qDebug() << QVariant(field).toDouble();
//            }
//            qDebug() << "<---";
            countLines ++ ;
        }

        qDebug() << "Count X / Y: " << Xvalues->size() << "/" << Yvalues->size() << " count: " << countLines;

        curve_file.close();

        //Criando entrada na tabela de arquivos
        TXT arquivoTXT;
        QFileInfo infofile = QFileInfo (fileName);
        arquivoTXT.name = infofile.baseName()+"."+infofile.suffix();
        arquivoTXT.X = Xvalues;
        arquivoTXT.Y = Yvalues;

        arquivoTXT.time = time;
        arquivoTXT.sigmaAxialDesv = sigmaAxialDesv;
        arquivoTXT.sigmaAxialTotal = sigmaAxialTotal;
        arquivoTXT.sigmaConf = sigmaConf;
        arquivoTXT.defAxial = defAxial;
        arquivoTXT.defLateral = defLateral;
        arquivoTXT.defVol = defVol;

        arquivoTXT.testType = testTypes(typeTest);

        int position = this->FilesList->size();
        this->FilesList->insert(position, arquivoTXT);
        //Criando entrada na listWidget
        QListWidgetItem *item;
        item = new QListWidgetItem();
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
        item->setCheckState(Qt::Unchecked);
        item->setText(this->FilesList->value(position).name);
        ui->listWidget->addItem(item);

    }
}

void MainWindow::on_actionExportPlot_triggered()
{
    currentPlot->exportPlot();
}
