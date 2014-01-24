#include "mainwindow.h"
#include "ui_mainwindow.h"
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
    // Test 3
    ui->Plot_3->setTitle("TEST 03");
    connect(ui->Plot_3->canvas_picker, SIGNAL(mouseLeftClicked(Plot *)),
            this, SLOT(clickedOnPlot(Plot *)));
    connect(ui->Plot_3->canvas_picker, SIGNAL(mouseDoubleClicked(Plot *)),
            this, SLOT(fullscreenOnPlot(Plot *)));
    // Test 4
    ui->Plot_4->setTitle("TEST 04");
    connect(ui->Plot_4->canvas_picker, SIGNAL(mouseLeftClicked(Plot *)),
            this, SLOT(clickedOnPlot(Plot *)));
    connect(ui->Plot_4->canvas_picker, SIGNAL(mouseDoubleClicked(Plot *)),
            this, SLOT(fullscreenOnPlot(Plot *)));

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
        int start_idx = newLabFile.Get_start_idx();
        int end_idx = newLabFile.Get_end_idx();
        int med_idx = DADOS.InsertLaboratoryData(newLabFile);
        int elastic_idx = newLabFile.Get_elastic_trans_idx();


        //Criando entrada na tabela de arquivos
        QFileInfo infofile = QFileInfo (fileName);
        QString fname = infofile.baseName()+"."+infofile.suffix();

        //Criando entrada na treeWidget
        QTreeWidgetItem *item;
        item = new QTreeWidgetItem();
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
        item->setCheckState(0, Qt::Unchecked);
        item->setText(0, fname);
        item->setData(0, Qt::UserRole, med_idx); //item "ID"
        ui->treeWidget->addTopLevelItem(item);

        //Criando entrada na comboBox de Medicoes
        ui->comboBoxMed->insertItem(ui->comboBoxMed->count() +1, fname, med_idx);

        //Setando valores minimo, maximo, inicial e final dos sliders e edits
        ui->start_idx_slider->setMaximum(end_idx);
        ui->start_idx_value->setMaxValue(end_idx);
        ui->end_idx_slider->setMaximum(end_idx);
        ui->end_idx_value->setMaxValue(end_idx);
        ui->elastic_trans_idx_slider->setMaximum(end_idx);
        ui->elastic_trans_idx_value->setMaxValue(end_idx);

        ui->start_idx_slider->setValue(start_idx);
        ui->start_idx_value->setValue(start_idx);
        ui->end_idx_slider->setValue(end_idx);
        ui->end_idx_value->setValue(end_idx);
        ui->elastic_trans_idx_slider->setValue(elastic_idx);
        ui->elastic_trans_idx_value->setValue(elastic_idx);
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
    QAction *aDeleteFile = myMenu.addAction("Unload file");

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

    if (aDeleteFile == selectedItem) {
//        int indexCurve = ui->treeWidget->indexAt(pos).row();
//        QTreeWidgetItem *tmp_item =  ui->treeWidget->topLevelItem( indexCurve );

//        indexCurve = -1;
//        indexCurve = tmp_item->data(0, Qt::UserRole).toInt();
//        qDebug() <<"INDEX CURVE!@%#$@!#$%$" <<indexCurve;
    }
}

void MainWindow::on_treeWidget_itemClicked(QTreeWidgetItem *item, int column)
{
    qDebug() << "CLICOU NO TREE WIDGET!!";
    int indexCurve = item->data(column, Qt::UserRole).toInt();
    qDebug() << "INDEX CURVE = "<<indexCurve;
    int checkStatus = item->checkState(column);
    qDebug() << "CHK STATUS = "<<checkStatus;

    if (checkStatus == 2) {
        qDebug() << "Criando curva";
        currentPlot->createCurve(indexCurve, checkStatus);
        //Recreate all envelopes in all plots, according to new curves
        if (ui->checkBox->isChecked())
            generateEnvelopeAllPlots();
    }
    else
    {
        qDebug() << "Apagando curva";
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
       int item_id = item->data(0,Qt::UserRole).toInt();
       item->setCheckState(0, Qt::CheckState(plotTmp->CurvesList.value(item_id).chk_status));
       for( int j = 0; j < item->childCount(); j++ )
       {
           QTreeWidgetItem *item_child = item->child(j);
           int item_child_id = item_child->data(0,Qt::UserRole).toInt();
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
        generateEnvelopeAllPlots();
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

void MainWindow::on_start_idx_slider_valueChanged(int value)
{
    QVariant testMed = ui->comboBoxMed->itemData(ui->comboBoxMed->currentIndex());
    if (testMed != QVariant::Invalid) {
        int indexCurve = testMed.toInt();
        DADOS.Set_Med_start_idx(indexCurve, value);

        ui->start_idx_value->setValue(value);
        setSymbAllPlots(value, indexCurve, Plot::startPoint);
    }
}

void MainWindow::on_end_idx_slider_valueChanged(int value)
{
    QVariant testMed = ui->comboBoxMed->itemData(ui->comboBoxMed->currentIndex());
    if (testMed != QVariant::Invalid) {
        int indexCurve = testMed.toInt();
        DADOS.Set_Med_end_idx(indexCurve, value);

        ui->end_idx_value->setValue(value);
        setSymbAllPlots(value, indexCurve, Plot::endPoint);
    }
}

void MainWindow::on_start_idx_value_valueChanged(double value)
{
    QVariant testMed = ui->comboBoxMed->itemData(ui->comboBoxMed->currentIndex());
    if (testMed != QVariant::Invalid) {
        int indexCurve = testMed.toInt();
        DADOS.Set_Med_start_idx(indexCurve, (int)value);

        ui->start_idx_slider->setValue((int)value);
        setSymbAllPlots((int)value, indexCurve, Plot::startPoint);
    }
}

void MainWindow::on_end_idx_value_valueChanged(double value)
{
    QVariant testMed = ui->comboBoxMed->itemData(ui->comboBoxMed->currentIndex());
    if (testMed != QVariant::Invalid) {
        int indexCurve = testMed.toInt();
        DADOS.Set_Med_end_idx(indexCurve, (int)value);

        ui->end_idx_slider->setValue((int)value);
        setSymbAllPlots((int)value, indexCurve, Plot::endPoint);
    }
}

void MainWindow::on_elastic_trans_idx_value_valueChanged(double value)
{
    QVariant testMed = ui->comboBoxMed->itemData(ui->comboBoxMed->currentIndex());
    if (testMed != QVariant::Invalid) {
        int indexCurve = testMed.toInt();
        DADOS.Set_Med_elastic_trans_idx(indexCurve, (int)value);

        ui->elastic_trans_idx_slider->setValue((int)value);
        setSymbAllPlots((int)value, indexCurve, Plot::elasticTransPoint);
    }
}

void MainWindow::on_elastic_trans_idx_slider_valueChanged(int value)
{
    QVariant testMed = ui->comboBoxMed->itemData(ui->comboBoxMed->currentIndex());
    if (testMed != QVariant::Invalid) {
        int indexCurve = testMed.toInt();
        DADOS.Set_Med_elastic_trans_idx(indexCurve, value);

        ui->elastic_trans_idx_value->setValue(value);
        setSymbAllPlots(value, indexCurve, Plot::elasticTransPoint);
    }
}

void MainWindow::setSymbAllPlots(int idx, int indexCurve, Plot::pointType ptnTyp) {
    ui->Plot_1->setSymbIndex(idx, indexCurve, ptnTyp);
    ui->Plot_2->setSymbIndex(idx, indexCurve, ptnTyp);
    ui->Plot_3->setSymbIndex(idx, indexCurve, ptnTyp);
    ui->Plot_4->setSymbIndex(idx, indexCurve, ptnTyp);
}

void MainWindow::on_runSimBtn_clicked(bool checked)
{
    try {
        ui->statusBar->clearMessage();

        QVariant testMed = ui->comboBoxMed->itemData(ui->comboBoxMed->currentIndex());
        if (testMed == QVariant::Invalid) {
            qDebug() << "PROBLEMA COM ID DA MEDICAO";
            return;
        }
        int indexCurve = testMed.toInt();

        qDebug() <<"aRunSimulation (" << indexCurve << ")";

        TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> sandlerObj;
        setParameters(sandlerObj);
//        DADOS.SetSandlerDimaggio(sandlerObj);

        TPBrStrainStressDataBase *basedata = DADOS.getObj(indexCurve);
        TPBrLaboratoryData *labdata = dynamic_cast<TPBrLaboratoryData *>(basedata);
        if(!labdata) DebugStop();

        DADOS.Set_Med_start_idx(indexCurve, ui->start_idx_value->value());
        DADOS.Set_Med_end_idx(indexCurve, ui->end_idx_value->value());
        DADOS.Set_Med_elastic_trans_idx(indexCurve, ui->elastic_trans_idx_value->value());

        qDebug() << "VAI SIMULAR: Start idx: " << ui->start_idx_value->value()
                 << " End idx: " << ui->end_idx_value->value()
                 << " Elastic idx: " << ui->elastic_trans_idx_value->value();

        int idx_sim = labdata->RunSimulation(sandlerObj);

        int idx_med = labdata->GlobalId();
        if(idx_med != indexCurve) DebugStop();

        qDebug() << "SIMULADO: Med idx: " << idx_med << " Sim idx: " << idx_sim;

        ///////////////////Criando entrada na treeWidget
        //Achando Medicao a qual essa simulacao pertence, para exibicao em forma de "tree"
        QTreeWidgetItem *item_parent;
        for(int i = 0; i < ui->treeWidget->topLevelItemCount(); i++)
        {
            QTreeWidgetItem *item_tmp = ui->treeWidget->topLevelItem(i);
            qDebug() << "ITEM_MASTER: " << item_tmp->text(0) << " id: " << item_tmp->data(0,Qt::UserRole).toInt();
            int item_id = item_tmp->data(0,Qt::UserRole).toInt();
            if (item_id == idx_med) {
                item_parent = item_tmp;
                qDebug() << "ITEM_PARENT: " << item_tmp->text(0) << " id: " << item_tmp->data(0,Qt::UserRole).toInt();
                break;
            }
        }
        QTreeWidgetItem *item;
        item = new QTreeWidgetItem(item_parent);
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
        item->setCheckState(0, Qt::Unchecked);
//        int count = DADOS.SizeLabData();
//        qDebug() << "COUNTM = " <<count << "!!!!!!!!!!!";
//        int counts = labdata->SizeSimData();
//        qDebug() << "COUNTS = " <<counts << "!!!!!!!!!!!";
        item->setText(0,ui->comboBoxSim->itemText(ui->comboBoxSim->currentIndex()));
        item->setData(0, Qt::UserRole, idx_sim); //item "ID"
        ui->treeWidget->addTopLevelItem(item);
        ui->treeWidget->expandAll();
        //update sim id in combobox sim
        ui->comboBoxSim->setItemData(ui->comboBoxSim->currentIndex(), idx_sim, Qt::UserRole);

        // lock parameters for this simulation
        //update lock / unlock simulation status button, stored in combo's Data role = Qt::UserRole+1
        ui->comboBoxSim->setItemData(ui->comboBoxSim->currentIndex(), true, Qt::UserRole+1);
        ui->lockParamsBtn->setChecked(true);
    }
    catch (...)
    {
        qDebug() << "SIMULACAO QUEBROU!!!";
        ui->statusBar->showMessage("SIMULACAO QUEBROU!!!");
    }
}

void MainWindow::on_comboBoxMed_currentIndexChanged(int index)
{
    int indexCurve = ui->comboBoxMed->itemData(index).toInt();
    int idx_max = DADOS.SizeMed(indexCurve) -1;
    int start_idx = DADOS.Get_Med_start_idx(indexCurve);
    int end_idx = DADOS.Get_Med_end_idx(indexCurve);

    //Setando valores minimo, maximo, inicial e final dos sliders e edits
    ui->start_idx_slider->setMaximum(idx_max);
    ui->end_idx_slider->setMaximum(idx_max);
    ui->start_idx_value->setMaxValue(idx_max);
    ui->end_idx_value->setMaxValue(idx_max);

    ui->start_idx_slider->setValue(start_idx);
    ui->start_idx_value->setValue(start_idx);
    ui->end_idx_slider->setValue(end_idx);
    ui->end_idx_value->setValue(end_idx);
}

void MainWindow::generateEnvelopeAllPlots() {
    REAL A, B, C;
    A = ui->A_counter->value();
    B = ui->B_counter->value();
    C = ui->C_counter->value();

    ui->Plot_1->Generate_Envelope(A, B, C);
    ui->Plot_2->Generate_Envelope(A, B, C);
    ui->Plot_3->Generate_Envelope(A, B, C);
    ui->Plot_4->Generate_Envelope(A, B, C);
}

void MainWindow::on_A_counter_valueChanged(double value)
{
    if (ui->checkBox->isChecked())
        generateEnvelopeAllPlots();
}

void MainWindow::on_B_counter_valueChanged(double value)
{
    if (ui->checkBox->isChecked())
        generateEnvelopeAllPlots();
}

void MainWindow::on_C_counter_valueChanged(double value)
{
    if (ui->checkBox->isChecked())
        generateEnvelopeAllPlots();
}

void MainWindow::on_comboBoxSim_currentIndexChanged(int index)
{
    //update lock / unlock simulation status button, stored in combo's Data role = Qt::UserRole+1
    ui->lockParamsBtn->setChecked(ui->comboBoxSim->itemData(ui->comboBoxSim->currentIndex(),Qt::UserRole+1).toBool());
}

void MainWindow::on_lockParamsBtn_toggled(bool checked)
{
    //lock / unlock simulation status, stored in Data role = Qt::UserRole+1
    bool locked = ui->comboBoxSim->itemData(ui->comboBoxSim->currentIndex(), Qt::UserRole+1).toBool();
    int id_sim = ui->comboBoxSim->itemData(ui->comboBoxSim->currentIndex(), Qt::UserRole).toInt();

    if (!checked and locked) {
        qDebug() << "Trying to unlock...";

        QMessageBox::StandardButton reply;
        reply = QMessageBox::question(this,
                                      "Warning",
                                      "Unlock it will invalidate (delete) all instances of this simulation, proceed?",
                                      QMessageBox::Yes|QMessageBox::No);
        if (reply == QMessageBox::Yes) {
          qDebug() << "Deletando as simulacoes plotadas de id = " << id_sim;
          //Apaga dos plots
          ui->Plot_1->deleteCurve(id_sim);
          ui->Plot_2->deleteCurve(id_sim);
          ui->Plot_3->deleteCurve(id_sim);
          ui->Plot_4->deleteCurve(id_sim);
          //Apaga da estrutura de controle de dados
          DADOS.DeleteGlobalId(id_sim);
          //Apaga do treewidget
          //Achando a simulacao na tree
          for(int i = 0; i < ui->treeWidget->topLevelItemCount(); i++)
          {
              QTreeWidgetItem *item_parent = ui->treeWidget->topLevelItem(i);
              qDebug() << "ITEMP: " << item_parent->text(0) << " id: " << item_parent->data(0,Qt::UserRole).toInt();
              for(int j = 0; j < item_parent->childCount(); j++) {
                  QTreeWidgetItem *item_tmp = item_parent->child(j);
                  int item_id = item_tmp->data(0,Qt::UserRole).toInt();
                  qDebug() << "ITEMC: " << item_tmp->text(0) << " id: " << item_id;
                  if (item_id == id_sim)
                    item_parent->removeChild(item_tmp);
              }
          }
        }
        else
        {
            //recheck the button and does not uncheck items
            ui->lockParamsBtn->setChecked(true);
            return;
        }

    }
    if (!locked and checked) {
        //uncheck the button and does not check items
        ui->lockParamsBtn->setChecked(false);
        return;
    }

    qDebug() << "checked: " << checked << " locked: " << locked;

    ui->A_counter->setDisabled(checked);
    ui->B_counter->setDisabled(checked);
    ui->C_counter->setDisabled(checked);
    ui->D_counter->setDisabled(checked);
    ui->R_counter->setDisabled(checked);
    ui->W_counter->setDisabled(checked);
    ui->poisson_counter->setDisabled(checked);
    ui->young_counter->setDisabled(checked);
    ui->start_idx_slider->setDisabled(checked);
    ui->start_idx_value->setDisabled(checked);
    ui->end_idx_slider->setDisabled(checked);
    ui->end_idx_value->setDisabled(checked);
    ui->elastic_trans_idx_slider->setDisabled(checked);
    ui->elastic_trans_idx_value->setDisabled(checked);
    ui->comboBoxMed->setDisabled(checked);
    ui->runSimBtn->setDisabled(checked);
    ui->identifyBtn->setDisabled(checked);

    //lock / unlock simulation status, stored in Data role = Qt::UserRole+1
    ui->comboBoxSim->setItemData(ui->comboBoxSim->currentIndex(), checked, Qt::UserRole+1);
}



void MainWindow::on_identifyBtn_clicked(bool checked)
{
    QVariant testMed = ui->comboBoxMed->itemData(ui->comboBoxMed->currentIndex());

    if (testMed != QVariant::Invalid) {

        int indexCurve = testMed.toInt();

        TPBrStrainStressDataBase *basedata = DADOS.getObj(indexCurve);
        TPBrLaboratoryData *labdata = dynamic_cast<TPBrLaboratoryData *>(basedata);

        REAL newYoung = ui->young_counter->value();
        REAL newPoisson = ui->poisson_counter->value();

        labdata->IdentifyElasticity(newYoung, newPoisson);

        ui->young_counter->setValue(newYoung);
        ui->poisson_counter->setValue(newPoisson);

    }
}
