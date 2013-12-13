#include "initialpointdock.h"
#include "ui_initialpointdock.h"

#include <QMessageBox>

initialpointdock::initialpointdock(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::initialpointdock)
{
    ui->setupUi(this);

    curveIndex = -1;
    isCut = false;

    connect (ui->sliderStart, SIGNAL(valueChanged(int)),
             this, SLOT(sliderStartValueChanged(int)));
    connect (ui->sliderEnd, SIGNAL(valueChanged(int)),
             this, SLOT(sliderEndValueChanged(int)));
}

initialpointdock::~initialpointdock()
{
    delete ui;
}

void initialpointdock::setGlobal_ID (int global_id) {
    this->curveIndex = global_id;

    TPBrLaboratoryData *symb_obj = dynamic_cast <TPBrLaboratoryData*> (DADOS.getObj(global_id));
    this->ui->sliderStart->setMinimum(0);
    this->ui->sliderStart->setMaximum(symb_obj->fSig_Ax.size()-1);
    this->ui->sliderEnd->setMinimum(0);
    this->ui->sliderEnd->setMaximum(symb_obj->fSig_Ax.size()-1);

    setSymbPoint(symb_obj->Get_start_idx(), global_id, Plot::startPoint);
    setSymbPoint(symb_obj->Get_end_idx(), global_id, Plot::endPoint);

    emit showSymbCurve (this->curveIndex);
}

//slot
void initialpointdock::sliderStartValueChanged ( int value ) {
    this->coordStartIndex = value;
    emit SymbPointChanged (this->coordStartIndex, this->curveIndex, Plot::startPoint);
}

//slot
void initialpointdock::sliderEndValueChanged ( int value ) {
    this->coordEndIndex = value;
    emit SymbPointChanged (this->coordEndIndex, this->curveIndex, Plot::endPoint);
}

//slot
void initialpointdock::setSymbPoint (int idx, int global_id, Plot::pointType typept) {
    if (this->curveIndex == curveIndex) {

//        TPBrLaboratoryData *symb_obj = dynamic_cast <TPBrLaboratoryData*> (DADOS.getObj(global_id));
        if (typept == Plot::startPoint) {
            this->coordStartIndex = idx;
            this->ui->sliderStart->setValue(idx);
//            this->ui->labelStartX->setText( QVariant (this->X.value(idx)).toString() );
//            this->ui->labelStartY->setText( QVariant (this->Y.value(idx)).toString() );
        }
        else
        {
            this->coordEndIndex = idx;
            this->ui->sliderEnd->setValue(idx);
//            this->ui->labelEndX->setText( QVariant (this->X.value(idx)).toString() );
//            this->ui->labelEndY->setText( QVariant (this->Y.value(idx)).toString() );
        }
    }
}


void initialpointdock::on_cutBtn_clicked()
{
    if (this->coordStartIndex >= this->coordEndIndex) {
        QMessageBox::warning(this, tr("Warning"),
                            tr("End coordinate cannot be before Start coordinate."),
                            QMessageBox::Ok);
        return;
    }
    else
    {
        emit hideSymbCurve (this->curveIndex);
        emit cutCurve(this->curveIndex, this->coordStartIndex, this->coordEndIndex);
        isCut = true;
        this->close();
    }
}

void initialpointdock::closeEvent(QCloseEvent *event)
{
//    if (!isCut) {
//        emit SymbPointChanged (0, this->curveIndex, Plot::startPoint);
//        emit SymbPointChanged (this->X.size()-1, this->curveIndex, Plot::endPoint);
//    }
    emit hideSymbCurve (this->curveIndex);
}
