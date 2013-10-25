#include "initialpointdock.h"
#include "ui_initialpointdock.h"

#include <QMessageBox>

initialpointdock::initialpointdock(QWidget *parent) :
    QDockWidget(parent),
    ui(new Ui::initialpointdock)
{
    ui->setupUi(this);

//    coordStartIndex = -1;
//    coordEndIndex = -1;
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

void initialpointdock::setIndexCurve (int idx) {
    this->curveIndex = idx;
    emit showSymbCurve (this->curveIndex);
}

void initialpointdock::setSlidersMaximum (int value) {
    this->ui->sliderStart->setMaximum(value);
    this->ui->sliderEnd->setMaximum(value);
}

void initialpointdock::setSlidersMinimum (int value) {
    this->ui->sliderStart->setMinimum(value);
    this->ui->sliderEnd->setMinimum(value);
}

void initialpointdock::setXdata (const QVector <double> &Xvec) {
    this->X = Xvec;
}

void initialpointdock::setYdata (const QVector <double> &Yvec){
    this->Y = Yvec;
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
void initialpointdock::setSymbPoint (int idx, int curveIndex, Plot::pointType typept) {
    if (this->curveIndex == curveIndex) {
        if (typept == Plot::startPoint) {
            this->coordStartIndex = idx;
            this->ui->sliderStart->setValue(idx);
            this->ui->labelStartX->setText( QVariant (this->X.value(idx)).toString() );
            this->ui->labelStartY->setText( QVariant (this->Y.value(idx)).toString() );
        }
        else
        {
            this->coordEndIndex = idx;
            this->ui->sliderEnd->setValue(idx);
            this->ui->labelEndX->setText( QVariant (this->X.value(idx)).toString() );
            this->ui->labelEndY->setText( QVariant (this->Y.value(idx)).toString() );
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
    if (!isCut) {
        emit SymbPointChanged (0, this->curveIndex, Plot::startPoint);
        emit SymbPointChanged (this->X.size()-1, this->curveIndex, Plot::endPoint);
    }
    emit hideSymbCurve (this->curveIndex);
}
