#ifndef PLOT_H
#define PLOT_H

#include <QMenu>
#include <QDebug>
#include "common.h"
#include <qwt_plot.h>
#include "qwt_plot_curve.h"
//class QwtPlotZoomer;
#include <qwt_plot_zoomer.h>
#include <QVector>
class QwtPlotPanner;
class QwtPlotGrid;
class CanvasPicker;
class QwtPlotCanvas;
#include "qwt_symbol.h"
#include "qwt_plot_marker.h"
#include "../PlasticityQT_CLI/TPBrStrainStressDataBase.h"
#include "../PlasticityQT_CLI/TPBrDataControl.h"
#include "../PlasticityQT_CLI/TPBrLaboratoryData.h"

class Zoomer: public QwtPlotZoomer
{
public:
    Zoomer(int xAxis, int yAxis, QwtPlotCanvas *canvas):
        QwtPlotZoomer(xAxis, yAxis, canvas)
    {
        setTrackerMode(QwtPicker::AlwaysOff);
        setRubberBand(QwtPicker::NoRubberBand);

        // RightButton: zoom out to full size

        setMousePattern(QwtEventPattern::MouseSelect2,
            Qt::RightButton, Qt::ControlModifier);
        setMousePattern(QwtEventPattern::MouseSelect3,
            Qt::RightButton);
    }
};

// This struct contains all QWT curves created for one Plot
class CURVE {
public:

    CURVE () {
        this->global_id = -1;
        this->chk_status = -1;
        this->initialPnt = -1;
        this->endPnt = -1;
        this->Xsmallest = 0;
        this->Xbiggest = 0;
        this->Ysmallest = 0;
        this->Ybiggest = 0;
        curve_ptr = NULL;
        curve_ptr2 = NULL;
        curve_ptr3 = NULL;
        symb_curve_ptr = NULL;
        mark_ptr1 = NULL;
        mark_ptr2 = NULL;
    }

    CURVE (int global_id) {
        this->global_id = global_id;
    }

    CURVE (const CURVE &copy) {
        curve_ptr = copy.curve_ptr;
        curve_ptr2 = copy.curve_ptr2;
        curve_ptr3 = copy.curve_ptr3;
        symb_curve_ptr = copy.symb_curve_ptr;
        mark_ptr1 = copy.mark_ptr1;
        mark_ptr2 = copy.mark_ptr2;
        chk_status = copy.chk_status;
        Xsmallest = copy.Xsmallest;
        Xbiggest = copy.Xbiggest;
        Ysmallest = copy.Ysmallest;
        Ybiggest = copy.Ybiggest;
    }

    //Operator =
    CURVE & operator=(const CURVE &copy) {
        if(this == &copy)
            return *this;
        curve_ptr = copy.curve_ptr;
        curve_ptr2 = copy.curve_ptr2;
        curve_ptr3 = copy.curve_ptr3;
        symb_curve_ptr = copy.symb_curve_ptr;
        mark_ptr1 = copy.mark_ptr1;
        mark_ptr2 = copy.mark_ptr2;
        chk_status = copy.chk_status;
        Xsmallest = copy.Xsmallest;
        Xbiggest = copy.Xbiggest;
        Ysmallest = copy.Ysmallest;
        Ybiggest = copy.Ybiggest;
        return *this;
    }

    //Destructor
    ~CURVE () {
    }

    int initialPnt;
    int endPnt;
    int global_id;
    double Xsmallest;
    double Xbiggest;
    double Ysmallest;
    double Ybiggest;
    int chk_status;
    QwtPlotCurve *curve_ptr;
    QwtPlotCurve *curve_ptr2;
    QwtPlotCurve *curve_ptr3;
//    QVector <double> X2;
//    QVector <double> Y2;
//    QVector <double> X3;
//    QVector <double> Y3;
    QwtPlotCurve *symb_curve_ptr;
    QwtPlotMarker *mark_ptr1;
    QwtPlotMarker *mark_ptr2;
    QwtSymbol *symb1;
    QwtSymbol *symb2;
};

class Plot: public QwtPlot
{
    Q_OBJECT

private slots:

public:
    enum pointType {startPoint = 0, endPoint = 1};

    Zoomer *zoomer;
    QwtPlotPicker *picker;
    QwtPlotPanner *panner;
    QwtPlotGrid   *grid;
    CanvasPicker  *canvas_picker;
    QHash <int, CURVE> CurvesList;
    QMenu *menuPlot;
//    QHash <int, QAction*> *ActionList;
    TPBrStrainStressDataBase::ECurveType curvetype;

    Plot( QWidget *parent = NULL );
    void setHighlighted (bool status);
    bool isHighlighted ();

    inline void createCurve (int global_id, int check_status){
//        qDebug() << "CREATE CURVE inicio!!!!!!";
//        int count = DADOS.SizeLabData();
//        qDebug() << "COUNTM = " <<count << "!!!!!!!!!!!";
//        int counts = DADOS.fMedicoes[global_id].fSimulacoes.size();
//        qDebug() << "COUNTS = " <<counts << "!!!!!!!!!!!";


        // VERIFICANDO SE A CURVA JA EXISTE PARA NAO RECRIA-LA
        if (this->CurvesList.contains(global_id)) return;

        CURVE new_curve(global_id);

        this->zoomer->setZoomBase();

        new_curve.curve_ptr = new QwtPlotCurve();
        new_curve.curve_ptr->setRenderHint( QwtPlotItem::RenderAntialiased );
        new_curve.curve_ptr->setPen( QPen( Qt::red ) );
        new_curve.curve_ptr->setLegendAttribute( QwtPlotCurve::LegendShowLine );
        new_curve.curve_ptr->setYAxis( QwtPlot::yLeft );
        new_curve.curve_ptr->attach( this );
        std::vector<REAL> X;
        std::vector<REAL> Y;

        TPBrStrainStressDataBase *labData = DADOS.getObj(global_id);
        labData->GeneratePlot(curvetype, X, Y);

        int start_idx = labData->Get_start_idx();
        int end_idx = labData->Get_end_idx();

        ////////////////////////////////////////////////////////////////////////
        qDebug() << "CREATE CURVE!!!!!!";
        int count2 = DADOS.SizeLabData();
        qDebug() << "COUNTM = " <<count2 << "!!!!!!!!!!!";
        //int counts2 = labData->SizeSimData();
        //qDebug() << "COUNTS = " <<counts2 << "!!!!!!!!!!!";
//        std::map<int,TPBrLaboratoryData>::iterator it;
//        qDebug() << "ID DAS MEDICOES:";
//        for(it = DADOS.fMedicoes.begin(); it != DADOS.fMedicoes.end(); ++it) {
//            qDebug() << (*it).first;
//        }
//        qDebug() << endl;
        DADOS.Print();
        ////////////////////////////////////////////////////////////////////////

        QVector <double> qVec_X = QVector <double>::fromStdVector(X);
        QVector <double> qVec_Y = QVector <double>::fromStdVector(Y);
        new_curve.curve_ptr->setSamples(qVec_X,qVec_Y);

        new_curve.curve_ptr2 = new QwtPlotCurve();
        new_curve.curve_ptr2->setRenderHint( QwtPlotItem::RenderAntialiased );
        new_curve.curve_ptr2->setPen( QPen( Qt::blue ) );
        new_curve.curve_ptr2->setLegendAttribute( QwtPlotCurve::LegendShowLine );
        new_curve.curve_ptr2->setYAxis( QwtPlot::yLeft );
        new_curve.curve_ptr2->attach( this );
        new_curve.curve_ptr2->hide();

        new_curve.curve_ptr3 = new QwtPlotCurve();
        new_curve.curve_ptr3->setRenderHint( QwtPlotItem::RenderAntialiased );
        new_curve.curve_ptr3->setPen( QPen( Qt::green ) );
        new_curve.curve_ptr3->setLegendAttribute( QwtPlotCurve::LegendShowLine );
        new_curve.curve_ptr3->setYAxis( QwtPlot::yLeft );
        new_curve.curve_ptr3->attach( this );
        new_curve.curve_ptr3->hide();

        double Xsmallest=numeric_limits<double>::max(), Xbiggest=numeric_limits<double>::min(),
               Ysmallest=numeric_limits<double>::max(), Ybiggest=numeric_limits<double>::min();

        for (int i = 0; i < X.size(); i++) {
            // scale
            if (X[i]<Xsmallest) Xsmallest = X[i];
            if (X[i]>Xbiggest)  Xbiggest = X[i];
            if (Y[i]<Ysmallest) Ysmallest = Y[i];
            if (Y[i]>Ybiggest) Ybiggest = Y[i];
        }

        new_curve.Xbiggest = Xbiggest;
        new_curve.Xsmallest = Xsmallest;
        new_curve.Ybiggest = Ybiggest;
        new_curve.Ysmallest = Ysmallest;

        new_curve.symb_curve_ptr = new QwtPlotCurve();
        new_curve.symb_curve_ptr->setRenderHint( QwtPlotItem::RenderAntialiased );
        new_curve.symb_curve_ptr->setPen( QPen( Qt::blue ) );
        new_curve.symb_curve_ptr->setLegendAttribute( QwtPlotCurve::LegendShowLine );
        new_curve.symb_curve_ptr->setYAxis( QwtPlot::yLeft );
        new_curve.symb_curve_ptr->setStyle(QwtPlotCurve::NoCurve);
        new_curve.symb_curve_ptr->attach( this );

        int sizes = 2;
        double Xs[2], Ys[2];
        Xs[0] = X[start_idx];
        Ys[0] = Y[start_idx];
        Xs[1] = X[end_idx];
        Ys[1] = Y[end_idx];
        qDebug() << "Size: " << X.size() << " start: " << X[start_idx] << " " << Y[start_idx] << endl;
        qDebug() << "Size: " << X.size() << " end: " << X[end_idx] << " " << Y[end_idx] ;
        new_curve.symb_curve_ptr->setSamples( Xs, Ys, sizes);
        new_curve.symb_curve_ptr->hide();

        // "Start" marker
        new_curve.symb1 = new QwtSymbol();
        new_curve.symb1->setBrush(QBrush(Qt::red, Qt::SolidPattern));
        new_curve.symb1->setStyle(QwtSymbol::Ellipse);
        new_curve.symb1->setSize(10);
        new_curve.mark_ptr1 = new QwtPlotMarker();
        new_curve.mark_ptr1->setSymbol(new_curve.symb1);
        new_curve.mark_ptr1->setXValue(Xs[0]);
        new_curve.mark_ptr1->setYValue(Ys[0]);
        new_curve.mark_ptr1->setLabelAlignment(Qt::Alignment(Qt::AlignLeft));
        new_curve.mark_ptr1->setLabel(QwtText("Start"));
        new_curve.mark_ptr1->attach(this);
        new_curve.mark_ptr1->hide();

        // End marker
        new_curve.symb2 = new QwtSymbol();
        new_curve.symb2->setBrush(QBrush(Qt::blue, Qt::SolidPattern));
        new_curve.symb2->setStyle(QwtSymbol::Ellipse);
        new_curve.symb2->setSize(10);
        new_curve.mark_ptr2 = new QwtPlotMarker ();
        new_curve.mark_ptr2->setSymbol(new_curve.symb2);
        new_curve.mark_ptr2->setXValue(Ys[1]);
        new_curve.mark_ptr2->setYValue(Ys[1]);
        new_curve.mark_ptr2->setLabelAlignment(Qt::Alignment(Qt::AlignLeft));
        new_curve.mark_ptr2->setLabel(QwtText("End"));
        new_curve.mark_ptr2->attach(this);
        new_curve.mark_ptr2->hide();

        new_curve.initialPnt = start_idx;
        new_curve.endPnt = end_idx;

        new_curve.chk_status = check_status;

        this->CurvesList.insert(global_id, new_curve);

        this->replot();
        this->AdjustScale();

        //naming axis
        QString DefaultTitleY = QChar (0x03C3);
        DefaultTitleY.append("ax (MPa)");
        this->setAxisTitle(QwtPlot::yLeft, DefaultTitleY);
        QString DefaultTitle1X = QChar (0x03B5);
        DefaultTitle1X.append("ax (%)");
        this->setAxisTitle(QwtPlot::xBottom, DefaultTitle1X);
    }

    inline void deleteCurve (int global_id){

        QwtPlotCurve *d_curve = this->CurvesList[global_id].curve_ptr;
        QwtPlotCurve *d_curve2 = this->CurvesList[global_id].curve_ptr2;
        QwtPlotCurve *d_curve3 = this->CurvesList[global_id].curve_ptr3;
        if ( d_curve != NULL) {
            d_curve->detach();
            delete d_curve;
            d_curve2->detach();
            delete d_curve2;
            d_curve3->detach();
            delete d_curve3;
        }
        QwtPlotCurve *d_curve_symbol = this->CurvesList[global_id].symb_curve_ptr;
        if ( d_curve_symbol != NULL ) {
            //removing markers
            QwtPlotMarker *d_plot_marker = NULL;
            d_plot_marker = this->CurvesList[global_id].mark_ptr1;
            d_plot_marker->detach();
            delete d_plot_marker;
            d_plot_marker = this->CurvesList[global_id].mark_ptr2;
            d_plot_marker->detach();
            delete d_plot_marker;
            //removing curve
            d_curve_symbol->detach();
            delete d_curve_symbol;
        }
        this->CurvesList.remove(global_id);

        this->replot();
        this->AdjustScale();

    }

    inline void AdjustScale () {

        double Xsmallest=numeric_limits<double>::max(), Xbiggest=numeric_limits<double>::min(),
               Ysmallest=numeric_limits<double>::max(), Ybiggest=numeric_limits<double>::min();

        foreach (int i, this->CurvesList.keys()) {

            if (this->CurvesList[i].Xsmallest<Xsmallest)
                Xsmallest = this->CurvesList.value(i).Xsmallest;

            if (this->CurvesList[i].Xbiggest>Xbiggest)
                Xbiggest = this->CurvesList[i].Xbiggest;

            if (this->CurvesList[i].Ysmallest<Ysmallest)
                Ysmallest = this->CurvesList[i].Ysmallest;

            if (this->CurvesList[i].Ybiggest>Ybiggest)
                Ybiggest = this->CurvesList[i].Ybiggest;
        }

                if (this->CurvesList.size() == 0)
                {
                    this->setAxisScale (QwtPlot::xBottom, 0, 10);
                    this->setAxisScale (QwtPlot::yLeft, 0, 10);
                }
                else {
                    this->setAxisScale (QwtPlot::xBottom, Xsmallest, Xbiggest);
                    this->setAxisScale (QwtPlot::yLeft, Ysmallest, Ybiggest);
                }

        this->replot();
    }


public Q_SLOTS:
    void mouseReleaseEvent( QMouseEvent * event );

    void AxisChanged_slot (QAction *action);
};

#endif // PLOT_H
