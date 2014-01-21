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
#include "TPBrStrainStressDataBase.h"
#include "TPBrDataControl.h"
#include "TPBrLaboratoryData.h"

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
        this->elastic_transPnt = -1;
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
        mark_ptr3 = NULL;
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
        mark_ptr3 = copy.mark_ptr3;
        chk_status = copy.chk_status;
        Xsmallest = copy.Xsmallest;
        Xbiggest = copy.Xbiggest;
        Ysmallest = copy.Ysmallest;
        Ybiggest = copy.Ybiggest;
        X1 = copy.X1;
        Y1 = copy.Y1;
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
        mark_ptr3 = copy.mark_ptr3;
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
    int elastic_transPnt;
    int global_id;
    double Xsmallest;
    double Xbiggest;
    double Ysmallest;
    double Ybiggest;
    int chk_status;
    QwtPlotCurve *curve_ptr;
    QwtPlotCurve *curve_ptr2;
    QwtPlotCurve *curve_ptr3;
    QVector <double> X1;
    QVector <double> Y1;
//    QVector <double> X2;
//    QVector <double> Y2;
//    QVector <double> X3;
//    QVector <double> Y3;
    QwtPlotCurve *symb_curve_ptr;
    QwtPlotMarker *mark_ptr1;
    QwtPlotMarker *mark_ptr2;
    QwtPlotMarker *mark_ptr3;
    QwtSymbol *symb1;
    QwtSymbol *symb2;
    QwtSymbol *symb3;
};

enum pointType {
        startPoint = 0,
        endPoint = 2,
        elasticTransPoint = 1
};

class Plot: public QwtPlot
{
    Q_OBJECT

private slots:

public:
    enum pointType {
        startPoint = 0,
        endPoint = 2,
        elasticTransPoint = 1
};

    Zoomer *zoomer;
    QwtPlotPicker *picker;
    QwtPlotPanner *panner;
    QwtPlotGrid   *grid;
    CanvasPicker  *canvas_picker;
    QHash <int, CURVE> CurvesList;
    QMenu *menuPlot;
    TPBrStrainStressDataBase::ECurveType curvetype;
    QwtPlotCurve *envelope_curve;
    int Xenv_max, Xenv_min;
    int Yenv_max, Yenv_min;
    int Env_status;

    Plot( QWidget *parent = NULL );
    void setHighlighted (bool status);
    bool isHighlighted ();

    // Move symbol to a given index of coords vector
    void setSymbIndex (int indexCoords, int global_id, pointType typept );

    void createCurve (int global_id, int check_status);

    void deleteCurve (int global_id);

    void Generate_Envelope(double A, double B, double C);

    void hide_envelope();

    void toggleSymb (int global_id, bool showsymb);

    void AdjustScale ();


public Q_SLOTS:
    void mouseReleaseEvent( QMouseEvent * event );

    void AxisChanged_slot (QAction *action);
};

#endif // PLOT_H
