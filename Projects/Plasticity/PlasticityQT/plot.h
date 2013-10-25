#ifndef PLOT_H
#define PLOT_H

#include <QMenu>
#include "common.h"
#include <qwt_plot.h>
class QwtPlotCurve;
//class QwtPlotZoomer;
#include <qwt_plot_zoomer.h>
class QwtPlotPanner;
class QwtPlotGrid;
class CanvasPicker;
class QwtPlotCanvas;
class QwtPlotMarker;

// This struct contains all QWT curves created for one Plot
class CURVE {
public:
    CURVE () : initialPnt (0), endPnt (0) {
    }

//    CURVE (QVector<double> &X, QVector<double> &Y, Plot *plot) {
//        setX(X);
//        setY(Y);
//        attach(plot);
//    }

    //Copy constructor
    CURVE (const CURVE &copy) {
        initialPnt = copy.initialPnt;
        endPnt = copy.endPnt;
        X = QVector<double> (copy.X);
        Y = QVector<double> (copy.Y);
    }

//    //Operator =
//    CURVE & operator=(const CURVE &copy) {
//        if(this == &copy)
//            return *this;
//    }

    //Destructor
    ~CURVE () {

    }

    void cutCurve (int initialPnt_n, int endPnt_n) {
        initialPnt = initialPnt_n;
        endPnt = endPnt_n;
    }

    int getInitialPoint () {
        return initialPnt;
    }
    int getEndPoint () {
        return endPnt;
    }
    void resetCurve () {
        initialPnt = 0;
        endPnt = X.size() -1;
    }
    int size () {
        return (endPnt - initialPnt) + 1;
    }

    QVector<double> getX () {
        return X.mid(initialPnt, endPnt-initialPnt +1);
    }
    QVector<double> getY () {
        return Y.mid(initialPnt, endPnt-initialPnt +1);
    }

    void setX (QVector<double> &newX) {
        X = newX;
        resetCurve();
        //re-adjusting scale
        double Xsmallest=numeric_limits<double>::max(), Xbiggest=numeric_limits<double>::min();
        for (int i = 0; i < X.size() ; i++) {
            // scale
            if (X.value(i) < Xsmallest) Xsmallest = X.value(i);
            if (X.value(i) > Xbiggest)  Xbiggest = X.value(i);
        }
    }
    void setY (QVector<double> &newY) {
        Y = newY;
        resetCurve();
        //re-adjusting scale
        double Ysmallest=numeric_limits<double>::max(), Ybiggest=numeric_limits<double>::min();
        for (int i = 0; i < Y.size() ; i++) {
            // scale
            if (Y.value(i) < Ysmallest) Ysmallest = Y.value(i);
            if (Y.value(i) > Ybiggest)  Ybiggest = Y.value(i);
        }
    }



private:
    int initialPnt;
    int endPnt;
    QVector <double> X;
    QVector <double> Y;
    double Xsmallest;
    double Xbiggest;
    double Ysmallest;
    double Ybiggest;
    int chk_status;
    QwtPlotCurve *curve_ptr;
    QwtPlotCurve *curve_ptr2;
    QwtPlotCurve *curve_ptr3;
    QVector <double> X2;
    QVector <double> Y2;
    QVector <double> X3;
    QVector <double> Y3;
    QwtPlotCurve *symb_curve_ptr;
    QwtPlotMarker *symb1;
    QwtPlotMarker *symb2;
};
//// This struct contains all QWT curves of symbols created for one Plot

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

    QHash <int, CURVE > CurvesList;
    testTypes testType;
    QMenu *menuPlot;

    QHash <int, QAction*> *ActionList;

    Plot( QWidget *parent = NULL );
    //Create curve
    void createCurve(int pos, const QVector <double> &X, const QVector <double> &Y, int checkStatus);
    //Remove curve
    void deleteCurve ( int pos );
    void AdjustScale ();
    // Get index of coords vector of a given curve
    int getSymbIndex ( int indexCurves, pointType typept );
    // Move symbol to a given index of coords vector
    void setSymbIndex ( int indexCoords, int indexCurves, pointType typept );
    void setHighlighted (bool status);
    bool isHighlighted ();
    void updateCurveOLD(int indexCurves, const QVector<double> &X, const QVector<double> &Y);
    void updateCurves();

public Q_SLOTS:
    void exportPlot();
    void showSymbCurve(int indexCurves);
    void hideSymbCurve(int indexCurves);


    void mouseReleaseEvent( QMouseEvent * event );

    void AxisChanged_slot (QAction *action);

signals:
    void AxisChanged_signal (Plot* plot_ptr, QString XCoord);
};

#endif // PLOT_H
