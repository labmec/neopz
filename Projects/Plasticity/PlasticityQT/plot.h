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
struct CURVE {int chk_status; QwtPlotCurve *curve_ptr;
              double Xsmallest; double Xbiggest; double Ysmallest; double Ybiggest;
              QVector <double> *X; QVector <double> *Y;
              QwtPlotCurve *symb_curve_ptr;
              int startIndexCoord; int endIndexCoord;
              QwtPlotMarker *startM; QwtPlotMarker *endM;};
//// This struct contains all QWT curves of symbols created for one Plot

class Zoomer: public QwtPlotZoomer
{
public:
    Zoomer(int xAxis, int yAxis, QwtPlotCanvas *canvas):
        QwtPlotZoomer(xAxis, yAxis, canvas)
    {
        setTrackerMode(QwtPicker::AlwaysOff);
        setRubberBand(QwtPicker::NoRubberBand);

        // RightButton: zoom out by 1
        // Ctrl+RightButton: zoom out to full size

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
    enum pointType {
            startPoint = 0,
            endPoint = 1
    };

    Zoomer *zoomer;
    QwtPlotPanner *panner;
    QwtPlotGrid   *grid;
    CanvasPicker  *canvas_picker;
//    QwtPlotPicker *picker;
    QHash <int, CURVE > *CurvesList;
    testTypes testType;
    QMenu *menuPlot;

    Plot( QWidget *parent = NULL );
    //Create curve
    void createCurve ( int pos, QVector <double>  *X, QVector <double>  *Y, int checkStatus);
    //Remove curve
    void deleteCurve ( int pos );
    //Check if curve exists
    int checkCurve ( int pos ) {
        //IMPLEMENTAR
    }
    void AdjustScale ();
    // Get index of coords vector of a given curve
    int getSymbIndex ( int indexCurves, pointType typept );
    // Get coords of coords vector of a given curve
    int getSymbCoords ( int indexCurves, pointType typept ) {
        //IMPLEMENTAR
    }
    // Move symbol to a given index of coords vector
    void setSymbIndex ( int indexCoords, int indexCurves, pointType typept );
    void setHighlighted (bool status);
    bool isHighlighted ();

public Q_SLOTS:
    void exportPlot();
    void showSymbCurve(int indexCurves);
    void hideSymbCurve(int indexCurves);
    void cutCurve(int indexCurves, int indexStartPoint, int indexEndPoint);
    void updateCurve(int indexCurves, int indexStartPoint, int indexEndPoint);

    void mouseReleaseEvent( QMouseEvent * event );

    void AxisChanged_slot (QAction *action);
//    void YTypeChanged (QAction *action);

signals:
    void AxisChanged_signal (Plot* plot_ptr, QString XCoord);
//    void YCoordsChanged (Plot* plot_ptr, QString YCoord);
};

#endif // PLOT_H
