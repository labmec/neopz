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

    //QHash <int, CURVE > CurvesList;
    //testTypes testType;
    QMenu *menuPlot;

    QHash <int, QAction*> *ActionList;

    Plot( QWidget *parent = NULL );
    void AdjustScale ();
    void setHighlighted (bool status);
    bool isHighlighted ();

public Q_SLOTS:
    void mouseReleaseEvent( QMouseEvent * event );

    void AxisChanged_slot (QAction *action);

signals:
    void AxisChanged_signal (Plot* plot_ptr, QString XCoord);
};

#endif // PLOT_H
