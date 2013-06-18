#ifndef CANVASPICKER_H
#define CANVASPICKER_H

#include <qobject.h>
#include "plot.h" //for point type enum

class QPoint;
class QMouseEvent;
class QCustomEvent;
class QwtPlotCurve;

class CanvasPicker: public QObject
{
    Q_OBJECT

public:
    CanvasPicker( Plot *plot );
    virtual bool eventFilter( QObject *, QEvent * );
    virtual bool event( QEvent * );

private:
    void select( const QPoint & );
    void move( const QPoint & );
    void release();
    Plot *plot();
    const Plot *plot() const;
    int d_selectedCurve;
    int d_selectedPoint;
    int d_selectedPointType;
    int pressed_flag;

signals:
    void SymbPointChanged ( int coordIndex, int curveIndex, Plot::pointType typept );
    void mouseLeftClicked( Plot * );
    void mouseDoubleClicked( Plot * );

public Q_SLOTS:
    void setSymbPoint ( int coordIndex, int curveIndex, Plot::pointType typept );
};

#endif // CANVASPICKER_H
