#include <qapplication.h>
#include <qevent.h>
#include <qwhatsthis.h>
#include <qpainter.h>
#include <qwt_plot.h>
#include <qwt_symbol.h>
#include <qwt_scale_map.h>
#include <qwt_plot_canvas.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_directpainter.h>
#include "canvaspicker.h"
#include "plot.h"

CanvasPicker::CanvasPicker( Plot *plot ):
    QObject( plot ),
    d_selectedCurve( /*NULL*/ -1 ),
    d_selectedPoint( -1 )
{
    QwtPlotCanvas *canvas = plot->canvas();

    canvas->installEventFilter( this );

    // We want the focus, but no focus rect. The
    // selected point will be highlighted instead.

    canvas->setFocusPolicy( Qt::StrongFocus );
#ifndef QT_NO_CURSOR
    canvas->setCursor( Qt::PointingHandCursor );
#endif
    canvas->setFocusIndicator( QwtPlotCanvas::ItemFocusIndicator );
    canvas->setFocus();

    pressed_flag = 0;
}

Plot *CanvasPicker::plot()
{
    return qobject_cast<Plot *>( parent() );
}

const Plot *CanvasPicker::plot() const
{
    return qobject_cast<const Plot *>( parent() );
}

bool CanvasPicker::event( QEvent *ev )
{
    if ( ev->type() == QEvent::User )
    {
        return true;
    }
    return QObject::event( ev );
}

bool CanvasPicker::eventFilter( QObject *object, QEvent *event )
{
    if ( plot() == NULL || object != plot()->canvas() )
        return false;

    switch( event->type() )
    {
        case QEvent::Paint:
        {
            QApplication::postEvent( this, new QEvent( QEvent::User ) );
            break;
        }
        case QEvent::MouseButtonPress:
        {
            QMouseEvent *ev = ( QMouseEvent * ) event;
            if(ev->button() == Qt::LeftButton) {
                select( ev->pos() );
                pressed_flag = 1;
            }

            return true;
        }
        case QEvent::MouseButtonRelease:
        {
            pressed_flag = 0;

            QMouseEvent *ev = ( QMouseEvent * ) event;
            if(ev->button() == Qt::LeftButton) {
                emit mouseLeftClicked(this->plot());
            }

            if(ev->button() == Qt::RightButton) {
                this->plot()->AdjustScale();
            }

            return true;
        }
        case QEvent::MouseButtonDblClick:
        {
            QMouseEvent *ev = ( QMouseEvent * ) event;
            if(ev->button() == Qt::LeftButton) {
                emit mouseDoubleClicked(this->plot());
            }
        }
        case QEvent::MouseMove:
        {
            QMouseEvent *ev = ( QMouseEvent * ) event;
            if(pressed_flag) {
                move( ev->pos() );
            }
            return true;
        }
    }

    return QObject::eventFilter( object, event );
}

// Select the point at a position. If there is no point
// deselect the selected point

void CanvasPicker::select( const QPoint &pos )
{
    double dist = 10e10, dists = 10e10;;
    int index = -1;
    int indexs = -1;
    int key = -1;

    foreach (int it, this->plot()->CurvesList.keys()){
        QwtPlotCurve *c=this->plot()->CurvesList.value(it).curve_ptr;
        QwtPlotCurve *cs=this->plot()->CurvesList.value(it).symb_curve_ptr;

        // if curves dont exist: abort
        if ( !c || !cs )
            continue;

        // check if curve is visible (if it is enabled to have its points selected)
        if ( cs->isVisible() == false )
            continue;

        double d, ds;
        int idx = c->closestPoint( pos, &d );
        int idxs = cs->closestPoint( pos, &ds );
        if ( d < dist )
        {
            index = idx;
            dist = d;
            key = it;
        }
        if ( ds < dists )
        {
            indexs = idxs;
            dists = ds;
            key = it;
        }
    }

//    if ( dist < 10 ) // 10 pixels tolerance
//    {
        d_selectedCurve = key;
        d_selectedPoint = index;
        d_selectedPointType = indexs;
//    }

//        qDebug() << key << " " << index << " " << indexs ;
}

// Move the selected point
void CanvasPicker::move( const QPoint &pos )
{
    if ( !d_selectedCurve < 0 )
        return;

    double dist = 10e10;
    int index = -1;

    QwtPlotCurve *c=this->plot()->CurvesList.value(d_selectedCurve).curve_ptr;

    // if curve does not exist: abort
    if ( !c )
        return;

    double d;
    int idx = c->closestPoint( pos, &d );
    if ( d < dist )
    {
            index = idx;
            dist = d;
    }

//    if ( d < 10 ) // 10 pixels tolerance
//    {
        d_selectedPoint = index;
//    }

//    this->plot()->setSymbIndex(d_selectedPoint, d_selectedCurve, Plot::pointType(d_selectedPointType));
//    emit SymbPointChanged (d_selectedPoint, d_selectedCurve, Plot::pointType(d_selectedPointType));
}

//slot
void CanvasPicker::setSymbPoint ( int coordIndex, int curveIndex, Plot::pointType typept ) {
//    this->d_selectedPoint = coordIndex;
//    this->d_selectedCurve = curveIndex;
//    this->d_selectedPointType = typept;

//    this->plot()->setSymbIndex(d_selectedPoint, d_selectedCurve, Plot::pointType(d_selectedPointType));
//    emit SymbPointChanged (d_selectedPoint, d_selectedCurve, typept);
}

