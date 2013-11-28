#include <QObject>
#include <QMouseEvent>
#include <QImageWriter>
#include <QFileDialog>
#include <limits>
#include <QTextDocument>
#include <QPrinter>
#include <QPainter>
using namespace std;

#include "plot.h"

#include <QMenu>
#include <qwt_symbol.h>
#include <qwt_plot_renderer.h>
#include <qwt_plot_marker.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_grid.h>
#include <qwt_plot.h>
#include <qwt_plot_zoomer.h>
#include <qwt_plot_panner.h>
#include <qwt_picker_machine.h>
#include <qwt_scale_engine.h>
#include "canvaspicker.h"

Plot::Plot( QWidget *parent ): QwtPlot (parent) {
    setCanvasBackground(QColor(Qt::white));
    setFont(QFont("System",-1,10,true));

    grid = new QwtPlotGrid();
    grid->enableXMin( true );
    grid->setMajPen( QPen( Qt::gray, 0, Qt::DashLine ) );
    grid->setMinPen( QPen( Qt::gray, 0 , Qt::DashLine ) );
    grid->attach( this );

    zoomer = new Zoomer(QwtPlot::xBottom, QwtPlot::yLeft, this->canvas());
    zoomer->setRubberBand(QwtPicker::RectRubberBand);
    zoomer->setRubberBandPen(QColor(Qt::green));
    zoomer->setTrackerMode(QwtPicker::ActiveOnly);
    zoomer->setTrackerPen(QColor(Qt::blue));


    //testType = testTypes(All);

    // BEGIN MENU
    QFont MenuBold;
    MenuBold.setBold(true);

    // Sigma ax x Epsilon ax
    QString Coords1 = QChar (0x03C3);
    Coords1.append("ax x ");
    Coords1.append(QChar (0x03B5));
    Coords1.append("ax");
    // Sigma rad x Epsilon rad
    QString Coords2 = QChar (0x03C3);
    Coords2.append("rad x ");
    Coords2.append(QChar (0x03B5));
    Coords2.append("rad");
    // Sigma vol x Epsilon vol
    QString Coords3 = QChar (0x03C3);
    Coords3.append("vol x ");
    Coords3.append(QChar (0x03B5));
    Coords3.append("vol");
    // Sigma ax x Epsilon ax Epsilon vol Epsilon rad
    QString Coords4 = QChar (0x03C3);
    Coords4.append("ax x ");
    Coords4.append(QChar (0x03B5));
    Coords4.append("ax ");
    Coords4.append(QChar (0x03B5));
    Coords4.append("vol ");
    Coords4.append(QChar (0x03B5));
    Coords4.append("rad");
    // I1 x sqrt J2
    QString Coords5 = "I1 x ";
    Coords5.append(QChar (0x221A));
    Coords5.append("J2");
    // Epsilon v x sqrt J2Epsilon
//    QString Coords6 = QChar (0x03B5);
//    Coords6.append("v x ");
//    Coords6.append(QChar (0x221A));
//    Coords6.append("J2");
//    Coords6.append(QChar (0x03B5));

    this->ActionList = new QHash <int, QAction*>;

    QMenu *mySubMenuAxis = new QMenu ("Set Axis");
    QActionGroup *groupSubmenuAxis = new QActionGroup (this);
    QAction *aSetAxis1 = mySubMenuAxis->addAction(Coords1);
    aSetAxis1->setCheckable(true);
    aSetAxis1->setChecked(true);
    aSetAxis1->setActionGroup(groupSubmenuAxis);
    ActionList->insert(0, aSetAxis1);
    QAction *aSetAxis2 = mySubMenuAxis->addAction(Coords2);
    aSetAxis2->setCheckable(true);
    aSetAxis2->setActionGroup(groupSubmenuAxis);
    ActionList->insert(1, aSetAxis2);
    QAction *aSetAxis3 = mySubMenuAxis->addAction(Coords3);
    aSetAxis3->setCheckable(true);
    aSetAxis3->setActionGroup(groupSubmenuAxis);
    ActionList->insert(2, aSetAxis3);
    QAction *aSetAxis4 = mySubMenuAxis->addAction(Coords4);
    aSetAxis4->setCheckable(true);
    aSetAxis4->setActionGroup(groupSubmenuAxis);
    ActionList->insert(3, aSetAxis4);
    QAction *aSetAxis5 = mySubMenuAxis->addAction(Coords5);
    aSetAxis5->setCheckable(true);
    aSetAxis5->setActionGroup(groupSubmenuAxis);
    ActionList->insert(4, aSetAxis5);
//    QAction *aSetAxis6 = mySubMenuAxis->addAction(Coords6);
//    aSetAxis6->setCheckable(true);
//    aSetAxis6->setActionGroup(groupSubmenuAxis);
//    ActionList->insert(5, aSetAxis6);

    menuPlot = new QMenu ();
    menuPlot->setFont(MenuBold);
    menuPlot->addMenu(mySubMenuAxis);

    connect (mySubMenuAxis, SIGNAL(triggered(QAction*)), this, SLOT(AxisChanged_slot (QAction*)));
    // END MENU

    // Picker with coords and green indication
    canvas_picker = new CanvasPicker( this );

    picker = new QwtPlotPicker(QwtPlot::xBottom, QwtPlot::yLeft,
                   QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
                   this->canvas());
    picker->setStateMachine(new QwtPickerDragPointMachine());
    picker->setRubberBandPen(QColor(Qt::green));
    picker->setRubberBand(QwtPicker::CrossRubberBand);
    picker->setTrackerPen(QColor(Qt::blue));

    // Mid-Button panner
    panner = new QwtPlotPanner(this->canvas());
    panner->setMouseButton(Qt::MidButton);

//    CurvesList = new QHash <int, CURVE > ;

    // scale divisions (Majors/Minors) and Engines
    this->setAxisMaxMajor( QwtPlot::xBottom, 5 );
    this->setAxisMaxMinor( QwtPlot::xBottom, 10 );
    this->setAxisMaxMajor( QwtPlot::yLeft, 5 );
    this->setAxisMaxMinor( QwtPlot::yLeft, 10 );
    this->setAxisScaleEngine( QwtPlot::xBottom, new QwtLinearScaleEngine );
    this->setAxisScaleEngine( QwtPlot::yLeft, new QwtLinearScaleEngine );

    this->setTitle("");
    this->setHighlighted(false);
}

void Plot::setHighlighted (bool status) {
    QwtText titlePlot = this->title();
    QFont fontTitle=titlePlot.font();
    fontTitle.setBold(status);
    titlePlot.setFont(fontTitle);
    this->setTitle(titlePlot);
}

bool Plot::isHighlighted () {
    QwtText titlePlot = this->title();
    QFont fontTitle=titlePlot.font();
    return fontTitle.bold();
}


void Plot::AdjustScale ()
{
//    double Xsmallest=numeric_limits<double>::max(), Xbiggest=numeric_limits<double>::min(),
//           Ysmallest=numeric_limits<double>::max(), Ybiggest=numeric_limits<double>::min();

//    foreach (int i, this->CurvesList.keys()) {
//        if (this->CurvesList.value(i).Xsmallest<Xsmallest)
//            Xsmallest = this->CurvesList.value(i).Xsmallest;

//        if (this->CurvesList.value(i).Xbiggest>Xbiggest)
//            Xbiggest = this->CurvesList.value(i).Xbiggest;

//        if (this->CurvesList.value(i).Ysmallest<Ysmallest)
//            Ysmallest = this->CurvesList.value(i).Ysmallest;

//        if (this->CurvesList.value(i).Ybiggest>Ybiggest)
//            Ybiggest = this->CurvesList.value(i).Ybiggest;
//    }

//    if (this->CurvesList.size() == 0)
//    {
//        // axis scale (values)
//        this->setAxisScale (QwtPlot::xBottom, 0, 10);
//        this->setAxisScale (QwtPlot::yLeft, 0, 10);
//    }
//    else {
//        this->setAxisScale (QwtPlot::xBottom, Xsmallest, Xbiggest);
//        this->setAxisScale (QwtPlot::yLeft, Ysmallest, Ybiggest);
//    }

//    this->replot();
}


//Click event = MousePress->MouseRelease
void Plot::mouseReleaseEvent( QMouseEvent * event ) {
    if(event->button() == Qt::RightButton) {
        QPoint globalPos = this->mapToGlobal(event->pos()); // for most widgets
        menuPlot->exec(globalPos);
    }
}

void Plot::AxisChanged_slot (QAction *action) {

}
