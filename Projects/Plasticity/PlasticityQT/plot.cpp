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


    testType = testTypes(All);

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

void Plot::createCurve(int pos, const QVector <double> &X, const QVector <double> &Y, int checkStatus)
{
    // VERIFICANDO SE A CURVA JA EXISTE PARA NAO RECRIA-LA
    if (!this->CurvesList.contains(pos)) return;

    CURVE new_curve (X, Y, this);

    this->zoomer->setZoomBase();

    this->CurvesList.insert(pos, new_curve);

    this->replot();
    this->AdjustScale();

    //naming axis
    QString DefaultTitleY = QChar (0x03C3);
    DefaultTitleY.append("ax (MPa)");
    this->setAxisTitle(QwtPlot::yLeft, DefaultTitleY);
    QString DefaultTitle1X = QChar (0x03B5);
    DefaultTitle1X.append("ax (%)");
    this->setAxisTitle(QwtPlot::xBottom, DefaultTitle1X);

    //refreshing axis menu
    QAction *Default_act = ActionList->value(0);
    Default_act->setChecked(true);

//    // VERIFICANDO SE A CURVA JA EXISTE PARA NAO RECRIA-LA
//    if (this->CurvesList[pos].curve_ptr) return;

//    QwtPlotCurve *d_curve, *d_curve2, *d_curve3;
//    d_curve = new QwtPlotCurve();
//    d_curve->setRenderHint( QwtPlotItem::RenderAntialiased );
//    d_curve->setPen( QPen( Qt::red ) );
//    d_curve->setLegendAttribute( QwtPlotCurve::LegendShowLine );
//    d_curve->setYAxis( QwtPlot::yLeft );
//    d_curve->attach( this );
//    d_curve->setSamples(X, Y);

//    d_curve2 = new QwtPlotCurve();
//    d_curve2->setRenderHint( QwtPlotItem::RenderAntialiased );
//    d_curve2->setPen( QPen( Qt::blue ) );
//    d_curve2->setLegendAttribute( QwtPlotCurve::LegendShowLine );
//    d_curve2->setYAxis( QwtPlot::yLeft );
//    d_curve2->attach( this );
//    d_curve2->hide();

//    d_curve3 = new QwtPlotCurve();
//    d_curve3->setRenderHint( QwtPlotItem::RenderAntialiased );
//    d_curve3->setPen( QPen( Qt::green ) );
//    d_curve3->setLegendAttribute( QwtPlotCurve::LegendShowLine );
//    d_curve3->setYAxis( QwtPlot::yLeft );
//    d_curve3->attach( this );
//    d_curve3->hide();

//    this->zoomer->setZoomBase();

//    double Xsmallest=numeric_limits<double>::max(), Xbiggest=numeric_limits<double>::min(),
//           Ysmallest=numeric_limits<double>::max(), Ybiggest=numeric_limits<double>::min();

//    for (int i = 0; i < X.size(); i++) {
//        // scale
//        if (X.value(i)<Xsmallest) Xsmallest = X.value(i);
//        if (X.value(i)>Xbiggest)  Xbiggest = X.value(i);
//        if (Y.value(i)<Ysmallest) Ysmallest = Y.value(i);
//        if (Y.value(i)>Ybiggest)  Ybiggest = Y.value(i);
//    }

//    QwtPlotCurve *d_curve_symbol;
//    d_curve_symbol = new QwtPlotCurve( );
//    d_curve_symbol->setRenderHint( QwtPlotItem::RenderAntialiased );
//    d_curve_symbol->setPen( QPen( Qt::blue ) );
//    d_curve_symbol->setLegendAttribute( QwtPlotCurve::LegendShowLine );
//    d_curve_symbol->setYAxis( QwtPlot::yLeft );
//    d_curve_symbol->setStyle(QwtPlotCurve::NoCurve);
//    d_curve_symbol->attach( this );

//    int sizes = 2;
//    double Xs[2], Ys[2];
//    Xs[0] = X.value(0);
//    Ys[0] = Y.value(0);
//    Xs[1] = X.value(X.size()-1);
//    Ys[1] = Y.value(Y.size()-1);
//    qDebug() << "Size: " << X.size() << " start: " << X.value(0) << " " << Y.value(0) << endl;
//    qDebug() << "Size: " << X.size() << " end: " << X.value(X.size()-1) << " " << Y.value(X.size()-1) ;
//    d_curve_symbol->setSamples( Xs, Ys, sizes);
//    d_curve_symbol->hide();

//    // "Start" marker
//    QwtSymbol *symb1 = new QwtSymbol();
//    symb1->setBrush(QBrush(Qt::red, Qt::SolidPattern));
//    symb1->setStyle(QwtSymbol::Ellipse);
//    symb1->setSize(10);
//    QwtPlotMarker *a = new QwtPlotMarker ();
//    a->setSymbol(symb1);
//    a->setXValue(Xs[0]);
//    a->setYValue(Ys[0]);
////    a->setLabelOrientation(Qt::Orientation(Qt::Vertical));
//    a->setLabelAlignment(Qt::Alignment(Qt::AlignLeft));
//    a->setLabel(QwtText("Start"));
//    a->attach(this);
//    a->hide();

//    // End marker
//    QwtSymbol *symb2 = new QwtSymbol();
//    symb2->setBrush(QBrush(Qt::blue, Qt::SolidPattern));
//    symb2->setStyle(QwtSymbol::Ellipse);
//    symb2->setSize(10);
//    QwtPlotMarker *b = new QwtPlotMarker ();
//    b->setSymbol(symb2);
//    b->setXValue(Ys[1]);
//    b->setYValue(Ys[1]);
//    b->setLabelAlignment(Qt::Alignment(Qt::AlignLeft));
//    b->setLabel(QwtText("End"));
//    b->attach(this);
//    b->hide();

//    CURVE curva;
//    curva.chk_status = checkStatus;
//    curva.curve_ptr = d_curve;
//    curva.curve_ptr2 = d_curve2;
//    curva.curve_ptr3 = d_curve3;
//    curva.Xbiggest = Xbiggest;
//    curva.Xsmallest = Xsmallest;
//    curva.Ybiggest = Ybiggest;
//    curva.Ysmallest = Ysmallest;
//    int initialPnt = 0;
//    int endPnt = X.size() -1;
//    curva.startIndexCoord = initialPnt;
//    curva.endIndexCoord = endPnt;
//    curva.symb_curve_ptr = d_curve_symbol;
//    curva.startM = a;
//    curva.endM = b;
//    curva.X = X;
//    curva.Y = Y;
//    this->CurvesList.insert(pos, curva);

//    this->replot();
//    this->AdjustScale();

//    //naming axis
//    QString DefaultTitleY = QChar (0x03C3);
//    DefaultTitleY.append("ax (MPa)");
//    this->setAxisTitle(QwtPlot::yLeft, DefaultTitleY);
//    QString DefaultTitle1X = QChar (0x03B5);
//    DefaultTitle1X.append("ax (%)");
//    this->setAxisTitle(QwtPlot::xBottom, DefaultTitle1X);

//    //refreshing axis menu
//    QAction *Default_act = ActionList->value(0);
//    Default_act->setChecked(true);
}

// Get index of coords vector of a given curve
int Plot::getSymbIndex ( int indexCurves, pointType typept ) {
    QwtPlotCurve *d_curve_symbol = this->CurvesList[indexCurves].symb_curve_ptr;

    //Check if curve exists
    if (!d_curve_symbol) {
        qDebug() << "Curve does not exist. Check getSymbIndex method. (returned value: -999)";
        return -999;
    }

    int returnIdx = -1;

    if (typept == startPoint)
        returnIdx = this->CurvesList[indexCurves].startIndexCoord;
    else
        returnIdx = this->CurvesList[indexCurves].endIndexCoord;

    qDebug() << "Plot::getSymbIndex: indexCurves = " << indexCurves << " return Pt idx: " << returnIdx;

    return returnIdx;
}

//Remove curve
void Plot::deleteCurve ( int pos ) {
    QwtPlotCurve *d_curve = this->CurvesList[pos].curve_ptr;
    QwtPlotCurve *d_curve2 = this->CurvesList[pos].curve_ptr2;
    QwtPlotCurve *d_curve3 = this->CurvesList[pos].curve_ptr3;
    if ( d_curve != NULL) {
        d_curve->detach();
        delete d_curve;
        d_curve2->detach();
        delete d_curve2;
        d_curve3->detach();
        delete d_curve3;
    }
    QwtPlotCurve *d_curve_symbol = this->CurvesList[pos].symb_curve_ptr;
    if ( d_curve_symbol != NULL ) {
        //removing markers
        QwtPlotMarker *d_plot_marker = NULL;
        d_plot_marker = this->CurvesList[pos].startM;
        d_plot_marker->detach();
        delete d_plot_marker;
        d_plot_marker = this->CurvesList[pos].endM;
        d_plot_marker->detach();
        delete d_plot_marker;
        //removing curve
        d_curve_symbol->detach();
        delete d_curve_symbol;
    }
    this->CurvesList.remove(pos);

    this->replot();
    this->AdjustScale();
}

void Plot::AdjustScale ()
{
    double Xsmallest=numeric_limits<double>::max(), Xbiggest=numeric_limits<double>::min(),
           Ysmallest=numeric_limits<double>::max(), Ybiggest=numeric_limits<double>::min();

    foreach (int i, this->CurvesList.keys()) {
        if (this->CurvesList.value(i).Xsmallest<Xsmallest)
            Xsmallest = this->CurvesList.value(i).Xsmallest;

        if (this->CurvesList.value(i).Xbiggest>Xbiggest)
            Xbiggest = this->CurvesList.value(i).Xbiggest;

        if (this->CurvesList.value(i).Ysmallest<Ysmallest)
            Ysmallest = this->CurvesList.value(i).Ysmallest;

        if (this->CurvesList.value(i).Ybiggest>Ybiggest)
            Ybiggest = this->CurvesList.value(i).Ybiggest;
    }

    if (this->CurvesList.size() == 0)
    {
        // axis scale (values)
        this->setAxisScale (QwtPlot::xBottom, 0, 10);
        this->setAxisScale (QwtPlot::yLeft, 0, 10);
    }
    else {
        this->setAxisScale (QwtPlot::xBottom, Xsmallest, Xbiggest);
        this->setAxisScale (QwtPlot::yLeft, Ysmallest, Ybiggest);
    }

    this->replot();
}

//slot
void Plot::setSymbIndex ( int indexCoords, int indexCurves, pointType typept ) {
    QwtPlotCurve *d_curve = this->CurvesList[indexCurves].curve_ptr;

    //Check if curve exists
    if (!d_curve) return;

    CURVE d_curve_tmp = this->CurvesList.take(indexCurves);

    //updating symbol position
    int sizes = 2;
    double Xs[2], Ys[2];
    Xs[0] = d_curve_tmp.symb_curve_ptr->sample(0).x();
    Ys[0] = d_curve_tmp.symb_curve_ptr->sample(0).y();
    Xs[1] = d_curve_tmp.symb_curve_ptr->sample(1).x();
    Ys[1] = d_curve_tmp.symb_curve_ptr->sample(1).y();

    d_curve_tmp.startM->setXValue(Xs[0]);
    d_curve_tmp.startM->setYValue(Ys[0]);
    d_curve_tmp.endM->setXValue(Xs[1]);
    d_curve_tmp.endM->setYValue(Ys[1]);

    if (typept == startPoint) {
        //updating curve
        Xs[0] = d_curve_tmp.curve_ptr->sample(indexCoords).x();
        Ys[0] = d_curve_tmp.curve_ptr->sample(indexCoords).y();
        d_curve_tmp.symb_curve_ptr->setSamples( Xs, Ys, sizes);
        //updading markers (labels)
        d_curve_tmp.startM->setXValue(Xs[0]);
        d_curve_tmp.startM->setYValue(Ys[0]);
        //updating indexCoords control
        d_curve_tmp.startIndexCoord = indexCoords;
    }
    else {
        //updating curve
        Xs[1] = d_curve_tmp.curve_ptr->sample(indexCoords).x();
        Ys[1] = d_curve_tmp.curve_ptr->sample(indexCoords).y();
        d_curve_tmp.symb_curve_ptr->setSamples( Xs, Ys, sizes);
        //updading markers (labels)
        d_curve_tmp.endM->setXValue(Xs[1]);
        d_curve_tmp.endM->setYValue(Ys[1]);
        //updating indexCoords control
        d_curve_tmp.endIndexCoord = indexCoords;
    }

    //inserting new value
    this->CurvesList.insert(indexCurves, d_curve_tmp);

    qDebug() << "setSymbIndex: Pt idx: " << indexCoords << " type " << typept;

    this->replot();
    //this->AdjustScale();
}

//slot
void Plot::showSymbCurve(int indexCurves) {
    this->CurvesList[indexCurves].startM->show();
    this->CurvesList[indexCurves].endM->show();
    this->CurvesList[indexCurves].symb_curve_ptr->show();
    this->replot();
}

//slot
void Plot::hideSymbCurve(int indexCurves) {
    this->CurvesList[indexCurves].startM->hide();
    this->CurvesList[indexCurves].endM->hide();
    this->CurvesList[indexCurves].symb_curve_ptr->hide();
    this->replot();
}

void Plot::updateCurves() {

    this->replot();
    this->AdjustScale();
}

//slot
void Plot::updateCurveOLD(int indexCurves, const QVector <double> &X, const QVector <double> &Y) {

//    QwtPlotCurve *d_curve = this->CurvesList[indexCurves].curve_ptr;
//    QwtPlotCurve *d_curve2 = this->CurvesList[indexCurves].curve_ptr2;
//    QwtPlotCurve *d_curve3 = this->CurvesList[indexCurves].curve_ptr3;

//    //Check if curve exists
//    if (!d_curve) return;

//    d_curve->setSamples (X, Y);

//    CurvesList[indexCurves].X = X;
//    CurvesList[indexCurves].Y = Y;

////    qDebug() << X;
////    qDebug() << Y;

//    // Reposiciona marcadores no inicio/fim da curva cortada
//    // Parametros indexStartPoint e indexEndPoint estao com valores antigos.
//    int initialPnt = 0;
//    int endPnt = X.size() - 1;
//    CurvesList[indexCurves].startIndexCoord = initialPnt;
//    CurvesList[indexCurves].endIndexCoord = endPnt;

//    this->setSymbIndex(initialPnt,indexCurves,startPoint);
//    this->setSymbIndex(endPnt,indexCurves,endPoint);


//    qDebug() << "New points (update): " << initialPnt << " " << endPnt;

//    //re-adjusting scale
//    double Xsmallest=numeric_limits<double>::max(), Xbiggest=numeric_limits<double>::min(),
//           Ysmallest=numeric_limits<double>::max(), Ybiggest=numeric_limits<double>::min();
//    for (int i = 0; i < X.size() ; i++) {
//        // scale
//        if (X.value(i) < Xsmallest) Xsmallest = X.value(i);
//        if (X.value(i) > Xbiggest)  Xbiggest = X.value(i);
//        if (Y.value(i) < Ysmallest) Ysmallest = Y.value(i);
//        if (Y.value(i) > Ybiggest)  Ybiggest = Y.value(i);
//    }

//    CurvesList[indexCurves].Xbiggest=Xbiggest;
//    CurvesList[indexCurves].Xsmallest=Xsmallest;
//    CurvesList[indexCurves].Ybiggest=Ybiggest;
//    CurvesList[indexCurves].Ysmallest=Ysmallest;

//    this->replot();
//    this->AdjustScale();
}

//slot
void Plot::exportPlot()
{
//#ifndef QT_NO_PRINTER
//    QString fileName = "export.pdf";
//#else
//    QString fileName = "export.png";
//#endif

//#ifndef QT_NO_FILEDIALOG
//    const QList<QByteArray> imageFormats =
//        QImageWriter::supportedImageFormats();

//    QStringList filter;
//    filter += "PNG Documents (*.png)";
//    filter += "PDF Documents (*.pdf)";
//#ifndef QWT_NO_SVG
//    filter += "SVG Documents (*.svg)";
//#endif
//    filter += "Postscript Documents (*.ps)";

//    if ( imageFormats.size() > 0 )
//    {
//        QString imageFilter( "Images (" );
//        for ( int i = 0; i < imageFormats.size(); i++ )
//        {
//            if ( i > 0 )
//                imageFilter += " ";
//            imageFilter += "*.";
//            imageFilter += imageFormats[i];
//        }
//        imageFilter += ")";

//        filter += imageFilter;
//    }

//    fileName = QFileDialog::getSaveFileName(
//        this, "Export File Name", fileName,
//        filter.join( ";;" ), NULL, QFileDialog::DontConfirmOverwrite );
//#endif

//    if ( !fileName.isEmpty() )
//    {
        QwtPlotRenderer renderer;

        // flags to make the document look like the widget
//        renderer.setDiscardFlag( QwtPlotRenderer::DiscardBackground, true );
//        renderer.setLayoutFlag( QwtPlotRenderer::KeepFrames, true );

        renderer.renderDocument( this, "export.png" /*fileName*/, QSizeF( 150, 100 ), 100 );
//    }

    QTextDocument doc;
    doc.setHtml("<h1>Relatorio</h1><br><br><img src='export.png'>");
    QPrinter printer;
    printer.setOutputFileName("file.pdf");
    printer.setOutputFormat(QPrinter::PdfFormat);
    doc.print(&printer);
    printer.newPage();

}

//Click event = MousePress->MouseRelease
void Plot::mouseReleaseEvent( QMouseEvent * event ) {
    if(event->button() == Qt::RightButton) {
        QPoint globalPos = this->mapToGlobal(event->pos()); // for most widgets
        menuPlot->exec(globalPos);
    }
}

void Plot::AxisChanged_slot (QAction *action) {

  qDebug() << "Axis Changed !!! " << action->text();
  emit AxisChanged_signal(this, action->text());

}
