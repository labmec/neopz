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

    //MENU
    QFont MenuBold;
    MenuBold.setBold(true);


    QMenu *mySubMenuX = new QMenu("Set X");
    QActionGroup* groupSubMenuX = new QActionGroup( this );
    QAction *aSetX1 = mySubMenuX->addAction("Tension");
    aSetX1->setCheckable(true);
    aSetX1->setActionGroup(groupSubMenuX);
    QAction *aSetX2 = mySubMenuX->addAction("Deformation");
    aSetX2->setCheckable(true);
    aSetX2->setChecked(true);
    aSetX2->setActionGroup(groupSubMenuX);

    QMenu *mySubMenuY = new QMenu("Set Y");
    QActionGroup* groupSubMenuY = new QActionGroup( this );
    QAction *aSetY1 = mySubMenuY->addAction("Time");
    aSetY1->setCheckable(true);
    aSetY1->setActionGroup(groupSubMenuY);
    QAction *aSetY2 = mySubMenuY->addAction("Y");
    aSetY2->setCheckable(true);
    aSetY2->setChecked(true);
    aSetY2->setActionGroup(groupSubMenuY);

    menuPlot = new QMenu ();
    menuPlot->setFont(MenuBold);
    menuPlot->addMenu(mySubMenuX);
    menuPlot->addMenu(mySubMenuY);

    connect (mySubMenuX, SIGNAL(triggered(QAction*)), this, SLOT(XTypeChanged (QAction*)));
    connect (mySubMenuY, SIGNAL(triggered(QAction*)), this, SLOT(YTypeChanged (QAction*)));

//    connect (aSetX2, SIGNAL(triggered()), this, SLOT(X2Type()));


    //END MENU


    canvas_picker = new CanvasPicker( this );

    //VERIFICAR O PQ NAO FUNCIONA
//    picker = new QwtPlotPicker(QwtPlot::xBottom, QwtPlot::yLeft,
//                   QwtPlotPicker::CrossRubberBand, QwtPicker::AlwaysOn,
//                   this->canvas());
//    picker->setStateMachine(new QwtPickerDragPointMachine());
//    picker->setRubberBandPen(QColor(Qt::green));
//    picker->setRubberBand(QwtPicker::CrossRubberBand);
//    picker->setTrackerPen(QColor(Qt::blue));


    panner = new QwtPlotPanner(this->canvas());
    panner->setMouseButton(Qt::MidButton);

    CurvesList = new QHash <int, CURVE > ;

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

void Plot::createCurve(int pos, QVector <double> *X, QVector <double> *Y, int checkStatus)
{
    // VERIFICANDO SE A CURVA JA EXISTE PARA NAO RECRIA-LA
    if (this->CurvesList->value(pos).curve_ptr) return;

    QwtPlotCurve *d_curve;
    d_curve = new QwtPlotCurve();
    d_curve->setRenderHint( QwtPlotItem::RenderAntialiased );
    d_curve->setPen( QPen( Qt::red ) );
    d_curve->setLegendAttribute( QwtPlotCurve::LegendShowLine );
    d_curve->setYAxis( QwtPlot::yLeft );
    d_curve->attach( this );
    d_curve->setSamples(*X, *Y);

    this->zoomer->setZoomBase();

    double Xsmallest=numeric_limits<double>::max(), Xbiggest=numeric_limits<double>::min(),
           Ysmallest=numeric_limits<double>::max(), Ybiggest=numeric_limits<double>::min();

    for (int i = 0; i < X->size(); i++) {
        // scale
        if (X->value(i)<Xsmallest) Xsmallest = X->value(i);
        if (X->value(i)>Xbiggest)  Xbiggest = X->value(i);
        if (Y->value(i)<Ysmallest) Ysmallest = Y->value(i);
        if (Y->value(i)>Ybiggest)  Ybiggest = Y->value(i);
    }

    QwtPlotCurve *d_curve_symbol;
    d_curve_symbol = new QwtPlotCurve( );
    d_curve_symbol->setRenderHint( QwtPlotItem::RenderAntialiased );
    d_curve_symbol->setPen( QPen( Qt::blue ) );
    d_curve_symbol->setLegendAttribute( QwtPlotCurve::LegendShowLine );
    d_curve_symbol->setYAxis( QwtPlot::yLeft );
    d_curve_symbol->setStyle(QwtPlotCurve::NoCurve);
    d_curve_symbol->attach( this );

    int sizes = 2;
    double Xs[2], Ys[2];
    Xs[0] = X->value(0);
    Ys[0] = Y->value(0);
    Xs[1] = X->value(X->size()-1);
    Ys[1] = Y->value(Y->size()-1);
    qDebug() << "Size: " << X->size() << " start: " << X->value(0) << " " << Y->value(0) << endl;
    qDebug() << "Size: " << X->size() << " end: " << X->value(X->size()-1) << " " << Y->value(X->size()-1) ;
    d_curve_symbol->setSamples( Xs, Ys, sizes);
    d_curve_symbol->hide();

    // "Start" marker
    QwtSymbol *symb1 = new QwtSymbol();
    symb1->setBrush(QBrush(Qt::red, Qt::SolidPattern));
    symb1->setStyle(QwtSymbol::Ellipse);
    symb1->setSize(10);
    QwtPlotMarker *a = new QwtPlotMarker ();
    a->setSymbol(symb1);
    a->setXValue(Xs[0]);
    a->setYValue(Ys[0]);
//    a->setLabelOrientation(Qt::Orientation(Qt::Vertical));
    a->setLabelAlignment(Qt::Alignment(Qt::AlignLeft));
    a->setLabel(QwtText("Start"));
    a->attach(this);
    a->hide();

    // End marker
    QwtSymbol *symb2 = new QwtSymbol();
    symb2->setBrush(QBrush(Qt::blue, Qt::SolidPattern));
    symb2->setStyle(QwtSymbol::Ellipse);
    symb2->setSize(10);
    QwtPlotMarker *b = new QwtPlotMarker ();
    b->setSymbol(symb2);
    b->setXValue(Ys[1]);
    b->setYValue(Ys[1]);
    b->setLabelAlignment(Qt::Alignment(Qt::AlignLeft));
    b->setLabel(QwtText("End"));
    b->attach(this);
    b->hide();

    CURVE curva;
    curva.chk_status = checkStatus;
    curva.curve_ptr = d_curve;
    curva.Xbiggest = Xbiggest;
    curva.Xsmallest = Xsmallest;
    curva.Ybiggest = Ybiggest;
    curva.Ysmallest = Ysmallest;
    curva.startIndexCoord = 0;
    curva.endIndexCoord = X->size()-1;
    curva.symb_curve_ptr = d_curve_symbol;
    curva.startM = a;
    curva.endM = b;
    curva.X = X;
    curva.Y = Y;
    this->CurvesList->insert(pos, curva);

    this->replot();
    this->AdjustScale();
}

// Get index of coords vector of a given curve
int Plot::getSymbIndex ( int indexCurves, pointType typept ) {
    QwtPlotCurve *d_curve_symbol = this->CurvesList->value(indexCurves).symb_curve_ptr;

    //Check if curve exists
    if (!d_curve_symbol) return -1;

    int returnIdx = -1;

    if (typept == startPoint)
        returnIdx = this->CurvesList->value(indexCurves).startIndexCoord;
    else
        returnIdx = this->CurvesList->value(indexCurves).endIndexCoord;

    qDebug() << "getSymbIndex: Pt idx: " << returnIdx;

    return returnIdx;
}

//Remove curve
void Plot::deleteCurve ( int pos ) {
    QwtPlotCurve *d_curve = this->CurvesList->value(pos).curve_ptr;
    if ( d_curve != NULL) {
        d_curve->detach();
        delete d_curve;
    }
    QwtPlotCurve *d_curve_symbol = this->CurvesList->value(pos).symb_curve_ptr;
    if ( d_curve_symbol != NULL ) {
        //removing markers
        QwtPlotMarker *d_plot_marker = NULL;
        d_plot_marker = this->CurvesList->value(pos).startM;
        d_plot_marker->detach();
        delete d_plot_marker;
        d_plot_marker = this->CurvesList->value(pos).endM;
        d_plot_marker->detach();
        delete d_plot_marker;
        //removing curve
        d_curve_symbol->detach();
        delete d_curve_symbol;
    }
    this->CurvesList->remove(pos);

    this->replot();
    this->AdjustScale();
}

void Plot::AdjustScale ()
{
    double Xsmallest=numeric_limits<double>::max(), Xbiggest=numeric_limits<double>::min(),
           Ysmallest=numeric_limits<double>::max(), Ybiggest=numeric_limits<double>::min();

    foreach (int i, this->CurvesList->keys()) {
        if (this->CurvesList->value(i).Xsmallest<Xsmallest)
            Xsmallest = this->CurvesList->value(i).Xsmallest;

        if (this->CurvesList->value(i).Xbiggest>Xbiggest)
            Xbiggest = this->CurvesList->value(i).Xbiggest;

        if (this->CurvesList->value(i).Ysmallest<Ysmallest)
            Ysmallest = this->CurvesList->value(i).Ysmallest;

        if (this->CurvesList->value(i).Ybiggest>Ybiggest)
            Ybiggest = this->CurvesList->value(i).Ybiggest;
    }

    if (this->CurvesList->size() == 0)
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
    QwtPlotCurve *d_curve = this->CurvesList->value(indexCurves).curve_ptr;

    //Check if curve exists
    if (!d_curve) return;

    CURVE d_curve_tmp = this->CurvesList->take(indexCurves);

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
    this->CurvesList->insert(indexCurves, d_curve_tmp);

    qDebug() << "setSymbIndex: Pt idx: " << indexCoords << " type " << typept;

    this->replot();
    //this->AdjustScale();
}

//slot
void Plot::showSymbCurve(int indexCurves) {
    this->CurvesList->value(indexCurves).startM->show();
    this->CurvesList->value(indexCurves).endM->show();
    this->CurvesList->value(indexCurves).symb_curve_ptr->show();
    this->replot();
}

//slot
void Plot::hideSymbCurve(int indexCurves) {
    this->CurvesList->value(indexCurves).startM->hide();
    this->CurvesList->value(indexCurves).endM->hide();
    this->CurvesList->value(indexCurves).symb_curve_ptr->hide();
    this->replot();
}

//slot
void Plot::cutCurve(int indexCurves, int indexStartPoint, int indexEndPoint) {
    QwtPlotCurve *d_curve = this->CurvesList->value(indexCurves).curve_ptr;

    //Check if curve exists
    if (!d_curve) return;

    CURVE d_curve_tmp = this->CurvesList->take(indexCurves);

    int size_tmp = d_curve_tmp.X->size();
    qDebug() << "Removendo de " << indexEndPoint +1 << " count " << (size_tmp - 1) - indexEndPoint ;
    qDebug() << "Removendo de 0 count " << indexStartPoint ;
    d_curve_tmp.X->remove(indexEndPoint+1,(size_tmp - 1) - indexEndPoint);
    d_curve_tmp.Y->remove(indexEndPoint+1,(size_tmp - 1) - indexEndPoint);
    d_curve_tmp.X->remove(0,indexStartPoint);
    d_curve_tmp.Y->remove(0,indexStartPoint);


    d_curve_tmp.startIndexCoord = 0;
    d_curve_tmp.endIndexCoord = d_curve_tmp.X->size()-1;

    qDebug() << "New points (cut): " << d_curve_tmp.startIndexCoord << " " << d_curve_tmp.endIndexCoord;

    //re-inserting new value
    this->CurvesList->insert(indexCurves, d_curve_tmp);

    // update curve with new data
    this->updateCurve(indexCurves, d_curve_tmp.startIndexCoord, d_curve_tmp.endIndexCoord);
}

//slot
void Plot::updateCurve(int indexCurves, int indexStartPoint, int indexEndPoint) {
    QwtPlotCurve *d_curve = this->CurvesList->value(indexCurves).curve_ptr;

    //Check if curve exists
    if (!d_curve) return;

    CURVE d_curve_tmp = this->CurvesList->take(indexCurves);

    d_curve->setSamples (*d_curve_tmp.X, *d_curve_tmp.Y);

    // Reposiciona marcadores no inicio/fim da curva cortada
    // Parametros indexStartPoint e indexEndPoint estao com valores antigos.
    d_curve_tmp.startIndexCoord = 0;
    d_curve_tmp.endIndexCoord = d_curve_tmp.X->size()-1;

    this->setSymbIndex(d_curve_tmp.startIndexCoord,indexCurves,startPoint);
    this->setSymbIndex(d_curve_tmp.endIndexCoord,indexCurves,endPoint);

    qDebug() << "New points (update): " << d_curve_tmp.startIndexCoord << " " << d_curve_tmp.endIndexCoord;

    //re-adjusting scale
    double Xsmallest=numeric_limits<double>::max(), Xbiggest=numeric_limits<double>::min(),
           Ysmallest=numeric_limits<double>::max(), Ybiggest=numeric_limits<double>::min();
    for (int i = 0; i < d_curve_tmp.X->size() ; i++) {
        // scale
        if (d_curve_tmp.X->value(i) < Xsmallest) Xsmallest = d_curve_tmp.X->value(i);
        if (d_curve_tmp.X->value(i) > Xbiggest)  Xbiggest = d_curve_tmp.X->value(i);
        if (d_curve_tmp.Y->value(i) < Ysmallest) Ysmallest = d_curve_tmp.Y->value(i);
        if (d_curve_tmp.Y->value(i) > Ybiggest)  Ybiggest = d_curve_tmp.Y->value(i);
    }

    d_curve_tmp.Xbiggest=Xbiggest;
    d_curve_tmp.Xsmallest=Xsmallest;
    d_curve_tmp.Ybiggest=Ybiggest;
    d_curve_tmp.Ysmallest=Ysmallest;

    //re-inserting new value
    this->CurvesList->insert(indexCurves, d_curve_tmp);

    this->replot();
    this->AdjustScale();
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

void Plot::XTypeChanged (QAction *action) {

    qDebug() << "X TypeChanged " << action->text()  ;

    if (action->text() == "Tension") {
        qDebug() << "usar coluna tensao";
    }


}

void Plot::YTypeChanged (QAction *action) {

    qDebug() << "Y TypeChanged " << action->text() ;
}
