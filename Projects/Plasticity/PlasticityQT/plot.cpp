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

    envelope_curve = new QwtPlotCurve();
    envelope_curve->setRenderHint( QwtPlotItem::RenderAntialiased );
    envelope_curve->setPen( QPen( Qt::black ) );
    envelope_curve->setLegendAttribute( QwtPlotCurve::LegendShowLine );
    envelope_curve->setYAxis( QwtPlot::yLeft );

    Xenv_max=0;
    Xenv_min=0;
    Yenv_max=0;
    Yenv_min=0;
    Env_status = -1;

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
    // sqrt J2 x I1
    QString Coords5 = QChar (0x221A);
    Coords5.append("J2 x I1");

    QMenu *mySubMenuAxis = new QMenu ("Set Axis");
    QActionGroup *groupSubmenuAxis = new QActionGroup (this);
    QAction *aSetAxis1 = mySubMenuAxis->addAction(Coords1);
    aSetAxis1->setCheckable(true);
    aSetAxis1->setChecked(true);
    aSetAxis1->setActionGroup(groupSubmenuAxis);
    aSetAxis1->setData(TPBrStrainStressDataBase::EEpsaxSigax);
    QAction *aSetAxis2 = mySubMenuAxis->addAction(Coords2);
    aSetAxis2->setCheckable(true);
    aSetAxis2->setActionGroup(groupSubmenuAxis);
    aSetAxis2->setData(TPBrStrainStressDataBase::EEpsrSigr);
    QAction *aSetAxis3 = mySubMenuAxis->addAction(Coords3);
    aSetAxis3->setCheckable(true);
    aSetAxis3->setActionGroup(groupSubmenuAxis);
    aSetAxis3->setData(TPBrStrainStressDataBase::EEpsvSigv);
    QAction *aSetAxis4 = mySubMenuAxis->addAction(Coords4);
    aSetAxis4->setCheckable(true);
    aSetAxis4->setActionGroup(groupSubmenuAxis);
    aSetAxis4->setData(TPBrStrainStressDataBase::EEpsaxEpsrEpsvSigax);
    QAction *aSetAxis5 = mySubMenuAxis->addAction(Coords5);
    aSetAxis5->setCheckable(true);
    aSetAxis5->setActionGroup(groupSubmenuAxis);
    aSetAxis5->setData(TPBrStrainStressDataBase::EI1SqJ2);

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

    // scale divisions (Majors/Minors) and Engines
    this->setAxisMaxMajor( QwtPlot::xBottom, 5 );
    this->setAxisMaxMinor( QwtPlot::xBottom, 10 );
    this->setAxisMaxMajor( QwtPlot::yLeft, 5 );
    this->setAxisMaxMinor( QwtPlot::yLeft, 10 );
    this->setAxisScaleEngine( QwtPlot::xBottom, new QwtLinearScaleEngine );
    this->setAxisScaleEngine( QwtPlot::yLeft, new QwtLinearScaleEngine );

    this->setTitle("");
    this->setHighlighted(false);

    this->curvetype = TPBrStrainStressDataBase::EEpsaxSigax;
    QString X_title;
    QString Y_title;
    Y_title.append(QChar (0x03C3));
    Y_title.append("ax (MPa)");
    this->setAxisTitle(QwtPlot::yLeft, Y_title);
    X_title.append(QChar (0x03B5));
    X_title.append("ax (%)");
    this->setAxisTitle(QwtPlot::xBottom, X_title);
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


//Click event = MousePress->MouseRelease
void Plot::mouseReleaseEvent( QMouseEvent * event ) {
    if(event->button() == Qt::RightButton) {
        QPoint globalPos = this->mapToGlobal(event->pos()); // for most widgets
        menuPlot->exec(globalPos);
    }
}

void Plot::AxisChanged_slot (QAction *action) {

    this->curvetype = TPBrStrainStressDataBase::ECurveType (action->data().toInt());

    if(this->curvetype == TPBrStrainStressDataBase::EI1SqJ2)
    {
      this->envelope_curve->show();
      Env_status = 1;
    }
    else
    {
      this->envelope_curve->hide();
      Env_status = 0;
    }
    // plotting new curve
    foreach (int i, this->CurvesList.keys())
    {
        int chkStatus = this->CurvesList[i].chk_status;
        deleteCurve(i);
        createCurve(i, chkStatus);
    }

    // renaming axis
    QString X_title;
    QString Y_title;
    switch (curvetype)
    {
        case TPBrStrainStressDataBase::EEpsaxSigax:
        {
            Y_title.append(QChar (0x03C3));
            Y_title.append("ax (MPa)");
            this->setAxisTitle(QwtPlot::yLeft, Y_title);
            X_title.append(QChar (0x03B5));
            X_title.append("ax (%)");
            this->setAxisTitle(QwtPlot::xBottom, X_title);
            break;
        }

        case TPBrStrainStressDataBase::EEpsrSigr:
        {
            Y_title.append(QChar (0x03C3));
            Y_title.append("rad (MPa)");
            this->setAxisTitle(QwtPlot::yLeft, Y_title);
            X_title.append(QChar (0x03B5));
            X_title.append("rad (%)");
            this->setAxisTitle(QwtPlot::xBottom, X_title);
            break;
        }

        case TPBrStrainStressDataBase::EEpsvSigv:
        {
            Y_title.append(QChar (0x03C3));
            Y_title.append("vol (MPa)");
            this->setAxisTitle(QwtPlot::yLeft, Y_title);
            X_title.append(QChar (0x03B5));
            X_title.append("vol (%)");
            this->setAxisTitle(QwtPlot::xBottom, X_title);
            break;
        }

        case TPBrStrainStressDataBase::EEpsaxEpsrEpsvSigax:
        {
            Y_title.append(QChar (0x03C3));
            Y_title.append("ax (MPa)");
            this->setAxisTitle(QwtPlot::yLeft, Y_title);
            X_title.append(QChar (0x03B5));
            X_title.append("ax ");
            X_title.append(QChar (0x03B5));
            X_title.append("vol ");
            X_title.append(QChar (0x03B5));
            X_title.append("rad (%)");
            this->setAxisTitle(QwtPlot::xBottom, X_title);
            break;
        }

        case TPBrStrainStressDataBase::EI1SqJ2:
        {
            Y_title.append(QChar (0x221A));
            Y_title.append("J2");
            this->setAxisTitle(QwtPlot::yLeft, Y_title);
            X_title.append("I1");
            this->setAxisTitle(QwtPlot::xBottom, X_title);
            break;
        }
    }
    AdjustScale();

}

//slot
void Plot::setSymbIndex ( int indexCoords, int global_id, pointType typept ) {

    //Verificando curva existe
    if (!this->CurvesList.contains(global_id)) return;

    QwtPlotCurve *d_curve = this->CurvesList[global_id].curve_ptr;

    //Check if curve exists
    if (!d_curve) return;

    CURVE d_curve_tmp = this->CurvesList[global_id];

    //updating symbol position
    int sizes = 3;
    double Xs[3], Ys[3];

    //start
    Xs[0] = d_curve_tmp.symb_curve_ptr->sample(0).x();
    Ys[0] = d_curve_tmp.symb_curve_ptr->sample(0).y();
    //elastic trans
    Xs[1] = d_curve_tmp.symb_curve_ptr->sample(1).x();
    Ys[1] = d_curve_tmp.symb_curve_ptr->sample(1).y();
    //end
    Xs[2] = d_curve_tmp.symb_curve_ptr->sample(2).x();
    Ys[2] = d_curve_tmp.symb_curve_ptr->sample(2).y();

    d_curve_tmp.mark_ptr1->setXValue(Xs[0]);
    d_curve_tmp.mark_ptr1->setYValue(Ys[0]);
    d_curve_tmp.mark_ptr2->setXValue(Xs[2]);
    d_curve_tmp.mark_ptr2->setYValue(Ys[2]);
    d_curve_tmp.mark_ptr3->setXValue(Xs[1]);
    d_curve_tmp.mark_ptr3->setYValue(Ys[1]);

    switch (typept) {
    case startPoint:
        //updating curve
        Xs[0] = d_curve_tmp.curve_ptr->sample(indexCoords).x();
        Ys[0] = d_curve_tmp.curve_ptr->sample(indexCoords).y();
        d_curve_tmp.symb_curve_ptr->setSamples( Xs, Ys, sizes);
        //updading markers (labels)
        d_curve_tmp.mark_ptr1->setXValue(Xs[0]);
        d_curve_tmp.mark_ptr1->setYValue(Ys[0]);
        //updating indexCoords control
        d_curve_tmp.initialPnt = indexCoords;
        break;
    case endPoint:
        //updating curve
        Xs[2] = d_curve_tmp.curve_ptr->sample(indexCoords).x();
        Ys[2] = d_curve_tmp.curve_ptr->sample(indexCoords).y();
        d_curve_tmp.symb_curve_ptr->setSamples( Xs, Ys, sizes);
        //updading markers (labels)
        d_curve_tmp.mark_ptr2->setXValue(Xs[2]);
        d_curve_tmp.mark_ptr2->setYValue(Ys[2]);
        //updating indexCoords control
        d_curve_tmp.endPnt = indexCoords;
        break;
    case elasticTransPoint:
        //updating curve
        Xs[1] = d_curve_tmp.curve_ptr->sample(indexCoords).x();
        Ys[1] = d_curve_tmp.curve_ptr->sample(indexCoords).y();
        d_curve_tmp.symb_curve_ptr->setSamples( Xs, Ys, sizes);
        //updading markers (labels)
        d_curve_tmp.mark_ptr3->setXValue(Xs[1]);
        d_curve_tmp.mark_ptr3->setYValue(Ys[1]);
        //updating indexCoords control
        d_curve_tmp.elastic_transPnt = indexCoords;
        break;
    }

    //inserting new value
    this->CurvesList.insert(global_id, d_curve_tmp);

    qDebug() << "setSymbIndex: Pt idx: " << indexCoords << " type " << typept;

//    this->replot();
    this->AdjustScale();
}

void Plot::createCurve (int global_id, int check_status){

    // VERIFICANDO SE A CURVA JA EXISTE PARA NAO RECRIA-LA
    if (this->CurvesList.contains(global_id)) return;

    CURVE new_curve(global_id);

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
    int elastic_trans_idx = labData->Get_elastic_trans_idx();
    int isMed = DADOS.isMed(global_id);

    QVector <double> qVec_X = QVector <double>::fromStdVector(X);
    QVector <double> qVec_Y = QVector <double>::fromStdVector(Y);
    new_curve.curve_ptr->setSamples(qVec_X,qVec_Y);

    new_curve.X1 = qVec_X;
    new_curve.Y1 = qVec_Y;

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

    this->zoomer->setZoomBase();

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

    int sizes = 3;
    double Xs[3], Ys[3];
    Xs[0] = X[start_idx];
    Ys[0] = Y[start_idx];
    Xs[2] = X[end_idx];
    Ys[2] = Y[end_idx];
    Xs[1] = X[elastic_trans_idx];
    Ys[1] = Y[elastic_trans_idx];
    qDebug() << "Size: " << X.size() << " start: " << start_idx << " " << X[start_idx] << " " << Y[start_idx];
    qDebug() << "Size: " << X.size() << " end: " << end_idx << " " << X[end_idx] << " " << Y[end_idx];
    qDebug() << "Size: " << X.size() << " elastic trans: " << elastic_trans_idx << " " << X[elastic_trans_idx] << " " << Y[elastic_trans_idx];
    new_curve.symb_curve_ptr->setSamples( Xs, Ys, sizes);
    if (!isMed) new_curve.symb_curve_ptr->hide();

    // "Start" marker
    new_curve.symb1 = new QwtSymbol();
    new_curve.symb1->setBrush(QBrush(Qt::red, Qt::SolidPattern));
    new_curve.symb1->setStyle(QwtSymbol::Ellipse);
    new_curve.symb1->setSize(7);//(10);
    new_curve.mark_ptr1 = new QwtPlotMarker();
    new_curve.mark_ptr1->setSymbol(new_curve.symb1);
    new_curve.mark_ptr1->setXValue(Xs[0]);
    new_curve.mark_ptr1->setYValue(Ys[0]);
    new_curve.mark_ptr1->setLabelAlignment(Qt::Alignment(Qt::AlignLeft));
    //new_curve.mark_ptr1->setLabel(QwtText("Start"));
    new_curve.mark_ptr1->attach(this);
    if (!isMed) new_curve.mark_ptr1->hide();

    // End marker
    new_curve.symb2 = new QwtSymbol();
    new_curve.symb2->setBrush(QBrush(Qt::blue, Qt::SolidPattern));
    new_curve.symb2->setStyle(QwtSymbol::Ellipse);
    new_curve.symb2->setSize(7);//(10);
    new_curve.mark_ptr2 = new QwtPlotMarker ();
    new_curve.mark_ptr2->setSymbol(new_curve.symb2);
    new_curve.mark_ptr2->setXValue(Xs[2]);
    new_curve.mark_ptr2->setYValue(Ys[2]);
    new_curve.mark_ptr2->setLabelAlignment(Qt::Alignment(Qt::AlignLeft));
    //new_curve.mark_ptr2->setLabel(QwtText("End"));
    new_curve.mark_ptr2->attach(this);
    if (!isMed) new_curve.mark_ptr2->hide();

    // ElasticTrans marker
    new_curve.symb3 = new QwtSymbol();
    new_curve.symb3->setBrush(QBrush(Qt::green, Qt::SolidPattern));
    new_curve.symb3->setStyle(QwtSymbol::Ellipse);
    new_curve.symb3->setSize(7);//(10);
    new_curve.mark_ptr3 = new QwtPlotMarker ();
    new_curve.mark_ptr3->setSymbol(new_curve.symb3);
    new_curve.mark_ptr3->setXValue(Xs[1]);
    new_curve.mark_ptr3->setYValue(Ys[1]);
    new_curve.mark_ptr3->setLabelAlignment(Qt::Alignment(Qt::AlignLeft));
    //new_curve.mark_ptr3->setLabel(QwtText("Elast"));
    new_curve.mark_ptr3->attach(this);
    if (!isMed) new_curve.mark_ptr3->hide();


    new_curve.initialPnt = start_idx;
    new_curve.endPnt = end_idx;
    new_curve.elastic_transPnt = elastic_trans_idx;

    new_curve.chk_status = check_status;

    this->CurvesList.insert(global_id, new_curve);

    //naming axis
    if (curvetype == TPBrStrainStressDataBase::EEpsaxSigax){
        QString DefaultTitleY = QChar (0x03C3);
        DefaultTitleY.append("ax (MPa)");
        this->setAxisTitle(QwtPlot::yLeft, DefaultTitleY);
        QString DefaultTitle1X = QChar (0x03B5);
        DefaultTitle1X.append("ax (%)");
        this->setAxisTitle(QwtPlot::xBottom, DefaultTitle1X);
    }

    //this->replot();
    this->AdjustScale();
}

void Plot::deleteCurve (int global_id){

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
        d_plot_marker = this->CurvesList[global_id].mark_ptr3;
        d_plot_marker->detach();
        delete d_plot_marker;
        //removing curve
        d_curve_symbol->detach();
        delete d_curve_symbol;
    }
    this->CurvesList.remove(global_id);

    //this->replot();
    this->AdjustScale();

}

void Plot::Generate_Envelope(double A, double B, double C){

    std::vector<REAL> X;
    std::vector<REAL> Y;

    TPBrStrainStressDataBase *labdata = DADOS.getObj(0);
    labdata->SetEnvelope(A, B, C);

    double Xsmallest=numeric_limits<double>::max(), Xbiggest=numeric_limits<double>::min();

    foreach (int i, this->CurvesList.keys()) {

        if (this->CurvesList[i].Xsmallest<Xsmallest)
            Xsmallest = this->CurvesList.value(i).Xsmallest;

        if (this->CurvesList[i].Xbiggest>Xbiggest)
            Xbiggest = this->CurvesList[i].Xbiggest;
    }

    labdata->GenerateEnvelope(X, Y, Xsmallest, Xbiggest);

    QVector <double> qVec_X = QVector <double>::fromStdVector(X);
    QVector <double> qVec_Y = QVector <double>::fromStdVector(Y);

    envelope_curve->setSamples(qVec_X, qVec_Y);
    envelope_curve->attach(this);
    this->zoomer->setZoomBase();
    //this->replot();
    Xsmallest=numeric_limits<double>::max(), Xbiggest=numeric_limits<double>::min();
    double Ysmallest=numeric_limits<double>::max(), Ybiggest=numeric_limits<double>::min();

    for (int i = 0; i < X.size(); i++) {
        // scale
        if (X[i]<Xsmallest) Xsmallest = X[i];
        if (X[i]>Xbiggest)  Xbiggest = X[i];
        if (Y[i]<Ysmallest) Ysmallest = Y[i];
        if (Y[i]>Ybiggest) Ybiggest = Y[i];
    }

    Xenv_max=Xbiggest;
    Xenv_min=Xsmallest;
    Yenv_max=Ybiggest;
    Yenv_min=Ysmallest;

    if(this->curvetype == TPBrStrainStressDataBase::EI1SqJ2)
    {
      this->envelope_curve->show();
      Env_status = 1;
    }
    else
    {
      this->envelope_curve->hide();
      Env_status = 0;

    }
    
    this->AdjustScale();
}

void Plot::hide_envelope()
{
    envelope_curve->hide();
    Env_status = 0;
    this->AdjustScale();
}

void Plot::AdjustScale () {

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

    if (Env_status==1)
    {
        if (Xenv_min<Xsmallest)
            Xsmallest = Xenv_min;
        if (Xenv_max>Xbiggest)
            Xbiggest = Xenv_max;
        if (Yenv_min<Ysmallest)
            Ysmallest = Yenv_min;
        if (Yenv_max>Ybiggest)
            Ybiggest = Yenv_max;
    }

    if ((this->CurvesList.size() == 0) && (Env_status != 1))
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

void Plot::toggleSymb (int global_id, bool showsymb) {

//    //Verificando curva existe
//    if (!this->CurvesList.contains(global_id)) return;

//    CURVE d_curve_tmp = this->CurvesList[global_id];

//    if (showsymb) {
//        d_curve_tmp.mark_ptr1->show();
//        d_curve_tmp.mark_ptr2->show();
//    }
//    else
//    {
//        d_curve_tmp.mark_ptr1->hide();
//        d_curve_tmp.mark_ptr2->hide();
//    }

    DebugStop();

}
