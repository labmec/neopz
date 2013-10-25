#ifndef INITIALPOINTDOCK_H
#define INITIALPOINTDOCK_H

#include <QDockWidget>
#include "plot.h" //for point type enum

namespace Ui {
class initialpointdock;
}

class initialpointdock : public QDockWidget
{
    Q_OBJECT
    
public:
    explicit initialpointdock(QWidget *parent = 0);
    ~initialpointdock();
    void setIndexCurve (int idx);
    void setSlidersMaximum (int value);
    void setSlidersMinimum (int value);
    void setXdata (const QVector <double> &Xvec);
    void setYdata (const QVector <double> &Yvec);
    
private:
    Ui::initialpointdock *ui;
    int coordStartIndex;
    int coordEndIndex;
    int curveIndex;
    bool isCut;

    QVector <double> X;
    QVector <double> Y;

signals:
    void mouseDoubleClicked ( int );
    void SymbPointChanged ( int coordIndex, int curveIndex, Plot::pointType typept );
    void showSymbCurve ( int indexCurves );
    void hideSymbCurve ( int indexCurves );
    void cutCurve(int indexCurves, int indexStartPoint, int indexEndPoint);

public Q_SLOTS:
    void sliderStartValueChanged ( int value );
    void sliderEndValueChanged ( int value );
    void setSymbPoint (int coordIndex, int curveIndex, Plot::pointType typept );

private slots:
    void on_cutBtn_clicked();
    void closeEvent(QCloseEvent *event);
};

#endif // INITIALPOINTDOCK_H
