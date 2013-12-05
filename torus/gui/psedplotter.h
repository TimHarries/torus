#ifndef PSEDPLOTTER_H
#define PSEDPLOTTER_H

#include <QWidget>
#include <QDir>
#include <QMessageBox>
#include <QColor>
#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_grid.h>
#include <qwt_symbol.h>
#include <qwt_legend.h>
#include <qwt_date_scale_engine.h>

namespace Ui {
class PSEDPlotter;
}

class PSEDPlotter : public QWidget
{
    Q_OBJECT
    
public:

    enum SEDPlotType {
        SEDCombinedPlot,
        SEDStellarDirectPlot,
        SEDStellarScatteredPlot,
        SEDThermalDirectPlot,
        SEDThermalScatteredPlot
    };

    enum SEDParseStates {
        ParseStart,
        ParseExpectNumber,
        ParseReadNumber,
        ParseReadComment
    };

    explicit PSEDPlotter(QString dir, QString filePrefix = "", QWidget *parent = 0);
    ~PSEDPlotter();

    void setWorkingDirectory(QString dir);
    void setFilePrefix(QString prefix);
    void refreshSEDFiles();
    void PlotGraphs(QString file = "");

    //PlotType
    void PlotSED(QString prefixName = "", SEDPlotType type = SEDCombinedPlot);
    QPolygonF ParseSEDDataFile(QString filename);
    QPointF ParseSEDLine(QByteArray line);
    
private:
    Ui::PSEDPlotter *ui;
    void closeEvent(QCloseEvent *e);
    void InitCurves();

    QString workingDir;
    QString namePrefix;

    //Plot Stuff
    QwtLegend *legend;

    //Curves
    QwtPlotCurve *combinedCurve;
    QwtPlotCurve *stellarDirectCurve;
    QwtPlotCurve *stellarScatteredCurve;
    QwtPlotCurve *thermalDirectCurve;
    QwtPlotCurve *thermalScatteredCurve;

signals:
    void isGoingToClose();
private slots:
    void on_lst_SEDList_itemSelectionChanged();

public slots:
    void checkboxUpdateGraph();
};

#endif // PSEDPLOTTER_H
