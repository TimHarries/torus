#include "psedplotter.h"
#include "ui_psedplotter.h"

PSEDPlotter::PSEDPlotter(QString dir, QString filePrefix, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::PSEDPlotter)
{
    ui->setupUi(this);
    this->setWorkingDirectory(dir);
    this->setFilePrefix(filePrefix);

    //Scan for SEDs in the current directory
    this->refreshSEDFiles();

    //Setup the Graph
    ui->grp_SEDGraph->setAxisTitle(2, "Wavelength (microns)");
    ui->grp_SEDGraph->setAxisTitle(0, "Flux (W/m^2)");
    ui->grp_SEDGraph->setAxisAutoScale(0, true);
    ui->grp_SEDGraph->setAxisAutoScale(2, true);
    ui->grp_SEDGraph->setAxisScaleEngine(QwtPlot::xBottom, new QwtLogScaleEngine());
    ui->grp_SEDGraph->setAxisScaleEngine(QwtPlot::yLeft, new QwtLogScaleEngine());
    ui->grp_SEDGraph->setAutoDelete(false);

    //Setup Legend
    this->legend = new QwtLegend(this);
    ui->grp_SEDGraph->insertLegend(this->legend, QwtPlot::TopLegend);
    this->legend->show();

    //Init our Curves
    this->InitCurves();

    //Connect the Signals of the Checkboxes to Trigger the Graph Replot
    connect(ui->chk_Combined, SIGNAL(clicked()), this, SLOT(checkboxUpdateGraph()));
    connect(ui->chk_StelDir, SIGNAL(clicked()), this, SLOT(checkboxUpdateGraph()));
    connect(ui->chk_StelScat, SIGNAL(clicked()), this, SLOT(checkboxUpdateGraph()));
    connect(ui->chk_ThermDir, SIGNAL(clicked()), this, SLOT(checkboxUpdateGraph()));
    connect(ui->chk_ThermScat, SIGNAL(clicked()), this, SLOT(checkboxUpdateGraph()));

}

PSEDPlotter::~PSEDPlotter()
{
    delete ui;
}

void PSEDPlotter::setWorkingDirectory(QString dir) {
    this->workingDir = dir;
}

void PSEDPlotter::setFilePrefix(QString prefix) {
    this->namePrefix = prefix;
}

void PSEDPlotter::refreshSEDFiles() {
    //First Clear Our Listbox
    ui->lst_SEDList->clear();

    QDir wrkDir(this->workingDir);
    QStringList fileNameTemplates;
    QStringList listEntries;
    QString filePattern = "_stellar_direct.dat|_stellar_scattered.dat|_thermal_direct.dat|_thermal_scattered.dat|.dat";
    fileNameTemplates << this->namePrefix + "_inc*.dat";

    wrkDir.setNameFilters(fileNameTemplates);
    foreach(QString file, wrkDir.entryList()) {
        if(!listEntries.contains(file.remove(QRegExp(filePattern)))) {
            listEntries << file.remove(QRegExp(filePattern));
        }
    }

    //Add to listbox
    ui->lst_SEDList->addItems(listEntries);
}

void PSEDPlotter::closeEvent(QCloseEvent *e) {
    Q_UNUSED(e);
    emit this->isGoingToClose();
}

void PSEDPlotter::InitCurves() {
    //Allocate the Curves
    this->combinedCurve = new QwtPlotCurve();
    this->stellarDirectCurve = new QwtPlotCurve();
    this->stellarScatteredCurve = new QwtPlotCurve();
    this->thermalDirectCurve = new QwtPlotCurve();
    this->thermalScatteredCurve = new QwtPlotCurve();

    //Combined Curve Setup
    combinedCurve->setPen(Qt::black, 0, Qt::DashLine);
    combinedCurve->setSymbol(new QwtSymbol(QwtSymbol::Cross, QBrush(Qt::black, Qt::SolidPattern), QPen(Qt::black, 0, Qt::DashLine), QSize(2, 2)));
    combinedCurve->setTitle("Combined");

    //Stellar Direct Curve Setup
    stellarDirectCurve->setPen(Qt::blue, 0, Qt::SolidLine);
    stellarDirectCurve->setSymbol(new QwtSymbol(QwtSymbol::Cross, QBrush(Qt::blue, Qt::SolidPattern), QPen(Qt::blue, 0, Qt::DashLine), QSize(2, 2)));
    stellarDirectCurve->setTitle("Stellar Direct");

    //Stellar Scattered Curve Setup
    stellarScatteredCurve->setPen(Qt::green, 0, Qt::SolidLine);
    stellarScatteredCurve->setSymbol(new QwtSymbol(QwtSymbol::Cross, QBrush(Qt::green, Qt::SolidPattern), QPen(Qt::green, 0, Qt::DashLine), QSize(2, 2)));
    stellarScatteredCurve->setTitle("Stellar Scattered");

    //Thermal Direct Curve Setup
    thermalDirectCurve->setPen(Qt::red, 0, Qt::SolidLine);
    thermalDirectCurve->setSymbol(new QwtSymbol(QwtSymbol::Cross, QBrush(Qt::red, Qt::SolidPattern), QPen(Qt::red, 0, Qt::DashLine), QSize(2, 2)));
    thermalDirectCurve->setTitle("Thermal Direct");

    //Thermal Scattered Curve Curve Setup
    thermalScatteredCurve->setPen(QColor(255, 123, 0), 0, Qt::SolidLine);
    thermalScatteredCurve->setSymbol(new QwtSymbol(QwtSymbol::Cross, QBrush(QColor(255, 123, 0), Qt::SolidPattern), QPen(QColor(255, 123, 0), 0, Qt::DashLine), QSize(2, 2)));
    thermalScatteredCurve->setTitle("Thermal Scattered");

}

void PSEDPlotter::on_lst_SEDList_itemSelectionChanged() {
    QString file = this->workingDir + "/" + ui->lst_SEDList->currentItem()->text();
    this->PlotGraphs(file);
}

void PSEDPlotter::checkboxUpdateGraph() {
    QString file = this->workingDir + "/" + ui->lst_SEDList->currentItem()->text();
    this->PlotGraphs(file);
}

//Graph Plotting Function
void PSEDPlotter::PlotGraphs(QString file) {
    ui->grp_SEDGraph->detachItems(QwtPlotItem::Rtti_PlotCurve, false);

    if(ui->chk_Combined->isChecked()) {
        this->PlotSED(file, SEDCombinedPlot);
    }

    if(ui->chk_StelDir->isChecked()) {
        this->PlotSED(file, SEDStellarDirectPlot);
    }

    if(ui->chk_StelScat->isChecked()) {
        this->PlotSED(file, SEDStellarScatteredPlot);
    }

    if(ui->chk_ThermDir->isChecked()) {
        this->PlotSED(file, SEDThermalDirectPlot);
    }

    if(ui->chk_ThermScat->isChecked()) {
        this->PlotSED(file, SEDThermalScatteredPlot);
    }
}

//File Plotting Methods
void PSEDPlotter::PlotSED(QString prefixName, SEDPlotType type) {
    QPolygonF data;

    switch (type) {
        case SEDCombinedPlot:
                data = ParseSEDDataFile(prefixName + ".dat");
                combinedCurve->setSamples(data);
                combinedCurve->attach(ui->grp_SEDGraph);
            break;

        case SEDStellarDirectPlot:
            data = ParseSEDDataFile(prefixName + "_stellar_direct.dat");
            stellarDirectCurve->setSamples(data);
            stellarDirectCurve->attach(ui->grp_SEDGraph);
            break;

        case SEDStellarScatteredPlot:
            data = ParseSEDDataFile(prefixName + "_stellar_scattered.dat");
            stellarScatteredCurve->setSamples(data);
            stellarScatteredCurve->attach(ui->grp_SEDGraph);
            break;

        case SEDThermalDirectPlot:
            data = ParseSEDDataFile(prefixName + "_thermal_direct.dat");
            thermalDirectCurve->setSamples(data);
            thermalDirectCurve->attach(ui->grp_SEDGraph);
            break;

        case SEDThermalScatteredPlot:
            data = ParseSEDDataFile(prefixName + "_thermal_scattered.dat");
            thermalScatteredCurve->setSamples(data);
            thermalScatteredCurve->attach(ui->grp_SEDGraph);
            break;

        default:

            break;
    }
    ui->grp_SEDGraph->replot();
}

QPolygonF PSEDPlotter::ParseSEDDataFile(QString filename) {
    QFile file(filename);
    QPolygonF retval;
    QPointF lineValue;
    QMessageBox msgBox;

    qDebug() << "Loading File: " << filename;

    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return QPolygonF(0);

    while(file.atEnd() == false) {
        QByteArray line = file.readLine();
        lineValue = this->ParseSEDLine(line);

        if(!lineValue.isNull() ) {

            //Prevent the Scale Breaking

            if(lineValue.y() != 0)
                retval << lineValue;
        }
    }

    return retval;
}

QPointF PSEDPlotter::ParseSEDLine(QByteArray line) {
    SEDParseStates state = ParseStart;
    QString buffer;
    int columnNo = 0;
    QPointF returnPoint;

    for(int i = 0; i < line.count(); i++) {
        char character = line.at(i);

        switch(state) {
             case ParseStart :
                switch(character) {
                    case ' ' :
                        state = ParseExpectNumber;
                    break;

                    case '#' :
                        state = ParseReadComment;
                    break;

                    default :
                        buffer.push_back(character);
                        state = ParseReadNumber;
                    break;
                }
             break;

            case ParseExpectNumber :
                switch(character) {
                    case '#' :
                        state = ParseReadComment;
                    break;

                    case ' ' :
                    break;

                    default :
                        buffer.push_back(character);
                        state = ParseReadNumber;
                    break;
                }
            break;

            case ParseReadNumber :
            switch(character) {
                case '#' :
                    state = ParseReadComment;
                break;

                case ' ' :
                    if(columnNo < 2) {
                        if(columnNo == 0) {
                            returnPoint.setX((buffer.toDouble()));
                        } else if(columnNo == 1) {
                            returnPoint.setY((buffer.toDouble()));
                        }
                    }
                    buffer.clear();
                    columnNo++;
                    state = ParseExpectNumber;
                break;

                default :
                    buffer.push_back(character);
                    state = ParseReadNumber;
                break;
            }
            break;

            case ParseReadComment :
                //Do Nothing
                buffer.clear();
            break;

            default :
            break;
        }
    }

    return returnPoint;
}
