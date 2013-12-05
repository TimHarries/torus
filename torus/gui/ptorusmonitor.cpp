#include "ptorusmonitor.h"
#include "ui_ptorusmonitor.h"

PTorusMonitor::PTorusMonitor(QString workingDir, QStringList expectedFiles, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::PTorusMonitor)
{
    ui->setupUi(this);

    //Setup the Directories
    this->workingDirectory = workingDir;
    this->expectedOutputs = expectedFiles;
    this->noExpectedFiles = this->expectedOutputs.count();
    this->watcher = new QFileSystemWatcher(this);
    this->watcher->addPath(workingDirectory);
    this->watcher->addPaths(QStringList() << workingDirectory + "/albedo.dat" << workingDirectory + "/convergence_lucy.dat" << workingDirectory + "/info_grid.dat" << workingDirectory + "/gridFromFitsFile.vtu" << workingDirectory + "/rho.vtu");
    connect(watcher, SIGNAL(directoryChanged(const QString &)), this, SLOT(processDirectoryChange(const QString & )));
    connect(watcher, SIGNAL(fileChanged(QString)), this, SLOT(processDirectoryChange(const QString & )));

    //Setup the Timer
    connect(&execTimer, SIGNAL(timeout()), this, SLOT(updateTimer()));
    this->runningTime = 0;

    //Setup the Performance Tab
    ui->plt_TimeHist->enableAxis(QwtPlot::yLeft, false);
    prevFile = "torus.gfortran";
    this->dustCalculationsComplete = false;

    //Setup the Convergence Tab
    ui->grp_convPlot->setAxisTitle(QwtPlot::xBottom, "Iteration");
    ui->grp_convPlot->setAxisTitle(QwtPlot::yLeft, "Mean dT");
    ui->grp_convPlot->setAxisScale(QwtPlot::xBottom, 0, 1, 1);
    this->convIterations = 0;
}

PTorusMonitor::~PTorusMonitor()
{
    delete ui;
}

void PTorusMonitor::closeEvent(QCloseEvent *e) {
    e->ignore();
    this->hide();
}

void PTorusMonitor::updateTimer() {
    QString labelText;

    this->runningTime++;
    labelText = timeString(0, false);
    ui->lbl_RunningTime->setText("Running Time: " + labelText);
}

void PTorusMonitor::addEvent(PMonitorEvent::TorusMonitorEvents eventType, QString description) {
    PMonitorEvent event;
    event.time = this->runningTime;
    event.type = eventType;
    this->eventList.append(event);
    this->refreshEvents();
    this->prevEvent = &this->eventList.last();
}

void PTorusMonitor::addEvent(QString file) {
    PMonitorEvent event;
    event.time = this->runningTime;
    bool saveEvent = false;

    //qDebug() << "Modified File: " << file;

    //These are the rules to predict when events are happening.
    //These are guess at very best.

    if(file == "info_grid.dat") {
        if(prevFile != "info_grid.dat"){
            event.type = PMonitorEvent::GridCreationComplete;
            saveEvent = true;

            //Update the Grid Info in the UI
            updateGridInfo();
        }

    } else if(file == "convergence_lucy.dat") {
        if(this->prevFile != "convergence_lucy.dat" || this->prevEvent->time < (this->runningTime - 3)) {
            event.type = PMonitorEvent::LucyIterationComplete;
            saveEvent = true;

            //Plot Lucy Iterations
            this->updateConvergencePlot();
        }

    } else if(file == "beforesmooth.vtu") {
        event.type = PMonitorEvent::GridSmoothingStarted;
        saveEvent = true;

    } else if(file == "gridFromFitsFile.vtu") {
        if(prevFile != "gridFromFitsFile.vtu") {
            event.type = PMonitorEvent::FitsReadComplete;
            saveEvent = true;
        }

    } else if(file == "albedo.dat") {
        if(prevFile != "rho.vtu" && !this->dustCalculationsComplete) {
            event.type = PMonitorEvent::DustCalculationsComplete;
            this->dustCalculationsComplete = true;
            saveEvent = true;
        }

    } else if(file == "aftersmooth.vtu") {
        event.type = PMonitorEvent::GridSmoothingComplete;
        saveEvent = true;

    }

    if(this->expectedOutputs.contains(file)) {
        this->expectedOutputs.removeOne(file);

        emit this->progressBarNeedsUpdate(this->noExpectedFiles - this->expectedOutputs.count(), 0, this->noExpectedFiles);

        if(file.contains(QRegExp(".dat"))) {
            event.type = PMonitorEvent::SEDOutput;
            event.description = "Writen to " + file;
        } else {
            event.type = PMonitorEvent::FITSImageOutput;
            event.description = "Writen to " + file;
        }

        saveEvent = true;
    }


    if(saveEvent) {
        this->eventList.append(event);
        this->refreshEvents();
        this->prevEvent = &this->eventList.last();
    }

    this->prevFile = file;
}

void PTorusMonitor::processDirectoryChange(QString dir) {
    Q_UNUSED(dir);
    QString file = getLastModifiedFile();
    this->watcher->addPath(workingDirectory + "/" + file);
    this->addEvent(file);
}

void PTorusMonitor::processFileChange(QString file) {
    this->addEvent(file);
}

QString PTorusMonitor::getLastModifiedFile() {
    QDir *workingDir = new QDir(this->workingDirectory);
    QString file = workingDir->entryList(QStringList() << "*", QDir::Files, QDir::Time).first();
    return file;
}

void PTorusMonitor::refreshEvents() {
    ui->tbl_eventTable->clear();
    ui->tbl_eventTable->setRowCount(0);
    unsigned long int prevTime;
    int row = 0;

    ui->tbl_eventTable->setHorizontalHeaderLabels(QStringList() << "Event" << "Description" << "Time");

    foreach(PMonitorEvent e, this->eventList) {
        ui->tbl_eventTable->setRowCount(ui->tbl_eventTable->rowCount() + 1);
        QTableWidgetItem *rowEvent = new QTableWidgetItem();
        QTableWidgetItem *rowDesc = new QTableWidgetItem();
        QTableWidgetItem *rowTime = new QTableWidgetItem();

        //Event Item
        switch(e.type) {
            case PMonitorEvent::TorusStart :
                rowEvent->setText("Torus Started");
                break;

            case PMonitorEvent::FitsReadComplete :
                rowEvent->setText("FITS File Read Complete");
                break;

            case PMonitorEvent::GridCreationComplete :
                rowEvent->setText("Grid Creation Complete");
                break;

            case PMonitorEvent::LucyIterationComplete :
                rowEvent->setText("Lucy Iteration Complete");
                break;

            case PMonitorEvent::LucyConverged :
                rowEvent->setText("Lucy Algorithm Converged");
                break;

            case PMonitorEvent::DustCalculationsComplete :
                rowEvent->setText("Dust Calculations Complete");
                break;

            case PMonitorEvent::FITSImageOutput :
                rowEvent->setText("FITS Image Written");
                break;

            case PMonitorEvent::SEDOutput :
                rowEvent->setText("SED File Written");
                break;

            case PMonitorEvent::TorusFinish :
                rowEvent->setText("Torus Execution Complete");
                break;

            default:
                rowEvent->setText("N/A");
                break;
        }
        ui->tbl_eventTable->setItem(row, 0, rowEvent);

        //Desc Item
        rowDesc->setText(e.description);
        ui->tbl_eventTable->setItem(row, 1, rowDesc);

        //Time Item
        rowTime->setText(timeString(e.time));
        ui->tbl_eventTable->setItem(row, 2, rowTime);

        row++;
    }
}

QString PTorusMonitor::timeString(unsigned long time, bool timeset) {
    int timeMins;
    int timeSecs;
    int timeHours;
    unsigned long int rem;
    QString timeString;

    if(timeset == false) {
        time = runningTime;
    }

    timeHours = floor(time / 3600);
    rem = time % 3600;
    timeMins = floor(rem / 60);
    timeSecs = rem % 60;

    timeString.sprintf("%02d:%02d:%02d", timeHours, timeMins, timeSecs);
    return timeString;
}

void PTorusMonitor::updateConvergencePlot() {
    QList<QList<double> > data;
    QwtPlotCurve *convCurve = new QwtPlotCurve();
    QPolygonF plotSamples;

    qDebug() << "Lucy Plot Updated";

    ui->grp_convPlot->detachItems(QwtPlotItem::Rtti_PlotCurve, true);
    data = this->parseConvergenceFile();

    for(int i = 0; i < data.count(); i++) {
        QPointF point;
        point.setX(i + 1);
        point.setY(data.at(i).at(1));
        plotSamples << point;
    }

    convCurve->setSamples(plotSamples);
    convCurve->setPen(Qt::red, 0, Qt::SolidLine);
    convCurve->setSymbol(new QwtSymbol(QwtSymbol::Cross, QBrush(Qt::black, Qt::SolidPattern), QPen(Qt::black, 0, Qt::SolidLine), QSize(10, 10)));
    convCurve->setTitle("Convergence Plot");
    convCurve->attach(ui->grp_convPlot);
    ui->grp_convPlot->replot();
    ui->grp_convPlot->setAxisScale(QwtPlot::xBottom, 1, plotSamples.count() + 2, 1);

    convIterations++;
    ui->lbl_noIterations->setText("Lucy Iterations: " + QString::number(convIterations));
}

QList<QList<double> > PTorusMonitor::parseConvergenceFile(QString filename) {
    QString loadFile = this->workingDirectory + "/" + filename;
    QFile file(loadFile);
    QList<QList<double> > retval;
    QList<double> lineValue;
    QMessageBox msgBox;

    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return QList<QList<double> >();

    while(file.atEnd() == false) {
        QByteArray line = file.readLine();
        lineValue = this->parseConvergenceLine(line);

        if(!lineValue.isEmpty()) {
                retval << lineValue;
        }
    }

    return retval;
}

QList<double> PTorusMonitor::parseConvergenceLine(QByteArray line) {
    DataParseStates state = ParseStart;
    QString buffer;
    QList<double> returnList;

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
                    returnList.append(buffer.toDouble());
                    buffer.clear();
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

    return returnList;
}

/* This section is for parsing the grid info file */
void PTorusMonitor::updateGridInfo() {
    //Parse the info_grid.dat file and get QMap back
    QMap<QString, QString> infoGridData;
    infoGridData = this->parseInfoGridFile();
    this->gridInfo = infoGridData;
    this->formatGridInfoStrings();
    qDebug("Grid Info Update");

    //Update the UI
    ui->lbl_geometryType->setText("Geometry: " + gridGetValue("geometry"));
    ui->lbl_maxDepth->setText("Max AMR Depth: " + gridGetValue("maxDepth"));
    ui->lbl_halfSmallestSubcell->setText("Half Smallest Subcell: " + gridGetValue("halfSmallestSubcell"));
    ui->lbl_noOctals->setText("No. of Octals: " + gridGetValue("nOctals"));
    ui->lbl_noVoxels->setText("No. of Voxels: " + gridGetValue("nVoxels"));
    ui->lbl_smoothingFactor->setText("Smoothing Factor: " + gridGetValue("smoothingFactor"));
    ui->lbl_gridCentre->setText("Grid Centre: " + gridGetValue("grid center"));
    ui->lbl_sizeLargestCell->setText("Size of Largest Cell: " + gridGetValue("Size of largest cell"));
}

QMap<QString, QString> PTorusMonitor::parseInfoGridFile(QString filename) {
    QString loadFile = this->workingDirectory + "/" + filename;
    QFile file(loadFile);
    QMap<QString, QString> retval;
    QMap<QString, QString> lineValue;
    QMessageBox msgBox;

    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        return QMap<QString, QString>();

    while(file.atEnd() == false) {
        QByteArray line = file.readLine();
        lineValue = parseGridInfoLine(line);

        if(!lineValue.isEmpty()) {
            retval.unite(lineValue);
        }
    }

    qDebug() << this->gridInfoUnits;
    return retval;
}

QMap<QString, QString> PTorusMonitor::parseGridInfoLine(QByteArray line) {
    GridInfoParseStates state = IGParseStart;
    QString buffer;
    QMap<QString, QString> returnval;
    QString key;
    QString value;

    for(int i = 0; i < line.count(); i++) {
        char character = line.at(i);

        if(character == '#') { state = IGParseComment; }

        switch(state) {
            case IGParseStart :
                switch(character) {

                    case ' ' :
                        state = IGParseExpectKey;
                    break;

                    default :
                        buffer.append(character);
                        state = IGParseExpectEquals;
                    break;
                }

            break;

            case IGParseExpectKey :
                switch(character) {

                    case ' ' :
                    break;

                    default :
                        buffer.append(character);
                        state = IGParseExpectEquals;
                    break;
                }
            break;

            case IGParseExpectEquals:
                switch(character) {

                    case '=' :
                        key = buffer;
                        buffer.clear();
                        state = IGParseExpectValue;
                    break;

                    default :
                        buffer.append(character);
                    break;
                }
            break;

            case IGParseExpectValue :
                switch(character) {

                    case '[' :
                        value = buffer;
                        buffer.clear();
                        returnval.insert(key.trimmed(), value.trimmed());
                        state = IGParseReadUnit;
                    break;

                    default :
                        buffer.append(character);
                    break;
                }

                if(i == line.count() - 1) {
                    value = buffer;
                    buffer.clear();
                    returnval.insert(key.trimmed(), value.trimmed());
                }
            break;

            case IGParseReadUnit :
                switch(character) {

                    case ']' :
                        this->gridInfoUnits.insert(key.trimmed(), buffer.trimmed());
                        buffer.clear();
                        state = IGParseLineEnd;
                    break;

                    default :
                        buffer.append(character);
                    break;
                }
            break;

            case IGParseComment :
            break;

            case IGParseLineEnd :
            break;

            default :
            break;
        }
    }

    return returnval;
}

QString PTorusMonitor::gridGetValue(QString key) {
    if(this->gridInfo.value(key) != "") {
        QString value = this->gridInfo.value(key);

        if(gridGetUnits(key) != "") {
            value.append(" (" + gridGetUnits(key) + ")");
        }

        return value;
    }
    return "";
}

QString PTorusMonitor::gridGetUnits(QString key) {
    if(gridInfoUnits.value(key) != "") {
        return gridInfoUnits.value(key);
    }
    return "";
}

void PTorusMonitor::formatGridInfoStrings() {

    //Reformat the Doubles for the UI
    QString newString;
    QRegExp reg("([0-9]*\\.[0-9]*)");
    QStringList cap;

    foreach(QString key, this->gridInfo.keys()) {
        int pos = 0;
        QString value = gridInfo.value(key);
        while(reg.indexIn(value, pos) != -1) {
            cap << reg.cap(1);
            pos += reg.matchedLength();
        }

        if(!cap.isEmpty()) {
            for(int i = 0; i < cap.count(); i++) {
                QString dub = cap.at(i);
                newString += QString::number(dub.toDouble(), 'g') + " ";
            }
            this->gridInfo.insert(key, newString);
        }
        newString.clear();
        cap.clear();
    }
}
