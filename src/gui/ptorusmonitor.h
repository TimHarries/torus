#ifndef PTORUSMONITOR_H
#define PTORUSMONITOR_H

#include <QWidget>
#include <QCloseEvent>
#include <QTimer>
#include <QList>
#include <QFileSystemWatcher>
#include <QDir>
#include <QFile>
#include <QDebug>
#include <QMessageBox>

#include <qwt_plot_curve.h>
#include <qwt_symbol.h>

#include <pmonitorevent.h>
#include <math.h>

namespace Ui {
class PTorusMonitor;
}

class PTorusMonitor : public QWidget
{
    Q_OBJECT
    
public:
    enum DataParseStates {
        ParseStart,
        ParseExpectNumber,
        ParseReadNumber,
        ParseReadComment
    };

    enum GridInfoParseStates {
        IGParseStart,
        IGParseExpectKey,
        IGParseExpectEquals,
        IGParseExpectValue,
        IGParseReadUnit,
        IGParseLineEnd,
        IGParseComment
    };

    explicit PTorusMonitor(QString workingDir = "", QStringList expectedFiles = QStringList(), QWidget *parent = 0);
    ~PTorusMonitor();
    QWidget *attachedWindow;
    QTimer execTimer;
    QStringList expectedOutputs;
    int convIterations;

    
private:
    Ui::PTorusMonitor *ui;

    QList<PMonitorEvent> eventList;
    void closeEvent(QCloseEvent *e);
    long unsigned int runningTime;
    QString workingDirectory;
    QString getLastModifiedFile();
    QString timeString(unsigned long int time, bool timeset = true);

    QFileSystemWatcher *watcher;
    PMonitorEvent *prevEvent;
    QString prevFile;

    //Some Flags for Event Logging
    bool dustCalculationsComplete;
    int noExpectedFiles;

    //Stuff to deal with Convergence
    void updateConvergencePlot();
    QList<QList<double> > parseConvergenceFile(QString filename = "convergence_lucy.dat");
    QList<double> parseConvergenceLine(QByteArray line);

    //Stuff to deal with grid info parsing.
    void updateGridInfo();
    QMap<QString, QString> gridInfoUnits;
    QMap<QString, QString> gridInfo;
    QMap<QString, QString> parseInfoGridFile(QString filename = "info_grid.dat");
    QMap<QString, QString> parseGridInfoLine(QByteArray line);
    QString gridGetValue(QString key);
    QString gridGetUnits(QString key);
    void formatGridInfoStrings();

private slots:
    void updateTimer();
    void processDirectoryChange(const QString dir);
    void processFileChange(const QString file);

public slots:
    void addEvent(PMonitorEvent::TorusMonitorEvents eventType, QString description = "");
    void addEvent(QString file);
    void refreshEvents();

signals:
    void progressBarNeedsUpdate(int progress, int min, int max);
};

#endif // PTORUSMONITOR_H
