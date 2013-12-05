#ifndef TORUSCONSOLE_H
#define TORUSCONSOLE_H

#include <QWidget>
#include "qprocess.h"
#include <QDir>
#include <QDebug>
#include <QCloseEvent>
#include <QMessageBox>
#include <QFileSystemWatcher>
#include <QTimer>
#include <QString>
#include <math.h>
#include <ptorusmonitor.h>

namespace Ui {
class torusconsole;
}

class torusconsole : public QWidget
{
    Q_OBJECT
    
public:
    enum RunModes {
        DefaultRunMode,
        CheckRunMode
    };

    explicit torusconsole(QWidget *parent = 0, QString torusPath = "", QString workingPath = "", bool writeToLog = true, RunModes RunMode = torusconsole::DefaultRunMode, QStringList expectedFiles = QStringList());
    ~torusconsole();
    
private:
    Ui::torusconsole *ui;
    QProcess *proc;
    bool isRunning;
    bool writeToLog;
    QFile *logFile;
    QTextStream *logFileOut;
    void closeEvent(QCloseEvent *event);
    int noExpectedFiles;
    QTimer execTimer;
    long unsigned int RunningTime;

    //Progress Bar Properties
    QFileSystemWatcher *update;
    QString workingDirectory;
    QStringList expectedFiles;
    QStringList completeFiles;
    QString getInternalBuild();

    //Monitor Window
    PTorusMonitor *monitor;

private slots:
    void slotDataOnStdout();
    void slotProcessFinish();
    void slotProcessStart();
    void slotProcessError();
    void slotFileUpdate(QString dir);
    void on_btn_closeProcess_clicked();
    void on_btn_Monitor_clicked();
    void slotProgressUpdate(int progress, int min, int max);
};

#endif // TORUSCONSOLE_H
