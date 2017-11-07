#include "torusconsole.h"
#include "ui_torusconsole.h"

torusconsole::torusconsole(QWidget *parent, QString torusPath, QString workingPath, bool writeToLog, RunModes RunMode, QStringList expectedFiles) :
    QWidget(parent),
    ui(new Ui::torusconsole)
{
    ui->setupUi(this);
    proc = new QProcess(this);
    this->writeToLog = false;
    this->isRunning = false;
    qDebug() << expectedFiles;

    connect(proc, SIGNAL(readyReadStandardOutput()), this, SLOT(slotDataOnStdout()));
    connect(proc, SIGNAL(readyReadStandardError()), this, SLOT(slotProcessError()));
    connect(proc, SIGNAL(error(QProcess::ProcessError)), this, SLOT(slotProcessError()));
    connect(proc, SIGNAL(started()),this, SLOT(slotProcessStart()));
    connect(proc, SIGNAL(finished(int)), this, SLOT(slotProcessFinish()));

    //Get the binary of TORUS we're going to use
    if(torusPath == "" || QFile(torusPath).exists() == false) {
        torusPath = this->getInternalBuild();
        ui->txt_conOut->append("<b>Using internal build of Torus...</b><br />");
        ui->txt_conOut->append("Torus Binary: " + torusPath + "<br />");

        if(torusPath == "") {
            ui->txt_conOut->append("No internal build of Torus Found. Please update the internal build from the <b>Run</b> menu.");
            ui->prg_Progress->setMaximum(1);
        }
    } else {
        ui->txt_conOut->append("<b>Using external build of Torus...</b><br />");
        ui->txt_conOut->append("Torus Binary: " + torusPath + "<br />");
    }

    QString command;
    QStringList arguments;
    QString exeFile = torusPath.split("/").last();
    QString exePath = torusPath.remove(QRegExp("/[a-zA-Z.]+$"));
    QString workingDirStr = workingPath;
    QProcessEnvironment env;

    //Get the directory
    QDir runningDir(exePath);
    QDir workingDir(workingDirStr);

    //Set Run Mode
    switch(RunMode) {
        case DefaultRunMode :
            command = runningDir.filePath(exeFile);
        break;

        case CheckRunMode :
            command = runningDir.filePath(exeFile);
            arguments << "check";
        break;

        default :
            command = runningDir.filePath(exeFile);
            arguments << "check";
        break;
    }

    //Setup the Scrollbar
    //this->update = new QFileSystemWatcher(this);
    //this->update->addPath(workingDir.path());
    this->expectedFiles = expectedFiles;
    this->noExpectedFiles = expectedFiles.count();
    this->workingDirectory = workingDir.path();
    //connect(update, SIGNAL(directoryChanged(const QString &)), this, SLOT(slotFileUpdate(const QString &)));

    //Setup Environment
    env.insert("TORUS_DATA", workingDir.path());

    qDebug() << "Setting Working Directory to: " << workingDir.path();
    proc->setWorkingDirectory(workingDir.path());
    qDebug() << "Running: " << command;
    proc->setProcessEnvironment(env);
    proc->start(command, arguments);
    this->setWindowTitle(workingDir.path() + "/parameters.dat");

    //Setup Log file
    if(writeToLog) {
        this->writeToLog = true;
        this->logFile = new QFile(workingDir.filePath("torus.log"));
        this->logFile->open(QFile::WriteOnly | QFile::Text);
        this->logFileOut = new QTextStream(this->logFile);
    }

    //Setup the Monitor
    this->monitor = new PTorusMonitor(this->workingDirectory, this->expectedFiles);
    this->monitor->attachedWindow = this;
    connect(monitor, SIGNAL(progressBarNeedsUpdate(int,int,int)), this, SLOT(slotProgressUpdate(int,int,int)));
}

torusconsole::~torusconsole()
{
    proc->kill();
    delete ui;
}

/* Process Event Handlers */
void torusconsole::slotDataOnStdout() {
    QString read = proc->readAllStandardOutput();
    ui->txt_conOut->append(read);
    if(this->writeToLog == true)
        *this->logFileOut << read;
}

void torusconsole::slotProcessStart() {
    this->isRunning = true;
    this->monitor->execTimer.start(1000);
    this->monitor->addEvent(PMonitorEvent::TorusStart);
}

void torusconsole::slotProcessFinish() {
    this->isRunning = false;
    this->monitor->execTimer.stop();
    this->monitor->addEvent(PMonitorEvent::TorusFinish);
    ui->prg_Progress->setMaximum(1);
    ui->prg_Progress->setValue(1);

    if(writeToLog)
        this->logFile->close();

    switch(proc->exitStatus()) {
        case QProcess::NormalExit :
            qDebug() << "Process Complete";
            break;

        case QProcess::CrashExit :
            qDebug() << "Process Crashed: " << proc->errorString();
            break;

        default :
            qDebug() << "Process Exited: " << proc->errorString();
            ui->txt_conOut->append(proc->readAllStandardError());
            break;
    }
}

void torusconsole::slotProcessError() {
    QString read = proc->readAllStandardError();
    ui->txt_conOut->append(read);
    if(writeToLog)
        *this->logFileOut << read;
}

void torusconsole::on_btn_closeProcess_clicked()
{
    if(isRunning) {
        proc->close();
        ui->txt_conOut->append("--- Torus has been quit ---");
        this->isRunning = false;

        if(writeToLog)
            this->logFile->close();
    }
}

void torusconsole::closeEvent(QCloseEvent *event) {
    if(this->isRunning) {
        QMessageBox box;
        box.setText("Are you sure you want to close this window?");
        box.setInformativeText("Closing this window will also quit the current Torus session.");
        box.setStandardButtons(QMessageBox::Yes | QMessageBox::Cancel);

        if(box.exec() == QMessageBox::Yes) {
            proc->kill();
            proc->waitForFinished();
            this->monitor->hide();
            delete this->monitor;
            event->accept();
        } else {
            event->ignore();
        }
    } else {
        this->monitor->hide();
        delete this->monitor;
        event->accept();
    }
}

void torusconsole::slotFileUpdate(QString dir) {
    Q_UNUSED(dir);

    QDir *workingDir = new QDir(this->workingDirectory);
    workingDir->setNameFilters(this->expectedFiles);
    if(workingDir->entryList().count() > 0) {
        //Stop the Progress Bar being Indeterminate
        if(ui->prg_Progress->maximum() == 0) {
            ui->prg_Progress->setMinimum(0);
            ui->prg_Progress->setMaximum(this->expectedFiles.count());
        }

        QString file = workingDir->entryList(this->expectedFiles, QDir::Files, QDir::Time).first();
        if(this->expectedFiles.contains(file)) {
            this->expectedFiles.removeOne(file);
            ui->prg_Progress->setValue(this->noExpectedFiles - this->expectedFiles.count());
        }
    }
}

void torusconsole::slotProgressUpdate(int progress, int min, int max) {
    ui->prg_Progress->setMinimum(min);
    ui->prg_Progress->setMaximum(max);
    ui->prg_Progress->setValue(progress);
}

QString torusconsole::getInternalBuild() {
    QDir bundleBinDir;
    QString bundleBinDirStr = QCoreApplication::applicationDirPath();
    bundleBinDir.cd(bundleBinDirStr);
    QString filePath = bundleBinDir.path() + "/torus.gfortran";

    if(QFile(filePath).exists()) {
        return filePath;
    }

    return "";
}

void torusconsole::on_btn_Monitor_clicked()
{
    this->monitor->show();
}
