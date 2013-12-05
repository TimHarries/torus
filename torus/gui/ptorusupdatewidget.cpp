#include "ptorusupdatewidget.h"
#include "ui_ptorusupdatewidget.h"

PTorusUpdateWidget::PTorusUpdateWidget(RunMode runMode, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PTorusUpdateWidget)
{
    ui->setupUi(this);

    this->buildProcess = 0x0;
    this->extractProcess = 0x0;
    this->makeCommand = "make";

    //Determine OS Compatibility
    #ifdef Q_OS_UNIX
        this->isUnixLikeSystem = true;
    #else
        this->isUnixLikeSystem = false;
    #endif

    QString file;
    QMessageBox message;
    switch(runMode) {
        case PTorusUpdateWidget::UpdateFromWeb :
            this->close();
            break;

        case PTorusUpdateWidget::UpdateFromBin :
            file = QFileDialog::getOpenFileName(this, "Select the Torus Executable", QDir::homePath(), "*.gfortran");
            if(file != "") {
                copyToAppBundle(file);
                ui->txt_UpdateLog->append("<b>Copying Torus Execuable: " + file + "</b>");
                ui->txt_UpdateLog->append("<b>Torus has been Successfuly Updated.</b>");
                ui->progressBar->setMaximum(1);
                ui->progressBar->setValue(1);
                ui->btn_Cancel->buttons().at(0)->setText("Close");

            } else {
                ui->txt_UpdateLog->append("<b>Torus update has been cancelled.</b>");
                ui->progressBar->setMaximum(1);
                ui->progressBar->setValue(0);
                ui->btn_Cancel->buttons().at(0)->setText("Close");
            }
            break;

        case PTorusUpdateWidget::UpdateFromTar :
            compileToBin();
            this->close();
            break;

        default :
            break;
    }
}

PTorusUpdateWidget::~PTorusUpdateWidget()
{
    delete ui;
}

void PTorusUpdateWidget::on_btn_Cancel_rejected()
{
    if(this->extractProcess != 0x0)
        this->extractProcess->kill();
    if(this->buildProcess != 0x0)
        this->buildProcess->kill();

    this->close();
}

void PTorusUpdateWidget::copyToAppBundle(QString torusExeFile) {
    QDir bundleBinDir;
    QString bundleBinDirStr = QCoreApplication::applicationDirPath();
    bundleBinDir.cd(bundleBinDirStr);
    if(QFile::exists(bundleBinDir.path() + "/torus.gfortran"))
        QFile::remove(bundleBinDir.path() + "/torus.gfortran");
    QFile::copy(torusExeFile, bundleBinDir.path() + "/torus.gfortran");
}

void PTorusUpdateWidget::compileToBin() {

    if(this->isUnixLikeSystem) {

        PCompileConfigureDialog compileConfigDialog;
        compileConfigDialog.exec();
        this->cfitsLibs = compileConfigDialog.getCfitsioDir();
        this->cfitsioEnabled = compileConfigDialog.getCfitsioEnabled();
        this->openMpEnabled = compileConfigDialog.getOpenMpEnabled();
        this->MpiEnabled = compileConfigDialog.getMpiEnabled();
        this->tarFilePath = compileConfigDialog.getSourceFile();
        this->wrkDirStr = compileConfigDialog.getBuildDir();
        this->buildDir = new QDir(this->wrkDirStr);
        this->copyToBundle = compileConfigDialog.getCopyToBundle();
        int result = compileConfigDialog.result();

        if(result == QDialog::Accepted) {

            //Do Possible Cleanup
            QFile::remove(wrkDirStr + "/torus.tar.gz");
            QDir cleanup(wrkDirStr + "/torus");
            cleanup.removeRecursively();

            //Assume that we have everything we need
            QFile tarFile(tarFilePath);
            tarFile.copy(buildDir->path() + "/torus.tar.gz");

            this->extractProcess = new QProcess(this);
            connect(extractProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(writeExtractToConsole()));
            connect(extractProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(extractFinished(int, QProcess::ExitStatus)));
            extractProcess->setWorkingDirectory(buildDir->path());
            extractProcess->start("tar", QStringList() << "zxvf" << buildDir->path() + "/torus.tar.gz");
            ui->txt_UpdateLog->append("Extracting...");
        } else {
            ui->txt_UpdateLog->append("<b>Torus build has been cancelled.</b>");
            ui->progressBar->setMaximum(1);
            ui->btn_Cancel->buttons().at(0)->setText("Close");
        }

    } else {
        QMessageBox box;
        box.setText("Sorry, In-App Compilation doesn't currently support your platform");
    }
}


void PTorusUpdateWidget::writeExtractToConsole() {
    QString read = extractProcess->readAllStandardOutput();
    ui->txt_UpdateLog->append(read);
}

void PTorusUpdateWidget::writeBuildToConsole() {
    QString read = buildProcess->readAllStandardOutput();
    ui->txt_UpdateLog->append(read);
}

void PTorusUpdateWidget::extractFinished(int status, QProcess::ExitStatus exitStatus) {
    Q_UNUSED(status);
    qDebug() << "Extract Finished";
    QProcessEnvironment env;
    ui->txt_UpdateLog->append(extractProcess->readAll());

    if(exitStatus == QProcess::NormalExit) {
        this->buildDir->cd("torus");

        //Compile Torus
        this->buildProcess = new QProcess();
        connect(buildProcess, SIGNAL(readyReadStandardOutput()), this, SLOT(writeBuildToConsole()));
        connect(buildProcess, SIGNAL(finished(int, QProcess::ExitStatus)), this, SLOT(buildFinished(int, QProcess::ExitStatus)));

        //Setup File System Watcher
        this->watcher = new QFileSystemWatcher();
        this->fileList = this->buildDir->entryList(QStringList() << "*.f90" << "*.F90", QDir::Files).replaceInStrings(QRegExp("(f90|F90)"), "o");
        watcher->addPath(this->buildDir->path());
        connect(watcher, SIGNAL(directoryChanged(QString)), this, SLOT(buildDirChanged(QString)));
        ui->progressBar->setMaximum(this->fileList.count() + 1);

        env.insert("PATH", "/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin");
        this->buildProcess->setWorkingDirectory(this->buildDir->path());
        this->buildProcess->setProcessEnvironment(env);
        QStringList args;

        if(this->cfitsioEnabled) {
            if(this->cfitsLibs != "") {
                args.append("CFITSIOLIBS=-lcfitsio -L" + this->cfitsLibs);
            } else {
                args.append("CFITSIOLIBS=-lcfitsio ");
            }
            args.append("cfitsio=yes");
        }

        if(this->openMpEnabled) {
            args.append("openmp=yes");
        } else {
            args.append("openmp=no");
        }

        if(this->MpiEnabled) {
            args.append("mpi=yes");
        } else {
            args.append("mpi=no");
        }

        buildProcess->start(this->makeCommand, args);
        ui->txt_UpdateLog->append("Building...");

    } else {
        QMessageBox box;
        box.setText("There was a problem extracting the Torus Tarball");
        box.setInformativeText(extractProcess->errorString());
        box.exec();
    }
}

void PTorusUpdateWidget::buildFinished(int status, QProcess::ExitStatus exitStatus) {
    qDebug() << "Build Finished";
    if(exitStatus == QProcess::NormalExit) {
        if(this->copyToBundle) {
            this->copyToAppBundle(this->buildDir->path() + "/torus.gfortran");
            ui->txt_UpdateLog->append("<b>Updated Internal Torus Build</b>");
        }
        ui->progressBar->setValue(ui->progressBar->maximum());

        if(status != 0) {
            ui->txt_UpdateLog->append("<b style=' color: red;'>There was a Problem Compiling Torus: </b><br />");
            ui->txt_UpdateLog->append("<p style='color: red;'>" + buildProcess->readAllStandardError() + "</p>");
            ui->btn_Cancel->buttons().at(0)->setText("Close");
            ui->lbl_Heading->setText("Update Failed");
        } else {
            ui->txt_UpdateLog->append("<b>Torus has been Successfully Updated</b>");
            ui->btn_Cancel->buttons().at(0)->setText("Close");
            ui->lbl_Heading->setText("Update Complete!");
        }
    } else {
        ui->txt_UpdateLog->append("<b>There was a Problem Compiling Torus</b><br />");
        ui->btn_Cancel->buttons().at(0)->setText("Close");
        ui->lbl_Heading->setText("Update Failed");
    }
}

void PTorusUpdateWidget::buildDirChanged(QString dir) {
    Q_UNUSED(dir);
    this->buildDir->setNameFilters(this->fileList);
    if(this->buildDir->entryList().count() > 0) {
        QString file = this->buildDir->entryList(this->fileList, QDir::Files, QDir::Time).first();
        if(this->fileList.contains(file)) {
            ui->progressBar->setValue(ui->progressBar->value() + 1);
            this->fileList.removeOne(file);
        }
    }
}

void PTorusUpdateWidget::buildFileChanged(QString file) {
    if(this->fileList.contains(QFileInfo(file).fileName())) {
        ui->progressBar->setValue(ui->progressBar->value() + 1);
        this->fileList.removeOne(QFileInfo(file).fileName());
    }
}
