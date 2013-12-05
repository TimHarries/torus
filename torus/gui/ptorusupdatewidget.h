#ifndef PTORUSUPDATEWIDGET_H
#define PTORUSUPDATEWIDGET_H

#include <QDialog>
#include <QTemporaryDir>
#include <QFileSystemWatcher>
#include <QFile>
#include <QDir>
#include <QProcess>
#include <QProcessEnvironment>
#include <QMessageBox>
#include <QFileDialog>
#include <QDebug>
#include <QFileSystemWatcher>
#include <QFormLayout>
#include <QCheckBox>
#include <QLineEdit>
#include <QPushButton>
#include <pcompileconfiguredialog.h>

namespace Ui {
class PTorusUpdateWidget;
}

class PTorusUpdateWidget : public QDialog
{
    Q_OBJECT
    
public:
    enum RunMode{
        UpdateFromTar,
        UpdateFromBin,
        UpdateFromWeb
    };

    explicit PTorusUpdateWidget(PTorusUpdateWidget::RunMode runMode, QWidget *parent = 0);
    ~PTorusUpdateWidget();

    void copyToAppBundle(QString torusExeFile);
    void compileToBin();
    
private slots:
    void on_btn_Cancel_rejected();

private:
    Ui::PTorusUpdateWidget *ui;
    bool isUnixLikeSystem;
    QProcess *buildProcess;
    QProcess *extractProcess;
    QTemporaryDir *tempDir;
    QDir *buildDir;
    QFileSystemWatcher *watcher;
    QStringList fileList;
    QStringList createdFiles;

    //Compile Properties
    QString cfitsLibs;
    QString makeCommand;
    QString wrkDirStr;
    QString tarFilePath;
    bool openMpEnabled;
    bool MpiEnabled;
    bool cfitsioEnabled;
    bool copyToBundle;

    void configureCompileConfigDialog(QDialog &dialog);

private slots:
    void writeExtractToConsole();
    void writeBuildToConsole();
    void extractFinished(int status, QProcess::ExitStatus exitStatus);
    void buildFinished(int status, QProcess::ExitStatus exitStatus);
    void buildDirChanged(QString dir);
    void buildFileChanged(QString file);
};

#endif // PTORUSUPDATEWIDGET_H
