#include "pcompileconfiguredialog.h"


PCompileConfigureDialog::PCompileConfigureDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PCompileConfigureDialog)
{
    ui->setupUI(this);
}

QString PCompileConfigureDialog::getCfitsioDir() {
    return ui->txt_cfitsioDir->text();
}

QString PCompileConfigureDialog::getMakeCommand() {
    return ui->txt_makeCommand->text();
}

QString PCompileConfigureDialog::getSourceFile() {
    return ui->txt_sourceFile->text();
}

QString PCompileConfigureDialog::getBuildDir() {
    return ui->txt_buildDir->text();
}

bool PCompileConfigureDialog::getOpenMpEnabled() {
    return ui->chk_openMpFlag->isChecked();
}

bool PCompileConfigureDialog::getMpiEnabled() {
    return ui->chk_mpiFlag->isChecked();
}

bool PCompileConfigureDialog::getCfitsioEnabled() {
    return ui->chk_cfitsioFlag->isChecked();
}

bool PCompileConfigureDialog::getCopyToBundle() {
    return ui->chk_copyToBundle->isChecked();
}

void PCompileConfigureDialog::openFitsLibDir() {
    QFileDialog fitsOpen;
    fitsOpen.setFileMode(QFileDialog::Directory);
    fitsOpen.setOptions(QFileDialog::ShowDirsOnly);
    fitsOpen.setWindowTitle("Select the FITS Library Directory");
    fitsOpen.exec();
    if(fitsOpen.result()) {
        ui->txt_cfitsioDir->setText(fitsOpen.selectedFiles()[0]);
    }
}

void PCompileConfigureDialog::openBuildDir() {
    QFileDialog buildDirOpen;
    buildDirOpen.setFileMode(QFileDialog::Directory);
    buildDirOpen.setOptions(QFileDialog::ShowDirsOnly);
    buildDirOpen.setWindowTitle("Select the Torus Build Directory");
    buildDirOpen.exec();
    if(buildDirOpen.result()) {
        ui->txt_buildDir->setText(buildDirOpen.selectedFiles()[0]);
    }
}

void PCompileConfigureDialog::openSourceFile() {
    QFileDialog sourceFileOpen;
    sourceFileOpen.setWindowTitle("Select the Torus Source File");
    sourceFileOpen.exec();
    if(sourceFileOpen.result()) {
        ui->txt_sourceFile->setText(sourceFileOpen.selectedFiles()[0]);
    }
}

void PCompileConfigureDialog::validateAccept() {
    qDebug() << "Accepted";
    QMessageBox *warning = new QMessageBox();
    if(ui->txt_buildDir->text() != "") {
        if(ui->txt_sourceFile->text() != "") {
            if(ui->txt_makeCommand->text() != "") {
                this->accept();
            } else {
                warning->setText("Please make sure that Make Command is not empty.");
                warning->exec();
            }
        } else {
            warning->setText("Please make sure that you've selected the Torus Source File.");
            warning->exec();
        }
    } else {
        warning->setText("Please make sure that you've selected a Torus build directory");
        warning->exec();
    }
}
