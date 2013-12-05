#ifndef PCOMPILECONFIGUREDIALOG_H
#define PCOMPILECONFIGUREDIALOG_H

#include <QDialog>
#include <QCheckBox>
#include <QLineEdit>
#include <QPushButton>
#include <QLabel>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QFileDialog>
#include <QDir>
#include <QGroupBox>
#include <QMessageBox>
#include <QDebug>

namespace Ui {
    class PCompileConfigureDialog;
}

class PCompileConfigureDialog : public QDialog
{
    Q_OBJECT
public:
    explicit PCompileConfigureDialog(QWidget *parent = 0);

    //Accessors
    QString getCfitsioDir();
    QString getMakeCommand();
    QString getAdvFlags();
    QString getBuildDir();
    QString getSourceFile();
    bool getOpenMpEnabled();
    bool getMpiEnabled();
    bool getCfitsioEnabled();
    bool getCopyToBundle();

    void toggleAdvancedGroup();

private:
    Ui::PCompileConfigureDialog *ui;
    
signals:
    
public slots:

private slots:
    void openFitsLibDir();
    void openBuildDir();
    void openSourceFile();
    void validateAccept();
};


/* Initialise UI */
class Ui_PCompileConfigureDialog {
public :
    QLabel *lbl_heading;
    QCheckBox *chk_openMpFlag;
    QCheckBox *chk_mpiFlag;
    QCheckBox *chk_cfitsioFlag;
    QLineEdit *txt_cfitsioDir;
    QDialogButtonBox *grp_buttons;
    QLineEdit *txt_makeCommand;
    QLineEdit *txt_advFlags;
    QFormLayout *layout;
    QHBoxLayout *cfisDirContainer;
    QVBoxLayout *masterLayout;
    QLineEdit *advFlags;
    QDialogButtonBox *buttons;
    QPushButton *openCfitsDir;
    QPushButton *compileButton;
    QPushButton *cancelButton;
    QLabel *advancedToggle;
    QFormLayout *advLayout;
    QGroupBox *advGroup;
    QLabel *heading;
    QHBoxLayout *lay_buildDir;
    QLineEdit *txt_buildDir;
    QPushButton *btn_openBuildDir;
    QHBoxLayout *lay_sourceFile;
    QLineEdit *txt_sourceFile;
    QPushButton *btn_openSourceFile;
    QLabel *lbl_buildFlags;
    QCheckBox *chk_copyToBundle;

    void setupUI(PCompileConfigureDialog *compileConfig) {
        layout = new QFormLayout();
        masterLayout = new QVBoxLayout();
        compileConfig->setFixedSize(compileConfig->width(), compileConfig->height());

        //Add Heading
        QFont headFont;
        headFont.setFamily(headFont.defaultFamily());
        headFont.setPointSize(16);
        heading = new QLabel();
        heading->setText("Please Configure your Torus Build");
        heading->setFont(headFont);
        layout->addRow(heading);

        //Torus Build Directory Code
        txt_buildDir = new QLineEdit();
        btn_openBuildDir = new QPushButton();
        lay_buildDir = new QHBoxLayout();
        btn_openBuildDir->setText("Select");
        txt_buildDir->setMinimumSize(txt_buildDir->minimumSizeHint());
        lay_buildDir->addWidget(txt_buildDir);
        lay_buildDir->addWidget(btn_openBuildDir);
        layout->addRow("Torus Build Directory", lay_buildDir);
        QObject::connect(btn_openBuildDir, SIGNAL(clicked()), compileConfig, SLOT(openBuildDir()));

        //Torus Tarball Code
        txt_sourceFile = new QLineEdit();
        btn_openSourceFile = new QPushButton();
        lay_sourceFile = new QHBoxLayout();
        btn_openSourceFile->setText("Select");
        txt_sourceFile->setMinimumSize(txt_sourceFile->minimumSizeHint());
        lay_sourceFile->addWidget(txt_sourceFile);
        lay_sourceFile->addWidget(btn_openSourceFile);
        layout->addRow("Torus Source File", lay_sourceFile);
        QObject::connect(btn_openSourceFile, SIGNAL(clicked()), compileConfig, SLOT(openSourceFile()));

        //CFITS Library Directory
        txt_cfitsioDir = new QLineEdit();
        openCfitsDir = new QPushButton();
        cfisDirContainer = new QHBoxLayout();

        openCfitsDir->setText("Select");
        txt_cfitsioDir->setMinimumSize(txt_cfitsioDir->minimumSizeHint());
        cfisDirContainer->addWidget(txt_cfitsioDir);
        cfisDirContainer->addWidget(openCfitsDir);
        layout->addRow("FITS Library Directory", cfisDirContainer);
        QObject::connect(openCfitsDir, SIGNAL(clicked()), compileConfig, SLOT(openFitsLibDir()));

        //Add Checkboxes
        chk_openMpFlag = new QCheckBox();
        chk_mpiFlag = new QCheckBox();
        chk_cfitsioFlag = new QCheckBox();
        lbl_buildFlags = new QLabel("Build Flags");
        lbl_buildFlags->setFont(headFont);
        layout->addRow(lbl_buildFlags);
        chk_openMpFlag->setText("Use OpenMP");
        layout->addRow(chk_openMpFlag);
        chk_mpiFlag->setText("Use MPI");
        layout->addRow(chk_mpiFlag);
        chk_cfitsioFlag->setText("Use FITS File Support");
        layout->addRow(chk_cfitsioFlag);

        //Setup the Buttons
        buttons = new QDialogButtonBox();
        buttons->addButton("Compile Torus", QDialogButtonBox::AcceptRole);
        buttons->addButton("Cancel", QDialogButtonBox::RejectRole);

        //Advanced Controls
        advancedToggle = new QLabel("Advanced Options");
        advancedToggle->setFont(headFont);
        layout->addRow(advancedToggle);
        masterLayout->addLayout(layout);
        advLayout = new QFormLayout();

        advGroup = new QGroupBox();
        txt_makeCommand = new QLineEdit();
        advFlags = new QLineEdit();
        chk_copyToBundle = new QCheckBox();
        txt_makeCommand->setText("make");
        advLayout->addRow("Make Command", txt_makeCommand);
        advLayout->addRow("Additional Flags", advFlags);
        advLayout->addRow("Update Internal Torus Build", chk_copyToBundle);
        advGroup->setLayout(advLayout);
        masterLayout->addWidget(advGroup);

        QObject::connect(buttons, SIGNAL(accepted()), compileConfig, SLOT(validateAccept()));
        QObject::connect(buttons, SIGNAL(rejected()), compileConfig, SLOT(reject()));
        masterLayout->addWidget(buttons);
        compileConfig->setLayout(masterLayout);
    }
};

namespace Ui {
    class PCompileConfigureDialog : public Ui_PCompileConfigureDialog {};
}

#endif // PCOMPILECONFIGUREDIALOG_H
