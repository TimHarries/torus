#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "pparameterdata.h"
#include <panglegen.h>
#include <unitconverter.h>
#include <torusconsole.h>
#include <QString>
#include <QFileDialog>
#include <QMessageBox>
#include <QDropEvent>
#include <QMimeData>
#include <QUrl>
#include <QUndoStack>
#include <QUndoView>
#include "math.h"
#include <QMapIterator>
#include <QQueue>
#include <QLayout>
#include <ptorusupdatewidget.h>
#include <paboutdialog.h>
#include <QClipboard>
#include <typeinfo>
#include <psedplotter.h>
#include <qmath.h>

#define PI 3.14159265

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    PParameterData params;

    void refreshUI(void);
    void resetUI(void);
    bool canCloseCurrentDocument();
    void updateModel();
    void updateGeometry();
    void checkForTorusExecutable(QString filePath);
    QString getFileDialogPath();
    void renderGeometryUI();
    
private slots:

    void on_btn_AddSource_clicked();
    void on_actionClose_triggered();
    void on_actionOpen_triggered();
    void on_actionSave_triggered();
    void on_actionSaveAs_triggered();
    void on_btn_RemoveSource_clicked();
    void on_btn_SaveSource_clicked();
    void on_lst_Sources_itemSelectionChanged();
    void on_tck_spectrum_toggled(bool checked);
    void on_tck_images_toggled(bool checked);
    void on_lst_Images_itemSelectionChanged();
    void on_actionNew_triggered();
    void on_btn_SaveImage_clicked();
    void on_btn_AddImage_clicked();
    void on_btn_RemoveImage_clicked();
    void on_btn_sedGenInc_clicked();
    void on_spn_sedAddInc_clicked();
    void on_spn_sedRemoveInc_clicked();
    void on_btn_sedGenPos_clicked();
    void on_spn_sedAddPos_clicked();
    void on_spn_sedRemovePos_clicked();
    void on_actionPaste_triggered();
    void on_actionCopy_triggered();
    void on_actionCut_triggered();
    void on_actionRun_Torus_triggered();
    void on_actionCheck_Paramaters_triggered();
    void on_btn_selectTorusExecutable_clicked();
    void on_actionUnit_Converter_triggered();
    void on_tck_dustPhysics_toggled();
    void on_btn_openGridInputFile_clicked();
    void on_btn_openGridOutputFile_clicked();

    //Custom Slots for Undo Framework
    void spinBoxChange(double i);
    void textBoxChange(QString text);

    void on_btn_AddPopulation_clicked();

    void on_btn_RemovePopulation_clicked();

    void on_btn_SavePopulation_clicked();

    void on_lst_DustPopulations_itemSelectionChanged();

    void on_lst_Geometries_itemSelectionChanged();

    void on_actionUpdate_Torus_from_Source_triggered();

    void on_actionUpdate_Torus_from_Executable_triggered();

    void on_actionAbout_triggered();
    void on_tck_readGrid_toggled(bool  checked);
    void on_tck_writeGrid_toggled(bool checked);

    //SED Plotter Slots
    void on_actionSED_Viewer_triggered();
    void on_sedPlot_Destroyed();

protected:
    void dropEvent(QDropEvent *event);
    void dragEnterEvent(QDragEnterEvent *event);

private:
    Ui::MainWindow *ui;
    QUndoStack *undoStack;
    QUndoView *undoView;
    QString torusExecutable;
    bool enableUndoFramework;
    bool documentChanged;
    void closeEvent(QCloseEvent *event);
    void setupUndoFramework();
    bool saveDocument();
    QStringList getExpectedOutputFiles();
    QGridLayout *currentLayout;
    QList<QWidget *> geomWidgets;

    PSEDPlotter *activeSEDPlotter;
};

#endif // MAINWINDOW_H
