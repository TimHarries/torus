#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    setAcceptDrops(true);
    QCoreApplication::setApplicationName("Parametor");
    setWindowTitle("Parametor");

    //Get Geometires and Populate
    QGroupBox *group = new QGroupBox;
    ui->lst_Geometries->addItems(this->params.geometries.keys());
    group->setLayout(ui->geomFormLayout);
    ui->scrollArea->setWidget(group);

    //Set MaximumValues
    ui->spn_ImageDistance->setMaximum(1e+20);
    ui->spn_GridSize->setMaximum(1e+20);

    //Setup the Atmospheres in the Combo Box (Might Move)
    ui->com_StarContFlux->addItem("Blackbody Curve", "blackbody");
    ui->com_StarContFlux->addItem("Kurucz Model Atmosphere", "kurucz");

    //Add Dust Types
    ui->com_DustGrainType->addItem("Draine and Lee Silicates", "sil_dl");
    ui->com_DustGrainType->addItem("Amorphous Carbon", "amc");
    ui->com_DustGrainType->addItem("Draine and Lee ortho graphite", "gr1_dl");
    ui->com_DustGrainType->addItem("Draine and Lee para graphite", "gr2_dl");
    ui->com_DustGrainType->addItem("Christophe Pinte interstellar grains", "pinteISM");
    ui->com_DustGrainType->addItem("Draine Silicates", "draine_sil");
    ui->com_DustGrainType->addItem("Amorphous Carbon (amc_zb)", "amc_zb");

    ui->com_gridDim->addItems(QStringList() << "1D" << "2D" << "3D");

    enableUndoFramework = false;
    this->setupUndoFramework();
    this->documentChanged = false;
    this->refreshUI();
    setWindowModified(false);

    //Setup the Pointers for Single Instance Windows
    activeSEDPlotter = 0x0;
}

MainWindow::~MainWindow()
{
    delete ui;
}

/*
 *Initialises the Undo Framework
 */
void MainWindow::setupUndoFramework() {
    if(this->enableUndoFramework) {
        //Undo Stack
        undoStack = new QUndoStack(this);
        undoView = new QUndoView(undoStack);
        undoView->setWindowTitle("History");
        undoView->setAttribute(Qt::WA_QuitOnClose, false);
        undoView->show();
    }

    //Append Undo / Redo Controls to Widgets
    QList<QDoubleSpinBox*> spins = this->findChildren<QDoubleSpinBox*>();
    foreach(QDoubleSpinBox *w, spins) {
        connect(w, SIGNAL(valueChanged(double)), this, SLOT(spinBoxChange(double)));
    }

    QList<QLineEdit*> texts = this->findChildren<QLineEdit*>();
    foreach(QLineEdit *w, texts) {
        connect(w, SIGNAL(textEdited(QString)), this, SLOT(textBoxChange(QString)));
    }
}

/*
 *Function responsible for pushing all values in the model to the UI
 */
void MainWindow::refreshUI(void) {
    //General Tab
    ui->txt_GridOutputFile->setText(this->params.stringParams.value("outputfile"));
    ui->txt_gridInputFile->setText(this->params.stringParams.value("inputfile"));
    ui->tck_images->setChecked(this->params.boolParams.value("image"));
    ui->tck_spectrum->setChecked(this->params.boolParams.value("spectrum"));
    ui->tck_isoScatter->setChecked(this->params.boolParams.value("iso_scatter"));
    ui->tck_writeGrid->setChecked(this->params.boolParams.value("writegrid"));
    ui->tck_readGrid->setChecked(this->params.boolParams.value("readgrid"));
    ui->tck_dustPhysics->setChecked(this->params.boolParams.value("dustphysics"));
    ui->txt_TorusExecutable->setText(this->torusExecutable);
    ui->spn_ImagePhotons->setValue(params.numParams.value("nphotons"));
    ui->spn_ImageDistance->setValue(params.numParams.value("distance"));
    ui->spn_AMRMinDepth->setValue(this->params.numParams.value("mindepthamr"));
    ui->spn_AMRMaxDepth->setValue(this->params.numParams.value("maxdepthamr"));
    ui->spn_GridSize->setValue(this->params.getNum("amrgridsize"));

    if(this->params.getBool("amr1d")) {
        ui->com_gridDim->setCurrentIndex(0);
    } else if(this->params.getBool("amr2d")) {
        ui->com_gridDim->setCurrentIndex(1);
    } else if (this->params.getBool("amr3d")) {
        ui->com_gridDim->setCurrentIndex(2);
    }

    //Geometries Tab
    for(int i = 0; i < ui->lst_Geometries->count(); i++) {
        QListWidgetItem *item = ui->lst_Geometries->item(i);
        if(item->text() == this->params.getString("geometry")) {
            ui->lst_Geometries->setCurrentItem(item);
        }
        if(ui->lst_Geometries->currentItem() == item) {
            this->renderGeometryUI();
        }
    }

    //Sources Tab
    int nsource = this->params.sources.count();
    int i;
    for(i = 1; i <= nsource; i++) {
        ui->lst_Sources->addItem(QString("Source " +  QString::number(i)));
    }

    //Dust Tab
    for(int i = 0; i < this->params.dusts.count(); i++) {
        ui->lst_DustPopulations->addItem(this->params.dusts.at(i).label);
    }

    //Spectrum Tab
    ui->chk_sedSIUnits->setChecked(this->params.getBool("sised"));
    ui->chk_sedWavLin->setChecked(this->params.getBool("sedwavlin"));
    ui->txt_sedNamePrefix->setText(this->params.getString("filename"));
    ui->spn_sedMinLambda->setValue(this->params.getNum("sedlammin"));
    ui->spn_sedMaxLambda->setValue(this->params.getNum("sedlammax"));

    /* Generate the Angles for the Angles Boxes */
    if(this->params.stringParams.contains("inclinations")) {
        QString incString = this->params.getString("inclinations");
        QStringList incs = incString.split(QRegExp(" "));

        QString itemStr;
        foreach(itemStr, incs) {
            QListWidgetItem *item = new QListWidgetItem(itemStr, ui->lst_sedIncList);
            item->setFlags(item->flags() | Qt::ItemIsEditable);
            ui->lst_sedPosList->addItem(item);
        }

    } else if (this->params.numParams.contains("firstinc") && this->params.numParams.contains("lastinc")) {
        double firstInc = this->params.getNum("firstinc");
        double lastInc = this->params.getNum("lastinc");
        int ninc = this->params.getNum("ninc");
        bool cosSpaced = this->params.getBool("sedcosspaced");

        QList<double> angles = this->params.getAngles(firstInc, lastInc, ninc, cosSpaced);
        double angle;

        foreach(angle, angles) {
            QListWidgetItem *item = new QListWidgetItem(QString::number(angle), ui->lst_sedIncList);
            item->setFlags(item->flags() | Qt::ItemIsEditable);
            ui->lst_sedIncList->addItem(item);
        }
    }

    if(this->params.stringParams.contains("posangs")) {
        QString incString = this->params.getString("posangs");
        QStringList incs = incString.split(QRegExp(" "));

        QString itemStr;
        foreach(itemStr, incs) {
            QListWidgetItem *item = new QListWidgetItem(itemStr, ui->lst_sedPosList);
            item->setFlags(item->flags() | Qt::ItemIsEditable);
            ui->lst_sedPosList->addItem(item);
        }

    } else if (this->params.numParams.contains("firstPA") && this->params.numParams.contains("lastPA")) {
        double firstPA = this->params.getNum("firstPA");
        double lastPA = this->params.getNum("lastPA");
        int ninc = this->params.getNum("ninc");
        bool cosSpaced = this->params.getBool("sedcosspaced");

        QList<double> angles = this->params.getAngles(firstPA, lastPA, ninc, cosSpaced);
        double angle;

        foreach(angle, angles) {
            QListWidgetItem *item = new QListWidgetItem(QString::number(angle), ui->lst_sedPosList);
            item->setFlags(item->flags() | Qt::ItemIsEditable);
            ui->lst_sedPosList->addItem(item);
        }
    }


    //Images Tab
    int nimage = this->params.numParams.value("nimage");
    for(i = 1; i <= nimage; i++) {
        ui->lst_Images->addItem(QString(this->params.images.at(i-1).name));
    }
}

/*
 *Responsible for clearing the UI
 */
void MainWindow::resetUI(void) {
    //General Tab
    ui->txt_GridOutputFile->clear();
    ui->txt_gridInputFile->clear();
    ui->tck_images->setChecked(false);
    ui->tck_spectrum->setChecked(false);
    ui->tck_isoScatter->setChecked(false);
    ui->tck_writeGrid->setChecked(false);
    ui->tck_readGrid->setChecked(false);
    ui->tck_dustPhysics->setChecked(false);
    ui->txt_TorusExecutable->clear();
    ui->tck_writeLogFile->setChecked(false);
    ui->spn_ImagePhotons->clear();
    ui->spn_ImageDistance->clear();
    ui->spn_AMRMinDepth->clear();
    ui->spn_AMRMaxDepth->clear();
    ui->spn_GridSize->clear();

    //Sources Tab
    ui->lst_Sources->clear();
    ui->txt_StartName->clear();
    ui->spn_StarMass->setValue(0);
    ui->spn_StarRadius->setValue(0);
    ui->spn_StarTemp->setValue(0);
    ui->spn_StarPosX->setValue(0);
    ui->spn_StarPosY->setValue(0);
    ui->spn_StarPosZ->setValue(0);

    //Dust Physics Tab
    ui->lst_DustPopulations->clear();
    ui->spn_DustMinSize->clear();
    ui->spn_DustMaxSize->clear();
    ui->spn_DustPowerLawIndex->clear();
    ui->txt_DustLabel->clear();
    ui->com_DustGrainType->setCurrentIndex(0);
    ui->spn_DustGrainFrac->clear();

    //SEDs Tab
    ui->chk_sedSIUnits->setChecked(false);
    ui->chk_sedWavLin->setChecked(false);
    ui->txt_sedNamePrefix->clear();
    ui->spn_sedMinLambda->clear();
    ui->spn_sedMaxLambda->clear();
    ui->lst_sedIncList->clear();
    ui->lst_sedPosList->clear();

    //Images Tab
    ui->lst_Images->clear();

    //This tells the program that it doesn't need to ask before closing
    this->documentChanged = false;
    setWindowModified(false);
}

/*
 *Tells the app whether it can close the current document
 */
bool MainWindow::canCloseCurrentDocument() {
    if(this->documentChanged) {
        bool response;
        QMessageBox saveDialog;
        saveDialog.setText("Save Changes?");
        saveDialog.setInformativeText("Would you like to save changes to the current document?");
        saveDialog.setStandardButtons(QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
        int ret = saveDialog.exec();

        switch (ret) {
        case QMessageBox::Save:
            response = this->saveDocument();
            return response;
            break;
        case QMessageBox::Discard :
            return true;
            break;
        case QMessageBox::Cancel :
            return false;
            break;
        default:
            break;
        }
    }

    return true;
}

/*
 *Saves all UI changes back to the model
 */
void MainWindow::updateModel() {
    //General Tab
    this->params.setString("outputfile", ui->txt_GridOutputFile->text());
    this->params.setString("inputfile", ui->txt_gridInputFile->text());
    this->params.setBool("image", ui->tck_images->isChecked());
    this->params.setBool("spectrum", ui->tck_spectrum->isChecked());
    this->params.setBool("iso_scatter", ui->tck_isoScatter->isChecked());
    this->params.setBool("writegrid", ui->tck_writeGrid->isChecked());
    this->params.setBool("readgrid", ui->tck_readGrid->isChecked());
    this->params.setBool("dustphysics", ui->tck_dustPhysics->isChecked());
    this->params.setNum("nphotons", ui->spn_ImagePhotons->value());
    this->params.setNum("distance", ui->spn_ImageDistance->value());
    this->params.setNum("mindepthamr", ui->spn_AMRMinDepth->value());
    this->params.setNum("maxdepthamr", ui->spn_AMRMaxDepth->value());
    this->params.setNum("amrgridsize", ui->spn_GridSize->value());

    switch(ui->com_gridDim->currentIndex()){
        case 0 :
            params.setBool("amr1d", true);
            params.setBool("amr2d", false);
            params.setBool("amr3d", false);
            break;
        case 1 :
            params.setBool("amr1d", false);
            params.setBool("amr2d", true);
            params.setBool("amr3d", false);
            break;
        default :
            params.setBool("amr1d", false);
            params.setBool("amr2d", false);
            params.setBool("amr3d", true);
            break;
    }

    //Geometry File
    //This is going to get complex, so we're going to hand off to a dedicated method
    this->updateGeometry();

    //Sources Tab
    //Nothing here as all done in the tab

    //SEDs Tab
    this->params.setBool("sised", ui->chk_sedSIUnits->isChecked());
    this->params.setBool("sedwavlin", ui->chk_sedWavLin->isChecked());
    this->params.setString("filename", ui->txt_sedNamePrefix->text());
    this->params.setNum("sedlammin", ui->spn_sedMinLambda->value());
    this->params.setNum("sedlammax", ui->spn_sedMaxLambda->value());

    QList<double> incs;
    QList<double> pas;
    //Work out the way we're going to represent angles here
    for(int i = 0; i < ui->lst_sedIncList->count(); i++){
        incs.append(ui->lst_sedIncList->item(i)->text().toDouble());
    }
    for(int i = 0; i < ui->lst_sedPosList->count(); i++){
        pas.append(ui->lst_sedPosList->item(i)->text().toDouble());
    }

    this->params.setAngles("inclinations", incs);
    this->params.setAngles("posangs", pas);
    if(incs.count() >= pas.count()) {
        this->params.setNum("ninc", incs.count());
    } else {
        this->params.setNum("ninc", pas.count());
    }

    //Images Tab
    //Nothing here as all done in the tab
}

/*
 *Checks for a Torus Execurable in the current working directory
 */
void MainWindow::checkForTorusExecutable(QString filePath) {
    QDir directory(QFileInfo(filePath).path());
    QStringList filters(QStringList() << "torus.gfortran" << "torus.exe" << "torus");

    directory.setNameFilters(filters);
    QFileInfoList fList = directory.entryInfoList();

    if(fList.count() > 0) {
        this->torusExecutable = fList.at(0).filePath();
        qDebug() << "Torus Executable: " << this->torusExecutable;
        ui->txt_TorusExecutable->setText(this->torusExecutable);
        ui->actionCheck_Paramaters->setEnabled(true);
    }
}

/*
 *Gets the initial path for the file dialogues
 */
QString MainWindow::getFileDialogPath() {
    QString path;

    if(this->params.openFile != "") {
        //Return the folder
        path = QFileInfo(this->params.openFile).dir().path();
        qDebug() << path;
    } else {
        path = QDir::homePath();
    }

    return path;
}

//-------------------Slots to Deal with Undo Framework---------------------------------
void MainWindow::spinBoxChange(double i) {
    Q_UNUSED(i);
    /*
    QUndoCommand *undo = new QSpinUndoCommand(ui->spn_ImageAspect, i);
    undo->setText("Value set to " + QString::number(i) + " on object ");
    this->undoStack->push(undo);
    */
    this->documentChanged = true;
    setWindowModified(true);
}

void MainWindow::textBoxChange(QString text) {
    Q_UNUSED(text);
    /*
    QUndoCommand *undo = new QUndoCommand();
    undo->setText("Text set to " + text);
    this->undoStack->push(undo);
    */
    this->documentChanged = true;
    setWindowModified(true);
}

bool MainWindow::saveDocument() {
    QString fileName;
    if(this->params.openFile == "") {
        fileName = QFileDialog::getSaveFileName(this, tr("Save Parameter File"), this->getFileDialogPath());
    } else {
        //Save the file straight
        fileName = this->params.openFile;
    }

    if(fileName != "") {
        this->updateModel();
        this->params.saveToFile(fileName);
        this->documentChanged = false;
        setWindowModified(false);
        this->setWindowFilePath(fileName);
        this->setWindowTitle("Parametor - " + fileName);
        this->checkForTorusExecutable(fileName);
        qDebug() << "Filename: " << fileName;
        return true;
    } else {
        return false;
    }
    return false;
}

//-------------------Event Handlers Follow---------------------------------

void MainWindow::on_actionClose_triggered() {
    qApp->activeWindow()->close();
}

void MainWindow::dragEnterEvent(QDragEnterEvent *event) {
    if (event->mimeData()->hasFormat("text/uri-list"))
        event->acceptProposedAction();
}

void MainWindow::dropEvent(QDropEvent *event) {
    const QMimeData *mimeData = event->mimeData();
    // check for our needed mime type, here a file or a list of files
       if (mimeData->hasUrls())
       {
           QStringList pathList;
           QList<QUrl> urlList = mimeData->urls();

           // extract the local paths of the files
           for (int i = 0; i < urlList.size() && i < 32; ++i)
           {
               pathList.append(urlList.at(i).toLocalFile());
           }

           // call a function to open the files
           QString fileName = pathList.last();
           if(fileName != "") {
               //todo: Check that the file is a parameter.dat file
               QStringList segments = fileName.split("/");
               QString file = segments.last();

               if(file.contains(".dat")) {
                   this->resetUI();
                   this->params.clear();
                   this->params.loadFile(pathList.last());
                   this->checkForTorusExecutable(pathList.last());
                   this->setWindowFilePath(fileName);
                   this->setWindowTitle("Parametor - " + fileName);
                   this->refreshUI();
               }
           }
       }
}


void MainWindow::on_actionOpen_triggered()
{   
    if(this->canCloseCurrentDocument()) {
        QString fileName = QFileDialog::getOpenFileName(this, tr("Open Parameter File"), this->getFileDialogPath());
        if(fileName != "") {
            //todo: Check that the file is a parameter.dat file
            QStringList segments = fileName.split("/");
            QString file = segments.last();

            if(file.contains(".dat")) {
                this->resetUI();
                this->params.clear();
                this->params.loadFile(fileName);
                this->checkForTorusExecutable(fileName);
                this->setWindowFilePath(fileName);
                this->setWindowTitle("Parametor - " + fileName);
                this->refreshUI();
            } else {
                QMessageBox msgBox;
                msgBox.setText("This doesn't appear to be a valid parameters file");
                msgBox.exec();
            }
        }
    }
}

void MainWindow::on_actionSave_triggered()
{
    this->saveDocument();
}

void MainWindow::on_actionSaveAs_triggered()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save Parameter File As..."), this->getFileDialogPath());

    if(fileName != "") {
        this->updateModel();
        this->params.saveToFile(fileName);
        this->checkForTorusExecutable(fileName);
        this->documentChanged = false;
        setWindowModified(false);
        this->setWindowFilePath(fileName);
        this->setWindowTitle("Parametor - " + fileName);
    }
}

void MainWindow::on_btn_AddSource_clicked()
{
    PSourceData source;

    //Testing Code - Remove later
    source.name = QString("Source " + QString::number(this->params.sources.count() + 1));
    this->params.sources.insert(params.sources.count(), source);
    ui->lst_Sources->addItem(source.name);

}

void MainWindow::on_btn_RemoveSource_clicked()
{
    //Code to remove object from sources array
    if(ui->lst_Sources->count() > 1 && ui->lst_Sources->currentIndex().isValid()) {
        int itid = ui->lst_Sources->currentIndex().row();
        QListWidgetItem *item = ui->lst_Sources->takeItem(itid);
        this->params.sources.removeAt(itid);
        delete item;
    }
}

void MainWindow::on_btn_SaveSource_clicked()
{
    if(ui->lst_Sources->currentIndex().isValid()) {
        int id = ui->lst_Sources->currentIndex().row();
        PSourceData *source = &this->params.sources[id];
        source->contflux = ui->com_StarContFlux->itemData(ui->com_StarContFlux->currentIndex()).toString();
        source->mass = ui->spn_StarMass->value();
        source->radius = ui->spn_StarRadius->value();
        source->name = ui->txt_StartName->text();
        source->teff = ui->spn_StarTemp->value();

        //Position Vector Update
        source->position[0] = ui->spn_StarPosX->value();
        source->position[1] = ui->spn_StarPosY->value();
        source->position[2] = ui->spn_StarPosZ->value();

        ui->lst_Sources->currentItem()->setText(source->name);
    } else {
        MainWindow::on_btn_AddSource_clicked();

        int id = 0;
        PSourceData *source = &this->params.sources[id];
        source->contflux = ui->com_StarContFlux->itemData(ui->com_StarContFlux->currentIndex()).toString();
        source->mass = ui->spn_StarMass->value();
        source->radius = ui->spn_StarRadius->value();
        source->name = ui->txt_StartName->text();
        source->teff = ui->spn_StarTemp->value();

        //Position Vector Update
        source->position[0] = ui->spn_StarPosX->value();
        source->position[1] = ui->spn_StarPosY->value();
        source->position[2] = ui->spn_StarPosZ->value();

        ui->lst_Sources->setCurrentRow(0);
    }
}

void MainWindow::on_lst_Sources_itemSelectionChanged()
{
    int sid = ui->lst_Sources->currentIndex().row();
    int index = -1;

    if(sid < this->params.sources.count()) {
        PSourceData source = this->params.sources.at(sid);
        ui->txt_StartName->setText(source.name);
        ui->spn_StarMass->setValue(source.mass);
        ui->spn_StarRadius->setValue(source.radius);
        ui->spn_StarTemp->setValue(source.teff);

        index = ui->com_StarContFlux->findData(source.contflux);
        if (index != -1) {
           ui->com_StarContFlux->setCurrentIndex(index);
        }

        //Set Vectors
        ui->spn_StarPosX->setValue(source.position.at(0));
        ui->spn_StarPosY->setValue(source.position.at(1));
        ui->spn_StarPosZ->setValue(source.position.at(2));
    }
}

void MainWindow::on_lst_Images_itemSelectionChanged()
{
    int sid = ui->lst_Images->currentIndex().row();

    if(sid < this->params.images.count()) {
        PImageData image = this->params.images.at(sid);
        ui->txt_ImageName->setText(image.name);
        ui->spn_ImagePixels->setValue(image.npixels);
        ui->spn_ImageInclination->setValue(image.inclination);
        ui->spn_ImagePA->setValue(image.positionangle);
        ui->spn_ImageWavelength->setValue(image.lambdaimage);
        ui->spn_ImageAspect->setValue(image.aspect);
        ui->spn_ImageCentreX->setValue(image.imagecentrex);
        ui->spn_ImageCentreY->setValue(image.imagecentrey);
    }
}

void MainWindow::on_tck_spectrum_toggled(bool checked)
{
    ui->tab_SED->setEnabled(checked);
    if(ui->tck_images->isChecked() || ui->tck_spectrum->isChecked()) {
        ui->spn_ImagePhotons->setEnabled(true);
        ui->spn_ImageDistance->setEnabled(true);
    } else {
        ui->spn_ImagePhotons->setEnabled(false);
        ui->spn_ImageDistance->setEnabled(false);
    }
}

void MainWindow::on_tck_images_toggled(bool checked)
{
    ui->tab_Image->setEnabled(checked);
    if(ui->tck_images->isChecked() || ui->tck_spectrum->isChecked()) {
        ui->spn_ImagePhotons->setEnabled(true);
        ui->spn_ImageDistance->setEnabled(true);
    } else {
        ui->spn_ImagePhotons->setEnabled(false);
        ui->spn_ImageDistance->setEnabled(false);
    }
}

void MainWindow::on_actionNew_triggered()
{
    if(this->canCloseCurrentDocument()) {
        this->params.clear();
        this->resetUI();
        this->refreshUI();
        this->documentChanged = false;
        setWindowModified(false);
    }
}

void MainWindow::on_btn_SaveImage_clicked()
{
    if(ui->lst_Images->currentIndex().isValid()) {
        int id = ui->lst_Images->currentIndex().row();
        PImageData *image = &this->params.images[id];
        image->name = ui->txt_ImageName->text();
        image->npixels = ui->spn_ImagePixels->value();
        image->inclination = ui->spn_ImageInclination->value();
        image->positionangle = ui->spn_ImagePA->value();
        image->lambdaimage = ui->spn_ImageWavelength->value();
        image->aspect = ui->spn_ImageAspect->value();

        //Position Vector Update
        image->imagecentrex = ui->spn_ImageCentreX->value();
        image->imagecentrey = ui->spn_ImageCentreY->value();

        ui->lst_Images->currentItem()->setText(image->name);
    } else {
        MainWindow::on_btn_AddImage_clicked();

        int id = 0;
        PImageData *image = &this->params.images[id];
        image->name = ui->txt_ImageName->text();
        image->nphotons = ui->spn_ImagePhotons->value();
        image->npixels = ui->spn_ImagePixels->value();
        image->inclination = ui->spn_ImageInclination->value();
        image->positionangle = ui->spn_ImagePA->value();
        image->lambdaimage = ui->spn_ImageWavelength->value();
        image->aspect = ui->spn_ImageAspect->value();

        //Position Vector Update
        image->imagecentrex = ui->spn_ImageCentreX->value();
        image->imagecentrey = ui->spn_ImageCentreY->value();

        ui->lst_Images->setCurrentRow(0);
    }
}

void MainWindow::on_btn_AddImage_clicked()
{
    PImageData image;

    //Testing Code - Remove later
    image.name = QString("Image " + QString::number(this->params.images.count() + 1));
    this->params.images.insert(params.images.count(), image);
    ui->lst_Images->addItem(image.name);
}

void MainWindow::on_btn_RemoveImage_clicked()
{
    //Code to remove object from images array
    if(ui->lst_Images->count() > 1 && ui->lst_Images->currentIndex().isValid()) {
        int itid = ui->lst_Images->currentIndex().row();
        QListWidgetItem *item = ui->lst_Images->takeItem(itid);
        this->params.images.removeAt(itid);
        delete item;
    }
}

void MainWindow::on_btn_sedGenInc_clicked()
{
    PAngleGen dialog;
    int inc;
    double minAngle;
    double maxAngle;
    bool cosSpaced;

    if(dialog.exec() == QDialog::Accepted) {
        dialog.getSettings(&inc, &minAngle, &maxAngle, &cosSpaced);
        QList<double> angles = this->params.getAngles(minAngle, maxAngle, inc, cosSpaced);

        double angle;
        foreach(angle, angles) {
            QListWidgetItem *item = new QListWidgetItem(QString::number(angle), ui->lst_sedIncList);
            item->setFlags(item->flags() | Qt::ItemIsEditable);
            ui->lst_sedIncList->addItem(item);
        }
    }
}

void MainWindow::on_spn_sedAddInc_clicked()
{
    QListWidgetItem *item = new QListWidgetItem(QString::number(0), ui->lst_sedIncList);
    item->setFlags(item->flags() | Qt::ItemIsEditable);
    ui->lst_sedIncList->addItem(item);
}

void MainWindow::on_spn_sedRemoveInc_clicked()
{
    if(ui->lst_sedIncList->currentIndex().isValid()) {
        QListWidgetItem *item;
        item = ui->lst_sedIncList->takeItem(ui->lst_sedIncList->currentIndex().row());
        delete item;
    }
}

void MainWindow::on_btn_sedGenPos_clicked()
{
    PAngleGen dialog;
    int inc;
    double minAngle;
    double maxAngle;
    bool cosSpaced;

    if(dialog.exec() == QDialog::Accepted) {
        dialog.getSettings(&inc, &minAngle, &maxAngle, &cosSpaced);
        QList<double> angles = this->params.getAngles(minAngle, maxAngle, inc, cosSpaced);

        double angle;
        foreach(angle, angles) {
            QListWidgetItem *item = new QListWidgetItem(QString::number(angle), ui->lst_sedPosList);
            item->setFlags(item->flags() | Qt::ItemIsEditable);
            ui->lst_sedPosList->addItem(item);
        }
    }
}

void MainWindow::on_spn_sedAddPos_clicked()
{
    QListWidgetItem *item = new QListWidgetItem(QString::number(0), ui->lst_sedPosList);
    item->setFlags(item->flags() | Qt::ItemIsEditable);
    ui->lst_sedPosList->addItem(item);
}

void MainWindow::on_spn_sedRemovePos_clicked()
{
    if(ui->lst_sedPosList->currentIndex().isValid()) {
        QListWidgetItem *item;
        item = ui->lst_sedPosList->takeItem(ui->lst_sedPosList->currentIndex().row());
        delete item;
    }
}
/*
 *------------------------------------------------------
 * Clipboard Manipulation
 *------------------------------------------------------
 */
void MainWindow::on_actionPaste_triggered()
{
    QWidget *object = this->focusWidget();
    QString widgetType = typeid(object).name();
    QClipboard *clip = QApplication::clipboard();

    if(widgetType == "QLineEdit") {
        ((QLineEdit*)object)->setText(clip->text());
    } else if(widgetType == "QDoubleSpinBox") {
        QString value = clip->text();
        if(value.toDouble()) {
            ((QDoubleSpinBox*)object)->setValue(value.toDouble());
        }
    }
}

void MainWindow::on_actionCopy_triggered()
{
}

void MainWindow::on_actionCut_triggered()
{
}

/*
 *------------------------------------------------------
 * Deals with Executing Torus
 *------------------------------------------------------
 */
void MainWindow::on_actionRun_Torus_triggered()
{
    this->torusExecutable = ui->txt_TorusExecutable->text();
    QString fileDir = QFileInfo(this->params.openFile).dir().canonicalPath();
    torusconsole *torcon = new torusconsole(0 , this->torusExecutable, fileDir, ui->tck_writeLogFile->isChecked(), torusconsole::DefaultRunMode, this->getExpectedOutputFiles());
    torcon->show();
}

void MainWindow::on_actionCheck_Paramaters_triggered()
{
    this->torusExecutable = ui->txt_TorusExecutable->text();
    QString fileDir = QFileInfo(this->params.openFile).dir().canonicalPath();
    torusconsole *torcon = new torusconsole(0 , this->torusExecutable, fileDir, ui->tck_writeLogFile->isChecked(), torusconsole::CheckRunMode, QStringList(""));
    torcon->show();
}

QStringList MainWindow::getExpectedOutputFiles() {
    QStringList files;
    QString workingDir = QFileInfo(this->params.openFile).dir().path();

    //Get Images
    if(this->params.getBool("image")) {
        for(int i = 1; i <= this->params.images.count(); i++) {
            files.append(this->params.getString("imagefile" + QString::number(i)));
        }
    }

    //GetSEDs
    if(this->params.getBool("spectrum")) {
        int posidx = 0;
        int incidx = 0;

        for(int i = 0; i < this->params.getNum("ninc"); i++) {
            QString inc;
            QString pos;

            if(i < ui->lst_sedIncList->count()) {
                inc = inc.sprintf("%03g", round(ui->lst_sedIncList->item(incidx)->text().toFloat()));
            } else { inc = inc.sprintf("%03d", 0); }

            if(i < ui->lst_sedPosList->count()) {
                pos = pos.sprintf("%03g", round(ui->lst_sedPosList->item(posidx)->text().toFloat()));
            } else { pos = pos.sprintf("%03d", 0); }

            qDebug() << inc << " - " << pos;

            if(pos != "000") {
                files.append(this->params.getString("filename") + "_inc" + inc + "_PA" + pos  + ".dat");
            } else {
                files.append(this->params.getString("filename") + "_inc" + inc  + ".dat");
            }

            posidx++;
            incidx++;
        }
    }

    return files;
}

/*
 *------------------------------------------------------
 * This is to deal with Torus Executable Selection
 *------------------------------------------------------
 */
void MainWindow::on_btn_selectTorusExecutable_clicked()
{
    QString fileName;
    fileName = QFileDialog::getOpenFileName(this, "Select the Torus Executable", this->getFileDialogPath());

    if(fileName != "") {
        ui->txt_TorusExecutable->setText(fileName);
        this->torusExecutable = fileName;
    }
}

void MainWindow::on_actionUnit_Converter_triggered()
{
    UnitConverter *unitc = new UnitConverter();
    unitc->show();
}

void MainWindow::closeEvent(QCloseEvent *event) {
    bool canClose = this->canCloseCurrentDocument();
    if(canClose) {
        event->accept();
    } else {
        event->ignore();
    }
}

void MainWindow::on_tck_dustPhysics_toggled()
{
    if(ui->tck_dustPhysics->isChecked()) {
        ui->tab_Dust->setEnabled(true);
    } else {
        ui->tab_Dust->setEnabled(false);
    }
}

void MainWindow::on_btn_openGridInputFile_clicked()
{
    QString fileName;
    fileName = QFileDialog::getOpenFileName(this, "Select the Grid Input File", this->getFileDialogPath());

    if(fileName != "") {
        ui->txt_gridInputFile->setText(fileName);
    }
}

void MainWindow::on_btn_openGridOutputFile_clicked()
{
    QString fileName;
    fileName = QFileDialog::getSaveFileName(this, "Select the Grid Output File", this->getFileDialogPath());

    if(fileName != "") {
        ui->txt_GridOutputFile->setText(fileName);
    }
}

/*
 *------------------------------------------------------
 * Event Handlers for Dust Physics UI
 *------------------------------------------------------
 */
void MainWindow::on_btn_AddPopulation_clicked()
{
    PDustData dust;

    //Testing Code - Remove later
    dust.label = QString("Dust Population " + QString::number(this->params.dusts.count() + 1));
    this->params.dusts.insert(params.dusts.count(), dust);
    ui->lst_DustPopulations->addItem(dust.label);
}

void MainWindow::on_btn_RemovePopulation_clicked()
{
    //Code to remove object from dusts array
    if(ui->lst_DustPopulations->count() > 1 && ui->lst_DustPopulations->currentIndex().isValid()) {
        int itid = ui->lst_DustPopulations->currentIndex().row();
        QListWidgetItem *item = ui->lst_DustPopulations->takeItem(itid);
        this->params.dusts.removeAt(itid);
        delete item;
    }
}

void MainWindow::on_btn_SavePopulation_clicked()
{
    if(ui->lst_DustPopulations->currentIndex().isValid()) {
        int id = ui->lst_DustPopulations->currentIndex().row();
        PDustData *dust = &this->params.dusts[id];
        dust->label = ui->txt_DustLabel->text();
        dust->amin = ui->spn_DustMinSize->value();
        dust->amax = ui->spn_DustMaxSize->value();
        dust->qdist = ui->spn_DustPowerLawIndex->value();
        dust->a0 = ui->spn_DustGrainScaleLength->value();
        dust->pdist = ui->spn_DustExpCutoffExp->value();
        dust->grainfrac = ui->spn_DustGrainFrac->value();
        dust->graintype = ui->com_DustGrainType->itemData(ui->com_DustGrainType->currentIndex()).toString();

        ui->lst_DustPopulations->currentItem()->setText(dust->label);
    } else {
        MainWindow::on_btn_AddImage_clicked();

        int id = 0;
        PDustData *dust = &this->params.dusts[id];
        dust->label = ui->txt_DustLabel->text();
        dust->amin = ui->spn_DustMinSize->value();
        dust->amax = ui->spn_DustMaxSize->value();
        dust->qdist = ui->spn_DustPowerLawIndex->value();
        dust->a0 = ui->spn_DustGrainScaleLength->value();
        dust->pdist = ui->spn_DustExpCutoffExp->value();
        dust->graintype = ui->com_DustGrainType->itemData(ui->com_DustGrainType->currentIndex()).toString();

        ui->lst_Images->setCurrentRow(0);
        ui->lst_DustPopulations->currentItem()->setText(dust->label);
    }
}

void MainWindow::on_lst_DustPopulations_itemSelectionChanged()
{
    int sid = ui->lst_DustPopulations->currentIndex().row();

    if(sid < this->params.dusts.count()) {
        PDustData dust = this->params.dusts.at(sid);
        ui->txt_DustLabel->setText(dust.label);
        ui->spn_DustMinSize->setValue(dust.amin);
        ui->spn_DustMaxSize->setValue(dust.amax);
        ui->spn_DustPowerLawIndex->setValue(dust.qdist);
        ui->spn_DustGrainScaleLength->setValue(dust.a0);
        ui->spn_DustExpCutoffExp->setValue(dust.pdist);
        ui->spn_DustGrainFrac->setValue(dust.grainfrac);

        int index = ui->com_DustGrainType->findData(dust.graintype);
        if (index != -1) {
           ui->com_DustGrainType->setCurrentIndex(index);
        }
    }
}

//Clever UI Stuff to Follow!
void MainWindow::on_lst_Geometries_itemSelectionChanged()
{
    this->renderGeometryUI();
}

void MainWindow::renderGeometryUI() {

    //Remove Current For
    this->geomWidgets.clear();
    while(QLayoutItem *deletionItem = ui->geomFormLayout->takeAt(0)) {
        deletionItem->widget()->hide();
        delete deletionItem->widget();
        delete deletionItem;
    }

    QString geomSelect = ui->lst_Geometries->item(ui->lst_Geometries->currentIndex().row())->text();
    QList<PParameterDefinition> geomDefs = this->params.geometries.value(geomSelect);

    if(!geomDefs.isEmpty()) {
        for(int i = 0; i < geomDefs.count(); i++) {
            PParameterDefinition def = geomDefs.at(i);

            QQueue<QLineEdit *> linedit;
            QQueue<QDoubleSpinBox *> spins;
            QQueue<QCheckBox *> checks;
            QQueue<QLabel *> labels;

            if(def.type == "string") {
                QLineEdit *w;
                linedit.enqueue(new QLineEdit());
                w = linedit.last();
                w->setObjectName(def.name);
                if(this->params.stringParams.contains(def.name)) {
                    w->setText(this->params.getString(def.name));
                } else {
                    w->setText(def.defaultValue);
                }

                ui->geomFormLayout->addRow(def.label, w);
                this->geomWidgets.append(w);

            } else if(def.type == "double") {
                QDoubleSpinBox *w;
                spins.enqueue(new QDoubleSpinBox());
                w = spins.last();
                w->setObjectName(def.name);
                w->setMinimum(def.min);
                w->setMaximum(def.max);
                w->setDecimals(def.decimals);
                if(this->params.numParams.contains(def.name)) {
                    w->setValue(this->params.getNum(def.name));
                } else {
                    w->setValue(def.defaultValue.toDouble());
                }

                ui->geomFormLayout->addRow(def.label, w);
                this->geomWidgets.append(w);

            } if(def.type == "int") {
                QDoubleSpinBox *w;
                spins.enqueue(new QDoubleSpinBox());
                w = spins.last();
                w->setObjectName(def.name);
                w->setMinimum(def.min);
                w->setMaximum(def.max);
                w->setDecimals(0);
                if(this->params.numParams.contains(def.name)) {
                    w->setValue(this->params.getNum(def.name));
                } else {
                    w->setValue(def.defaultValue.toInt());
                }

                ui->geomFormLayout->addRow(def.label, w);
                this->geomWidgets.append(w);

            } if(def.type == "bool") {
                QCheckBox *w;
                checks.enqueue(new QCheckBox());
                w = checks.last();
                w->setObjectName(def.name);
                if(this->params.boolParams.contains(def.name)) {
                    w->setChecked(this->params.getBool(def.name));
                } else {
                    w->setChecked(def.defaultValue.toInt());
                }
                ui->geomFormLayout->addRow(def.label, w);
                this->geomWidgets.append(w);
            } if(def.type == "heading") {
                QLabel *w;
                QFont font;
                font.setFamily(font.defaultFamily());
                font.setPointSize(16);

                labels.enqueue(new QLabel());
                w = labels.last();
                w->setText(def.label);
                w->setFont(font);
                ui->geomFormLayout->addRow(w);
                this->geomWidgets.append(w);
            }
        }
    }
}


void MainWindow::updateGeometry() {
    //Remove Existing Params
    this->params.geomParams.clear();

    //Write the Geometry Variable from the selected geometry
    this->params.setString("geometry", this->ui->lst_Geometries->currentItem()->text());

    //Loop through the UI Elements and Save to Geometry
    for(int i = 0; i < this->geomWidgets.count(); i++) {
        QWidget *obj = geomWidgets.at(i);

        if(QString(obj->metaObject()->className()) == "QLineEdit") {
            this->params.geomParams.insert(obj->objectName(), ((QLineEdit*)obj)->text());

        } else if(QString(obj->metaObject()->className()) == "QDoubleSpinBox") {
            this->params.geomParams.insert(obj->objectName(), QString::number(((QDoubleSpinBox*)obj)->value(), 'g', 5));

        } else if(QString(obj->metaObject()->className()) == "QCheckBox") {
            this->params.geomParams.insert(obj->objectName(), this->params.boolToStr(((QCheckBox*)obj)->isChecked()));
        }
    }
}

void MainWindow::on_actionUpdate_Torus_from_Source_triggered()
{
    PTorusUpdateWidget *w = new PTorusUpdateWidget(PTorusUpdateWidget::UpdateFromTar);
    w->exec();
}


void MainWindow::on_actionUpdate_Torus_from_Executable_triggered()
{
    PTorusUpdateWidget *w = new PTorusUpdateWidget(PTorusUpdateWidget::UpdateFromBin);
    w->exec();
}

void MainWindow::on_actionAbout_triggered()
{
    PAboutDialog *box = new PAboutDialog();
    box->exec();
}

void MainWindow::on_tck_readGrid_toggled(bool checked)
{
    ui->txt_gridInputFile->setEnabled(checked);
    ui->btn_openGridInputFile->setEnabled(checked);
}


void MainWindow::on_tck_writeGrid_toggled(bool checked)
{
    ui->txt_GridOutputFile->setEnabled(checked);
    ui->btn_openGridOutputFile->setEnabled(checked);
}

void MainWindow::on_actionSED_Viewer_triggered()
{
    //Get the init directory
    QString workingDir;
    if(this->params.openFile != "") {
        workingDir = QFileInfo(this->params.openFile).dir().absolutePath();
    } else {
        workingDir = qApp->applicationDirPath();
    }

    if(this->activeSEDPlotter == 0x0) {
        this->activeSEDPlotter = new PSEDPlotter(workingDir, this->ui->txt_sedNamePrefix->text());
        activeSEDPlotter->show();
        connect(activeSEDPlotter, SIGNAL(isGoingToClose(void)), this, SLOT(on_sedPlot_Destroyed(void)));
    } else {
        activeSEDPlotter->raise();
        activeSEDPlotter->activateWindow();
        activeSEDPlotter->setWorkingDirectory(workingDir);
        activeSEDPlotter->setFilePrefix(ui->txt_sedNamePrefix->text());
        activeSEDPlotter->refreshSEDFiles();
    }
}

void MainWindow::on_sedPlot_Destroyed() {
    delete this->activeSEDPlotter;
    this->activeSEDPlotter = 0x0;
}
