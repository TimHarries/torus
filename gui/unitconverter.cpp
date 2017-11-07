#include "unitconverter.h"
#include "ui_unitconverter.h"

UnitConverter::UnitConverter(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::UnitConverter)
{
    ui->setupUi(this);

    //Populate Combo Boxes
    QStringList distanceUnits;
    QStringList wavelengthUnits;

    distanceUnits << "Metres (m)" << "Centimetres (cm)" << "Parsecs (pc)" << "Lightyears (ly)" << "Torus Units (10^10cm)";
    wavelengthUnits << "Metres (m)" << "Centimetres (cm)" << "Microns (µm)" << "Nanometres (nm)" << "Angstroms (Å)";
    ui->com_DistanceFrom->addItems(distanceUnits);
    ui->com_DistanceTo->addItems(distanceUnits);
    ui->com_WavelengthFrom->addItems(wavelengthUnits);
    ui->com_WavelengthTo->addItems(wavelengthUnits);

    //Initialise Units Arrays
    this->distanceUnits = new QList<double>();
    this->wavelengthUnits = new QList<double>();
    *this->distanceUnits << 1E0 << 1E-2 << 3.08567758E16 << 9.4605284E15 << 1E8;
    *this->wavelengthUnits << 1E0 << 1E-2 << 1E-6 << 1E-9 << 1E-10;
}

UnitConverter::~UnitConverter()
{
    delete ui;
}

void UnitConverter::on_btn_Close_clicked()
{
    this->close();
}

void UnitConverter::on_btn_Convert_clicked()
{
    double fromUnit;
    double toUnit;

    switch(ui->tabWidget->currentIndex()) {
        case 0 :
            fromUnit = this->distanceUnits->at(ui->com_DistanceFrom->currentIndex());
            toUnit = this->distanceUnits->at(ui->com_DistanceTo->currentIndex());
            ui->txt_DistanceTo->setText(QString::number(ui->txt_DistanceFrom->text().toDouble() * (fromUnit / toUnit)));
            break;

        case 1:
            fromUnit = this->wavelengthUnits->at(ui->com_WavelengthFrom->currentIndex());
            toUnit = this->wavelengthUnits->at(ui->com_WavelengthTo->currentIndex());
            ui->txt_WavelengthTo->setText(QString::number(ui->txt_WavelengthFrom->text().toDouble() * (fromUnit / toUnit)));
            break;

         default :
            break;
    }
}
