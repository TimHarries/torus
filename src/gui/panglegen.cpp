#include "panglegen.h"
#include "ui_panglegen.h"
#include <QMap>

PAngleGen::PAngleGen(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PAngleGen)
{
    ui->setupUi(this);
}

PAngleGen::~PAngleGen()
{
    delete ui;
}

void PAngleGen::getSettings(int *inc, double *minAngle, double *maxAngle, bool *cosSpaced) {
    *inc = ui->spn_AngleIncs->value();
    *minAngle = ui->spn_minAngle->value();
    *maxAngle = ui->spn_MaxAngle->value();
    *cosSpaced = ui->chk_cosSpaced->checkState();
}
