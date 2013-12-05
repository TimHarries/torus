#include "paboutdialog.h"
#include "ui_paboutdialog.h"
#include <QGraphicsPixmapItem>
#include <QGraphicsScene>

PAboutDialog::PAboutDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PAboutDialog)
{
    ui->setupUi(this);
    this->setFixedSize(this->width(), this->height());
}

PAboutDialog::~PAboutDialog()
{
    delete ui;
}
