#ifndef UNITCONVERTER_H
#define UNITCONVERTER_H

#include <QWidget>

namespace Ui {
class UnitConverter;
}

class UnitConverter : public QWidget
{
    Q_OBJECT
    
public:
    explicit UnitConverter(QWidget *parent = 0);
    ~UnitConverter();
    QList<double> *distanceUnits;
    QList<double> *wavelengthUnits;
    
private slots:
    void on_btn_Close_clicked();
    void on_btn_Convert_clicked();

private:
    Ui::UnitConverter *ui;
};

#endif // UNITCONVERTER_H
