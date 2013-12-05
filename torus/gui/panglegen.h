#ifndef PANGLEGEN_H
#define PANGLEGEN_H

#include <QDialog>

namespace Ui {
class PAngleGen;
}

class PAngleGen : public QDialog
{
    Q_OBJECT
    
public:
    explicit PAngleGen(QWidget *parent = 0);
    ~PAngleGen();

    void getSettings(int *inc, double *minAngle, double *maxAngle, bool *cosSpaced);
    
private:
    Ui::PAngleGen *ui;
};

#endif // PANGLEGEN_H
