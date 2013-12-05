#ifndef PABOUTDIALOG_H
#define PABOUTDIALOG_H

#include <QDialog>

namespace Ui {
class PAboutDialog;
}

class PAboutDialog : public QDialog
{
    Q_OBJECT
    
public:
    explicit PAboutDialog(QWidget *parent = 0);
    ~PAboutDialog();
    
private:
    Ui::PAboutDialog *ui;
};

#endif // PABOUTDIALOG_H
