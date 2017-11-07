#ifndef PSOURCEDATA_H
#define PSOURCEDATA_H
#include <QString>
#include <QVector>

class PSourceData
{
public:
    PSourceData();

    //Properties
    QString name;
    bool pointsource;
    double radius;
    double teff;
    double mass;
    QString contflux;
    double probsource;
    QVector<double> position;
    QVector<double> velocity;

    //Position
    double z;
    double y;
    double x;

    //Velocity
    double vx;
    double vy;
    double vz;

};

#endif // PSOURCEDATA_H
