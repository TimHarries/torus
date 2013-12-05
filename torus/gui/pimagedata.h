#ifndef PIMAGEDATA_H
#define PIMAGEDATA_H
#include <QString>

class PImageData
{
public:
    PImageData();

    //Variables
    QString name;
    int nphotons;
    int npixels;
    double inclination;
    double positionangle;
    double lambdaimage;
    double aspect;
    bool polimage;
    QString imageaxisunits;

    double imagecentrex;
    double imagecentrey;
};

#endif // PIMAGEDATA_H
