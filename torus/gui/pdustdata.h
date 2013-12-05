#ifndef PDUSTDUST_H
#define PDUSTDUST_H
#include <QString>

class PDustData
{
public:
    PDustData();

    QString label;
    double amin;
    double amax;
    double qdist;
    double grainfrac;
    double a0;
    double pdist;
    QString graintype;
};

#endif // PDUSTDUST_H
