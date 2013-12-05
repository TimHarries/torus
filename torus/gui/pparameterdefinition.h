#ifndef PPARAMETERDEFINITION_H
#define PPARAMETERDEFINITION_H

#include <QString>

class PParameterDefinition
{
public:
    PParameterDefinition();

    QString name;
    QString defaultValue;
    QString type;
    QString value;
    bool isRequired;
    bool canGroup;
    QString label;

    double min;
    double max;
    int decimals;
};

#endif // PPARAMETERDEFINITION_H
