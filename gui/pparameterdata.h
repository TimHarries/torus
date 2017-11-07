#ifndef PPARAMETERDATA_H
#define PPARAMETERDATA_H

#include <psourcedata.h>
#include <pimagedata.h>
#include <pdustdata.h>
#include <QDebug>
#include <QList>
#include <QMap>
#include <QFile>
#include <QMessageBox>
#include <QMapIterator>
#include <QXmlStreamReader>
#include <QMessageBox>
#include <QResource>
#include <pparameterdefinition.h>
#include <QVector>
#include <QApplication>
#include <QHash>

#include "math.h"
#define PI 3.14159265

class PParameterData : QWidget
{
public:

    enum FileSource {
        InternalFile,
        ExternalFile
    };

    enum ParseState {
        Start,
        PropertyRead,
        ValueRead
    };

    //Maps for holding variables
    QMap<QString, bool> boolParams;
    QMap<QString, QString> stringParams;
    QMap<QString, double> numParams;
    QMap<QString, QString> rawParams;
    QMap<QString, QString> geomParams;

    //Arrays
    QList<PImageData> images;
    QList<PSourceData> sources;
    QList<PDustData> dusts;

    QString openFile;

    PParameterData();
    ~PParameterData();

    int loadFile(QString filename);
    int saveFile(QString filename);
    void loadDefinitions(PParameterData::FileSource source = InternalFile);
    void loadGeometries(PParameterData::FileSource source = InternalFile);
    void clear(void);
    bool hasKey(QString key);

    double getNum(QString key);
    bool getBool(QString key);
    QString getString(QString key);
    QVector<double> getVector(QString key);

    void setNum(QString key, double value);
    void setBool(QString key, bool value);
    void setString(QString key, QString value);
    void setVector(QString key, double x, double y, double z);

    void getDefault(QString key, QString &value);
    void getDefault(QString key, double &value );
    void getDefault(QString key, bool &value);

    QList<double> getAngles(double minAng, double maxAng, int incs, bool cosSpaced);
    void setAngles(QString key, QList<double> values);

    QString renderParameterString();
    void writeComment(QString text);
    void writeTuple(QString key, QString value);
    void writeLine();
    void writeGroup(QString groupName);
    void renderImageParams();
    void renderSourceParams();
    void renderDustParams();
    void renderGeometryParams();
    QString boolToStr(bool value);
    QString doubleToStr(double value);
    QString vecToStr(QVector<double> value);
    void saveToFile(QString fileName);

    QMap<QString, QMap<QString, PParameterDefinition> > paramDefinitions;
    QMap<QString, QList<PParameterDefinition> > geometries;

private:
    QString parameterString;
    QStringList definedParameters;

    void parseLine(QByteArray line);
    void refineData(void);
};

#endif // PPARAMETERDATA_H
