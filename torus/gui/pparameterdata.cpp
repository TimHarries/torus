#include "pparameterdata.h"

PParameterData::PParameterData()
{
    this->openFile = "";
    //First we load the internal files compiled into the executable then we get the external ones
    loadDefinitions();
    loadDefinitions(ExternalFile);
    loadGeometries();
    loadGeometries(ExternalFile);
}

PParameterData::~PParameterData() {
    this->openFile = "";
    this->definedParameters.clear();
    this->geometries.clear();
}

void PParameterData::loadDefinitions(PParameterData::FileSource source) {
    //Load the definitions files
    QXmlStreamReader *xmlReader;
    QString externalFile = QApplication::applicationDirPath() + "/torus-input-definitions.xml";
    QFile xmlFile;
    if(source == PParameterData::ExternalFile && QFile(externalFile).exists()) {
        xmlFile.setFileName(externalFile);
    } else {
        xmlFile.setFileName(":/torus-input-definitions.xml");
    }

    if(!xmlFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QMessageBox::critical(this, "Torus Definition File Error", "Couldn't load the Torus Inputs Definitions File.", QMessageBox::Ok);
        return;
    }
    xmlReader = new QXmlStreamReader(&xmlFile);

    //Parse the XML into definitions
    bool inGroup = false;
    QString groupName = "";
    QMap<QString, PParameterDefinition> group;

    while(!xmlReader->atEnd()) {
        //Read the Next Element
        QXmlStreamReader::TokenType token = xmlReader->readNext();
        if(token == QXmlStreamReader::StartDocument) {
            continue;
        }

        if(xmlReader->hasError()) {
            qDebug() << xmlReader->errorString();
        }

        //Read the Element
        if(token == QXmlStreamReader::StartElement) {
            if(xmlReader->name() == "group") {
                if(!inGroup) {
                    inGroup = true;
                    groupName = xmlReader->attributes().value("name").toString();
                    if(!this->paramDefinitions.contains(groupName)) {
                        this->paramDefinitions.insert(groupName, QMap<QString, PParameterDefinition>());
                    }
                }
            } else if(xmlReader->name() == "input") {
                QXmlStreamAttributes att = xmlReader->attributes();
                PParameterDefinition def;
                def.name = att.value("name").toString();
                def.type = att.value("type").toString();
                def.defaultValue = att.value("default").toString();
                def.canGroup = att.value("canGroup").toString().toInt();
                def.label = att.value("label").toString();
                if(att.value("isRequired").toString() == "true") {
                    def.isRequired = true;
                } else {
                    def.isRequired = false;
                }

                //Optional Parameters
                if(att.value("maximum") != "") {

                }

                if(att.value("minimum") != "") {

                }

                if(att.value("decimals") != "") {

                }

                this->paramDefinitions[groupName].insert(def.name, def);
            }
        }

        if(token == QXmlStreamReader::EndElement) {
            if(xmlReader->name() == "group" && inGroup == true) {
                inGroup = false;
            }
        }

        if(token == QXmlStreamReader::EndDocument) {
            continue;
        }
    }
}

void PParameterData::loadGeometries(PParameterData::FileSource source) {
    //Load the definitions files
    QXmlStreamReader *xmlReader;
    QFile xmlFile;
    QString externalFile = QApplication::applicationDirPath() + "/geometries.xml";

    if(source == PParameterData::ExternalFile && QFile(externalFile).exists()) {
        xmlFile.setFileName(externalFile);
    } else {
        xmlFile.setFileName(":/geometries.xml");
    }

    if(!xmlFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QMessageBox::critical(this, "Torus Geometry File Error", "Couldn't load the Torus Geometry Definitions File.", QMessageBox::Ok);
        return;
    }
    xmlReader = new QXmlStreamReader(&xmlFile);

    //Parse the XML into definitions
    bool inGroup = false;
    QString groupName = "";
    QMap<QString, PParameterDefinition> group;

    while(!xmlReader->atEnd()) {
        //Read the Next Element
        QXmlStreamReader::TokenType token = xmlReader->readNext();
        if(token == QXmlStreamReader::StartDocument) {
            continue;
        }

        if(xmlReader->hasError()) {
            qDebug() << xmlReader->errorString();
        }

        //Read the Element
        if(token == QXmlStreamReader::StartElement) {
            if(xmlReader->name() == "geometry") {
                if(!inGroup) {
                    inGroup = true;
                    groupName = xmlReader->attributes().value("name").toString();
                    if(!this->geometries.contains(groupName)) {
                        this->geometries.insert(groupName, QList<PParameterDefinition>());
                    }
                }
            } else if(xmlReader->name() == "input") {
                QXmlStreamAttributes att = xmlReader->attributes();
                PParameterDefinition def;
                def.name = att.value("name").toString();
                def.type = att.value("type").toString();
                def.defaultValue = att.value("default").toString();
                def.canGroup = att.value("canGroup").toString().toInt();
                def.label = att.value("label").toString();
                if(att.value("isRequired").toString() == "true") {
                    def.isRequired = true;
                } else {
                    def.isRequired = false;
                }

                //Check to see if Geometries already exist and replace if needed
                bool exists = false;
                int replIndex;
                if(!geometries.value(groupName).isEmpty() && def.name != "") {
                    for(int i = 0; i < geometries.value(groupName).count(); i++) {
                        if(geometries.value(groupName).at(i).name == def.name) {
                            exists = true;
                            replIndex = i;
                            break;
                        }
                    }
                } else {
                    exists = false;
                }

                if(exists) {
                    geometries[groupName].replace(replIndex, def);
                } else {
                    geometries[groupName].append(def);
                }

            }
        }

        if(token == QXmlStreamReader::EndElement) {
            if(xmlReader->name() == "geometry" && inGroup == true) {
                inGroup = false;
                group.clear();
            }
        }

        if(token == QXmlStreamReader::EndDocument) {
            continue;
        }
    }
}


/*
 *Responsible for loading a parameters.dat file into memory
 */
int PParameterData::loadFile(QString filename) {
    QFile file(filename);
    QMessageBox msgBox;

    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
          return 0;

    this->openFile = filename;
    while(file.atEnd() == false) {
        QByteArray line = file.readLine();
        this->parseLine(line);
    }

    this->refineData();

    return 0;
}

/*
 *Responsible for saving a parameters.dat file
 */
int PParameterData::saveFile(QString filename) {
    QFile file(filename);
    return 0;
}

/*
 *Basic finite state machine for parsing in the parameters file
 *The States are 1: Property Name Read, 2: Value Read
 *It has the ability to ignore comments
 */
/*
void PParameterData::parseLine(QByteArray line) {
    char chr;
    bool nextLine = false;
    int state = 1;
    bool isVector = false;
    QString property("");
    QString value("");

    for (int i = 0; i < line.count(); i++) {
        chr = line.at(i);

        //Check to see if Character a comment
        if(chr == '!') {
            nextLine = true;
        }

        //Checks to see if not commented
        if(!nextLine) {
            switch(state) {
                //Property Read State
                case 1:
                    switch (chr) {

                        case ' ' :
                            state = 2;
                            break;

                        default:
                            property.append(QString(chr));
                            break;
                    }
                break;

                //Value read State
                case 2 :
                    switch(chr) {

                        case 10 :
                            isVector = false;
                            break;

                        case '.' :
                            value.append(QString(chr));
                            isVector = true;
                            break;

                        case ' ' :
                            if(isVector) {
                                value.append(QString(chr));
                            }
                            isVector = false;
                            break;

                        default :
                            value.append(QString(chr));
                            isVector = false;
                            break;
                    }
                 break;
             }
        }
    }

    //Write the value to params var for parsing later
    if(property.count() != 0 && value.count() != 0) {
        this->rawParams.insert(property, value);
    }
} */

void PParameterData::parseLine(QByteArray line) {
    char chr;
    bool nextLine = false;
    ParseState state = Start;
    bool isVector = false;
    QString property("");
    QString value("");
    QString tempBuffer;

    for (int i = 0; i < line.count(); i++) {
        chr = line.at(i);

        //Check to see if Character a comment
        if(chr == '!') {
            tempBuffer.clear();
            nextLine = true;
        }

        //Checks to see if not commented
        if(!nextLine) {
            switch(state) {
                //Property Read State
            case Start:
                switch (chr) {

                    case ' ' :
                        state = ValueRead;
                        break;

                    default:
                        property.append(QString(chr));
                        state = PropertyRead;
                        continue;
                        break;
                }

                case PropertyRead:
                    switch (chr) {

                        case ' ' :
                            state = ValueRead;
                            break;

                        default:
                            property.append(QString(chr));
                            break;
                    }
                break;

                //Value read State
                case ValueRead :
                    switch(chr) {

                        case 10 :
                            isVector = false;
                            break;

                        case ' ' :
                            tempBuffer.append(QString(chr));
                            break;

                        default :
                            if(!tempBuffer.isEmpty()) {
                                value.append(tempBuffer);
                                tempBuffer.clear();
                            }
                            value.append(QString(chr));
                            isVector = false;
                            break;
                    }
                 break;
             }
        }
    }

    //Write the value to params var for parsing later
    if(property.count() != 0 && value.count() != 0) {
        this->rawParams.insert(property, value);
    }
}

/*
 *Refines the raw data read in to the independent arrays
 */
void PParameterData::refineData(void) {
    QMapIterator<QString, QString> params(this->rawParams);

    while(params.hasNext()) {
        params.next();

        qDebug() << params.key() << ": " << params.value();

        if(params.value() == "T" || params.value() == "F") {
            //This is a bool
            bool boolVal;

            if(params.value() == "T") {
                boolVal = true;
            } else {
                boolVal = false;
            }

            this->boolParams.insert(params.key(), boolVal);

        } else if(params.value().toDouble() || params.value() == "0.0" || params.value() == "0") {
            //This is a double
            this->numParams.insert(params.key(), params.value().toDouble());
        } else {
            //This is a string
            this->stringParams.insert(params.key(), params.value());
        }
    }

    //Getting Image Parameters into a QList
    int nimage = 0;
    if(hasKey("nimage")) {
        nimage = this->getNum("nimage");
    } else {
        //Count from variables names
        while(this->hasKey("filename" + QString::number(nimage + 1))) {
            nimage++;
        }
    }

    for(int i = 1; i <= nimage; i++) {
        PImageData image;
        image.name = this->getString("imagefile" + QString::number(i));
        image.npixels = this->getNum("npixels" + QString::number(i));
        image.inclination = this->getNum("inclination" + QString::number(i));
        image.positionangle = this->getNum("positionangle" + QString::number(i));
        image.lambdaimage = this->getNum("lambdaimage" + QString::number(i));
        image.imagecentrex = this->getNum("imagecentrex" + QString::number(i));
        image.imagecentrey = this->getNum("imagecentrey" + QString::number(i));

        this->images.append(image);
    }

    //Getting Source Parameters into a QList
    int nsource = 0;
    if(hasKey("nsource")) {
        nsource = this->getNum("nsource");
    } else {
        //Count from variables names
        while(this->hasKey("radius" + QString::number(nsource + 1)) || this->hasKey("mass" + QString::number(nsource + 1)) || this->hasKey("teff" + QString::number(nsource + 1)) || this->hasKey("contflux" + QString::number(nsource + 1))) {
            nsource++;
        }
    }

     for(int i = 1; i <= nsource; i++) {
        PSourceData source;
        source.name = QString("Source " + QString::number(i));
        source.pointsource = this->getBool("pointsource" + QString::number(i));
        source.pointsource = this->getNum("radius" + QString::number(i));
        source.teff = this->getNum("teff" + QString::number(i));
        source.mass = this->getNum("mass" + QString::number(i));
        source.radius = this->getNum("radius" + QString::number(i));
        source.contflux = this->getString("contflux" + QString::number(i));
        source.probsource = this->getNum("prodsource" + QString::number(i));
        source.position = this->getVector("sourcepos" + QString::number(i));
        source.velocity = this->getVector("velocity" + QString::number(i));

        //Append to the Sources List
        this->sources.append(source);
     }

     //Load Load in Dust Params
     int ndusttype = 0;
     if(hasKey("ndusttype")) {
         ndusttype = this->getNum("ndusttype");
     } else {
         //Count from variables names
         while(this->hasKey("amin" + QString::number(ndusttype + 1)) || this->hasKey("amax" + QString::number(ndusttype + 1)) || this->hasKey("grainfrac" + QString::number(ndusttype + 1)) || this->hasKey("graintype" + QString::number(ndusttype + 1))) {
             ndusttype++;
         }
     }

     for(int i = 1; i <= ndusttype; i++) {
         PDustData dust;
         dust.label = "Dust Population " + QString::number(i);
         dust.graintype = this->getString("graintype" + QString::number(i));
         dust.amin = this->getNum("amin" + QString::number(i));
         dust.amax = this->getNum("amax" + QString::number(i));
         dust.qdist = this->getNum("qdist" + QString::number(i));
         this->dusts.append(dust);
         i++;
     }

     //Clear the raw list to save memory
     this->rawParams.clear();
}


/*
 *---------------Safe Accessors ----------------------
 */
double PParameterData::getNum(QString key) {
    if(this->numParams.contains(key)) {
        return this->numParams.value(key);
    } else if(key.contains(QRegExp("[0-9]+$"))) {
        key.remove(QRegExp("[0-9]+$"));
        if(this->numParams.contains(key)) {
            return this->numParams.value(key);
        }
    }

    double retval;
    key.remove(QRegExp("[0-9]+$"));
    this->getDefault(key, retval);
    return retval;
}

bool PParameterData::getBool(QString key) {
    if(this->boolParams.contains(key)) {
        return this->boolParams.value(key);
    } else if(key.contains(QRegExp("[0-9]+$"))) {
        key.remove(QRegExp("[0-9]+$"));
        if(this->boolParams.contains(key)) {
            return this->boolParams.value(key);
        }
    }

    bool retval;
    key.remove(QRegExp("[0-9]+$"));
    this->getDefault(key, retval);
    return retval;
}

QString PParameterData::getString(QString key) {
    if(this->stringParams.contains(key)) {
        return this->stringParams.value(key);
    }else if(key.contains(QRegExp("[0-9]+$"))) {
        key.remove(QRegExp("[0-9]+$"));
        if(this->stringParams.contains(key)) {
            return this->stringParams.value(key);
        }
    }

    QString retval;
    this->getDefault(key, retval);
    return retval;
}

QVector<double> PParameterData::getVector(QString key) {
    QVector<double> returnvec(3, 0x0);
    QString valString;
    QStringList split;

    if(this->stringParams.contains(key)){
        valString = this->stringParams.value(key);
        split = valString.split(". ");

        if(split.count() == 3) {
            returnvec[0] = split.at(0).toDouble();
            returnvec[1] = split.at(1).toDouble();
            returnvec[2] = split.at(2).toDouble();
        }

    } else if(key.contains(QRegExp("[0-9]+$"))) {
        key.remove(QRegExp("//[0-9]+$//"));

        if(this->stringParams.contains(key)) {
            valString = this->stringParams.value(key);
            split = valString.split(".");

            if(split.count() == 3) {
                returnvec[0] = split.at(0).toDouble();
                returnvec[1] = split.at(1).toDouble();
                returnvec[2] = split.at(2).toDouble();
            }
        }
    }

    return returnvec;
}


/*
 *---------------Safe Mutators ----------------------
 */
void PParameterData::setNum(QString key, double value) {
    this->numParams.insert(key, value);
}

void PParameterData::setBool(QString key, bool value) {
    if(value == true) {
        this->boolParams.insert(key, true);
    } else {
        this->boolParams.insert(key, false);
    }
}

void PParameterData::setString(QString key, QString value) {
    this->stringParams.insert(key, value);
}

void PParameterData::setVector(QString key, double x, double y, double z) {
    QString valueString = QString(QString::number(x) + " " + QString::number(y) + " " + QString::number(z));
    this->stringParams.insert(key, valueString);
}


void PParameterData::clear(void) {
    this->sources.clear();
    this->images.clear();
    this->dusts.clear();

    this->boolParams.clear();
    this->numParams.clear();
    this->stringParams.clear();
    this->rawParams.clear();
}

QList<double> PParameterData::getAngles(double minAngle, double maxAngle, int inc, bool cosSpaced = false) {
    QList<double> retList;
    if(cosSpaced) {
        //Don't forget that Math in C uses rads
        double minAngleRad = minAngle / (180 / PI);
        double maxAngleRad = maxAngle / (180 / PI);
        double dcosinc = (cos(maxAngleRad) - cos(minAngleRad)) / (inc - 1);
        double value = minAngle;

        for(int i = 0; i < inc; i++) {
            value = acos(cos(minAngleRad) + dcosinc * i) * 180 / PI;
            retList.append(value);
        }

    } else {
        double increment = (maxAngle - minAngle) / (inc - 1);
        double value = minAngle;

        for(int i = 0; i < inc; i++) {
            value = minAngle + (increment * i);
            retList.append(value);
        }
    }

    return retList;
}

void PParameterData::setAngles(QString key, QList<double> values) {
    if(!values.isEmpty()) {
        QString angles = QString::number(values.takeFirst());
        foreach(double d, values) {
            angles += " " + QString::number(d);
        }
        this->stringParams.insert(key, angles);
    }
}

/*-------------------------------------------------------------------------
 *Everything beyond this point is to deal with rendering the paramters file
 *-------------------------------------------------------------------------
 */

void PParameterData::saveToFile(QString fileName) {
    QFile file(fileName);
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    out << this->renderParameterString();
    file.close();

    //Update the open file string with the new file
    this->openFile = fileName;
}

QString PParameterData::renderParameterString() {
    parameterString.clear();
    definedParameters.clear();
   //Assuming at this point that the model has been updated from the UI

    //General Settings
    writeComment("Torus v2 Parameters File");

    definedParameters << "amr1d" << "amr2d" << "amr3d";
    if(this->getBool("amr1d") == true) {
        writeTuple("amr1d", boolToStr(getBool("amr1d")));
    } else if(this->getBool("amr2d") == true) {
        writeTuple("amr2d", boolToStr(getBool("amr2d")));
    } else if(this->getBool("amr3d") == true) {
        writeTuple("amr3d", boolToStr(getBool("amr3d")));
    }

    //Set Geometries (diabled)
    /*
    if(this->hasKey("inputfile")) {
        writeTuple("geometry", "fitsfile");
    } */

    writeGroup("misc");

    //Geometry
    writeLine();
    writeComment("----------Geometry----------");
    renderGeometryParams();

    //Sources
    writeLine();
    writeComment("----------Sources----------");
    renderSourceParams();
    writeGroup("source");

    writeComment("----------Hydro----------");
    writeGroup("hydro");

    writeComment("----------Dust Properties----------");
    renderDustParams();
    writeGroup("dust");

    //SEDs
    writeComment("----------SED Parameters----------");
    writeTuple("spectrum", boolToStr(getBool("spectrum")));

    writeGroup("sed");
    writeLine();

    //Images
    writeComment("----------Images----------");
    writeTuple("image", boolToStr(getBool("image")));

    renderImageParams();
    writeGroup("image");

    return parameterString;
}

void PParameterData::writeComment(QString text) {
    parameterString.append("!" + text + "\n");
}

void PParameterData::writeTuple(QString key, QString value) {
    parameterString.append(key + " " + value + "\n");
    definedParameters.append(key.remove(QRegExp("[0-9]*$")));
}

void PParameterData::writeLine() {
    parameterString.append("\n");
}

void PParameterData::writeGroup(QString groupName) {
    QMap<QString, PParameterDefinition> params;
    if(paramDefinitions.contains(groupName)) {
        params = paramDefinitions.value(groupName);

        QMapIterator<QString, PParameterDefinition> parit(params);
        while(parit.hasNext()) {
            parit.next();
            PParameterDefinition def = parit.value();

            if((def.isRequired == true || this->hasKey(def.name)) && !definedParameters.contains(def.name)) {
                if(this->hasKey(def.name)) {

                    if(def.type == "string") {
                        writeTuple(def.name, this->getString(def.name));
                    } else if(def.type == "double") {
                        writeTuple(def.name, QString::number(this->getNum(def.name)));
                    } else if(def.type == "int") {
                        writeTuple(def.name, QString::number((int)this->getNum(def.name)));
                    } else if(def.type == "bool") {
                        writeTuple(def.name, boolToStr(this->getBool(def.name)));
                    }

                } else {
                    if(def.type == "string") {
                        writeTuple(def.name, def.defaultValue);
                    } else if(def.type == "double") {
                        writeTuple(def.name, def.defaultValue);
                    } else if(def.type == "int") {
                        writeTuple(def.name, def.defaultValue);
                    } else if(def.type == "bool") {
                        writeTuple(def.name, boolToStr(def.defaultValue.toInt()));
                    }
                }
            }
        }
        writeLine();
    }
}

QString PParameterData::boolToStr(bool value) {
    QString retval;
    if(value == true) {
        retval = "T";
    } else {
        retval = "F";
    }
    return retval;
}

QString PParameterData::doubleToStr(double value) {
    int exp = 0;
    while(value > 10) {
        value = value / pow(10, exp);
        exp++;
    }
    return QString::number(value) + "d" + QString::number(exp, 10, 0);
}

QString PParameterData::vecToStr(QVector<double> value) {
    QString retval = "";
    if(value.count() > 0) {
        retval = QString::number(value.at(0));
        value.remove(0);
        for(int i = 0; i < value.count(); i++) {
            retval += " " + QString::number(value.at(i));
        }
    }

    return retval;
}

void PParameterData::renderImageParams() {
    QListIterator<PImageData> *imgs = new QListIterator<PImageData>(this->images);
    PImageData img;
    int i = 0;

    writeTuple("nimage", QString::number(this->images.count()));

    if(this->images.count() > 1) {
        while(imgs->hasNext()) {
            img = imgs->next();
            i++;

            writeTuple("imagefile" + QString::number(i), img.name);
            writeTuple("positionangle" + QString::number(i), QString::number(img.positionangle, 'f', 3));
            writeTuple("inclination" + QString::number(i), QString::number(img.inclination, 'f', 3));
            writeTuple("lambdaimage" + QString::number(i), QString::number(img.lambdaimage, 'f', 3));
            writeTuple("npixels" + QString::number(i), QString::number(img.npixels));

            writeLine();
        }
    } else {
        while(imgs->hasNext()) {
            img = imgs->next();
            i++;

            writeTuple("imagefile", img.name);
            writeTuple("positionangle", QString::number(img.positionangle, 'f', 3));
            writeTuple("inclination", QString::number(img.inclination, 'f', 3));
            writeTuple("lambdaimage", QString::number(img.lambdaimage, 'f', 3));
            writeTuple("npixels", QString::number(img.npixels));

            writeLine();
        }
    }
}

void PParameterData::renderSourceParams() {
    QListIterator<PSourceData> *sources = new QListIterator<PSourceData>(this->sources);
    PSourceData source;
    int i = 0;

    writeTuple("nsource", QString::number(this->sources.count()));

    while(sources->hasNext()) {
        source = sources->next();
        i++;

        writeTuple("radius" + QString::number(i), QString::number(source.radius));
        writeTuple("teff" + QString::number(i), QString::number(source.teff));
        writeTuple("contflux" + QString::number(i), source.contflux);
        writeTuple("mass" + QString::number(i), QString::number(source.mass));
        writeTuple("sourcepos" + QString::number(i), vecToStr(source.position));

        writeLine();
    }
}

void PParameterData::renderDustParams() {
    QListIterator<PDustData> *dustList = new QListIterator<PDustData>(this->dusts);
    PDustData dust;
    int i = 0;

    writeTuple("ndusttype", QString::number(this->dusts.count()));

    while(dustList->hasNext()) {
        dust = dustList->next();
        i++;

        writeTuple("amin" + QString::number(i), QString::number(dust.amin));
        writeTuple("amax" + QString::number(i), QString::number(dust.amax));
        writeTuple("qdist" + QString::number(i), QString::number(dust.qdist));
        writeTuple("grainFrac" + QString::number(i), QString::number(dust.grainfrac));
        writeTuple("graintype" + QString::number(i), dust.graintype);
        writeTuple("a0" + QString::number(i), QString::number(dust.a0));
        writeTuple("pdist" + QString::number(i), QString::number(dust.pdist));

        writeLine();
    }
}

void PParameterData::renderGeometryParams() {
    foreach(QString key, this->geomParams.keys()) {
        QString value = this->geomParams.value(key);
        writeTuple(key, value);
    }
}

bool PParameterData::hasKey(QString key) {
    if(stringParams.contains(key)) {
        return true;
    } else if(numParams.contains(key)) {
        return true;
    } else if(boolParams.contains(key)) {
        return true;
    }
    return false;
}

/* Default Accessors */
void PParameterData::getDefault(QString key, QString &value) {

    foreach(QString groupKey, this->paramDefinitions.keys()) {
        foreach(QString defKey, this->paramDefinitions.value(groupKey).keys()) {
            if(defKey == key) {
                PParameterDefinition def = this->paramDefinitions.value(groupKey).value(defKey);
                value = def.defaultValue;
            }
        }
    }
}

void PParameterData::getDefault(QString key, double &value) {
    value = 0;
    foreach(QString groupKey, this->paramDefinitions.keys()) {
        foreach(QString defKey, this->paramDefinitions.value(groupKey).keys()) {
            if(defKey == key) {
                PParameterDefinition def = this->paramDefinitions.value(groupKey).value(defKey);
                value = def.defaultValue.toDouble();
            }
        }
    }

}

void PParameterData::getDefault(QString key, bool &value) {
    value = false;

    foreach(QString groupKey, this->paramDefinitions.keys()) {
        foreach(QString defKey, this->paramDefinitions.value(groupKey).keys()) {
            if(defKey == key) {
                PParameterDefinition def = this->paramDefinitions.value(groupKey).value(defKey);
                value = def.defaultValue.toInt();
            }
        }
    }
}
