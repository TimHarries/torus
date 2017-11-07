#include "pimagedata.h"

PImageData::PImageData()
{
    //Set the defaults for the values
    this->name = "";
    this->nphotons = 10000;
    this->npixels = 200;
    this->inclination = 0;
    this->positionangle = 0;
    this->lambdaimage = 6562.8;
    this->aspect = 1.0;
    this->polimage = false;
    this->imageaxisunits = "AU";
    this->imagecentrex = 0;
    this->imagecentrey = 0;
}
