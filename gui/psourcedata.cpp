#include "psourcedata.h"

PSourceData::PSourceData()
{
    //Properties
    this->name = "";
    this->pointsource = false;
    this->radius = 1;
    this->teff = 1;
    this->mass = 1;
    this->contflux = "blackbody";
    this->probsource = 1;

    this->position.resize(3);
    this->velocity.resize(3);
}
