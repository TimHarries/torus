#ifndef PMONITOREVENT_H
#define PMONITOREVENT_H

#include <QString>

class PMonitorEvent
{
public:
    PMonitorEvent();

    enum TorusMonitorEvents {
        TorusStart,
        FitsReadComplete,
        GridCreationComplete,
        DustCalculationsComplete,
        GridSmoothingStarted,
        GridSmoothingComplete,
        LucyIterationComplete,
        LucyConverged,
        SEDOutput,
        FITSImageOutput,
        TorusFinish
    };

    enum TorusMonitorEventType {
        Torus,
        Lucy,
        SED,
        Image
    };

    TorusMonitorEvents type;
    long unsigned int time;
    QString description;
};

#endif // PMONITOREVENT_H
