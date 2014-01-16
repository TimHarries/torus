#-------------------------------------------------
#
# Project created by QtCreator 2013-06-11T20:48:11
#
#-------------------------------------------------
QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Parametor
TEMPLATE = app
RC_FILE = parametor.rc
QT += svg

unix {
    LIBS += -L/usr/local/qwt-6.1.1-svn/lib -lqwt
    INCLUDEPATH += /usr/local/qwt-6.1.1-svn/include
}

macx {
    #APP_BUNDLE_FILES.files = torus_input_definitions.xml
    #APP_BUNDLE_FILES.path = Contents/Resources
    #QMAKE_BUNDLE_DATA += APP_BUNDLE_FILES
    QMAKE_INFO_PLIST = Info.plist
    ICON = ParametorIcon.icns

    #Qwt libraries for Mac OS X
    LIBS += -L/usr/local/qwt-6.1.0/lib -lqwt
    INCLUDEPATH += /Users/smorrell/build/qwt/qwt-6.1.0/src
    CONFIG += qwt
}

win32 {
    CONFIG(debug, debug|release) {
        LIBS += -LC:\Qwt-6.1.0\lib -lqwtd
    } else {
        LIBS += -LC:\Qwt-6.1.0\lib -lqwt
    }
    INCLUDEPATH += C:\Qwt-6.1.0\include
    CONFIG += qwt
}

SOURCES += main.cpp\
        mainwindow.cpp \
        pparameterdata.cpp \
    psourcedata.cpp \
    pimagedata.cpp \
    panglegen.cpp \
    torusconsole.cpp \
    unitconverter.cpp \
    pparameterdefinition.cpp \
    pdustdata.cpp \
    ptorusupdatewidget.cpp \
    paboutdialog.cpp \
    pcompileconfiguredialog.cpp \
    ptorusmonitor.cpp \
    psedplotter.cpp \
    pmonitorevent.cpp

HEADERS  += mainwindow.h \
    	pparameterdata.h \
    psourcedata.h \
    pimagedata.h \
    panglegen.h \
    torusconsole.h \
    unitconverter.h \
    pparameterdefinition.h \
    pdustdata.h \
    ptorusupdatewidget.h \
    paboutdialog.h \
    pcompileconfiguredialog.h \
    ptorusmonitor.h \
    psedplotter.h \
    pmonitorevent.h

FORMS    += mainwindow.ui \
    panglegen.ui \
    torusconsole.ui \
    unitconverter.ui \
    ptorusupdatewidget.ui \
    paboutdialog.ui \
    ptorusmonitor.ui \
    psedplotter.ui

RESOURCES += \
    Resources.qrc

