QT += core
QT -= gui

CONFIG += c++11

TARGET = TestModulesQuick
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    tests.cpp \
    mylib.cpp \
    myarraylib.cpp \
    ttable.cpp \
    kynemparam.cpp \
    mymathlib.cpp \
    array2dsizes.cpp

HEADERS += \
    tests.h \
    mylib.h \
    myarraylib.h \
    ttable.h \
    kynemparam.h \
    mymathlib.h \
    array2dsizes.h
