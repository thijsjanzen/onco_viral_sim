CONFIG += console
CONFIG += c++17
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS_RELEASE += -O3

SOURCES += \
        ../Simulation/node.cpp \
        ../Simulation/setup.cpp \
        ../Simulation/simulation.cpp \
        main.cpp

HEADERS += \
    ../Simulation/node.hpp \
    ../Simulation/parameters.hpp \
    ../Simulation/random_thijs.hpp \
    ../Simulation/rndutils.hpp \
    ../Simulation/simulation.hpp \
    config_parser.h

DISTFILES += \
    config.ini
