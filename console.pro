CONFIG += console
CONFIG += c++17
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS_RELEASE += -O3

SOURCES += \
        ./Simulation/analysis.cpp \
        ./Simulation/node.cpp \
        ./Simulation/setup.cpp \
        ./Simulation/simulation.cpp \
        ./console/main.cpp

HEADERS += \
    ./Simulation/analysis.hpp \
    ./Simulation/node.hpp \
    ./Simulation/parameters.hpp \
    ./Simulation/random_thijs.hpp \
    ./Simulation/rndutils.hpp \
    ./Simulation/simulation.hpp \
    ./Simulation/voronoi.hpp \
    config_parser.h

DISTFILES += \
    config.ini

# In release mode, turn on profiling
CONFIG(release, debug|release) {

  DEFINES += NDEBUG

  # gprof
  QMAKE_CXXFLAGS += -pg
  QMAKE_LFLAGS += -pg
}

TARGET = CONSOLE.app

# Boost.Test
LIBS += -lboost_unit_test_framework
LIBS += -L"/usr/local/Cellar/boost/1.70.0/lib"
INCLUDEPATH += "/usr/local/Cellar/boost/1.70.0/include"
