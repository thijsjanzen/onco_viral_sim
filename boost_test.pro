TEMPLATE = app
CONFIG -= qt
CONFIG -= app_bundle
CONFIG += console

CONFIG += c++17
QMAKE_CXXFLAGS += -std=c++17
# QMAKE_CXXFLAGS += "-Wno-old-style-cast"

BOOST_INCLUDE_DIR="/usr/local/Cellar/boost/1.75.0_2/include"
!isEmpty(BOOST_INCLUDE_DIR): INCLUDEPATH *= $${BOOST_INCLUDE_DIR}

isEmpty(BOOST_INCLUDE_DIR): {
    message("BOOST_INCLUDE_DIR is not set, assuming Boost can be found automatically in your system")
}

SOURCES += \
    ./Simulation/analysis.cpp \
    ./Simulation/node.cpp \
    ./Simulation/setup.cpp \
    ./Simulation/simulation.cpp \
    ./boost_test_files/main.cpp

HEADERS += \
    ./Simulation/analysis.hpp \
    ./Simulation/node.hpp \
    ./Simulation/parameters.hpp \
    ./Simulation/random_thijs.hpp \
    ./Simulation/rndutils.hpp \
    ./Simulation/simulation.hpp \
    ./Simulation/voronoi.hpp

CONFIG += debug_and_release

TARGET=BOOST_TEST.app

QMAKE_CXXFLAGS += --coverage
QMAKE_LFLAGS += --coverage

QMAKE_POST_LINK = rm -f "*.gcda"

CONFIG(debug, debug|release) {
  # gcov
  QMAKE_CXXFLAGS += -fprofile-arcs -ftest-coverage
}

# LIBS += -lboost_unit_test_framework
