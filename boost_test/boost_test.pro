TEMPLATE = app
CONFIG -= qt
CONFIG -= app_bundle
CONFIG += console
CONFIG += c++17


BOOST_INCLUDE_DIR="/usr/local/Cellar/boost/1.70.0/include"
!isEmpty(BOOST_INCLUDE_DIR): INCLUDEPATH *= $${BOOST_INCLUDE_DIR}

isEmpty(BOOST_INCLUDE_DIR): {
    message("BOOST_INCLUDE_DIR is not set, assuming Boost can be found automatically in your system")
}

SOURCES += \
    ../Simulation/analysis.cpp \
    ../Simulation/node.cpp \
    ../Simulation/setup.cpp \
    ../Simulation/simulation.cpp \
    main.cpp

HEADERS += \
    ../Simulation/analysis.hpp \
    ../Simulation/node.hpp \
    ../Simulation/parameters.hpp \
    ../Simulation/random_thijs.hpp \
    ../Simulation/rndutils.hpp \
    ../Simulation/simulation.hpp \
    ../Simulation/voronoi.hpp

CONFIG += debug_and_release

CONFIG(debug, debug|release) {

  # gcov
  QMAKE_CXXFLAGS += -fprofile-arcs -ftest-coverage
  LIBS += -lgcov
}


#  LIBS += -lboost_unit_test_framework
