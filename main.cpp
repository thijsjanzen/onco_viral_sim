//
//  main.cpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//

#include <iostream>
#include "simulation.hpp"

int main(int argc, const char* argv[]) {
  try {
    Param Params;
    Params.readFromIni("config.ini");
    model::Simulation Simulation(Params);

    Simulation.initialize_network();

    Simulation.run();
  } catch (const std::exception& err) {
    std::cerr << err.what() << '\n';
  } catch (const char* err) {
    std::cerr << err << '\n';
  } catch (...) {
    std::cerr << "unknown exception\n";
  }
  std::cerr << "bailing out.\n";
  return -1;
}

