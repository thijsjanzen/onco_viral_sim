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

  std::cout << "Welcome to this Oncolytic Virus Simulator\n";
  std::cout << "You are using v0.2\n";

  try {
    Param all_parameters;
 //   all_parameters.readFromIni("config.ini");
    simulation Simulation(all_parameters);

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

