//
//  parameters.h
//
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//

#ifndef parameters_h
#define parameters_h

#include "./config_parser.h"

struct Param {

  int maximum_time;

  float birth_normal;
  float birth_cancer;
  float birth_infected;
  float death_normal;
  float death_cancer;
  float death_infected;

  int infection_routine;
  float time_adding_virus;

  int num_cells;


  void readFromIni(const std::string& file_name) {
        ConfigFile from_config(file_name);

        maximum_time                = from_config.getValueOfKey<int>("maximum_time");

    }

    Param(const Param& other) {
        maximum_time = other.maximum_time;

    }

  Param() {
	  // default values, these don't necessarily make sense!
    maximum_time = 75; // default from Berg et al 2019
    time_adding_virus = 7; // default from Berg et al 2019
    num_cells = 100;

    float birth_normal = 0.5f;
    float birth_cancer = 1.0f;
    float birth_infected = 1.2f;
    float death_normal = 0.2f;
    float death_cancer = 0.1f;
    float death_infected = 0.1f;

    int infection_routine = 0;


  }


};

#endif /* parameters_h */
