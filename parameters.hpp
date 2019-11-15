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
  float time_adding_cancer;

  int initial_number_cancer_cells;

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
    maximum_time = 500; // default from Berg et al 2019
    time_adding_virus = 7; // default from Berg et al 2019
    time_adding_cancer = 10000; // normally the population is first allowed to stabilize

    num_cells = 100*100;
    initial_number_cancer_cells = 10;

    birth_normal = 0.5f;
    birth_cancer = 1.0f;
    birth_infected = 1.2f;
    death_normal = 0.2f;
    death_cancer = 0.1f;
    death_infected = 0.1f;

    infection_routine = 0;




  }


};

#endif /* parameters_h */
