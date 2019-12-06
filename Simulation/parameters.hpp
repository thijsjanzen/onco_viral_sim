//
//  parameters.h
//
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//

#ifndef parameters_h
#define parameters_h

#include <cstdio>
#include <cmath>

enum infection_routine {random_infection, center_infection};
enum start_type {full, grow};

struct Param {

  int maximum_time;

  float birth_normal;
  float birth_cancer;
  float birth_infected;
  float death_normal;
  float death_cancer;
  float death_infected;

  float birth_cancer_resistant;
  float death_cancer_resistant;
  float freq_resistant;


  int time_adding_virus;
  int time_adding_cancer;

  int initial_number_cancer_cells;
  int initial_number_normal_cells;

  float percent_infected;
  infection_routine infection_type;

  float prob_normal_infection;

  float distance_infection_upon_death;
  float prob_infection_upon_death;



  start_type start_setup;


    Param(const Param& other) {
        maximum_time = other.maximum_time;

    }

  Param() {
	  // default values, these don't necessarily make sense!
    maximum_time = 500; // default from Berg et al 2019
    time_adding_virus = 3000; // default from Berg et al 2019
    time_adding_cancer = 1000; // normally the population is first allowed to stabilize

    initial_number_cancer_cells = 5;

    birth_normal = 0.5f;
    birth_cancer = 1.0f;
    birth_infected = 1.2f;
    death_normal = 0.2f;
    death_cancer = 0.1f;
    death_infected = 0.1f;

    death_cancer_resistant = 0.1f;
    birth_cancer_resistant = 0.9f;

    infection_type = random_infection;
    percent_infected = 0.1f;

    prob_normal_infection = 0.01f;
    freq_resistant = 0.01f;
    distance_infection_upon_death = 1;
    prob_infection_upon_death = 0.f;
  }


};

#endif /* parameters_h */
