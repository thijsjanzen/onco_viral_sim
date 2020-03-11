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
enum start_type {full, grow, converge};

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

  size_t initial_number_cancer_cells;
  size_t initial_number_normal_cells;

  float percent_infected;
  infection_routine infection_type;

  float prob_normal_infection;

  float distance_infection_upon_death;
  float prob_infection_upon_death;

  size_t sq_num_cells;
  size_t sq_num_pixels;
  bool use_voronoi_grid;

  start_type start_setup;

  Param(const Param& other) {
      maximum_time = other.maximum_time;
      time_adding_virus = other.time_adding_virus; // default from Berg et al 2019
      time_adding_cancer = other.time_adding_cancer; // normally the population is first allowed to stabilize

      initial_number_cancer_cells = other.initial_number_cancer_cells;
      initial_number_normal_cells = other.initial_number_normal_cells;

      birth_normal = other.birth_normal;
      birth_cancer = other.birth_cancer;
      birth_infected = other.birth_infected;
      death_normal = other.death_normal;
      death_cancer = other.death_cancer;
      death_infected = other.death_infected;

      death_cancer_resistant = other.death_cancer_resistant;
      birth_cancer_resistant = other.birth_cancer_resistant;

      infection_type = other.infection_type;
      percent_infected = other.percent_infected;

      prob_normal_infection = other.prob_normal_infection;
      freq_resistant = other.freq_resistant;
      distance_infection_upon_death = other.distance_infection_upon_death;
      prob_infection_upon_death = other.prob_infection_upon_death;

      start_setup = other.start_setup;

      use_voronoi_grid = other.use_voronoi_grid;

      sq_num_cells = other.sq_num_cells;
      sq_num_pixels = other.sq_num_pixels;
    }


  Param() {
    // default values, these don't necessarily make sense!
    maximum_time = 500; // default from Berg et al 2019
    time_adding_virus = 3000; // default from Berg et al 2019
    time_adding_cancer = 1000; // normally the population is first allowed to stabilize

    initial_number_cancer_cells = 500;
    initial_number_normal_cells = 1000;

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

    prob_normal_infection = 0.0f;
    freq_resistant = 0.0f;
    distance_infection_upon_death = 1;
    prob_infection_upon_death = 0.f;

    start_setup = converge;

    sq_num_cells = 100;
    sq_num_pixels = 100;

    use_voronoi_grid = false;
  }
};

#endif /* parameters_h */
