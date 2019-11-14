//
//  simulation.hpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//

#ifndef simulation_hpp
#define simulation_hpp

#include <stdio.h>
#include <vector>
#include "parameters.hpp"
#include "random_thijs.hpp"
#include <array>
#include <cmath>
#include "node.hpp"

class simulation {
public:
  simulation(const Param& param);
  void run();
  void initialize_network();

private:

  int num_nodes;
  int sq_size;
  Param parameters;
  rnd_t rndgen;

  std::vector< node > world;
  std::vector< float > growth_prob_normal;
  std::vector< float > growth_prob_cancer;
  std::vector< float > growth_prob_infected;

  std::vector< float > death_prob_normal;
  std::vector< float > death_prob_cancer;
  std::vector< float > death_prob_infected;

  float max_growth_prob_normal;
  float max_growth_prob_cancer;
  float max_growth_prob_infected;

  std::array<float, 4> num_cell_types;

  void update_cell_types(); // this is very heavy, should not be ran often
  void update_rates(std::array< float, 6>& rates);
  void update_total_growth_prob();

  size_t pick_event(const std::array< float, 6>& rates, float sum);
  void update_cell_counts(int event);
  void do_event(int event);
  void update_growth_probabilities();
  void implement_death(const cell_type& parent);
  void implement_growth(const cell_type& parent);
};



#endif /* simulation_hpp */
