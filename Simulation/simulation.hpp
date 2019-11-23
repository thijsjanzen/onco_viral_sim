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
#include <array>
#include <cmath>
#include "node.hpp"
#include "rndutils.hpp"
#include "random_thijs.hpp"

class simulation {
public:
  simulation(const Param& param);
  void run();
  void initialize_network();

  float t;
  size_t sq_size; // size of one side of the rectangle of the world
  void update_one_step();

  std::vector< node > world;
  std::vector<int> get_cell_numbers();

  rnd_t rndgen;


private:

  int num_nodes;

  Param parameters;

  std::vector< rndutils::mutable_discrete_distribution<int, rndutils::all_zero_policy_throw > > growth_prob_rnd;
  std::vector< rndutils::mutable_discrete_distribution<int, rndutils::all_zero_policy_throw > > death_prob_rnd;

  std::vector< std::vector< float >> growth_probs;
  std::vector< std::vector< float >> death_probs;


  std::array<float, 3> num_cell_types;

  std::array< float, 6> rates;

  void count_cell_types(); // this is very heavy, should not be ran often
  void update_rates(std::array< float, 6>& rates);
  void update_total_growth_prob();

  size_t pick_event(const std::array< float, 6>& rates, float sum);
  void do_event(size_t event);
  void update_growth_probabilities();
  void implement_death(const cell_type& parent);
  void implement_growth(const cell_type& parent);
  void update_growth_prob(size_t pos);

  void add_cells(cell_type focal_cell_type);

  void print_to_file(float t);

  void update_death_prob_vectors(const cell_type& parent, int pos);
  void update_growth_prob_vectors(int pos);
};



#endif /* simulation_hpp */
