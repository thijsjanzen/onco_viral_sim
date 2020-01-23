//
//  simulation.hpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//

#ifndef simulation_hpp
#define simulation_hpp

#include <cstdio>
#include <vector>
#include <array>
#include <cmath>

#include "parameters.hpp"
#include "rndutils.hpp"
#include "random_thijs.hpp"
#include "node.hpp"

class simulation {
public:
  simulation(const Param& param);
  void run();
  void initialize_network();

  float t;
  void update_one_step();

  std::vector< node > world;
  rnd_t rndgen;

  std::array< binned_distribution, 4 > growth_prob;
  std::array< binned_distribution, 4 > death_prob;

  void add_infected(); // has to be public to allow for interaction by pressing
                       // a button

  void set_infection_type(const infection_routine& infect_type);
  void set_percent_infected(const float& percent);
  std::array<int, 5> num_cell_types;
  std::array<int, 5> count_cell_types();

private:
  int num_cells;
  int sq_size;

  Param parameters;

  std::array< float, 8> rates;

  std::vector<double> long_distance_infection_probability;



  void update_rates();
  void update_total_growth_prob();

  size_t pick_event(const std::array< float, 8>& rates, float sum);
  void do_event(size_t event);
  void update_growth_probabilities();
  void implement_death(const cell_type& parent);
  void implement_growth(const cell_type& parent);
  void update_cell(size_t pos);
  void update_growth_prob(size_t pos);
  void update_death_prob(size_t pos);
  void update_death_prob(size_t pos,
                         cell_type old_type,
                         cell_type new_type);

  void add_cells(const cell_type& focal_cell_type);

  void print_to_file(float t);

  void update_death_cdf(const cell_type& parent, size_t pos);
  void update_growth_cdf(size_t pos);
  void update_death_cdf_all();


  void infect_random();
  void infect_center();
  void infect_all_cancer();

  void initialize_full();

  void infect_long_distance(size_t pos);

  void setup_voronoi();


  void change_cell_type(const size_t& pos, const cell_type& focal_cell_type);
  void ask_infect_neighbours(int depth, float p, size_t pos);
};



#endif /* simulation_hpp */
