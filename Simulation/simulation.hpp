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
#include "parameters.hpp"
#include <array>
#include <cmath>
#include "node.hpp"
#include "rndutils.hpp"
#include "random_thijs.hpp"


class simulation {
public:
  simulation();
  simulation(const Param& param);
  void run();
  void initialize_network();

  float t;
  void update_one_step();

  std::vector< node > world;
  std::vector<int> get_cell_numbers();

  rnd_t rndgen;

  const static size_t sq_size = 1000;

<<<<<<< Updated upstream
  const static size_t num_cells = sq_size * sq_size;

  std::array< binned_distribution<sq_size>, 4 > growth_prob_rnd;
  std::array< binned_distribution<sq_size>, 4 > death_prob_rnd;

  std::vector< std::vector< float >> growth_probs;
  std::vector< std::vector< float >> death_probs;
  void set_infection_type(const infection_routine& infect_type) ;
=======
  void add_infected(); // has to be public to allow for interaction by pressing
                       // a button
>>>>>>> Stashed changes

  void set_percent_infected(const float& percent);
    void add_infected();
private:

  int num_nodes;

  Param parameters;

  std::array<float, 4> num_cell_types;

  std::array< float, 8> rates;

  std::vector<double> long_distance_infection_probability;

  void count_cell_types(); // this is very heavy, should not be ran often
  void update_rates();
  void update_total_growth_prob();

  size_t pick_event(const std::array< float, 8>& rates, float sum);
  void do_event(size_t event);
  void update_growth_probabilities();
  void implement_death(const cell_type& parent);
  void implement_growth(const cell_type& parent);
  void update_growth_prob(size_t pos);
  void set_growth_prob(size_t pos);

  void update_death_prob(size_t pos);
  void update_death_prob(size_t pos, cell_type old_type);

  void add_cells(const cell_type& focal_cell_type);

  void print_to_file(float t);


  void infect_random();
  void infect_center();
  void infect_all_cancer();

  void initialize_full();

  void infect_long_distance(size_t pos);


//  void change_cell_type(const size_t& pos, const cell_type& focal_cell_type);
  void ask_infect_neighbours(int depth, float p, size_t pos);


};



#endif /* simulation_hpp */
