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
  void initialize_network(std::vector< std::vector< voronoi_point > >& all_polys,
                          grid_type used_grid_type);

  float t;
  size_t num_cells;
  size_t sq_size;

  void update_one_step();

  std::vector< node > world;
  rnd_t rndgen;

  std::array< binned_distribution, 4 > growth_prob;
  std::array< binned_distribution, 4 > death_prob;

  void set_percent_infected(float);
  void set_infection_type(infection_routine infect_routine);
  void add_infected(infection_routine infect_type,
                    float fraction); // has to be public to allow for interaction by pressing
                       // a button

  std::array<size_t, 5> count_cell_types() const;
  void check_cell_type_counts();
  std::array<size_t, 5> get_count_cell_types() const;

  void add_cells(const cell_type& focal_cell_type);

  void obtain_equilibrium();

  std::array<size_t, 5> num_cell_types;

  Param get_parameters();
  infection_routine get_infection_type() const;
  float get_percent_infected() const;
  void set_start_setup(start_type new_type);

 // size_t find_center(const cell_type& focal_cell_type) const;
  size_t find_central_cell(const cell_type& focal_cell_type) const;
  size_t find_central_cell(const std::vector< size_t >& positions) const;

  void test_change_cell_type(const size_t& pos, const cell_type& new_cell_type);
  void test_event(size_t event);
  void test_update_rates();
  void test_ask_infect_neighbours(size_t depth, size_t pos);
  void test_increase_t_cell_concentration(size_t pos);
  void test_diffuse();
  void test_infect_periphery(float frac);
  void test_infect_random(float frac);
  void test_infect_center(float frac);
  void test_infect_center_largest(float frac);
  void test_infect_all_cancer();
  void test_infect_long_distance(size_t pos);

  float get_rates(size_t event);
  size_t test_pick_event(const std::array<float, 8>& v, float s);

  float calc_max_t_cell_rate();

private:
  Param parameters;

  float total_t_cell_concentration;
  void diffuse();
  void increase_t_cell_concentration(size_t pos);

  std::array< float, 8> rates;

  std::vector<double> long_distance_infection_probability;

  void update_rates();
  void update_total_growth_prob();

  size_t pick_event(const std::array< float, 8>& rates, float sum);
  void do_event(size_t event);
  void implement_death(const cell_type& parent);
  void implement_growth(const cell_type& parent);
  void update_growth_prob(size_t pos);
  void update_death_prob(size_t pos);
  void update_death_prob_cancer(float t_cell_rate,
                                            size_t pos);

  void update_death_prob(size_t pos,
                         cell_type old_type,
                         cell_type new_type);

  void print_to_file(float t);

  void update_death_cdf(const cell_type& parent, size_t pos);
  void update_growth_cdf(size_t pos);
  void update_death_cdf_all();


  void infect_random(float fraction);
  void infect_center(float fraction);
  void infect_center_largest(float fraction);
  void infect_periphery(float fraction);
  void infect_all_cancer();

  void initialize_full();

  void infect_long_distance(size_t pos);

  void setup_voronoi(std::vector< std::vector< voronoi_point > >& all_polys,
                     grid_type used_grid_type);

  void change_cell_type(const size_t& pos, const cell_type& new_cell_type);

  void ask_infect_neighbours(size_t depth, size_t pos);

  void update_count(cell_type old_type, cell_type new_type);
  // float calc_t_cell_death_rate(float concentration);
};



#endif /* simulation_hpp */
