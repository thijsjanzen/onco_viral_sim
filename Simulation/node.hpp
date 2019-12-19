//
//  node.hpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//

#ifndef node_hpp
#define node_hpp

#include <cstdio>
#include <vector>
#include <memory>

enum cell_type {normal, cancer, infected, resistant, empty};

typedef struct node node;

struct node {
  node();
  node(size_t p, float norm_infection_rate);
  node(const node& other);

  void operator=(const node& other);


  size_t pos;
  size_t x_;
  size_t y_;

  float prob_normal_infected;
  std::vector< node* > neighbors;
  std::vector< float > freq_type_neighbours;
  std::vector< float > prob_of_growth;


  void set_coordinates(size_t row_size);
  void update_neighbors(std::vector< node >& world,
                        size_t world_size);

  void update_freq_neighbours();

  std::vector< cell_type > return_neighbor_types();

  std::vector< float> prob_of_growth_;
  std::vector< float > neighbor_freqs;

 // std::vector< float> calc_prob_of_growth();
  void calc_prob_of_growth();
  float freq_type_neighbors(const cell_type& ref_type);

  void update_neighbor_freq(const cell_type& old_type,
                            const cell_type& new_type);
  void set_all_neighbor_freqs();

  float neighbor_freq; // frequency of one neighbor, e.g. 1 / num neighbors


  void die();

  void set_node_type(const cell_type& new_type);
  cell_type get_node_type() const;
  void initialize_node_type(const cell_type& new_type);
private:
    cell_type node_type;
};


#endif /* node_hpp */
