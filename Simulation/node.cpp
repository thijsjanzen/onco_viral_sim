//
//  node.cpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//
#include <assert.h>
#include "node.hpp"
#include <iostream>
#include <array>

node::node() {
  node_type = empty;
}

node::node(size_t p) : pos(p) {
  node_type = empty;
}

void node::update_neighbors(std::vector< node >& world,
                            size_t world_size) {

  for(int i = -1; i <= 1; i++) {
    for(int j = -1; j <= 1; j++) {
      if(abs(i - j) == 1) {
        int other_x = static_cast<int>(x_) + i;
        int other_y = static_cast<int>(y_) + j;
        int other_pos = other_x * static_cast<int>(world_size) + other_y;
        if(other_pos >= 0 && other_pos < static_cast<int>(world.size())) {
          node* neighbor = &world[static_cast<size_t>(other_pos)];
          neighbors.push_back(neighbor);
        }
      }
    }
  }
}

void node::update_neighbor_types() {
  neighbor_types = std::vector< cell_type>(neighbors.size());
  size_t j = 0;
  for(auto it = neighbors.begin(); it != neighbors.end(); ++it) {
    cell_type other = (*it)->node_type;
    neighbor_types[j] = other;
    j++;
  }
}

std::array<float, 3> node::calc_prob_of_growth() {
  update_neighbor_types(); // this may be removed later perhaps
  std::array<float, 3> prob_of_growth = {0.f, 0.f, 0.f};

  if(node_type == normal) {
    // no other type can grow into this
  }
  if(node_type == cancer) {
    // infected nodes can grow into this
    prob_of_growth[infected] = freq_type_neighbours(infected);
  }
  if(node_type == empty) {
    prob_of_growth[normal] = freq_type_neighbours(normal);
    prob_of_growth[cancer] = freq_type_neighbours(cancer);
  }
  return prob_of_growth;
}

float node::freq_type_neighbours(const cell_type& ref_type) {
  int count = 0;
  for(auto i : neighbors) {
    if(i->node_type == ref_type) count++;
  }
  return 1.f * count / neighbors.size();
}

void node::set_coordinates(size_t row_size) {
    x_ = pos / row_size;
    y_ = pos % row_size;
}



