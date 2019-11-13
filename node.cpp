//
//  node.cpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//
#include <assert.h>
#include "node.hpp"

node::node() {
  node_type = empty;
}

void node::set_coordinates(int x_
                           , int y_) {
  x = x_;
  y = y_;
}

void node::update_neigbors(const std::vector<std::vector< node >>& world) {
  size_t num_cells = world.size();
  for(int xi = x-1; xi <= x+1; ++xi) {
    for(int yi = y-1; yi <= y+1; ++yi) {
      if(xi != x && yi != y) {
        int other_x = xi % num_cells;
        int other_y = yi % num_cells;
        std::shared_ptr<node> neighbor = std::make_shared<node>(world[other_x][other_y]);
        neighbors.push_back(neighbor);
      }
    }
  }
  if(neighbors.size() != 4) {
    std::cout << "the number of neighbors is incorrect\n";
  }
}

void node::update_neighbor_types() {
  neighbor_types = std::vector< cell_type>(neighbors.size());
  int j = 0;
  for(auto i : neighbors) {
    neighbor_types[j] = i->node_type;
    j++;
  }
}


std::vector<cell_type> node::return_neighbor_types() {
  return neighbor_types;
}
