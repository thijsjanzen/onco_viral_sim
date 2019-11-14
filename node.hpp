//
//  node.hpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//

#ifndef node_hpp
#define node_hpp

#include <stdio.h>
#include <vector>
#include <memory>

enum cell_type {normal, cancer, infected, empty};

typedef struct node node;

struct node {
  node();

  cell_type node_type;

  std::vector< std::shared_ptr< node > > neighbors;
  std::vector< cell_type > neighbor_types;

  void set_coordinates(int x_, int y_);
  void update_neighbors(const std::vector< node >& world,
                       int i,
                       int world_size);
  std::vector< cell_type > return_neighbor_types();
  void update_neighbor_types();

  std::array<float, 3> calc_prob_of_growth();
  float freq_type_neighbours(const cell_type& ref_type);

  void die();
};


#endif /* node_hpp */
