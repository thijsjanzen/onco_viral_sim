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

enum cell_type {normal, cancer, infected, empty};

typedef struct node node;

struct node {
  node();
  node(size_t p);

  cell_type node_type;
  size_t pos;
  size_t x_;
  size_t y_;

  std::vector< node* > neighbors;

  void set_coordinates(size_t row_size);
  void update_neighbors(std::vector< node >& world,
                        size_t world_size);
  std::vector< cell_type > return_neighbor_types();

  std::array<float, 3> calc_prob_of_growth();
  float freq_type_neighbours(const cell_type& ref_type);

  void die();
};


#endif /* node_hpp */
