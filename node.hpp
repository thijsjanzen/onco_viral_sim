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
  int x;
  int y;
  std::vector< std::shared_ptr< node > > neighbors;
  std::vector< cell_type > neighbor_types;

  void set_coordinates(int x_, int y_);
  void update_neigbors(const std::vector< std::vector< node >>& world);
  std::vector< cell_type > return_neighbor_types();
  void update_neighbor_types();
};


#endif /* node_hpp */
