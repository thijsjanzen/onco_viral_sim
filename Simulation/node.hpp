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
  ~node();
  node(const node& other);

  void operator=(const node& other);

  cell_type node_type;
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


  void update_prob_of_growth();

  void die();
};


#endif /* node_hpp */
