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

struct plot_edge {
    float x0_, y0_, x1_, y1_;
    int left;
    int right;
    plot_edge(float x0, float y0, float x1, float y1, int l, int r) :
        x0_(x0), y0_(y0), x1_(x1), y1_(y1), left(l), right(r) {
    }
};


struct node {
  node();
  node(size_t p, float norm_infection_rate);
  node(size_t p, float norm_infected, float x_, float y_);

  cell_type node_type;
  size_t pos;
  float x_;
  float y_;

  float prob_normal_infected;
  std::vector< node* > neighbors;

  std::vector< plot_edge > edges; // x0, y0, x1, y1
  void sort_edges();

  void set_coordinates(size_t row_size);
  void update_neighbors(std::vector< node >& world,
                        size_t world_size);

  std::vector< cell_type > return_neighbor_types();

  void add_neighbor(std::vector< node >& world,
                    size_t other_pos);

  std::array<float, 4> calc_prob_of_growth();
  float freq_type_neighbours(const cell_type& ref_type);

  void die();
};

void sort_edges(std::vector< plot_edge>& edges);

#endif /* node_hpp */
