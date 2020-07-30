//
//  node.cpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//
#include <cassert>
#include <iostream>
#include <array>
#include <cmath>
#include <array>
#include <fstream>

#include "node.hpp"
#include "random_thijs.hpp"

node::node() {
  node_type = empty;
  pos = 0;
  x_ = 0;
  y_ = 0;
  inv_num_neighbors = 0.f;
  prob_normal_infected = 0.f;
  t_cell_concentration = 0.f;
}

void node::add_neighbor(std::vector< node >& world,
                        size_t other_pos) {

    if(!neighbors.empty()) {
        for(auto i : neighbors) {
            if(i->pos == other_pos) {
                return; // if the neighbor is already added, don't add again.
            }
        }
    }

    node* neighbor = &world[static_cast<size_t>(other_pos)];
    neighbors.push_back(neighbor);
}


bool within_limits(int x, int y, int z, int limit) {
  if(x < 0) return false;
  if(y < 0) return false;
  if(z < 0) return false;
  if(x >= limit) return false;
  if(y >= limit) return false;
  if(z >= limit) return false;

  return true;
}

void node::update_neighbors(std::vector< node >& world,
                            size_t world_size) {

    static int relative_points[6][3] = { {-1,  0,  0},
                                         { 1,  0,  0},
                                         { 0,  1,  0},
                                         { 0, -1,  0},
                                         { 0,  0,  1},
                                         { 0,  0, -1} };

    for(int i = 0; i < 6; ++i) {
        int other_x = static_cast<int>(x_) + relative_points[i][0];
        int other_y = static_cast<int>(y_) + relative_points[i][1];
        int other_z = static_cast<int>(z_) + relative_points[i][2];
        if(within_limits(other_x, other_y, other_z, world_size)) {
            int other_pos = other_y + (other_x + other_z * static_cast<int>(world_size))
                                * static_cast<int>(world_size);

            if(other_pos >= 0 && other_pos < static_cast<int>(world.size())) {
                node* neighbor = &world[static_cast<size_t>(other_pos)];

                assert(static_cast<int>(neighbor->x_) == other_x);
                assert(static_cast<int>(neighbor->y_) == other_y);
                assert(static_cast<int>(neighbor->z_) == other_z);
                neighbors.push_back(neighbor);
            }
        }
    }
    inv_num_neighbors = 1.f / neighbors.size();
}

std::array<float, 4> node::calc_prob_of_growth() const {
  std::array<float, 4> prob_of_growth = {0.f, 0.f, 0.f, 0.f};

  if(node_type == normal) {
    // infected can grow into this, but at lower frequency:
    prob_of_growth[infected]  =  prob_normal_infected * freq_type_neighbours(infected);
  }
  if(node_type == cancer) {
    // infected nodes can grow into this
    prob_of_growth[infected]  = freq_type_neighbours(infected);
  }
  if(node_type == empty) {
    prob_of_growth[normal]    = freq_type_neighbours(normal);
    prob_of_growth[cancer]    = freq_type_neighbours(cancer);
    prob_of_growth[resistant] = freq_type_neighbours(resistant);
  }

  return prob_of_growth;
}

float node::freq_type_neighbours(const cell_type& ref_type) const {
  int count = 0;
  for(auto i : neighbors) {
    if(i->get_cell_type() == ref_type) count++;
  }
  return static_cast<float>(count * inv_num_neighbors);
}

void node::set_coordinates(size_t row_size) {
  z_ = pos / (row_size * row_size);
  size_t pos2 = pos - z_ * row_size * row_size;
  x_ = pos2 / row_size;
  y_ = pos2 % row_size;
}

void node::add_t_cell(float amount) {
  t_cell_concentration += amount;
}

void invert_edges(std::vector< voronoi_edge>& edges, size_t pos) {
    // check for inverted edges.
    for(auto& i : edges) {
        if(i.right != pos) {
            size_t tmp = i.right;
            i.right = i.left;
            i.left = tmp;

            auto tmp_2 = i.start;
            i.start = i.end;
            i.end = tmp_2;
        }
    }
}

std::vector< double > calc_dist(const std::vector< voronoi_edge >& v,
                                const voronoi_edge& point) {
  std::vector< double > output(v.size());
  for(size_t i = 0; i < v.size(); ++i) {
    double dx = point.end.x_ - v[i].start.x_;
    double dy = point.end.y_ - v[i].start.y_;
    output[i] = dx * dx + dy * dy;
  }
  return output;
}

std::vector< voronoi_point> clean_edges(const std::vector< voronoi_edge >& input_edges,
                                        size_t pos) {
    // the goal is to connect all edges, and then provide
    // all the outer points
    std::vector< voronoi_edge > edges = input_edges;
    invert_edges(edges, pos);

    std::vector< voronoi_edge > new_edges;

    voronoi_edge focal_edge = edges.back();
    new_edges.push_back(focal_edge);

    while(!edges.empty()) {
       std::vector< double > distances = calc_dist(edges, focal_edge);
       auto min_val = std::min_element(distances.begin(), distances.end());
       auto match = std::distance(distances.begin(), min_val);

       focal_edge = edges[match];
       edges[match] = edges.back();
       edges.pop_back();
       new_edges.push_back(focal_edge);
    }
    // allright, we have connected them now
    // now we collect the starting points
    std::vector< voronoi_point > outer_points;
    for(auto i : new_edges) {
        voronoi_point temp_point = i.start;
        outer_points.push_back(temp_point);
    }
    return outer_points;
}

std::vector< size_t > node::get_cancer_neighbours() const {
  std::vector< size_t > output;
  for(auto i : neighbors) {
    if(i->get_cell_type() == cancer) {
      output.push_back(i->pos);
    }
  }
  return output;
}

float node::calc_t_cell_added_death_rate(float t_cell_rate,
                                         float t_cell_density_scaler) const {
  if (t_cell_concentration < 1e-5f) {
      return 0.f;
   }

  float added_t_cell_death_rate = expf(t_cell_rate * t_cell_concentration);

  float mult = 1.0f - t_cell_density_scaler *
                      freq_type_neighbours(cancer);
  if(mult < 0.f) mult = 0.f;
  float output = mult * added_t_cell_death_rate;
  if (std::isinf(output)) {
      output = 1e6;
  }
  return output;
}
