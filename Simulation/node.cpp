//
//  node.cpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//
#include <cassert>
#include "node.hpp"
#include <iostream>
#include <array>
#include <cmath>
#include <array>
#include "random_thijs.hpp"

node::node() {
  node_type = empty;
  prob_of_growth_ = std::vector<float>(4, 0.f);
  neighbor_freqs = std::vector<float>(5, 0.f);
  neighbor_freq = 0.25;
}

node::node(size_t p, float norm_infected) :
    pos(p),
    prob_normal_infected(norm_infected) {
    node_type = empty;
    prob_of_growth_ = std::vector<float>(4, 0.f);
    neighbor_freqs = std::vector<float>(5, 0.f);
    neighbor_freq = 0.25;
}

node::~node() {
    std::cout << "destructor\n";
    neighbors.clear();
}

node::node(const node& other) {
    node_type = other.node_type;
    pos = other.pos;
    x_ = other.x_;
    y_ = other.y_;

    prob_normal_infected = other.prob_normal_infected;
    freq_type_neighbours = other.freq_type_neighbours;
    prob_of_growth = other.prob_of_growth;

    neighbors = other.neighbors;
}

void node::operator=(const node& other) {
    node_type = other.node_type;
    pos = other.pos;
    x_ = other.x_;
    y_ = other.y_;

    prob_normal_infected = other.prob_normal_infected;
    freq_type_neighbours = other.freq_type_neighbours;
    prob_of_growth = other.prob_of_growth;

    neighbors = other.neighbors;
}

void node::update_neighbors(std::vector< node >& world,
                            size_t world_size) {


    static int relative_points[4][2] = { {-1, 0},
                                        {1, 0},
                                        {0, 1},
                                        {0, -1} };

    for(int i = 0; i < 4; ++i) {
        int other_x = static_cast<int>(x_) + relative_points[i][0];
        int other_y = static_cast<int>(y_) + relative_points[i][1];
        if(other_x >= 0 &&
           other_y >= 0 &&
           other_x < static_cast<int>(world_size) &&
           other_y < static_cast<int>(world_size)) {
            int other_pos = other_y + other_x * static_cast<int>(world_size);
            if(other_pos >= 0 && other_pos < static_cast<int>(world.size())) {
                node* neighbor = &world[static_cast<size_t>(other_pos)];

                assert(static_cast<int>(neighbor->x_) == other_x);
                assert(static_cast<int>(neighbor->y_) == other_y);
                neighbors.push_back(neighbor);
            }
        }
    }
    neighbor_freq = 1.0f / neighbors.size();
}

<<<<<<< Updated upstream
void node::update_freq_neighbours() {
    freq_type_neighbours = {0.f, 0.f, 0.f, 0.f};
    float one = 0.25;
    for(auto i : neighbors) {
        freq_type_neighbours[i->node_type] += one;
    }
}


void node::update_prob_of_growth() {
  prob_of_growth = {0.f, 0.f, 0.f, 0.f};

  update_freq_neighbours();

  if(node_type == normal) {
    // infected can grow into this, but at lower frequency:
    prob_of_growth[infected]  =  prob_normal_infected * freq_type_neighbours[infected];
  }
  if(node_type == cancer) {
    // infected nodes can grow into this
    prob_of_growth[infected]  = freq_type_neighbours[infected];
  }
  if(node_type == empty) {
    prob_of_growth[normal]    = freq_type_neighbours[normal];
    prob_of_growth[cancer]    = freq_type_neighbours[cancer];
    prob_of_growth[resistant] = freq_type_neighbours[resistant];
=======
void node::calc_prob_of_growth() {
  std::fill(prob_of_growth_.begin(), prob_of_growth_.end(), 0.f);

  if(node_type == normal) {
    // infected can grow into this, but at lower frequency:
    prob_of_growth_[infected]  =  prob_normal_infected * neighbor_freqs[infected];
  }
  if(node_type == cancer) {
    // infected nodes can grow into this
    prob_of_growth_[infected]  = neighbor_freqs[infected];
  }
  if(node_type == empty) {
    prob_of_growth_[normal]    = neighbor_freqs[normal];
    prob_of_growth_[cancer]    = neighbor_freqs[cancer];
    prob_of_growth_[resistant] = neighbor_freqs[resistant];
  }

  return;
}

float node::freq_type_neighbors(const cell_type& ref_type) {
  int count = 0;
  for(auto i : neighbors) {
    if(i->node_type == ref_type) count++;
>>>>>>> Stashed changes
  }
}

void node::set_coordinates(size_t row_size) {
    x_ = pos / row_size;
    y_ = pos % row_size;
}

void node::update_neighbor_freq(const cell_type& old_type,
                                const cell_type& new_type)  {
    if(new_type == old_type) return;

    neighbor_freqs[old_type] -= neighbor_freq;
    neighbor_freqs[new_type] += neighbor_freq;

    assert(neighbor_freqs[old_type] > -1e-6f);
    assert((neighbor_freqs[0] + neighbor_freqs[1] + neighbor_freqs[2] +
           neighbor_freqs[3] + neighbor_freqs[4]) >= (1.f - 1e-6f));

    if(neighbor_freqs[old_type] < 0.f)  neighbor_freqs[old_type] = 0.f;
}

void node::set_node_type(const cell_type& new_type){
    cell_type old_type = node_type;
    node_type = new_type;
    for(auto i: neighbors) {
        i->update_neighbor_freq(old_type, new_type);
    }
}

void node::initialize_node_type(const cell_type& new_type) {
    node_type = new_type;
}



cell_type node::get_node_type() const {
    return node_type;
}


void node::set_all_neighbor_freqs() {
   // std::fill(neighbor_freqs.begin(), neighbor_freqs.end(), 0);
    neighbor_freqs.assign(5, 0.f); // 5 types, although #5 is never used
    for(auto i : neighbors) {
        neighbor_freqs[ i->node_type]+= neighbor_freq;
    }
    assert((neighbor_freqs[0] + neighbor_freqs[1] + neighbor_freqs[2] +
           neighbor_freqs[3] + neighbor_freqs[4]) == 1.f);
}


