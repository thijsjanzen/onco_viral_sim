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
}

node::node(size_t p, float norm_infected) :
    pos(p),
    prob_normal_infected(norm_infected) {
    node_type = empty;
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
}

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
  }
}

void node::set_coordinates(size_t row_size) {
    x_ = pos / row_size;
    y_ = pos % row_size;
}



