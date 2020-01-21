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
#include <cmath>

enum cell_type {normal, cancer, infected, resistant, empty};

typedef struct node node;

struct voronoi_point {
    float x_, y_;

    voronoi_point() {
        x_ = 0.f; y_ = 0.f;
    }

    voronoi_point(const voronoi_point& other) {
        x_ = other.x_;
        y_ = other.y_;
    }

    voronoi_point(float x, float y) : x_(x), y_(y) {}





    bool operator=(const voronoi_point& other) {
        x_ = other.x_;
        y_ = other.y_;
    }

    bool operator==(const voronoi_point& other) const {
        // explicitly, this doesn't keep track of left and right!
        if(fabs(x_ - other.x_) > 1e-4f) return false;
        if(fabs(y_ - other.y_) > 1e-4f) return false;
        return true;
    }

    bool operator<(const voronoi_point& other) const {
        if(fabs(x_ - other.x_) < 1e-6f) return y_ < other.y_;
        return x_ < other.x_;
    }

    bool operator!=(const voronoi_point& other) const {
        return !(*this == other);
    }
};

struct voronoi_edge {

    voronoi_edge(voronoi_point s, voronoi_point e,
                 int l, int r) : start(s), end(e), left(l), right(r) {}

    voronoi_point start;
    voronoi_point end;
    int left;
    int right;

    bool operator<(const voronoi_edge& other) const {
        if(start.x_ == other.start.x_) return start.y_ < other.start.y_;
        return start.x_ < other.start.x_;
    }

    bool operator==(const voronoi_edge& other) const {
        // explicitly, this doesn't keep track of left and right!
        if(fabs(start.x_ - other.start.x_) > 1e-6f) return false;
        if(fabs(start.y_ - other.start.y_) > 1e-6f) return false;
        if(fabs(end.x_ - other.end.x_) > 1e-6f) return false;
        if(fabs(end.y_ - other.end.y_) > 1e-6f) return false;
        return true;
    }

    bool operator!=(const voronoi_edge& other) const {
        return !(*this == other);
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

  std::vector< voronoi_edge > edges;
  std::vector< voronoi_point > outer_points;
  void clean_edges();
  void invert_edges();
  void check_distances(float max_val);
  float calc_distance(const node& other);

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

#endif /* node_hpp */
