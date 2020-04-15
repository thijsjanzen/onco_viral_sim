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
#include <iostream>

enum cell_type {normal, cancer, infected, resistant, empty, max_num};

typedef struct node node;

struct voronoi_point {
    double x_, y_;

    voronoi_point() {
        x_ = 0.0; y_ = 0.0;
    }

    voronoi_point& operator=(const voronoi_point& other) {
      x_ = other.x_;
      y_ = other.y_;
      return *this;
    }

    voronoi_point(const voronoi_point& other) {
        x_ = other.x_;
        y_ = other.y_;
    }

    voronoi_point(double x, double y) : x_(x), y_(y) {}

    bool operator==(const voronoi_point& other) const {
        // explicitly, this doesn't keep track of left and right!
      //  if(x_ != other.x_) return false;
      //  if(y_ != other.y_) return false;
      if(abs(x_ - other.x_) > 1e-3) return false;
      if(abs(y_ - other.y_) > 1e-3) return false;

      return true;
    }

    bool operator<(const voronoi_point& other) const {
       // if(fabs(x_ - other.x_) < 1e-4f) return y_ < other.y_;
       // return x_ < other.x_;
       if(x_ == other.x_) return y_ < other.y_;
       return x_ < other.x_;
    }

    bool operator!=(const voronoi_point& other) const {
        return !(*this == other);
    }
};

struct voronoi_edge {
  voronoi_edge() {
  }

    voronoi_edge(voronoi_point s, voronoi_point e,
                 size_t l, size_t r) : start(s), end(e), left(l), right(r) {}

    voronoi_point start;
    voronoi_point end;
    size_t left;
    size_t right;

    voronoi_edge& operator=(const voronoi_edge& other) {
      left  = other.left;
      right = other.right;
      start = other.start;
      end   = other.end;
      return *this;
    }

    voronoi_edge(const voronoi_edge& other) {
      left  = other.left;
      right = other.right;
      start = other.start;
      end   = other.end;
    }

    bool operator<(const voronoi_edge& other) const {
      return start < other.start;
    }

    bool operator==(const voronoi_edge& other) const {
     if(start != other.start) return false;
     if(end   != other.end)   return false;
     return true;
    }

    bool operator!=(const voronoi_edge& other) const {
        return !(*this == other);
    }

    bool check() const {
      if(std::isnan(start.x_)) {
          std::cout << "start.x_ isnan\n";
          return false;
      }
      if(std::isnan(start.y_)) {
          std::cout << "start.y_ isnan\n";
          return false;
      }
      if(std::isnan(end.x_)) {
          std::cout << "end.x_ isnan\n";
          return false;
      }
      if(std::isnan(end.y_)) {
          std::cout << "end.y_ isnan\n";
          return false;
      }
      return true;
    }

    double calc_dist() {
      double x_x = (start.x_ - end.x_) * (start.x_ - end.x_);
      double y_y = (start.y_ - end.y_) * (start.y_ - end.y_);
      return(std::sqrt(x_x + y_y));
    }
};

std::vector< voronoi_point> clean_edges(const std::vector< voronoi_edge >& input_edges,
                                        size_t pos);
void invert_edges(std::vector< voronoi_edge>& edges, size_t pos);

struct node {
  node(node&&) = delete;
  const node& operator=(node&&) = delete;
  node(const node&) = delete;
  const node& operator=(const node&) = delete;

  node();

  size_t pos;
  float x_;
  float y_;

  float inv_num_neighbors;
  float prob_normal_infected;
  float t_cell_concentration;

  std::vector< node* > neighbors;

  void check_distances(float max_val);
  float calc_distance(const node& other);

  void set_coordinates(size_t row_size);
  void update_neighbors(std::vector< node >& world,
                        size_t world_size);

  void add_t_cell(float amount);


  std::vector< cell_type > return_neighbor_types();

  void add_neighbor(std::vector< node >& world,
                    size_t other_pos);

  std::array<float, 4> calc_prob_of_growth() const;
  float freq_type_neighbours(const cell_type& ref_type) const;

  void die();

  cell_type get_cell_type() const {
    return node_type;
  }

  void set_cell_type(cell_type new_type) {
    node_type = new_type;
  }

private:
    cell_type node_type;
};

#endif /* node_hpp */
