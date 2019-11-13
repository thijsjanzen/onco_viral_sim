//
//  simulation.cpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//

#include "simulation.hpp"

using namespace model;

void Simulation::run() {
  // do stuff

  std::array< float, 6> rates;
  float t = 0.f;
  while(t < parameters.maximum_time) {
    update_rates(rates);
    float lambda = std::accumulate(rates.begin(), rates.end(), 0.0);
    float dt = rndgen.Expon(lambda);
    int event = pick_event(rates, lambda);
    do_event(event);
    t += dt;
  }
}

void Simulation::initialize_network() {
  // we first make a regular network
  // we will do voronoi later...
  std::vector< node > temp(parameters.num_cells);
  world.resize(parameters.num_cells, temp);

  // now we have to loop!
  for(int i = 0; i < parameters.num_cells; ++i) {
    for(int j = 0; j < parameters.num_cells; ++j) {
      world[i][j].set_coordinates(i,j);
      world[i][j].update_neigbors(world);
      world[i][j].update_neighbor_types();
    }
  }

  update_cell_types();
}

void Simulation::update_cell_types() {
  for(auto i : world) {
    for(auto j : i) {
      num_cell_types[j.node_type]++;
    }
  }
}

void Simulation::do_event(int event) {

  switch(event) {
    case 0: {
      // birth normal
    }
    case 1: {
      // death normal
    }
    case 2: {
      // birth cancer
    }
    case 3: {
      // death cancer
    }
    case 4: {
      // birth infection
    }
    case 5: {
      // death infection
    }
    default: {
      // do nothing
    }
  }

  update_cell_counts(event);

  return;
}

void Simulation::update_rates(std::array< float, 6>& rates) {
  rates[0] = parameters.birth_normal * num_cell_types[0];
  rates[1] = parameters.death_normal * num_cell_types[0];

  rates[2] = parameters.birth_cancer * num_cell_types[1];
  rates[3] =  parameters.death_cancer * num_cell_types[1];

  rates[4] = parameters.birth_infected * num_cell_types[2];
  rates[5] = parameters.death_infected * num_cell_types[2];
}

int Simulation::pick_event(const std::array< float, 6>& rates, float sum) {
  float r = rndgen.uniform() * sum;
  for(size_t i = 0; i < rates.size(); ++i) {
    r -= rates[i];
    if(r <= 0) {
      return i;
    }
  }
  return -1;
}

void Simulation::update_cell_counts(int event) {
  switch(event) {
    case 0: {
      num_cell_types[0]++; // birth normal
    }
    case 1: {
      num_cell_types[0]--;  // death normal
      num_cell_types[3]++;
    }
    case 2: {
      num_cell_types[1]++;  // birth cancer
    }
    case 3: {
      num_cell_types[1]--;  // death cancer
      num_cell_types[3]++;
    }
    case 4: {
      num_cell_types[2]++;  // birth infection
      num_cell_types[1]--;
    }
    case 5: {
      num_cell_types[2]--;
      num_cell_types[3]++;
        // death infection
    }
    default: {
        // do nothing
    }
  }
}

Simulation::Simulation(const Param& param) {
  parameters = param;
  num_cell_types.resize(4, 0);
}
