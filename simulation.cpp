//
//  simulation.cpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//

#include <algorithm>
#include "simulation.hpp"




void simulation::run() {
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

void simulation::initialize_network() {
  // we first make a regular network
  // we will do voronoi later...
  
  for(int i = 0; i < world.size(); ++i) {
    world[i].update_neighbors(world, i, parameters.num_cells);
    world[i].update_neighbor_types();
  }

  update_cell_types();
  update_growth_probabilities();
}

void simulation::update_cell_types() {
  for(auto i : world) {
      num_cell_types[i.node_type]++;
  }
}

void simulation::update_growth_probabilities() {
  int pos = 0;

  max_growth_prob_infected = -1.f;
  max_growth_prob_normal   = -1.f;
  max_growth_prob_cancer   = -1.f;

  for(auto i : world) {
    std::array<float, 3> probs = i.calc_prob_of_growth();
    growth_prob_normal[pos]   = probs[0];
    growth_prob_cancer[pos]   = probs[1];
    growth_prob_infected[pos] = probs[2];

    if(probs[0] > max_growth_prob_normal)   max_growth_prob_normal   = probs[0];
    if(probs[1] > max_growth_prob_cancer)   max_growth_prob_cancer   = probs[1];
    if(probs[2] > max_growth_prob_infected) max_growth_prob_infected = probs[2];

    switch(i.node_type) {
      case normal:
        death_prob_normal[pos] = 1.f;
        break;
      case cancer:
        death_prob_cancer[pos] = 1.f;
        break;
      case infected:
        death_prob_infected[pos] = 1.f;
        break;
      case empty:
        break;
    }

    pos++;
  }
}


void simulation::do_event(int event) {

  switch(event) {
    case 0: {
      implement_growth(normal);       // birth normal
    }
    case 1: {
      implement_death(normal);        // death normal
      break;
    }
    case 2: {
      implement_growth(cancer);       // birth cancer
      break;
    }
    case 3: {
      implement_death(cancer);        // death cancer
      break;
    }
    case 4: {
      implement_growth(infected);     // birth infection
      break;
    }
    case 5: {
      implement_death(infected);     // death infection
      break;
    }
    default: {
      // do nothing
      break;
    }
  }

  update_cell_counts(event);

  return;
}

void simulation::update_rates(std::array< float, 6>& rates) {
  rates[0] = parameters.birth_normal   * std::accumulate(growth_prob_normal.begin(), growth_prob_normal.end(), 0.0);
  rates[1] = parameters.death_normal   * num_cell_types[0];

  rates[2] = parameters.birth_cancer   * std::accumulate(growth_prob_cancer.begin(), growth_prob_cancer.end(), 0.0);
  rates[3] = parameters.death_cancer   * num_cell_types[1];

  rates[4] = parameters.birth_infected * std::accumulate(growth_prob_infected.begin(), growth_prob_infected.end(), 0.0);
  rates[5] = parameters.death_infected * num_cell_types[2];
}

size_t simulation::pick_event(const std::array< float, 6>& rates, float sum) {
  float r = rndgen.uniform() * sum;
  for(size_t i = 0; i < rates.size(); ++i) {
    r -= rates[i];
    if(r <= 0) {
      return i;
    }
  }
  return -1;
}


void simulation::update_cell_counts(int event) {
  switch(event) {
    case 0: {
      num_cell_types[0]++; // birth normal
    }
    case 1: {
      num_cell_types[0]--;  // death normal
    }
    case 2: {
      num_cell_types[1]++;  // birth cancer
    }
    case 3: {
      num_cell_types[1]--;  // death cancer
    }
    case 4: {
      num_cell_types[2]++;  // birth infection
      num_cell_types[1]--;
    }
    case 5: {
      num_cell_types[2]--;
        // death infection
    }
    default: {
        // do nothing
    }
  }
}


void simulation::implement_death(const cell_type& parent) {
  int pos;
  switch(parent) {
    case normal:
      pos = rndgen.draw_from_dist(death_prob_normal.begin(), death_prob_normal.end(), 1.f);
      break;
    case cancer:
      pos = rndgen.draw_from_dist(death_prob_cancer.begin(), death_prob_cancer.end(), 1.f);
      break;
    case infected:
      pos = rndgen.draw_from_dist(death_prob_infected.begin(), death_prob_infected.end(), 1.f);
      break;
    case empty:
      std::cout << "ERROR! empty node is dying\n";
      pos = -1;
      break;
  }

  world[pos].die();
  death_prob_normal[pos] = 0.f;
  std::array<float, 3> probs = world[pos].calc_prob_of_growth();
  growth_prob_normal[pos]   = probs[0];
  growth_prob_cancer[pos]   = probs[1];
  growth_prob_infected[pos] = probs[2];

  // TODO: update growth probabilities of neighbors.
}

void simulation::implement_growth(const cell_type& parent) {
  int pos;
  switch(parent) {
    case normal:
      pos = rndgen.draw_from_dist(growth_prob_normal.begin(), growth_prob_normal.end(), max_growth_prob_normal);
      break;
    case cancer:
      pos = rndgen.draw_from_dist(growth_prob_cancer.begin(), growth_prob_cancer.end(), max_growth_prob_cancer);
      break;
    case infected:
      pos = rndgen.draw_from_dist(growth_prob_infected.begin(), growth_prob_infected.end(), max_growth_prob_infected);
      break;
    case empty:
      std::cout << "ERROR! empty node is growing\n";
      pos = -1;
      break;
  }

  // now we have a node picked, and it will become the new type
  world[pos].node_type = parent;
  death_prob_normal[pos] = 1.f;
  std::array<float, 3> probs = world[pos].calc_prob_of_growth();
  growth_prob_normal[pos]   = probs[0];
  growth_prob_cancer[pos]   = probs[1];
  growth_prob_infected[pos] = probs[2];

  for(auto i : world[pos].neighbors) {
    (*i).update_neighbor_types();
      // TODO: update growth probabilities of neighbors.
  }

}



simulation::simulation(const Param& param) {
  parameters = param;

  world.resize(parameters.num_cells);
  growth_prob_normal.resize(  parameters.num_cells, 0.f);
  growth_prob_cancer.resize(  parameters.num_cells, 0.f);
  growth_prob_infected.resize(parameters.num_cells, 0.f);

  death_prob_normal.resize(  parameters.num_cells, 0.f);
  death_prob_cancer.resize(  parameters.num_cells, 0.f);
  death_prob_infected.resize(parameters.num_cells, 0.f);

  max_growth_prob_cancer = 1.f;
  max_growth_prob_normal = 1.f;
  max_growth_prob_infected = 1.f;
}
