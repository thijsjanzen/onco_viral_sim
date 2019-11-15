//
//  simulation.cpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//

#include <algorithm>
#include "simulation.hpp"
#include <chrono>



void simulation::run() {
  // do stuff

  std::array< float, 6> rates;
  float t = 0.f;
  auto start = std::chrono::system_clock::now();
  while(t < parameters.maximum_time) {
    update_rates(rates);
    float lambda = std::accumulate(rates.begin(), rates.end(), 0.0);
    float dt = rndgen.Expon(lambda);
    int event = pick_event(rates, lambda);
    do_event(event);

    if(t < parameters.time_adding_cancer && t+dt >= parameters.time_adding_cancer) {
      add_cells(cancer);
    }

    if((int)t < (int)(t+dt)) {
      print_to_file(t);
      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;

      std::cout << t << "\t" << elapsed_seconds.count() << "\n";

      start = std::chrono::system_clock::now();
    }

    t += dt;

    if(num_cell_types[0] == 0 ) {
      std::cout << "all normal cells are dead\n";
      return;
    }
  }
}

void simulation::initialize_network() {
  // we first make a regular network
  // we will do voronoi later...

  for(size_t i = 0; i < world.size(); ++i) {
    world[i].update_neighbors(world, i, sq_size);
  }

  add_cells(normal);

  for(auto it = world.begin(); it != world.end(); ++it) {
    (*it).update_neighbor_types();
  }

  count_cell_types();
  update_growth_probabilities();
}

void simulation::count_cell_types() {
  num_cell_types = {0, 0, 0};
  for(auto i : world) {
      num_cell_types[i.node_type]++;
  }
}

void simulation::update_growth_probabilities() {

  max_growth_prob = {-1.f, -1.f, -1.f};

  for(auto i : world) {
    std::array<float, 3> probs = i.calc_prob_of_growth();

    int pos = i.pos;

    growth_probs[normal][pos] = probs[normal];
    growth_probs[cancer][pos] = probs[cancer];
    growth_probs[infected][pos] = probs[infected];

    if(probs[normal] > max_growth_prob[normal])   max_growth_prob[normal]   = probs[normal];
    if(probs[cancer] > max_growth_prob[cancer])   max_growth_prob[cancer]   = probs[cancer];
    if(probs[infected] > max_growth_prob[infected]) max_growth_prob[infected] = probs[infected];

    switch(i.node_type) {
      case normal:
        death_probs[normal][pos] = 1.f;
        break;
      case cancer:
        death_probs[cancer][pos] = 1.f;
        break;
      case infected:
        death_probs[infected][pos] = 1.f;
        break;
      case empty:
        break;
    }
  }
  return;
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
  rates[0] = parameters.birth_normal   * std::accumulate(growth_probs[normal].begin(), growth_probs[normal].end(), 0.0);
  rates[1] = parameters.death_normal   * num_cell_types[normal];

  rates[2] = parameters.birth_cancer   * std::accumulate(growth_probs[cancer].begin(), growth_probs[cancer].end(), 0.0);
  rates[3] = parameters.death_cancer   * num_cell_types[cancer];

  rates[4] = parameters.birth_infected * std::accumulate(growth_probs[infected].begin(), growth_probs[infected].end(), 0.0);
  rates[5] = parameters.death_infected * num_cell_types[infected];
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
      num_cell_types[normal]++;    // birth normal
      break;
    }
    case 1: {
      num_cell_types[normal]--;    // death normal
      break;
    }
    case 2: {
      num_cell_types[cancer]++;    // birth cancer
      break;
    }
    case 3: {
      num_cell_types[cancer]--;    // death cancer
      break;
    }
    case 4: {
      num_cell_types[infected]++;  // birth infection
      num_cell_types[cancer]--;
      break;
    }
    case 5: {
      num_cell_types[infected]--;  // death infection
      break;
    }
    default: {
        // do nothing
      break;
    }
  }
}


void simulation::implement_death(const cell_type& parent) {
  int position_of_dying_cell = -1;
  switch(parent) {
    case normal:
      position_of_dying_cell = rndgen.draw_from_dist(death_probs[normal].begin(), death_probs[normal].end(), 1.f);
      break;
    case cancer:
      position_of_dying_cell = rndgen.draw_from_dist(death_probs[cancer].begin(), death_probs[cancer].end(), 1.f);
      break;
    case infected:
      position_of_dying_cell = rndgen.draw_from_dist(death_probs[infected].begin(), death_probs[infected].end(), 1.f);
      break;
    case empty:
      std::cout << "ERROR! empty node is dying\n";
      position_of_dying_cell = -1;
      break;
  }

  world[position_of_dying_cell].node_type = empty;
  update_growth_prob(position_of_dying_cell);
  for(auto i : world[position_of_dying_cell].neighbors) {
    update_growth_prob(i->pos);
  }
}

void simulation::implement_growth(const cell_type& parent) {
  int position_of_grown_cell = -1;
  switch(parent) {
    case normal:
      position_of_grown_cell = rndgen.draw_from_dist(growth_probs[normal].begin(), growth_probs[normal].end(), max_growth_prob[normal]);
      break;
    case cancer:
      position_of_grown_cell = rndgen.draw_from_dist(growth_probs[cancer].begin(), growth_probs[cancer].end(), max_growth_prob[cancer]);
      break;
    case infected:
      position_of_grown_cell = rndgen.draw_from_dist(growth_probs[infected].begin(), growth_probs[infected].end(), max_growth_prob[infected]);
      break;
    case empty:
      std::cout << "ERROR! empty node is growing\n";
      position_of_grown_cell = -1;
      break;
  }

  world[position_of_grown_cell].node_type = parent;
  update_growth_prob(position_of_grown_cell);
  for(auto i : world[position_of_grown_cell].neighbors) {
    update_growth_prob(i->pos);
  }
}

void simulation::add_cells(cell_type focal_cell_type) {
    // pick center node:
  int x = sq_size / 2;
  int y = sq_size / 2;
  int focal_pos = x * sq_size + y;
  std::vector<int> cells_turned(1, focal_pos);
  world[focal_pos].node_type = focal_cell_type;
  int counter = 0;
  int total_new_cells = 0;
  int max_number_of_cells = parameters.initial_number_cancer_cells;
  if(focal_cell_type == normal) max_number_of_cells = parameters.num_cells * 0.1;
  while(total_new_cells < max_number_of_cells) {
    focal_pos = cells_turned[counter];
    for(size_t i = 0; i < world[focal_pos].neighbors.size(); ++i) {
      int other_pos = world[focal_pos].neighbors[i]->pos;
      if(world[other_pos].node_type != focal_cell_type) {
        total_new_cells++;
        world[other_pos].node_type = focal_cell_type;
        cells_turned.push_back(other_pos);
      }
    }
    counter++;
  }

  for(auto i : cells_turned) {
    update_growth_prob(i);
    for(auto j : world[i].neighbors) {
      update_growth_prob(j->pos);
    }
  }

}



void simulation::update_growth_prob(int pos) {
  std::array<float, 3> probs = world[pos].calc_prob_of_growth();
  for(int i = 0; i < 3; ++i) {
    if(growth_probs[i][pos] == max_growth_prob[i] && probs[i] < max_growth_prob[i]) {
      growth_probs[i][pos] = probs[i];
      auto max_val = std::max_element(growth_probs[i].begin(), growth_probs[i].end());
      max_growth_prob[i] = *max_val;
    } else {
      growth_probs[i][pos] = probs[i];
    }
  }
}

simulation::simulation(const Param& param) {
  parameters = param;

  sq_size = sqrt(parameters.num_cells);

  world.resize(parameters.num_cells);
  for(size_t i = 0; i < world.size(); ++i) {
    world[i].pos = i;
  }
  std::vector< float > temp(parameters.num_cells, 0.f);

  growth_probs.resize(3, temp);
  death_probs.resize(3, temp);

  max_growth_prob = {0.f, 0.f, 0.f};
}

void simulation::print_to_file(float t) {
  std::string file_name = "world_file_" + std::to_string((int)t) + ".txt";
  std::ofstream outfile(file_name);
  for(size_t i = 0; i < world.size(); ++i) {
    if(i % sq_size == 0) outfile << "\n";
    char output = '0';
    switch(world[i].node_type)  {
      case normal:
        output = 'N'; break;
      case cancer:
        output = 'C'; break;
      case infected:
        output = 'I'; break;
      case empty:
        output = '0'; break;
    }
    outfile << output << " ";

  }
  outfile.close();
}
