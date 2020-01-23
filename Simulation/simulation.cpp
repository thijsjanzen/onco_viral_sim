//
//  simulation.cpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>

#include "simulation.hpp"

void simulation::update_one_step() {
    update_rates();
    float lambda = std::accumulate(rates.begin(), rates.end(), 0.0f);

    float dt = rndgen.Expon(lambda);
    while(std::isinf(dt)) {
        static int counter = 0;
        dt = rndgen.Expon(lambda);
        counter++;
        if(counter > 100) {
            std::cout << "could not generate non infinite dt\n";
            exit(1);
        }
    }

    size_t event = pick_event(rates, lambda);
    do_event(event);

    if(parameters.start_setup == grow) {
        if(t < parameters.time_adding_cancer && t+dt >= parameters.time_adding_cancer) {
          add_cells(cancer);
        }

        if(t < parameters.time_adding_virus && t+dt >= parameters.time_adding_virus) {
          add_infected();
        }
    }

    t += dt;
}

void simulation::do_event(size_t event) {

  switch(event) {
    case 0: {
      implement_growth(normal);       // birth normal
      break;
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

    case 6: {
        implement_growth(resistant);     // birth infection
        break;
    }
    case 7: {
        implement_death(resistant);     // death infection
        break;
    }

    default: {
      // do nothing
      break;
    }
  }
}

void simulation::update_cell(size_t pos) {
  update_death_prob(pos);
  update_growth_prob(pos);
  for(auto i : world[pos].neighbors) {
    update_growth_prob(i->pos);
  }
}

void simulation::implement_death(const cell_type& parent) {
  size_t position_of_dying_cell = death_prob[parent].draw_explicit(rndgen);

  cell_type previous_type = world[position_of_dying_cell].node_type;

  world[position_of_dying_cell].node_type = empty;

  //update_cell(position_of_dying_cell);
  update_death_prob(position_of_dying_cell, previous_type, empty);
  update_growth_prob(position_of_dying_cell);
  for(auto i : world[position_of_dying_cell].neighbors) {
    update_growth_prob(i->pos);
  }

  if(parent == infected) {
      infect_long_distance(position_of_dying_cell);
  }
}

void simulation::implement_growth(const cell_type& parent) {

   // drawing position of growth is slow/bottleneck:
  size_t position_of_grown_cell = growth_prob[parent].draw_explicit(rndgen);

  cell_type previous_type = world[position_of_grown_cell].node_type;

  world[position_of_grown_cell].node_type = parent;

  if(parent == cancer) {
      if(rndgen.uniform() < parameters.freq_resistant) {
          world[position_of_grown_cell].node_type = resistant;
      }
  }

  update_death_prob(position_of_grown_cell, previous_type, parent);

  update_growth_prob(position_of_grown_cell);
  for(auto i : world[position_of_grown_cell].neighbors) {
    update_growth_prob(i->pos);
  }
}

void simulation::ask_infect_neighbours(int depth, float p, size_t pos) {
    if(p < 1e-6f) return;
    if(depth > 1) {
        depth--;
        for(auto& n : world[pos].neighbors) {
            ask_infect_neighbours(depth, p, n->pos);
        }

    } else {
        for(auto& n : world[pos].neighbors) {
            if(n->node_type == cancer) {
                if(rndgen.uniform() < p) {
                    n->node_type = infected;
                    update_death_prob(n->pos);
                    update_growth_prob(n->pos);
                    for(auto i : n->neighbors) {
                      update_growth_prob(i->pos);
                    }
                }
            }
        }
    }
}

void simulation::infect_long_distance(size_t pos) {
   for(size_t i = 1; i < long_distance_infection_probability.size(); ++i) {
       ask_infect_neighbours(i, long_distance_infection_probability[i], pos);
   }
}

void simulation::update_growth_prob(size_t pos) {
  std::array<float, 4> probs = world[pos].calc_prob_of_growth();
  for (size_t i = 0; i < growth_prob.size(); ++i) {
      growth_prob[i].update_entry(pos, probs[i]);
  }
}

void simulation::update_death_prob(size_t pos) { 
    for (size_t i = 0; i < 4; ++i) {
        float new_val = 0.f;
        if (i == world[pos].node_type) new_val = 1.f;

        death_prob[i].update_entry(pos, new_val);
    }
}

// specific death update
void simulation::update_death_prob(size_t pos,
                                   cell_type old_type,
                                   cell_type new_type) {
    if(old_type < death_prob.size()) death_prob[old_type].update_entry(pos, 0.f);
    if(new_type < death_prob.size()) death_prob[new_type].update_entry(pos, 1.f);
}


void simulation::update_rates() {
  rates[0] = parameters.birth_normal   * growth_prob[normal].get_total_sum();
  rates[1] = parameters.death_normal   * death_prob[normal].get_total_sum();

  rates[2] = parameters.birth_cancer   * growth_prob[cancer].get_total_sum();
  rates[3] = parameters.death_cancer   * death_prob[cancer].get_total_sum();

  rates[4] = parameters.birth_infected * growth_prob[infected].get_total_sum();
  rates[5] = parameters.death_infected * death_prob[infected].get_total_sum();

  rates[6] = parameters.birth_cancer_resistant * growth_prob[resistant].get_total_sum();
  rates[7] = parameters.death_cancer_resistant * death_prob[resistant].get_total_sum();
}

size_t simulation::pick_event(const std::array< float, 8>& rates, float sum) {
  float r = rndgen.uniform() * sum;
  for(size_t i = 0; i < rates.size(); ++i) {
    r -= rates[i];
    if(r <= 0) {
      return i;
    }
  }
  return 0;
}

std::array<int, 5> simulation::count_cell_types() {
  std::array<int, 5> total_num_cell_types = {0, 0, 0, 0, 0};
  for(auto i : world) {
      total_num_cell_types[i.node_type]++;
  }
  num_cell_types = total_num_cell_types;
  return total_num_cell_types;
}
