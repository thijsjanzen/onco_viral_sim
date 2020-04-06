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
            exit(1);
        }
    }

    size_t event = pick_event(rates, lambda);
    do_event(event);

    if(parameters.start_setup == grow) {
        if(t < parameters.time_adding_cancer &&
           t+dt >= parameters.time_adding_cancer) {
          add_cells(cancer);
        }

        if(t < parameters.time_adding_virus &&
           t+dt >= parameters.time_adding_virus) {
          add_infected();
        }
    }

    // check if the distributions are not getting close to zero,
    // in that case, some numerical irregularities might pop up
    // check only every hour:
    int delta_t = static_cast<int>(t + dt) - static_cast<int>(t);
    if (delta_t > 0) {
      for(size_t i = 0; i < 4; ++i) {
        if (num_cell_types[i] < 100 &&
            growth_prob[i].get_total_sum() != 0.f) {
          growth_prob[i].update_all();
        }
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

void simulation::implement_death(const cell_type& parent) {
  size_t position_of_dying_cell = death_prob[parent].draw_explicit(rndgen);

  change_cell_type(position_of_dying_cell, empty);

  if(parent == infected && parameters.prob_infection_upon_death > 0.f) {
      infect_long_distance(position_of_dying_cell);
  }
}

void simulation::implement_growth(const cell_type& parent) {

   // drawing position of growth is slow/bottleneck:
  size_t position_of_grown_cell = growth_prob[parent].draw_explicit(rndgen);

  cell_type new_type = parent;

  if(parent == cancer) {
      if(rndgen.uniform() < parameters.freq_resistant) {
          new_type = resistant;
      }
  }

  change_cell_type(position_of_grown_cell, new_type);
}

void simulation::ask_infect_neighbours(size_t depth, float p, size_t pos) {
    if(p < 1e-6f) return;
    if(std::isnan(p)) return;

    if(depth > 1) {
        depth--;
        for(const auto& n : world[pos].neighbors) {
            ask_infect_neighbours(depth, p, n->pos);
        }

    } else {
        for(auto& n : world[pos].neighbors) {
            if(n->get_cell_type() == cancer) {
                if(rndgen.uniform() < p) {
                    change_cell_type(n->pos, infected);
                }
            }
            if(n->get_cell_type() == normal) {
                 if(rndgen.uniform() < p) {
                     if(rndgen.uniform() < n->prob_normal_infected) {
                         change_cell_type(n->pos, infected);
                     }
                 }
            }
        }
    }
}

void simulation::infect_long_distance(size_t pos) {
   for(size_t i = 1; i < long_distance_infection_probability.size(); ++i) {
       ask_infect_neighbours(i,
                             static_cast<float>(long_distance_infection_probability[i]),
                             pos);
   }
}


void simulation::update_growth_prob(size_t pos) {
  std::array<float, 4> probs = world[pos].calc_prob_of_growth();
  for (size_t i = 0; i < 4; ++i) {
      growth_prob[i].update_entry(pos, probs[i]);
  }
}

void simulation::update_death_prob(size_t pos) { 
    for (size_t i = 0; i < 4; ++i) {
        float new_val = 0.f;
        if (i == world[pos].get_cell_type()) new_val = 1.f;

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
  assert(rates[1] == num_cell_types[normal] * parameters.death_normal);

  rates[2] = parameters.birth_cancer   * growth_prob[cancer].get_total_sum();
  rates[3] = parameters.death_cancer   * death_prob[cancer].get_total_sum();
  assert(rates[3] == num_cell_types[cancer] * parameters.death_cancer);


  rates[4] = parameters.birth_infected * growth_prob[infected].get_total_sum();
  rates[5] = parameters.death_infected * death_prob[infected].get_total_sum();
  assert(rates[5] == num_cell_types[infected] * parameters.death_infected);


  rates[6] = parameters.birth_cancer_resistant * growth_prob[resistant].get_total_sum();
  rates[7] = parameters.death_cancer_resistant * death_prob[resistant].get_total_sum();
  assert(rates[7] == num_cell_types[resistant] * parameters.death_cancer_resistant);

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

void simulation::change_cell_type(const size_t& pos,
                                  const cell_type& new_cell_type) {
  cell_type previous_type = world[pos].get_cell_type();
  world[pos].set_cell_type(new_cell_type);

  update_death_prob(pos, previous_type, new_cell_type);

  update_growth_prob(pos);
  for(const auto& i : world[pos].neighbors) {
    update_growth_prob(i->pos);
  }
  update_count(previous_type, new_cell_type);
}

void simulation::obtain_equilibrium() {

  std::vector< float > densities(10, 0);
  float prev_t = t;
  std::array<size_t, 5> cell_counts = num_cell_types;
  int count = 0;
  while(t < parameters.maximum_time) {
      update_one_step();

      if(static_cast<int>(t) - static_cast<int>(prev_t) == 10) {
           cell_counts = num_cell_types;
           auto density_normal = 1.f * cell_counts[normal] / (sq_size * sq_size);
           densities[count % 10] = cell_counts[normal];
           count++;
           if(count / 10 > 1) {
               float sum_first_half = 0.f;
               float sum_second_half = 0.f;
               for(size_t i = 0; i < 5; ++i) {
                   sum_first_half  += densities[i ];
                   sum_second_half += densities[i + 5];
               }
               std::cout << t << "\t" << sum_first_half * 0.2f
                         << "\t"      << sum_second_half * 0.2f << "\n";

               std::ofstream logfile("equilibrium_log.txt", std::ios::app);
               logfile << t << "\t" << sum_first_half * 0.2f
                             << "\t"      << sum_second_half * 0.2f << "\n";
               logfile.close();

               if(sum_first_half >= sum_second_half && density_normal > 0.5) {
                   break;
               }
           }
           prev_t = t;
      }
  }
  return;
}

void simulation::update_count(cell_type old_type, cell_type new_type) {
  num_cell_types[old_type]--;
  num_cell_types[new_type]++;
}

std::array<size_t, 5> simulation::count_cell_types() const {
  std::array<size_t, 5> total_num_cell_types = {0, 0, 0, 0, 0};
  for(const auto& i : world) {
      total_num_cell_types[ i.get_cell_type()]++;
  }
  return total_num_cell_types;
}

std::array<size_t, 5> simulation::get_count_cell_types() const {
  return num_cell_types;
}



void simulation::set_percent_infected(float percent_infected) {
  parameters.percent_infected = percent_infected;
}
void simulation::set_infection_type(infection_routine infect_routine) {
  parameters.infection_type = infect_routine;
}



/*
void simulation::check_cell_type_counts() {
  std::array<size_t, 5> manual_count = count_cell_types();
  std::array<size_t, 5> book_keep_count = num_cell_types;

  for(int i = 0; i < 5; ++i) {
    assert(manual_count[i] == book_keep_count[i]);
  }
  return;
}
*/


