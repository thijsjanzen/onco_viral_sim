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
#include <fstream>
#include "rndutils.hpp"
#include <iostream>

void simulation::update_one_step() {
    update_rates();
    float lambda = std::accumulate(rates.begin(), rates.end(), 0.0f);
    float dt = rndgen.Expon(lambda);
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



void simulation::run() {
  t = 0.f;
  while(t < parameters.maximum_time) {
     update_one_step();
  }
}

void simulation::initialize_network() {
  // we first make a regular network
  // we will do voronoi later...

  for(auto& i : world) {
    i.update_neighbors(world, sq_size);
    i.set_node_type(empty);
  }

  if(parameters.start_setup == grow) {

    add_cells(normal);

    for(auto& i : world) {
<<<<<<< Updated upstream
        update_growth_prob(i.pos);
=======
        i.set_all_neighbor_freqs();
        i.calc_prob_of_growth();
        set_growth_prob(i.pos);
>>>>>>> Stashed changes
        update_death_prob(i.pos);
    }

    for(size_t i = 0; i < 4; ++i) {
        growth_prob_rnd[i] = binned_distribution<sq_size>(growth_probs[i].begin(), growth_probs[i].end());

        death_prob_rnd[i] = binned_distribution<sq_size>(death_probs[i].begin(), death_probs[i].end());
    }
  }

  if(parameters.start_setup == full) {
    for(size_t i = 0; i < 4; ++i) {
         growth_prob_rnd[i] = binned_distribution<sq_size>();
         death_prob_rnd[i] = binned_distribution<sq_size>();
    }
    initialize_full();
  }
}

void simulation::count_cell_types() {
  num_cell_types = {0, 0, 0};
  for(auto i : world) {
      num_cell_types[i.get_node_type()]++;
  }
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

void simulation::update_rates() {
  rates[0] = parameters.birth_normal   * growth_prob_rnd[normal].get_total_sum();
  rates[1] = parameters.death_normal   * death_prob_rnd[normal].get_total_sum();

  rates[2] = parameters.birth_cancer   * growth_prob_rnd[cancer].get_total_sum();
  rates[3] = parameters.death_cancer   * death_prob_rnd[cancer].get_total_sum();

  rates[4] = parameters.birth_infected * growth_prob_rnd[infected].get_total_sum();
  rates[5] = parameters.death_infected * death_prob_rnd[infected].get_total_sum();

  rates[6] = parameters.birth_cancer_resistant * growth_prob_rnd[resistant].get_total_sum();
  rates[7] = parameters.death_cancer_resistant * death_prob_rnd[resistant].get_total_sum();
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

void simulation::implement_death(const cell_type& parent) {
  size_t position_of_dying_cell = 0;
  switch(parent) {
    case normal:
      position_of_dying_cell = static_cast<size_t>(death_prob_rnd[normal].draw_explicit(death_probs[normal].begin(),
                                                                                        death_probs[normal].end(),
                                                                                        rndgen)); //static_cast<size_t>(death_prob_rnd[normal](rndgen.rndgen_));
      break;
    case cancer:
      position_of_dying_cell = static_cast<size_t>(death_prob_rnd[cancer].draw_explicit(death_probs[cancer].begin(),
                                                                                        death_probs[cancer].end(),
                                                                                        rndgen));//static_cast<size_t>(death_prob_rnd[cancer](rndgen.rndgen_));
      break;
    case infected:
      position_of_dying_cell = static_cast<size_t>(death_prob_rnd[infected].draw_explicit(death_probs[infected].begin(),
                                                                                          death_probs[infected].end(),
                                                                                            rndgen));
      break;
    case resistant:
      position_of_dying_cell = static_cast<size_t>(death_prob_rnd[resistant].draw_explicit(death_probs[resistant].begin(),
                                                                                          death_probs[resistant].end(),
                                                                                            rndgen));
      break;
    case empty:
      position_of_dying_cell = 0;
      break;
  }

  cell_type previous_type = world[position_of_dying_cell].get_node_type();
  world[position_of_dying_cell].set_node_type(empty);
  update_death_prob(position_of_dying_cell, previous_type);

  update_growth_prob(position_of_dying_cell);

  for(auto i : world[position_of_dying_cell].neighbors) {
    update_growth_prob(i->pos);
  }
  
  if(parent == infected) {
      infect_long_distance(position_of_dying_cell);
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
<<<<<<< Updated upstream
            if(n->node_type != infected) {
=======
            if(n->get_node_type() == cancer) {
>>>>>>> Stashed changes
                if(rndgen.uniform() < p) {
                    n->set_node_type(infected);
                    update_death_prob(n->pos, cancer);
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


void simulation::implement_growth(const cell_type& parent) {

   // drawing position of growth is slow/bottleneck:
  size_t position_of_grown_cell = 0;
  switch(parent) {
    case normal:
      position_of_grown_cell = static_cast<size_t>(growth_prob_rnd[normal].draw_explicit(growth_probs[normal].begin(),
                                                                     growth_probs[normal].end(),
                                                                     rndgen));
      break;
    case cancer:
      position_of_grown_cell = static_cast<size_t>(growth_prob_rnd[cancer].draw_explicit(growth_probs[cancer].begin(),
                                                                     growth_probs[cancer].end(),
                                                                     rndgen));
      break;
    case infected:
      position_of_grown_cell = static_cast<size_t>(growth_prob_rnd[infected].draw_explicit(growth_probs[infected].begin(),
                                                                        growth_probs[infected].end(),
                                                                        rndgen));
      break;
    case resistant:
      position_of_grown_cell = static_cast<size_t>(growth_prob_rnd[resistant].draw_explicit(growth_probs[resistant].begin(),
                                                                      growth_probs[resistant].end(),
                                                                      rndgen));
        break;
    case empty:
      position_of_grown_cell = 0;
      break;
  }

  cell_type old_type = world[position_of_grown_cell].get_node_type();
  world[position_of_grown_cell].set_node_type(parent);
  if(parent == cancer) {
      if(rndgen.uniform() < parameters.freq_resistant) {
          world[position_of_grown_cell].set_node_type(resistant);
      }
  }

  // update growth probability of cell
  update_growth_prob(position_of_grown_cell);
  update_death_prob(position_of_grown_cell, old_type);

  // update growth probabilities of the neighbors:
  for(auto i : world[position_of_grown_cell].neighbors) {
    update_growth_prob(i->pos);
  }
}

void simulation::update_death_cdf(const cell_type& parent, size_t pos) {
    death_prob_rnd[parent].mutate(death_probs[parent].begin(),
                                  death_probs[parent].end(),
                                  pos);
}

void simulation::update_death_cdf_all() {
    for(size_t i = 0; i < death_prob_rnd.size(); ++i) {
        death_prob_rnd[i].mutate_all(death_probs[i].begin(),
                                     death_probs[i].end());
    }
}

void simulation::update_growth_cdf(size_t pos) {
    for(size_t i = 0; i < growth_prob_rnd.size(); ++i) {
        growth_prob_rnd[i].mutate(growth_probs[i].begin(),
                                  growth_probs[i].end(),
                                  pos);
    }
}



void simulation::add_cells(const cell_type& focal_cell_type) {
    // pick center node:
  size_t x = sq_size / 2;
  size_t y = sq_size / 2;
  size_t focal_pos = x * sq_size + y;
  std::vector<size_t> cells_turned(1, focal_pos);
  std::vector<cell_type> cell_types_turned(1, world[focal_pos].get_node_type());

  auto max_number_of_cells = static_cast<size_t>(parameters.initial_number_cancer_cells);
  if (focal_cell_type == normal) max_number_of_cells = parameters.initial_number_normal_cells;
  if (max_number_of_cells > world.size()) max_number_of_cells = world.size();
  size_t counter = 0;
  while (cells_turned.size() < max_number_of_cells) {
    focal_pos = cells_turned[counter];
    for (size_t i = 0; i < world[focal_pos].neighbors.size(); ++i) {
      size_t other_pos = world[focal_pos].neighbors[i]->pos;
      if (world[other_pos].get_node_type() != focal_cell_type) {
        cell_types_turned.push_back(world[other_pos].get_node_type());
        //change_cell_type(other_pos, focal_cell_type);
        world[other_pos].set_node_type(focal_cell_type);
        cells_turned.push_back(other_pos);

        if (cells_turned.size() >= max_number_of_cells) break;
      }
    }
    counter++;
  }

  for (size_t i = 0; i < cells_turned.size(); ++i) {
    world[ cells_turned[i] ].set_all_neighbor_freqs();
    update_growth_prob(cells_turned[i]);
    update_death_prob(cells_turned[i], cell_types_turned[i]);
    for (auto& j : world[cells_turned[i]].neighbors) {
      j->set_all_neighbor_freqs();
      update_growth_prob(j->pos);
    }
  }
}

void simulation::update_growth_prob(size_t pos) {
<<<<<<< Updated upstream
  world[pos].update_prob_of_growth();
  for (size_t i = 0; i < growth_prob_rnd.size(); ++i) {
      growth_prob_rnd[i].update_entry(growth_probs[i].begin(), growth_probs[i].end(),
                                      pos, world[pos].prob_of_growth[i]);
      growth_probs[i][pos] = world[pos].prob_of_growth[i];
=======
  std::vector<float> old_probs = world[pos].prob_of_growth_;
  world[pos].calc_prob_of_growth();
  for (size_t i = 0; i < growth_prob_rnd.size(); ++i) {
      if(old_probs[i] != world[pos].prob_of_growth_[i]) {
          growth_prob_rnd[i].update_entry(growth_probs[i].begin(), growth_probs[i].end(),
                                          pos, world[pos].prob_of_growth_[i]);
          growth_probs[i][pos] = world[pos].prob_of_growth_[i];
      }
  }
}

void simulation::set_growth_prob(size_t pos) {
  world[pos].calc_prob_of_growth();
  for (size_t i = 0; i < growth_prob_rnd.size(); ++i) {
      growth_prob_rnd[i].update_entry(growth_probs[i].begin(), growth_probs[i].end(),
                                      pos, world[pos].prob_of_growth_[i]);
      growth_probs[i][pos] = world[pos].prob_of_growth_[i];
>>>>>>> Stashed changes
  }
}




void simulation::update_death_prob(size_t pos, cell_type old_type) {

    if(old_type != empty) {
        death_prob_rnd[old_type].update_entry(death_probs[old_type].begin(),
                                               death_probs[old_type].end(),
                                               pos,
                                               0.f);

        death_probs[old_type][pos] = 0.f;
    }

    cell_type new_type = world[pos].get_node_type();
    if(new_type != empty) {
        death_prob_rnd[new_type].update_entry(death_probs[new_type].begin(),
                                           death_probs[new_type].end(),
                                           pos,
                                           1.f);

        death_probs[new_type][pos] = 1.f;
    }
}

void simulation::update_death_prob(size_t pos) {
    // only update if there was a change in node type:
    if(world[pos].node_type != empty &&
            death_probs[ world[pos].node_type ][pos] != 0.f) return;

    // loop over all and change where needed:
    for (size_t i = 0; i < 4; ++i) {
        float new_val = 0.f;
        if (i == world[pos].get_node_type()) new_val = 1.f;

        death_prob_rnd[i].update_entry(death_probs[i].begin(),
                                       death_probs[i].end(),
                                       pos,
                                       new_val);
        death_probs[i][pos] = new_val;
    }
}


<<<<<<< Updated upstream
=======
simulation::simulation() {
    parameters = Param(); // use default parameters
    world.resize(num_cells);

    for (size_t i = 0; i < world.size(); ++i) {
      world[i].pos = i;
      world[i].set_coordinates(sq_size);
      world[i].prob_normal_infected = parameters.prob_normal_infection;
    }

    growth_probs.resize(4, std::vector< float >(num_cells, 0.f));
    death_probs.resize( 4, std::vector< float >(num_cells, 0.f));

    long_distance_infection_probability = std::vector<double>(1, 0);
    float lambda = 1.0f / parameters.distance_infection_upon_death;
    for(size_t d = 1; d < sq_size; ++d) {
        double local_prob = parameters.prob_infection_upon_death * lambda * expf(-lambda * d);
        long_distance_infection_probability.push_back(local_prob);
        if(local_prob < 1e-3) break;
    }
}

>>>>>>> Stashed changes
simulation::simulation(const Param& param) {
  parameters = param;

  world.resize(num_cells);

  for (size_t i = 0; i < world.size(); ++i) {
    world[i].pos = i;
    world[i].set_coordinates(sq_size);
    world[i].prob_normal_infected = parameters.prob_normal_infection;
  }

  growth_probs.resize(4, std::vector< float >(num_cells, 0.f));
  death_probs.resize( 4, std::vector< float >(num_cells, 0.f));

  long_distance_infection_probability = std::vector<double>(1, 0);
  float lambda = 1.0f / parameters.distance_infection_upon_death;
  for(size_t d = 1; d < sq_size; ++d) {
      double local_prob = parameters.prob_infection_upon_death * lambda * expf(-lambda * d);
      long_distance_infection_probability.push_back(local_prob);
      if(local_prob < 1e-3) break;
  }

}

// initialization routines:
void simulation::update_growth_probabilities() {

<<<<<<< Updated upstream
  for(auto i : world) {
    i.update_prob_of_growth();

    size_t pos = i.pos;

    growth_probs[normal][pos] = i.prob_of_growth[normal];
    growth_probs[cancer][pos] = i.prob_of_growth[cancer];
    growth_probs[infected][pos] = i.prob_of_growth[infected];
    growth_probs[resistant][pos] = i.prob_of_growth[resistant];
=======
  for(auto& i : world) {
    i.set_all_neighbor_freqs();
    i.calc_prob_of_growth();

    size_t pos = i.pos;

    growth_probs[normal][pos]       = i.prob_of_growth_[normal];
    growth_probs[cancer][pos]       = i.prob_of_growth_[cancer];
    growth_probs[infected][pos]     = i.prob_of_growth_[infected];
    growth_probs[resistant][pos]    = i.prob_of_growth_[resistant];
>>>>>>> Stashed changes

    switch(i.get_node_type()) {
      case normal:
        death_probs[normal][pos] = 1.f;
        break;
      case cancer:
        death_probs[cancer][pos] = 1.f;
        break;
      case infected:
        death_probs[infected][pos] = 1.f;
        break;
      case resistant:
        death_probs[resistant][pos] = 1.f;
        break;
      case empty:
        break;
    }
  }
}

std::vector<int> simulation::get_cell_numbers() {
    std::vector<int> output(4,0);
    output[normal]      =    static_cast<int>(std::accumulate(death_probs[normal].begin(),    death_probs[normal].end(), 0));
    output[cancer]      =    static_cast<int>(std::accumulate(death_probs[cancer].begin(),    death_probs[cancer].end(), 0));
    output[infected]    =    static_cast<int>(std::accumulate(death_probs[infected].begin(),  death_probs[infected].end(), 0));
    output[resistant]   =    static_cast<int>(std::accumulate(death_probs[resistant].begin(), death_probs[resistant].end(), 0));
    return output;
}

void simulation::infect_random() {
    // count number of cancer cells
    int num_cancer_cells = static_cast<int>(std::accumulate(death_probs[cancer].begin(), death_probs[cancer].end(), 0));

    int infected_cells = 0;
    int to_be_infected = static_cast<int>(parameters.percent_infected * num_cancer_cells);
    if(to_be_infected == 0) return;
    while(infected_cells < to_be_infected && num_cancer_cells > 0) {
        size_t position_of_grown_cell = static_cast<size_t>(death_prob_rnd[cancer](rndgen.rndgen_));

<<<<<<< Updated upstream
        world[position_of_grown_cell].node_type = infected;
=======
    std::vector< size_t > cancer_pos(num_cancer_cells);
    int j = 0;
    for(auto i : world) {
        if(i.get_node_type() == cancer) {
            cancer_pos[j] = static_cast<size_t>(i.pos);
            j++;
        }
    }
>>>>>>> Stashed changes

        // update growth probability of cell
        update_growth_prob(position_of_grown_cell);

<<<<<<< Updated upstream
        // update growth probabilities of the neighbors:
        for(auto i : world[position_of_grown_cell].neighbors) {
          update_growth_prob(i->pos);
        }
=======
        world[position_of_infected_cell].set_node_type(infected);
>>>>>>> Stashed changes

        // update death probability of new cell:
        update_death_prob(position_of_grown_cell);
        num_cancer_cells--; // easy count for now.
        infected_cells++;
    }

<<<<<<< Updated upstream
    update_growth_cdf(0); // update all vectors from the start.
    update_death_cdf_all();
=======
    // hard, full update
    for(auto i : world) {
        i.set_all_neighbor_freqs();
        update_growth_prob(i.pos);
        update_death_prob(i.pos);
    }
>>>>>>> Stashed changes
}

void simulation::infect_center() {
    int num_cancer_cells = static_cast<int>(std::accumulate(death_probs[cancer].begin(), death_probs[cancer].end(), 0));

    size_t to_be_infected = static_cast<size_t>(parameters.percent_infected * num_cancer_cells);
    if(to_be_infected == 0) return;

    float avg_x = 0;
    float avg_y = 0;
    for(auto i : world) {
        if(i.get_node_type() == cancer) {
            avg_x += i.x_;
            avg_y += i.y_;
        }
    }

    size_t start_x = static_cast<size_t>(avg_x / num_cancer_cells);
    size_t start_y = static_cast<size_t>(avg_y / num_cancer_cells);

    //now find starting cell to infect
    size_t starting_pos = start_x * sq_size + start_y;
    while(world[starting_pos].get_node_type() != cancer) {
        for(size_t i = 0; i < world[starting_pos].neighbors.size(); ++i) {
            if( world[starting_pos].neighbors[i]->get_node_type() == cancer) {
                starting_pos = world[starting_pos].neighbors[i]->pos;
                break;
            }
        }
        int rand_index = rndgen.random_number(world[starting_pos].neighbors.size());
        starting_pos = world[starting_pos].neighbors[rand_index]->pos;
    }

    size_t focal_pos = starting_pos;
    std::vector<size_t> cells_turned(1, focal_pos);
    world[focal_pos].set_node_type(infected);
    size_t counter = 0;
    while(cells_turned.size() < to_be_infected && num_cancer_cells > 0) {
      focal_pos = cells_turned[counter];
      for(size_t i = 0; i < world[focal_pos].neighbors.size(); ++i) {
        size_t other_pos = world[focal_pos].neighbors[i]->pos;
        if(world[other_pos].get_node_type() == cancer) {
           world[other_pos].set_node_type(infected);
           cells_turned.push_back(other_pos);
           num_cancer_cells--;

           if(cells_turned.size() >= to_be_infected) break;
        }
      }
      counter++;
      if(counter > cells_turned.size()) break;
    }

    for(auto i : world) {
        i.set_all_neighbor_freqs();
        update_growth_prob(i.pos);
        update_death_prob(i.pos);
    }
}

void simulation::infect_all_cancer() {
    for(auto& i : world) {
        if(i.get_node_type() == cancer) {
            i.set_node_type(infected);
        }
    }

    for(auto i : world) {
        update_growth_prob(i.pos);
        update_death_prob(i.pos);
    }
}


void simulation::add_infected() {

    if(parameters.percent_infected == 1.0f) {
        infect_all_cancer();
        return;
    }

    switch(parameters.infection_type) {
        case random_infection:
            infect_random();
            break;
        case center_infection:
            infect_center();
            break;
    }
}

void simulation::initialize_full() {

    for(auto& i : world) {
        i.initialize_node_type(normal);
        i.prob_normal_infected = parameters.prob_normal_infection;
    }
    // still have to update neighbor frequencies, because not all neighbors
    // were initialized in the previous loop
    for(auto& i : world) {
        i.set_all_neighbor_freqs();
    }

    parameters.initial_number_cancer_cells = static_cast<int>(0.1f * world.size());
    add_cells(cancer);

    // just for safety, do a full scan:
    for(auto& i : world) {
        i.set_all_neighbor_freqs();
        set_growth_prob(i.pos);
        update_death_prob(i.pos);
    }

    parameters.percent_infected = 0.1f;
    infect_center();

    // and update again, just for added safety:
    for(auto& i : world) {
        update_growth_prob(i.pos);
        update_death_prob(i.pos);
    }

    for(size_t i = 0; i < growth_probs.size(); ++i) {
        growth_prob_rnd[i].update_row_cdf(growth_probs[i].begin());
        death_prob_rnd[i].update_row_cdf(death_probs[i].begin());
    }
}
