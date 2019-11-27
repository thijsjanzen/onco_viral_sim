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
    update_rates(rates);
    float lambda = std::accumulate(rates.begin(), rates.end(), 0.0f);
    float dt = rndgen.Expon(lambda);
    size_t event = pick_event(rates, lambda);
    do_event(event);
   // verify_bounds();

    if(t < parameters.time_adding_cancer && t+dt >= parameters.time_adding_cancer) {
      add_cells(cancer);
    }

    if(t < parameters.time_adding_virus && t+dt >= parameters.time_adding_virus) {
      add_infected();
    }



    t += dt;
}



void simulation::run() {
  // do stuff

  t = 0.f;
  while(t < parameters.maximum_time) {
     update_one_step();
  }
}

void simulation::initialize_network() {
  // we first make a regular network
  // we will do voronoi later...

 // for(size_t i = 0; i < world.size(); ++i) {
  for(auto& i : world) {
    i.update_neighbors(world, sq_size);
  }

//  verify_bounds();

  add_cells(normal);

  for(auto& i : world) {
    i.update_neighbor_types();
  }
 // verify_bounds();

  count_cell_types();
  update_growth_probabilities();

  for(size_t i = 0; i < 3; ++i) {
    auto local_birth_dist = rndutils::mutable_discrete_distribution<int, rndutils::all_zero_policy_nothing>{};
    local_birth_dist.mutate(growth_probs[i].begin(), growth_probs[i].end());
    growth_prob_rnd.push_back(local_birth_dist);

    auto local_death_dist = rndutils::mutable_discrete_distribution<int, rndutils::all_zero_policy_nothing>{};
    local_death_dist.mutate(death_probs[i].begin(), death_probs[i].end());
    death_prob_rnd.push_back(local_death_dist);
  }
 // verify_bounds();
}

void simulation::count_cell_types() {
  num_cell_types = {0, 0, 0};
  for(auto i : world) {
      num_cell_types[i.node_type]++;
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

    default: {
      // do nothing
      break;
    }
  }
}

void simulation::update_rates(std::array< float, 6>& rates) {
  rates[0] = parameters.birth_normal   * static_cast<float>(growth_prob_rnd[normal].cdf().back()); //sum_growth_prob[normal]; //std::accumulate(growth_probs[normal].begin(), growth_probs[normal].end(), 0.0);
  rates[1] = parameters.death_normal   * static_cast<float>(death_prob_rnd[normal].cdf().back()); //sum_death_prob[normal];

  rates[2] = parameters.birth_cancer   * static_cast<float>(growth_prob_rnd[cancer].cdf().back()); //sum_growth_prob[cancer]; //std::accumulate(growth_probs[cancer].begin(), growth_probs[cancer].end(), 0.0);
  rates[3] = parameters.death_cancer   * static_cast<float>(death_prob_rnd[cancer].cdf().back()); //sum_death_prob[cancer];

  rates[4] = parameters.birth_infected * static_cast<float>(growth_prob_rnd[infected].cdf().back());  //std::accumulate(growth_probs[infected].begin(), growth_probs[infected].end(), 0.0);
  rates[5] = parameters.death_infected * static_cast<float>(death_prob_rnd[infected].cdf().back());
}

size_t simulation::pick_event(const std::array< float, 6>& rates, float sum) {
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
      position_of_dying_cell = static_cast<size_t>(death_prob_rnd[normal](rndgen.rndgen_));        //rndgen.draw_from_dist(death_probs[normal].begin(), death_probs[normal].end(), 1.f);
      break;
    case cancer:
      position_of_dying_cell = static_cast<size_t>(death_prob_rnd[cancer](rndgen.rndgen_)); //rndgen.draw_from_dist(death_probs[cancer].begin(), death_probs[cancer].end(), 1.f);
      break;
    case infected:
      position_of_dying_cell = static_cast<size_t>(death_prob_rnd[infected](rndgen.rndgen_)); //rndgen.draw_from_dist(death_probs[infected].begin(), death_probs[infected].end(), 1.f);
      break;
    case empty:
      //  std::cout << "ERROR! empty node is dying\n";
      position_of_dying_cell = 0;
      break;
  }

  world[position_of_dying_cell].node_type = empty;
  death_probs[parent][position_of_dying_cell] = 0.f;
  update_death_prob_vectors(parent, static_cast<long>(position_of_dying_cell));
  //sum_death_prob[parent]--; // update death probabilities

  size_t min_pos = world.size();
  update_growth_prob(position_of_dying_cell);
  for(auto i : world[position_of_dying_cell].neighbors) {
    update_growth_prob(i->pos);
    if(i->pos < min_pos) min_pos = i->pos;
  }
  update_growth_prob_vectors(static_cast<long>(min_pos));
}



void simulation::implement_growth(const cell_type& parent) {
  size_t position_of_grown_cell = 0;
  switch(parent) {
    case normal:
      position_of_grown_cell = static_cast<size_t>(growth_prob_rnd[normal](rndgen.rndgen_));
      break;
    case cancer:
      position_of_grown_cell = static_cast<size_t>(growth_prob_rnd[cancer](rndgen.rndgen_));
      break;
    case infected:
      position_of_grown_cell = static_cast<size_t>(growth_prob_rnd[infected](rndgen.rndgen_));
      break;
    case empty:
      // std::cout << "ERROR! empty node is growing\n";
      position_of_grown_cell = 0;
      break;
  }

  size_t min_pos = position_of_grown_cell;

  world[position_of_grown_cell].node_type = parent;

  // update growth probability of cell
  update_growth_prob(position_of_grown_cell);

  // update growth probabilities of the neighbors:
  for(auto i : world[position_of_grown_cell].neighbors) {
    update_growth_prob(i->pos);
    if(i->pos < min_pos) min_pos = i->pos;
  }

  update_growth_prob_vectors(static_cast<long>(min_pos));

  // update death probability of new cell:
  death_probs[parent][position_of_grown_cell] = 1.f;
  update_death_prob_vectors(parent, position_of_grown_cell);
}

void simulation::update_death_prob_vectors(const cell_type& parent, long pos) {
    death_prob_rnd[parent].mutate(death_probs[parent].begin(),
                                  death_probs[parent].end());

    //verify_death_probs();
}

void simulation::update_growth_prob_vectors(long pos) {

    for(size_t i = 0; i < 3; ++i) {
        growth_prob_rnd[i].mutate(growth_probs[i].begin(),
                                  growth_probs[i].end());

       // verify_growth_probs();

    }
}

/*void check_zeros(const std::vector<float>& p,
                 const std::vector<double>& cdf) {

    assert(p.size() == cdf.size());
    for(size_t i = 1; i < p.size(); ++i) {
        if(p[i] == 0) {
            //assert(cdf[i] == cdf[i-1]);
            if(cdf[i] != cdf[i-1]) {
                int a = 5;
            }
        }
    }
    return;
}*/

/*void simulation::verify_growth_probs() {
     for(size_t i = 0; i < 3; ++i) {
         check_zeros(growth_probs[i], growth_prob_rnd[i].cdf());
     }
}

void simulation::verify_death_probs() {
     for(size_t i = 0; i < 3; ++i) {
         check_zeros(death_probs[i], death_prob_rnd[i].cdf());
     }
}*/


void simulation::add_cells(cell_type focal_cell_type) {
    // pick center node:
  size_t x = sq_size / 2;
  size_t y = sq_size / 2;
  size_t focal_pos = x * sq_size + y;
  std::vector<size_t> cells_turned(1, focal_pos);
  world[focal_pos].node_type = focal_cell_type;
  auto max_number_of_cells = static_cast<size_t>(parameters.initial_number_cancer_cells);
  if(focal_cell_type == normal) max_number_of_cells = static_cast<size_t>(parameters.num_cells * 0.1);
  size_t counter = 0;
  while(cells_turned.size() < max_number_of_cells) {
    focal_pos = cells_turned[counter];
    for(size_t i = 0; i < world[focal_pos].neighbors.size(); ++i) {
      size_t other_pos = world[focal_pos].neighbors[i]->pos;
      if(world[other_pos].node_type != focal_cell_type) {
         world[other_pos].node_type = focal_cell_type;
        cells_turned.push_back(other_pos);
        if(cells_turned.size() >= max_number_of_cells) break;
      }
    }
    counter++;
  }

  for(auto& i : cells_turned) {
    update_growth_prob(i);
    for(auto& j : world[i].neighbors) {
      update_growth_prob(j->pos);
    }
  }
}

void simulation::update_growth_prob(size_t pos) {
  std::array<float, 3> probs = world[pos].calc_prob_of_growth();

  for(size_t i = 0; i < 3; ++i) {
      growth_probs[i][pos] = probs[i];
  }
}

simulation::simulation(const Param& param) {
  parameters = param;

  sq_size = static_cast<size_t>(sqrt(parameters.num_cells));

  world.resize(parameters.num_cells);

  for(size_t i = 0; i < world.size(); ++i) {
    world[i].pos = i;
    world[i].set_coordinates(sq_size);
  }

  growth_probs.resize(3, std::vector< float >(parameters.num_cells, 0.f));
  death_probs.resize(3, std::vector< float >(parameters.num_cells, 0.f));

  //max_growth_prob = {0.f, 0.f, 0.f};
}

void simulation::print_to_file(float t) {
  std::string file_name = "world_file_" + std::to_string(static_cast<int>(t)) + ".txt";
  std::ofstream outfile(file_name.c_str());
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

// initialization routines:
void simulation::update_growth_probabilities() {

  for(auto i : world) {
    std::array<float, 3> probs = i.calc_prob_of_growth();

    size_t pos = i.pos;

    growth_probs[normal][pos] = probs[normal];
    growth_probs[cancer][pos] = probs[cancer];
    growth_probs[infected][pos] = probs[infected];

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
}

std::vector<int> simulation::get_cell_numbers() {
    std::vector<int> output(3,0);
    output[normal] =    static_cast<int>(std::accumulate(death_probs[normal].begin(), death_probs[normal].end(), 0));
    output[cancer] =    static_cast<int>(std::accumulate(death_probs[cancer].begin(), death_probs[cancer].end(), 0));
    output[infected] =  static_cast<int>(std::accumulate(death_probs[infected].begin(), death_probs[infected].end(), 0));
    return output;
}


void simulation::infect_random() {
    // count number of cancer cells
    int num_cancer_cells = static_cast<int>(std::accumulate(death_probs[cancer].begin(), death_probs[cancer].end(), 0));

    int infected_cells = 0;
    int to_be_infected = static_cast<int>(parameters.percent_infected * num_cancer_cells);
    while(infected_cells < to_be_infected && num_cancer_cells > 0) {
        size_t position_of_grown_cell = static_cast<size_t>(death_prob_rnd[cancer](rndgen.rndgen_));

           world[position_of_grown_cell].node_type = infected;

        // update growth probability of cell
        update_growth_prob(position_of_grown_cell);

        // update growth probabilities of the neighbors:
        for(auto i : world[position_of_grown_cell].neighbors) {
          update_growth_prob(i->pos);
        }

        // update death probability of new cell:
        death_probs[infected][position_of_grown_cell] = 1.f;
        death_probs[cancer][position_of_grown_cell] = 0.f;
        update_death_prob_vectors(cancer, 0);
        num_cancer_cells--; // easy count for now.
        infected_cells++;
    }

    // superfluous check
    std::array<int, 3> num_cells = {0, 0, 0};
    for(auto i : world) {
        num_cells[i.node_type]++;
    }
    assert(num_cells[2] == infected_cells);


    update_growth_prob_vectors(0); // update all vectors from the start.
    update_death_prob_vectors(normal, 0);
    update_death_prob_vectors(cancer, 0);
    update_death_prob_vectors(infected, 0);

}

void simulation::infect_center() {
    int num_cancer_cells = static_cast<int>(std::accumulate(death_probs[cancer].begin(), death_probs[cancer].end(), 0));

    int to_be_infected = static_cast<int>(parameters.percent_infected * num_cancer_cells);


    // TODO: calculate the center of the cancer!

    size_t x = sq_size / 2;
    size_t y = sq_size / 2;
    size_t focal_pos = x * sq_size + y;
    std::vector<size_t> cells_turned(1, focal_pos);
    world[focal_pos].node_type = infected;
    int counter = 0;
    while(cells_turned.size() < to_be_infected && num_cancer_cells > 0) {
      focal_pos = cells_turned[counter];
      for(size_t i = 0; i < world[focal_pos].neighbors.size(); ++i) {
        size_t other_pos = world[focal_pos].neighbors[i]->pos;
        if(world[other_pos].node_type == cancer) {
           world[other_pos].node_type = infected;
           cells_turned.push_back(other_pos);
           num_cancer_cells--;

           if(cells_turned.size() >= to_be_infected) break;
        }
      }
      counter++;
    }

    for(auto& i : cells_turned) {
      update_growth_prob(i);
      for(auto& j : world[i].neighbors) {
        update_growth_prob(j->pos);
      }
    }
    update_growth_prob_vectors(0); // update all vectors from the start.
    update_death_prob_vectors(normal, 0);
    update_death_prob_vectors(infected, 0);
    update_death_prob_vectors(cancer, 0);
}


void simulation::add_infected() {
    switch(parameters.infection_type) {
        case random_infection:
            infect_random();
            break;
        case center_infection:
            infect_center();
            break;
        case multinode:
            // do something
            break;
        case perimeter:
            // do something
            break;
    }

    int num_infected = static_cast<int>(std::accumulate(death_probs[infected].begin(), death_probs[infected].end(), 0));
    assert(num_infected > 0);
}




