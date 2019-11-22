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

#include "rndutils.hpp"


void simulation::update_one_step() {
    update_rates(rates);
    float lambda = std::accumulate(rates.begin(), rates.end(), 0.0f);
    float dt = rndgen.Expon(lambda);
    size_t event = pick_event(rates, lambda);
    do_event(event);

    if(t < parameters.time_adding_cancer && t+dt >= parameters.time_adding_cancer) {
      add_cells(cancer);
    }

   /* if((int)t < (int)(t+dt)) {
   //   auto end = std::chrono::system_clock::now();
   //   std::chrono::duration<double> elapsed_seconds = end - start;

   //   std::cout << (int)t << "\t";
      for(int i = 0; i < 3; ++i) {
        std::cout << sum_death_prob[i] << "\t"; // these are the numbers of cells of each type
      }


 //     std::cout << elapsed_seconds.count() << "\n";

  //    start = std::chrono::system_clock::now();
    }*/

    t += dt;

    if(num_cell_types[0] < 1.f) {
      std::cout << "all normal cells are dead\n";
      return;
    }

}



void simulation::run() {
  // do stuff

  t = 0.f;
  auto start = std::chrono::system_clock::now();
  while(t < parameters.maximum_time) {
    update_rates(rates);
    float lambda = std::accumulate(rates.begin(), rates.end(), 0.f);
    float dt = rndgen.Expon(lambda);
    size_t event = pick_event(rates, lambda);
    do_event(event);

    if(t < parameters.time_adding_cancer && t+dt >= parameters.time_adding_cancer) {
      add_cells(cancer);
    }

    if(static_cast<int>(t) < static_cast<int>(t+dt)) {
    //  print_to_file(t);
      auto end = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_seconds = end - start;

      std::cout << static_cast<int>(t) << "\t";
      for(size_t i = 0; i < 3; ++i) {
        std::cout << sum_death_prob[i] << "\t"; // these are the numbers of cells of each type
      }


      std::cout << elapsed_seconds.count() << "\n";

      start = std::chrono::system_clock::now();
    }

    t += dt;

    if(num_cell_types[0] < 1 ) {
      std::cout << "all normal cells are dead\n";
      return;
    }
  }
}

void simulation::initialize_network() {
  // we first make a regular network
  // we will do voronoi later...

  for(size_t i = 0; i < world.size(); ++i) {
    world[i].update_neighbors(world, sq_size);
  }

  add_cells(normal);

  for(auto it = world.begin(); it != world.end(); ++it) {
    (*it).update_neighbor_types();
  }

  count_cell_types();
  update_growth_probabilities();

  for(size_t i = 0; i < 3; ++i) {
    sum_growth_prob[i] = std::accumulate(growth_probs[i].begin(), growth_probs[i].end(), 0.f);
    sum_death_prob[i] = std::accumulate(death_probs[i].begin(), death_probs[i].end(), 0.f);
  }
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
  return;
}

void simulation::update_rates(std::array< float, 6>& rates) {
  rates[0] = parameters.birth_normal   * sum_growth_prob[normal]; //std::accumulate(growth_probs[normal].begin(), growth_probs[normal].end(), 0.0);
  rates[1] = parameters.death_normal   * sum_death_prob[normal];

  rates[2] = parameters.birth_cancer   * sum_growth_prob[cancer]; //std::accumulate(growth_probs[cancer].begin(), growth_probs[cancer].end(), 0.0);
  rates[3] = parameters.death_cancer   * sum_death_prob[cancer];

  rates[4] = parameters.birth_infected * sum_growth_prob[infected]; //std::accumulate(growth_probs[infected].begin(), growth_probs[infected].end(), 0.0);
  rates[5] = parameters.death_infected * sum_death_prob[infected];

  bool check_sum = false;
  if(check_sum) {
    for(size_t i = 0; i < 3; ++i) {
      float sum_growth_prob_exact = std::accumulate(growth_probs[i].begin(), growth_probs[i].end(), 0.f);
      float sum_growth_prob_indirect = sum_growth_prob[i];

      if(fabs(sum_growth_prob_exact - sum_growth_prob_indirect) > 1e-5f) {
        std::cout << "growth prob incorrect! " << i << "\t";
        std::cout << sum_growth_prob_exact << "\t" << sum_growth_prob_indirect << "\t";
        std::cout << sum_growth_prob_exact - sum_growth_prob_indirect << "\n";
      }

      /*float sum_death_prob_exact = std::accumulate(death_probs[i].begin(), death_probs[i].end(), 0.0);
      float sum_death_prob_indirect = sum_death_prob[i];
      float num_cell_types_counted = num_cell_types[i];
      if(sum_death_prob_exact != sum_death_prob_indirect) {
        std::cout << "death prob incorrect! " << i << "\t";
        std::cout << sum_death_prob_exact << "\t" << sum_death_prob_indirect << "\t" << num_cell_types_counted << "\t";
        std::cout << sum_death_prob_exact - sum_death_prob_indirect << "\n";
      }*/
    }
  }
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
      position_of_dying_cell = 0;
      break;
  }

  world[position_of_dying_cell].node_type = empty;
  death_probs[parent][position_of_dying_cell] = 0.f;
  sum_death_prob[parent]--; // update death probabilities

  update_growth_prob(position_of_dying_cell);
  for(auto i : world[position_of_dying_cell].neighbors) {
    update_growth_prob(i->pos);
  }
}

void simulation::implement_growth(const cell_type& parent) {
  size_t position_of_grown_cell = 0;
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
      position_of_grown_cell = 0;
      break;
  }

  world[position_of_grown_cell].node_type = parent;

  // update growth probability of cell
  update_growth_prob(position_of_grown_cell);

  // update growth probabilities of the neighbors:
  for(auto i : world[position_of_grown_cell].neighbors) {
    update_growth_prob(i->pos);
  }

  // update death probability of new cell:
  death_probs[parent][position_of_grown_cell] = 1.f;
  sum_death_prob[parent]++;
}

void simulation::add_cells(cell_type focal_cell_type) {
    // pick center node:
  if(focal_cell_type == cancer) {
      int a  = 5;
  }
  size_t x = sq_size / 2;
  size_t y = sq_size / 2;
  size_t focal_pos = x * sq_size + y;
  std::vector<size_t> cells_turned(1, focal_pos);
  world[focal_pos].node_type = focal_cell_type;
  size_t max_number_of_cells = parameters.initial_number_cancer_cells;
  if(focal_cell_type == normal) max_number_of_cells = static_cast<int>(parameters.num_cells * 0.1);
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

  for(auto i = cells_turned.begin(); i != cells_turned.end(); ++i) {
    update_growth_prob((*i));
    for(auto j = world[(*i)].neighbors.begin(); j != world[(*i)].neighbors.end(); ++j) {
      update_growth_prob((*j)->pos);
    }
  }

  for(int i = 0; i < 3; ++i) {
      sum_death_prob[i] = std::accumulate(death_probs[i].begin(), death_probs[i].end(), 0.f);
      sum_growth_prob[i] = std::accumulate(growth_probs[i].begin(), growth_probs[i].end(), 0.f);
      auto m = std::max_element(growth_probs[i].begin(), growth_probs[i].end());
      max_growth_prob[i] = *m;
  }


}

void simulation::update_growth_prob(size_t pos) {
  std::array<float, 3> probs = world[pos].calc_prob_of_growth();

  for(size_t i = 0; i < 3; ++i) {

    float new_prob =  probs[i] ;
    float old_prob =  growth_probs[i][pos];

    float diff_prob = new_prob - old_prob;
    sum_growth_prob[i] += diff_prob;

    growth_probs[i][pos] = new_prob;

    if(old_prob == max_growth_prob[i] && new_prob < max_growth_prob[i]) {
      // re-calculate the maximum:
      auto max_val = std::max_element(growth_probs[i].begin(), growth_probs[i].end());
      max_growth_prob[i] = *max_val;
    }
  }
}

simulation::simulation(const Param& param) {
  rns = rndutils::make_random_engine<>();
  parameters = param;

  sq_size = static_cast<size_t>(sqrt(parameters.num_cells));

  world.resize(parameters.num_cells);

  for(size_t i = 0; i < world.size(); ++i) {
    world[i].pos = i;
    world[i].set_coordinates(sq_size);
  }

  growth_probs.resize(3, std::vector< float >(parameters.num_cells, 0.f));
  death_probs.resize(3, std::vector< float >(parameters.num_cells, 0.f));

  max_growth_prob = {0.f, 0.f, 0.f};
}

void simulation::print_to_file(float t) {
  std::string file_name = "world_file_" + std::to_string(static_cast<int>(t)) + ".txt";
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

// initialization routines:
void simulation::update_growth_probabilities() {

  max_growth_prob = {-1.f, -1.f, -1.f};

  for(auto i : world) {
    std::array<float, 3> probs = i.calc_prob_of_growth();

    size_t pos = i.pos;

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

std::vector<int> simulation::get_cell_numbers() {
    std::vector<int> output(3,0);
    output[0] = static_cast<int>(sum_death_prob[0]);
    output[1] = static_cast<int>(sum_death_prob[1]);
    output[2] = static_cast<int>(sum_death_prob[2]);
    return output;
}

