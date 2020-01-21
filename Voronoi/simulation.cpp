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

#include "voronoi.hpp"

#include "simulation.hpp"
#include "rndutils.hpp"


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



void simulation::run() {
  t = 0.f;
  while(t < parameters.maximum_time) {
     update_one_step();
  }
}

void simulation::initialize_network() {
  // we first make a regular network
  // we will do voronoi later...

  if(parameters.use_voronoi_grid) {
      setup_voronoi();
  } else {
      for(auto& i : world) {
        i.update_neighbors(world, sq_size);
        i.node_type = empty;
      }
  }

  if(parameters.start_setup == grow) {

    add_cells(normal);

    for(auto i : world) {
        update_growth_prob(i.pos);
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
  num_cell_types = {0, 0, 0, 0, 0};
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
  size_t position_of_dying_cell = static_cast<size_t>(death_prob_rnd[parent].draw_explicit(death_probs[parent].begin(),
                                                                                           death_probs[parent].end(),
                                                                                           rndgen));

  world[position_of_dying_cell].node_type = empty;
  update_death_prob(position_of_dying_cell);

  size_t min_pos = world.size();
  update_growth_prob(position_of_dying_cell);
  for(auto i : world[position_of_dying_cell].neighbors) {
    update_growth_prob(i->pos);
    if(i->pos < min_pos) min_pos = i->pos;
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


void simulation::implement_growth(const cell_type& parent) {

   // drawing position of growth is slow/bottleneck:
  size_t position_of_grown_cell = static_cast<size_t>(growth_prob_rnd[parent].draw_explicit(growth_probs[parent].begin(),
                                                                                            growth_probs[parent].end(),
                                                                                            rndgen));
  world[position_of_grown_cell].node_type = parent;

  // verify neighbors are of the right type
  int num_of_right_type = 0;
  for(auto i : world[position_of_grown_cell].neighbors) {
      if(i->node_type == parent) num_of_right_type++;
  }
  if(num_of_right_type == 0) {
      std::cout << "ABNORMAL BIRTH!\n";
      exit(0);
  }

  if(parent == cancer) {
      if(rndgen.uniform() < parameters.freq_resistant) {
          world[position_of_grown_cell].node_type = resistant;
      }
  }

  // update growth probability of cell
  update_growth_prob(position_of_grown_cell);
  update_death_prob(position_of_grown_cell);

  // update growth probabilities of the neighbors:
  for(auto i : world[position_of_grown_cell].neighbors) {
    update_growth_prob(i->pos);
  }
}

size_t find_central_cell(const std::vector< node > world, const cell_type& focal_cell_type) {
    // first calculate average x and y of cell type:
    float x = 0.f;
    float y = 0.f;
    int counter = 0;
    for(auto i : world) {
        if(i.node_type == focal_cell_type) {
            x += i.x_;
            y += i.y_;
            counter++;
        }
    }
    x *= 1.0f / counter;
    y *= 1.0f / counter;

    std::vector< float > dist(world.size(), 1e9);
    for(size_t i = 0; i < world.size(); ++i) {
        if(world[i].node_type == focal_cell_type) {
            dist[i] = (world[i].x_ - x) * (world[i].x_ - x) + (world[i].y_ - y) * (world[i].y_ - y);
        }
    }

    auto min = std::min_element(dist.begin(), dist.end());
    return(std::distance(dist.begin(), min));
}


void simulation::add_cells(const cell_type& focal_cell_type) {
    // pick center node:
//  size_t x = sq_size / 2;
 // size_t y = sq_size / 2;
 // size_t focal_pos = x * sq_size + y;\

  cell_type to_be_replaced = empty;
  if (focal_cell_type == cancer  ) to_be_replaced = normal;
  if (focal_cell_type == infected) to_be_replaced = cancer;

  int focal_pos = find_central_cell(world, to_be_replaced);


  std::vector<size_t> cells_turned(1, focal_pos);

  world[focal_pos].node_type = focal_cell_type;

  auto max_number_of_cells = static_cast<size_t>(parameters.initial_number_cancer_cells);
  if (focal_cell_type == normal) max_number_of_cells = parameters.initial_number_normal_cells;
  if (max_number_of_cells > world.size()) max_number_of_cells = world.size();
  size_t counter = 0;
  while (cells_turned.size() < max_number_of_cells) {
    focal_pos = cells_turned[counter];
    for (size_t i = 0; i < world[focal_pos].neighbors.size(); ++i) {
      size_t other_pos = world[focal_pos].neighbors[i]->pos;
      if (world[other_pos].node_type != focal_cell_type) {

        world[other_pos].node_type = focal_cell_type;
        cells_turned.push_back(other_pos);
        if (cells_turned.size() >= max_number_of_cells) break;
      }
    }
    counter++;
  }

  for (auto& i : cells_turned) {
    update_growth_prob(i);
    update_death_prob(i);
    for (auto& j : world[i].neighbors) {
      update_growth_prob(j->pos);
    }
  }
}

void simulation::update_growth_prob(size_t pos) {
  std::array<float, 4> probs = world[pos].calc_prob_of_growth();
  for (size_t i = 0; i < growth_prob_rnd.size(); ++i) {
      growth_probs[i][pos] = probs[i];

      growth_prob_rnd[i].update_entry(growth_probs[i].begin(), growth_probs[i].end(),
                                      pos, probs[i]);

  }
}

void simulation::update_death_prob(size_t pos) {
    for (size_t i = 0; i < 4; ++i) {
        float new_val = 0.f;
        if (i == world[pos].node_type) new_val = 1.f;
        death_probs[i][pos] = new_val;
        death_prob_rnd[i].update_entry(death_probs[i].begin(),
                                       death_probs[i].end(),
                                       pos,
                                       new_val);

    }
}


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

  for(auto i : world) {
    std::array<float, 4> probs = i.calc_prob_of_growth();

    size_t pos = i.pos;

    growth_probs[normal][pos] = probs[normal];
    growth_probs[cancer][pos] = probs[cancer];
    growth_probs[infected][pos] = probs[infected];
    growth_probs[resistant][pos] = probs[resistant];

    if(i.node_type != empty) {
        death_probs[i.node_type][pos] = 1.f;
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
    //size_t num_cancer_cells = static_cast<size_t>(std::accumulate(death_probs[cancer].begin(), death_probs[cancer].end(), 0));

    count_cell_types();
    int num_cancer_cells = num_cell_types[cancer];

    int infected_cells = 0;
    int to_be_infected = static_cast<int>(parameters.percent_infected * num_cancer_cells);
    if(to_be_infected == 0) return;

    std::vector< size_t > cancer_pos(num_cancer_cells);
    int j = 0;
    for(auto i : world) {
        if(i.node_type == cancer) {
            cancer_pos[j] = static_cast<size_t>(i.pos);
            j++;
        }
    }

    while(infected_cells < to_be_infected && num_cancer_cells > 0) {
        //size_t position_of_grown_cell = static_cast<size_t>(death_prob_rnd[cancer](rndgen.rndgen_));
        size_t index = static_cast<size_t>(rndgen.random_number(cancer_pos.size()));
        size_t position_of_infected_cell = cancer_pos[ index ];

        world[position_of_infected_cell].node_type = infected;

        num_cancer_cells--; // easy count for now.
        infected_cells++;
        cancer_pos[index] = cancer_pos.back();
        cancer_pos.pop_back();
    }

    for(auto i : world) {
        update_growth_prob(i.pos);
        update_death_prob(i.pos);
    }
}

void simulation::infect_center() {
    count_cell_types();
    int num_cancer_cells = num_cell_types[cancer];

    size_t to_be_infected = static_cast<size_t>(parameters.percent_infected * num_cancer_cells);
    if(to_be_infected == 0) return;

    //now find starting cell to infect
    size_t starting_pos = find_central_cell(world, cancer);
    while(world[starting_pos].node_type != cancer) {
        for(size_t i = 0; i < world[starting_pos].neighbors.size(); ++i) {
            if( world[starting_pos].neighbors[i]->node_type == cancer) {
                starting_pos = world[starting_pos].neighbors[i]->pos;
                break;
            }
        }
        int rand_index = rndgen.random_number(world[starting_pos].neighbors.size());
        starting_pos = world[starting_pos].neighbors[rand_index]->pos;
    }

    size_t focal_pos = starting_pos;
    std::vector<size_t> cells_turned(1, focal_pos);
    world[focal_pos].node_type = infected;
    size_t counter = 0;
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
      if(counter > cells_turned.size()) break;
    }

    for(auto i : world) {
        update_growth_prob(i.pos);
        update_death_prob(i.pos);
    }
}

void simulation::infect_all_cancer() {
    for(auto& i : world) {
        if(i.node_type == cancer) {
               i.node_type = infected;
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
        i.node_type = normal;
        i.prob_normal_infected = parameters.prob_normal_infection;
    }

    parameters.initial_number_cancer_cells = static_cast<int>(0.1f * world.size());
    add_cells(cancer);

    // just for safety, do a full scan:
    for(auto& i : world) {
        update_growth_prob(i.pos);
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

void simulation::set_infection_type(const infection_routine& infect_type) {
    parameters.infection_type = infect_type;
}

void simulation::set_percent_infected(const float& percent) {
    parameters.percent_infected = percent;
}

void simulation::setup_voronoi() {
   using namespace cinekine;

   voronoi::Sites sites;
   world.clear();

   std::vector< voronoi_point > v(num_cells);

   for(size_t i = 0; i < num_cells; ++i) {
    //  float x = i / sq_size;
    //   float y = i % sq_size;

     //  if((i % sq_size) % 2 == 0) x += 0.5;

      float x = rndgen.uniform() * sq_size;
      float y = rndgen.uniform() * sq_size;

      v[i] = voronoi_point(x, y);
   }

   for(auto i : v) {
      voronoi::Vertex temp_vertex(i.x_, i.y_);
      sites.push_back(temp_vertex);
   }

   voronoi::Graph graph = voronoi::build(std::move(sites), sq_size, sq_size);

   world.resize(num_cells);

   for(auto cell : graph.cells()) {

       size_t site_index = static_cast<size_t>(cell.site);

       voronoi::Site focal_site = graph.sites()[site_index];

       node temp_node(site_index, parameters.prob_normal_infection,
                      focal_site.x, focal_site.y);

       for(auto edge : cell.halfEdges) {
           voronoi::Edge focal_edge = graph.edges()[edge.edge];
           voronoi_point start(focal_edge.p0.x, focal_edge.p0.y);
           voronoi_point end(  focal_edge.p1.x, focal_edge.p1.y);
           voronoi_edge local_edge(start, end, focal_edge.leftSite, focal_edge.rightSite);
           temp_node.edges.push_back(local_edge);
       }
       world[site_index] = temp_node;
   }


   for(auto& i : world) {
       for(auto edge : i.edges) {
           int left = edge.left;
           int right = edge.right;

           if(left >= 0 && right >= 0 &&
              left < world.size() && right < world.size()){
               world[left].add_neighbor(world, right);
               world[right].add_neighbor(world, left);
           }
       }
   }

   for(auto& i : world) {
       i.clean_edges();
       update_growth_prob(i.pos);
       update_death_prob(i.pos);
   }

   for(size_t i = 0; i < growth_probs.size(); ++i) {
       growth_prob_rnd[i].update_row_cdf(growth_probs[i].begin());
       death_prob_rnd[i].update_row_cdf(death_probs[i].begin());
   }


   /// nice posterior check:
   ///
   float avg_distance = (1.0f * sq_size / sqrtf(num_cells)) * (1.0f * sq_size / sqrtf(num_cells)) * 2;
   for(auto i : world) {
      i.check_distances(avg_distance * 2.0f);
   }
}
