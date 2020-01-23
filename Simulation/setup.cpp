#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>

#include "voronoi.hpp"

#include "simulation.hpp"

simulation::simulation(const Param& param) :
world(param.sq_num_cells * param.sq_num_cells)
{
  parameters = param;

  sq_size = parameters.sq_num_cells;
  num_cells = sq_size * sq_size;

  for (size_t i = 0; i < world.size(); ++i) {
    world[i].pos = i;
    world[i].set_coordinates(sq_size);
    world[i].prob_normal_infected = parameters.prob_normal_infection;
  }

  binned_distribution temp(sq_size, num_cells);
  for(size_t i = 0; i < 4; ++i) {
      growth_prob[i] = temp;
      death_prob[i] = temp;
  }

  long_distance_infection_probability = std::vector<double>(1, 0);
  double lambda = static_cast<double>(1.0 / parameters.distance_infection_upon_death);
  for(size_t d = 1; d < sq_size; ++d) {
      double local_prob = static_cast<double>(parameters.prob_infection_upon_death) *
                          lambda * exp(-lambda * d);
      long_distance_infection_probability.push_back(local_prob);
      if(local_prob < 1e-3) break;
  }
}


size_t find_central_cell(const std::vector< node > world, const cell_type& focal_cell_type) {
    // first calculate average x and y of cell type:
    float x = 0.f;
    float y = 0.f;
    int counter = 0;
    for(const auto& i : world) {
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


void simulation::initialize_network() {
   // initialize default.
  for(size_t i = 0; i < 4; ++i) {
       growth_prob[i] = binned_distribution(sq_size, num_cells);
       death_prob[i] = binned_distribution(sq_size, num_cells);
  }

  for(auto& i : world) {
       i.prob_normal_infected = parameters.prob_normal_infection;
  }

  if(parameters.use_voronoi_grid == false) {
    for(auto& i : world) {
      i.update_neighbors(world, sq_size);
      i.node_type = empty;
    }
  }
  if(parameters.use_voronoi_grid == true) {
      setup_voronoi();
  }

  if(parameters.start_setup == grow) {
    add_cells(normal);
    for(auto& i : world) {
        update_growth_prob(i.pos);
        update_death_prob(i.pos);
    }
  }

  if(parameters.start_setup == full) {
    initialize_full();
  }

  for(size_t i = 0; i < growth_prob.size(); ++i) {
      growth_prob[i].update_all();
      death_prob[i].update_all();
  }
}

// initialization routines:
void simulation::update_growth_probabilities() {

  for(auto& i : world) {
    std::array<float, 4> probs = i.calc_prob_of_growth();

    size_t pos = i.pos;

    growth_prob[normal].update_entry(pos, probs[normal]);
    growth_prob[cancer].update_entry(pos, probs[cancer]);
    growth_prob[infected].update_entry(pos, probs[infected]);
    growth_prob[resistant].update_entry(pos, probs[resistant]);

    if(i.node_type != empty) {
        death_prob[i.node_type].update_entry(pos, 1.f);
    }

  }
}

void simulation::add_cells(const cell_type& focal_cell_type) {

  // pick center node:
  cell_type to_be_replaced = empty;
  if (focal_cell_type == cancer  ) to_be_replaced = normal;
  if (focal_cell_type == infected) to_be_replaced = cancer;

  size_t focal_pos = find_central_cell(world, to_be_replaced);

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

  for(auto& i : world) {
      update_growth_prob(i.pos);
      update_death_prob(i.pos);
  }

}

void simulation::infect_random() {
    count_cell_types();
    int num_cancer_cells = num_cell_types[cancer];

    int infected_cells = 0;
    int to_be_infected = static_cast<int>(parameters.percent_infected * num_cancer_cells);
    if(to_be_infected == 0) return;

    std::vector< size_t > cancer_pos(num_cancer_cells);
    int j = 0;
    for(const auto& i : world) {
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

    for(const auto& i : world) {
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
        size_t rand_index = rndgen.random_number(world[starting_pos].neighbors.size());
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

    for(auto i : cells_turned) {
       update_cell(i);
    }
}

void simulation::infect_all_cancer() {
    for(auto& i : world) {
        if(i.node_type == cancer) {
               i.node_type = infected;
        }
    }

    for(const auto& i : world) {
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
}

void simulation::setup_voronoi() {
   using namespace cinekine;

   voronoi::Sites sites;
   world.clear();

   std::vector< voronoi_point > v(num_cells);

   for(size_t i = 0; i < num_cells; ++i) {
      float x = rndgen.uniform() * sq_size;
      float y = rndgen.uniform() * sq_size;

      v[i] = voronoi_point(x, y);
   }

   for(auto i : v) {
      voronoi::Vertex temp_vertex(i.x_, i.y_);
      sites.push_back(temp_vertex);
   }

   voronoi::Graph graph = voronoi::build(std::move(sites), sq_size, sq_size);

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
}



void simulation::set_infection_type(const infection_routine& infect_type) {
    parameters.infection_type = infect_type;
}

void simulation::set_percent_infected(const float& percent) {
    parameters.percent_infected = percent;
}
