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

  num_cell_types = {0, 0, 0, 0, num_cells}; // all cells are empty

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


size_t simulation::find_central_cell(const cell_type& focal_cell_type) {
    // first calculate average x and y of cell type:
    float x = 0.f;
    float y = 0.f;
    int counter = 0;
    for(const auto& i : world) {
        if(i.get_cell_type() == focal_cell_type) {
            x += i.x_;
            y += i.y_;
            counter++;
        }
    }
    x *= 1.0f / counter;
    y *= 1.0f / counter;

    std::vector< float > dist(world.size(), 1e9);
    for(size_t i = 0; i < world.size(); ++i) {
        if(world[i].get_cell_type() == focal_cell_type) {
            dist[i] = (world[i].x_ - x) * (world[i].x_ - x) + (world[i].y_ - y) * (world[i].y_ - y);
        }
    }

    auto min = std::min_element(dist.begin(), dist.end());
    return(std::distance(dist.begin(), min));
}


void simulation::initialize_network(std::vector< std::vector< voronoi_point > >& all_polys) {
   // initialize default.
  std::cout << "Initializing network\n";
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
      change_cell_type(i.pos, empty);
    }
  }
  if(parameters.use_voronoi_grid == true) {
      std::cout << "setting up Voronoi grid\n";
      setup_voronoi(all_polys);
      std::cout << "Done setting up Voronoi grid\n";
  }

  if(parameters.start_setup == grow || parameters.start_setup == converge) {
    std::cout << "adding initial cells\n";
    add_cells(normal);
    //for(auto& i : world) {
    for(size_t i = 0; i < num_cells; ++i) {
        update_growth_prob(i);
        update_death_prob(i);
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

void simulation::add_cells(const cell_type& focal_cell_type) {

  // pick center node:
  cell_type to_be_replaced = empty;
  if (focal_cell_type == cancer  ) to_be_replaced = normal;
  if (focal_cell_type == infected) to_be_replaced = cancer;

  if(num_cell_types[to_be_replaced] == 0) {
      return;
  }

  size_t focal_pos = find_central_cell(to_be_replaced);

  std::vector<size_t> cells_turned(1, focal_pos);

  change_cell_type(focal_pos, focal_cell_type);

  auto max_number_of_cells = static_cast<size_t>(parameters.initial_number_cancer_cells);

  if (focal_cell_type == normal) max_number_of_cells = parameters.initial_number_normal_cells;

  if (max_number_of_cells > world.size()) max_number_of_cells = world.size();

  if(parameters.start_setup == converge) {
      max_number_of_cells = world.size() * 0.9f;
  }

  size_t counter = 0;
  while (cells_turned.size() < max_number_of_cells) {
    focal_pos = cells_turned[counter];
    for (size_t i = 0; i < world[focal_pos].neighbors.size(); ++i) {
      size_t other_pos = world[focal_pos].neighbors[i]->pos;
      if (world[other_pos].get_cell_type() != focal_cell_type) {
        change_cell_type(other_pos, focal_cell_type);
        cells_turned.push_back(other_pos);
        if (cells_turned.size() >= max_number_of_cells) break;
      }
    }
    counter++;
  }

  for(size_t i = 0; i < num_cells; ++i) {
      update_growth_prob(i);
      update_death_prob(i);
  }
}

void simulation::infect_random() {
    int num_cancer_cells = num_cell_types[cancer];

    int infected_cells = 0;
    int to_be_infected = static_cast<int>(parameters.percent_infected * num_cancer_cells);
    if(to_be_infected == 0) return;

    std::vector< size_t > cancer_pos(num_cancer_cells);
    int j = 0;
    for(const auto& i : world) {
        if(i.get_cell_type() == cancer) {
            cancer_pos[j] = static_cast<size_t>(i.pos);
            j++;
        }
    }

    while(infected_cells < to_be_infected && num_cancer_cells > 0) {
        //size_t position_of_grown_cell = static_cast<size_t>(death_prob_rnd[cancer](rndgen.rndgen_));
        size_t index = static_cast<size_t>(rndgen.random_number(cancer_pos.size()));
        size_t position_of_infected_cell = cancer_pos[ index ];

        change_cell_type(position_of_infected_cell, infected);

        num_cancer_cells--; // easy count for now.
        infected_cells++;
        cancer_pos[index] = cancer_pos.back();
        cancer_pos.pop_back();
    }

    for(size_t i = 0; i < num_cells; ++i) {
        update_growth_prob(i);
        update_death_prob(i);
    }
}

void simulation::infect_center() {
    int num_cancer_cells = num_cell_types[cancer];

    size_t to_be_infected = static_cast<size_t>(parameters.percent_infected * num_cancer_cells);
    if(to_be_infected == 0) return;

    //now find starting cell to infect
    size_t focal_pos = find_central_cell(cancer);
    std::vector<size_t> cells_turned(1, focal_pos);

    change_cell_type(focal_pos, infected);

    size_t counter = 0;
    while(cells_turned.size() < to_be_infected && num_cancer_cells > 0) {
      focal_pos = cells_turned[counter];
      for(size_t i = 0; i < world[focal_pos].neighbors.size(); ++i) {
        size_t other_pos = world[focal_pos].neighbors[i]->pos;
        if(world[other_pos].get_cell_type() == cancer) {

           change_cell_type(other_pos, infected);
           cells_turned.push_back(other_pos);
           num_cancer_cells--;

           if(cells_turned.size() >= to_be_infected) break;
        }
      }
      counter++;
      if(counter > cells_turned.size()) break;
    }

    for(const auto& i : cells_turned) {
       update_growth_prob(i);
       update_death_prob(i);
    }
}

void simulation::infect_all_cancer() {
    for(auto& i : world) {
        if(i.get_cell_type() == cancer) {
              change_cell_type(i.pos, infected);
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
       change_cell_type(i.pos, normal);
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

void simulation::setup_voronoi(std::vector< std::vector< voronoi_point > >& all_polys) {
   using namespace cinekine;

   voronoi::Sites sites;

   std::cout << "Generating centre points\n";
   std::vector< voronoi_point > v(num_cells);

   for(size_t i = 0; i < num_cells; ++i) {
      float x = rndgen.uniform() * sq_size;
      float y = rndgen.uniform() * sq_size;

      v[i] = voronoi_point(x, y);
   }

   std::cout << "convering centre points to vertices\n";
   for(auto i : v) {
      voronoi::Vertex temp_vertex(i.x_, i.y_);
      sites.push_back(temp_vertex);
   }

   std::cout << "creating voronoi graph\n";
   voronoi::Graph graph = voronoi::build(std::move(sites), sq_size, sq_size);

   std::vector< std::vector< voronoi_edge > > all_edges(world.size());

   std::cout << "Ready to build world\n";
   std::cout << "collecting all edges\n";
   for(const auto& cell : graph.cells()) {

       size_t site_index = static_cast<size_t>(cell.site);

       voronoi::Site focal_site = graph.sites()[site_index];

       world[site_index].x_ = focal_site.x;
       world[site_index].y_ = focal_site.y;

       for(const auto& edge : cell.halfEdges) {
           voronoi::Edge focal_edge = graph.edges()[edge.edge];
           voronoi_point start(focal_edge.p0.x, focal_edge.p0.y);
           voronoi_point end(  focal_edge.p1.x, focal_edge.p1.y);

           voronoi_edge local_edge(start, end, focal_edge.leftSite, focal_edge.rightSite);
           if(!local_edge.check()) {
            std::cout << site_index << "\n";
           }

           all_edges[site_index].push_back(local_edge);
       }
   }

   std::cout << "implementing all edges\n";
   for(auto i : all_edges) {
       for(auto edge : i) {
           size_t left  = edge.left;
           size_t right = edge.right;

           if(left < world.size() && right < world.size()) {
               world[left].add_neighbor(world, right);
               world[right].add_neighbor(world, left);
           }
       }
   }

   std::cout << "clean edges and add polies for plotting\n";
   for(size_t i = 0; i < num_cells; ++i) {
       std::vector< voronoi_point > poly = clean_edges(all_edges[i], i);
       all_polys.push_back(poly);
   }

   std::cout << "update neighbor information\n";
   for(size_t i = 0; i < num_cells; ++i) {
       world[i].inv_num_neighbors = 1.f / world[i].neighbors.size();
       update_growth_prob(i);
       update_death_prob(i);
   }
}
