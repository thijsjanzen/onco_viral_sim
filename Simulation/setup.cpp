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

  rndgen.set_seed(parameters.seed);
  std::ofstream seed_file("seed.txt", std::ios::app);
  seed_file << parameters.seed;
  seed_file.close();


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
  total_t_cell_concentration = 0.f;
  double lambda = static_cast<double>(1.0 / parameters.distance_infection_upon_death);
  for(size_t d = 1; d < sq_size; ++d) {
      double local_prob = static_cast<double>(parameters.prob_infection_upon_death) *
                          lambda * exp(-lambda * d);
      long_distance_infection_probability.push_back(local_prob);
    //  if(local_prob < 1e-5) break;
  }
}


size_t simulation::find_central_cell(const cell_type& focal_cell_type) const {
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

size_t simulation::find_central_cell(const std::vector< size_t >& positions) const {
    // first calculate average x and y of cell type:
    float x = 0.f;
    float y = 0.f;
    int counter = 0;
    for(const auto& i : positions) {
            x += world[i].x_;
            y += world[i].y_;
            counter++;
    }

    x *= 1.0f / counter;
    y *= 1.0f / counter;

    std::vector< float > dist(positions.size(), 1e9);
    for (size_t i = 0; i < positions.size(); ++i) {
      size_t pos = positions[i];
      dist[i] = (world[pos].x_ - x) * (world[pos].x_ - x) +
                 (world[pos].y_ - y) * (world[pos].y_ - y);
    }

    auto min = std::min_element(dist.begin(), dist.end());
    return(std::distance(dist.begin(), min));
}



void simulation::initialize_network(std::vector< std::vector< voronoi_point > >& all_polys) {
   // initialize default.
//   std::cout << "Initializing network\n";
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

  if(parameters.start_setup == converge && focal_cell_type == normal) {
      max_number_of_cells = world.size() * 0.95f;
  }

  size_t counter = 0;
  while (cells_turned.size() < max_number_of_cells &&
         counter < cells_turned.size()) {
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

void simulation::infect_random(float fraction) {
    int num_cancer_cells = num_cell_types[cancer];

    int infected_cells = 0;
    int to_be_infected = static_cast<int>(fraction * num_cancer_cells);
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

void remove_entries(std::vector< size_t >& source,
                    const std::vector< size_t >& to_remove) {
  for(auto i : to_remove) {
    for(size_t j = 0; j < source.size(); ++j)  {
       if(source[j] == i) {
           source[j] = source.back();
           source.pop_back();
           break;
       }
    }
  }
}


void simulation::infect_center_largest(float fraction) {
  size_t num_cancer_cells = num_cell_types[cancer];
  if (num_cancer_cells == 0) return;
  // here, we inject cancer in the largest tumour mass.
  // thus, we first have to collect all tumour cells, and group them.

  std::vector< size_t > cancer_cell_pos(num_cancer_cells);
  size_t cnt = 0;
  size_t j = 0;
  for (const auto& i : world) {
      if (i.get_cell_type() == cancer) {
          cancer_cell_pos[cnt] = j;
          cnt++;
      }
      j++;
  }

  std::vector< std::vector< size_t > > clusters;
  std::vector< size_t > cluster_sizes;
  while(!cancer_cell_pos.empty()) {
    std::vector< size_t > cluster;
    cluster.push_back(cancer_cell_pos[0]);

    for(size_t i = 0; i < cluster.size(); ++i) {
      std::vector< size_t > neighbours = world[cluster[i]].get_cancer_neighbours();
      for(size_t j = 0; j < neighbours.size(); ++j) {

          // this is not the fastest method, but it is very reliable.
          // probably a method using an unsorted set could work much faster.
          if(std::find(cluster.begin(), cluster.end(), neighbours[j]) ==
             cluster.end()) {
              cluster.push_back(neighbours[j]);
          }
      }
    }
    // now we have to remove these entries
    remove_entries(cancer_cell_pos, cluster);
    cluster_sizes.push_back(cluster.size());
    clusters.push_back(cluster);
  }

  auto m = std::max_element(cluster_sizes.begin(), cluster_sizes.end());

  size_t largest_cluster = 0;
  for (size_t i = 0; i < cluster_sizes.size(); ++i) {
      if (cluster_sizes[i] == *m) {
          largest_cluster = i;
          break;
       }
  }
  size_t to_be_infected = static_cast<size_t>(fraction * *m);
  if(to_be_infected == 0) return;

  if(*m < to_be_infected) {
      to_be_infected = *m;
  }

  size_t index_central_cell = find_central_cell(clusters[largest_cluster]);

  size_t focal_pos = clusters[largest_cluster][index_central_cell];
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

void simulation::infect_periphery(float fraction) {
  size_t num_cancer_cells = num_cell_types[cancer];
  if (num_cancer_cells == 0) return;
  // here, we inject cancer in the largest tumour mass.
  // thus, we first have to collect all tumour cells, and group them.

  std::vector< size_t > cancer_cell_pos(num_cancer_cells);
  size_t cnt = 0;
  size_t j = 0;
  for (const auto& i : world) {
      if (i.get_cell_type() == cancer) {
          cancer_cell_pos[cnt] = j;
          cnt++;
      }
      j++;
  }

  std::vector< std::vector< size_t > > clusters;
  std::vector< size_t > cluster_sizes;
  while(!cancer_cell_pos.empty()) {
    std::vector< size_t > cluster;
    cluster.push_back(cancer_cell_pos[0]);

    for(size_t i = 0; i < cluster.size(); ++i) {
      std::vector< size_t > neighbours = world[cluster[i]].get_cancer_neighbours();
      for(size_t j = 0; j < neighbours.size(); ++j) {

          // this is not the fastest method, but it is very reliable.
          // probably a method using an unsorted set could work much faster.
          if(std::find(cluster.begin(), cluster.end(), neighbours[j]) ==
             cluster.end()) {
              cluster.push_back(neighbours[j]);
          }
      }
    }
    // now we have to remove these entries
    remove_entries(cancer_cell_pos, cluster);
    cluster_sizes.push_back(cluster.size());
    clusters.push_back(cluster);
  }

  auto m = std::max_element(cluster_sizes.begin(), cluster_sizes.end());

  size_t largest_cluster = 0;
  for (size_t i = 0; i < cluster_sizes.size(); ++i) {
      if (cluster_sizes[i] == *m) {
          largest_cluster = i;
          break;
       }
  }

  // now we have the largest cluster.
  // let's now determine the distance to the center of the cluster.

  // first, we determine the center:
  float x = 0.f;
  float y = 0.f;
  int counter = 0;
  for(const auto& i : clusters[largest_cluster]) {
          x += world[i].x_;
          y += world[i].y_;
          counter++;
  }

  x *= 1.0f / counter;
  y *= 1.0f / counter;

  struct peri_cell {
    size_t pos;
    float freq_cancer;
    float dist_to_center;
  };

  std::vector< peri_cell > periphery;

  for (size_t i = 0; i < clusters[largest_cluster].size(); ++i) {
     size_t pos = clusters[largest_cluster][i];

     float freq_cancer = world[pos].freq_type_neighbours(cancer);

     if(freq_cancer < 1.0f &&
        world[pos].get_cell_type() == cancer) {
       peri_cell add;
       add.dist_to_center = (world[pos].x_ - x) * (world[pos].x_ - x) +
                            (world[pos].y_ - y) * (world[pos].y_ - y);
       add.pos = pos;
       add.freq_cancer = freq_cancer;
       periphery.push_back(add);
      }
    }

  size_t to_be_infected = static_cast<size_t>(fraction * periphery.size());
  if(to_be_infected == 0) return;

  std::sort(periphery.begin(), periphery.end(),
            [](auto const& a, auto const& b) {
              if (a.dist_to_center == b.dist_to_center) {
                 return a.freq_cancer < b.freq_cancer;
              }
              return a.dist_to_center > b.dist_to_center;});


  std::vector< size_t > cells_turned;
  for (size_t i = 0; i < to_be_infected; ++i) {
     size_t focal_pos = periphery[i].pos;
     if (world[focal_pos].get_cell_type() == cancer) {
       change_cell_type(focal_pos, infected);
       cells_turned.push_back(focal_pos);
     }
  }

  for(const auto& i : cells_turned) {
     update_growth_prob(i);
     update_death_prob(i);
  }
}


void simulation::infect_center(float fraction) {
 // infect_center_largest();
    int num_cancer_cells = num_cell_types[cancer];

    size_t to_be_infected = static_cast<size_t>(fraction * num_cancer_cells);
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


void simulation::add_infected(infection_routine infect_type,
                              float fraction) {

    if(fraction == 1.0f &&
       infect_type == random_infection) {
        infect_all_cancer();
        return;
    }

    switch(infect_type) {
        case random_infection:
            infect_random(fraction);
            break;
        case center_infection:
            infect_center_largest(fraction);
            break;
        case periphery_infection:
            infect_periphery(fraction);
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

    infect_center(0.1f);

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
           // uncomment to debug:
           // if(!local_edge.check()) {
           //  std::cout << site_index << "\n";
           // }

           if(local_edge.calc_dist() > 1e-2) {
             all_edges[site_index].push_back(local_edge);
           }
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
