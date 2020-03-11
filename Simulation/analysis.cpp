#include "analysis.hpp"




std::string do_analysis(Param all_parameters) {

  // simulation procedure was as follows:
  // A completely normal population was simulated until the population stabilized
  // Five hundred cancer nodes were then inserted into the center of the network
  // Seven days later, five percent of cancer cells at the center of the tumor were infected.
  // Simulations were allowed to run until either the tumor populations died out or a time of 1,000 days was reached, which ever occurred first

  all_parameters.initial_number_cancer_cells = 500;
  all_parameters.percent_infected = 0.05f;
  all_parameters.time_adding_virus = 7 * 24; // 7 days of 24 hours after adding cancer
  all_parameters.maximum_time = 1000 * 24; // 1000 days of 24 hours.

  std::cout << "initializing grid\n";
  simulation Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);
  std::cout << "starting simulation\n";
  Simulation.t = 0.f;
  float prev_t = Simulation.t;
  std::array<size_t, 5> cell_counts;

  if(all_parameters.start_setup == converge) {
      std::cout << "simulating until having reached equilibrium with normal cells\n";
      Simulation.obtain_equilibrium(); // this obtains a fully grown grid, with normal cells
  }

  Simulation.t = 0.f;
  bool cancer_added = false;
  bool virus_added = false;

  auto prev_timepoint = std::chrono::steady_clock::now();
  auto start_t = prev_timepoint;

  while(Simulation.t < all_parameters.maximum_time) {
          Simulation.update_one_step();

          if(prev_t < all_parameters.time_adding_cancer &&
             Simulation.t >= all_parameters.time_adding_cancer &&
             cancer_added == false) {
              std::cout << "adding cancer!\n";
              Simulation.add_cells(cancer);
              Simulation.t = 0.f; // reset time
              cancer_added = true;
          }
          if(prev_t < all_parameters.time_adding_virus &&
             Simulation.t >= all_parameters.time_adding_virus &&
             virus_added == false &&
             cancer_added  == true) {
              std::cout << "adding virus!\n";
             Simulation.add_cells(infected);
             Simulation.t = 0.f;
             virus_added = true;
          }

          auto next_t = std::chrono::steady_clock::now();

          auto diff_t = std::chrono::duration_cast<std::chrono::seconds>(next_t - prev_timepoint).count();

          if(diff_t > 3) {
            cell_counts = Simulation.num_cell_types;

            auto total_t = std::chrono::duration_cast<std::chrono::seconds>(next_t - start_t).count();

              std::cout << static_cast<int>(Simulation.t) << "\t" << cell_counts[normal] << "\t"
                        << cell_counts[cancer] << "\t" << cell_counts[infected] << "\t" <<
                           "total time spent: " << total_t << " seconds\n";

              std::ofstream logfile("logfile.txt", std::ios::app);
              logfile << static_cast<int>(Simulation.t) << "\t" << cell_counts[normal] << "\t"
                                 << cell_counts[cancer] << "\t" << cell_counts[infected] << "\t" <<
                                    "total time spent: " << total_t << " seconds\n";
              logfile.close();

             // cell_counts = Simulation.count_cell_types();

              for(int i = 0; i < 5; ++i) {
                if(Simulation.num_cell_types[i] != cell_counts[i]) {
                    std::cout << "wrong count!\n";
                  }
              }

              if(cell_counts[cancer] < 1 && cancer_added == true && cell_counts[infected] < 1) {
                  std::cout << "Simulation stopped because the cancer is gone\n";
                  break; // stop if cancer is extinct
              }
              if(cell_counts[normal] < 1 && cell_counts[infected] < 1 && virus_added == true) {
                  std::cout << "Simulation stopped because cancer has taken over\n";
                  break; // stop if normal and virus are extinct
              }

              prev_timepoint = next_t;
            }
          //  if(Simulation.t > 100) update_freq = 10;
           // if(Simulation.t > 1000) update_freq = 100;
           // if(Simulation.t > 5000) update_freq = 1000;

          prev_t = Simulation.t;
  }

  // we want to know:
  // A: tumor eradication, e.g. only normal cells remain
  // B: tumor victory, e.g. only tumor cells remain
  // C: co-existence of the three populations
  cell_counts = Simulation.count_cell_types();
  std::string outcome = get_outcome(cell_counts);

  return outcome;
}


std::string get_outcome(const std::array<size_t, 5>& cell_counts) {
  std::array<float, 4> freq = {0.f, 0.f, 0.f, 0.f};
  freq[normal] = cell_counts[normal];
  freq[cancer] = cell_counts[cancer];
  freq[infected] = cell_counts[infected];
  float total_num_cells = cell_counts[normal] + cell_counts[cancer] + cell_counts[infected];
  // we want to know:
  // A: tumor eradication, e.g. only normal cells remain
  // B: tumor victory, e.g. only tumor cells remain
  // C: co-existence of the three populations
  for(size_t i = 0; i < 3; ++i) freq[i] *= 1.0f / total_num_cells;

  if(freq[normal] >= (1-1e-6f) && freq[cancer] <= 1e-6f && freq[infected] <= 1e-6f) {
      return "A";
  }

  if(freq[normal] <= 1e-6f     && freq[cancer] >= (1-1e-6f) && freq[infected] <= 1e-6f) {
      return "B";
  }

  return "C";
}


void read_parameters_from_ini(Param& p, const std::string file_name) {

  ConfigFile from_config(file_name);

  p.maximum_time = from_config.getValueOfKey<int>("maximum_time");
  p.time_adding_cancer = from_config.getValueOfKey<int>("time_adding_cancer");
  p.time_adding_virus = from_config.getValueOfKey<int>("time_adding_virus");

  p.initial_number_cancer_cells = from_config.getValueOfKey<size_t>("initial_number_cancer_cells");
  p.initial_number_normal_cells = from_config.getValueOfKey<size_t>("initial_number_normal_cells");

  p.birth_normal = from_config.getValueOfKey<float>("birth_normal");
  p.death_normal = from_config.getValueOfKey<float>("death_normal");

  p.birth_cancer = from_config.getValueOfKey<float>("birth_cancer");
  p.death_cancer = from_config.getValueOfKey<float>("death_cancer");

  p.birth_infected = from_config.getValueOfKey<float>("birth_infected");
  p.death_infected = from_config.getValueOfKey<float>("death_infected");

  p.birth_cancer_resistant = from_config.getValueOfKey<float>("birth_resistant");
  p.death_cancer_resistant = from_config.getValueOfKey<float>("death_resistant");

  p.percent_infected = from_config.getValueOfKey<float>("percent_infected");
  p.prob_normal_infection = from_config.getValueOfKey<float>("prob_normal_infection");
  p.freq_resistant = from_config.getValueOfKey<float>("freq_resistant");

  p.distance_infection_upon_death = from_config.getValueOfKey<float>("distance_infection_upon_death");
  p.prob_infection_upon_death = from_config.getValueOfKey<float>("prob_infection_upon_death");

  p.sq_num_cells = from_config.getValueOfKey<size_t>("sq_num_cells");

  p.infection_type = random_infection;
  auto infection_string = from_config.getValueOfKey<std::string>("infection_type");
  if(infection_string == "Random")
    p.infection_type = random_infection;
  if(infection_string == "Center")
    p.infection_type = center_infection;

  auto start_string = from_config.getValueOfKey<std::string>("start_type");
  if(start_string == "Grow")
    p.start_setup = grow;
  if(start_string == "Full")
    p.start_setup = full;

  auto grid_string = from_config.getValueOfKey<std::string>("grid_type");
  if(grid_string == "regular")
    p.use_voronoi_grid = false;
  if(grid_string == "voronoi")
    p.use_voronoi_grid = true;

  return;
}
