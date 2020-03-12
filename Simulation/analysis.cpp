#include "analysis.hpp"
#include <fstream>

std::string do_analysis(Param all_parameters) {

  // simulation procedure was as follows:
  // A completely normal population was simulated until the population stabilized
  // Five hundred cancer nodes were then inserted into the center of the network
  // Seven days later, five percent of cancer cells at the center of the tumor were infected.
  // Simulations were allowed to run until either the tumor populations died out or a time of 1,000 days was reached, which ever occurred first

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

  float prev_cast_t = static_cast<int>(Simulation.t);

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
             Simulation.add_infected();
             Simulation.t = 0.f;
             virus_added = true;
          }

          auto next_t = std::chrono::steady_clock::now();

          auto diff_t = std::chrono::duration_cast<std::chrono::seconds>(next_t - prev_timepoint).count();

          int cast_t = static_cast<int>(Simulation.t);

          if((cast_t - prev_cast_t) >= 100) {
              prev_cast_t = cast_t;

            cell_counts = Simulation.num_cell_types;

            auto total_t = std::chrono::duration_cast<std::chrono::seconds>(next_t - start_t).count();

              std::cout << static_cast<int>(Simulation.t) << "\t" <<
                           Simulation.num_cell_types[normal] << "\t" <<
                           Simulation.num_cell_types[cancer] << "\t" <<
                           Simulation.num_cell_types[infected] << "\t" <<
                           Simulation.num_cell_types[resistant] << "\t" <<
                           "total time spent: " << total_t << " seconds\n";

              std::ofstream logfile("logfile.txt", std::ios::app);
              logfile << static_cast<int>(Simulation.t) << "\t" <<
                         Simulation.num_cell_types[normal] << "\t" <<
                         Simulation.num_cell_types[cancer] << "\t" <<
                         Simulation.num_cell_types[infected] << "\t" <<
                         Simulation.num_cell_types[resistant] << "\t" <<
                         "total time spent: " << total_t << " seconds\n";
              logfile.close();

              if(all_parameters.start_setup == grow ||
                 all_parameters.start_setup == converge) {

                if(Simulation.num_cell_types[cancer] < 1 &&
                   virus_added == true) {
                    std::cout << "cancer eradicated\n";
                    break;
                }
                if(Simulation.num_cell_types[normal] < 1) {
                    std::cout << "normal tissue gone\n";
                    break;
                }
              }
              if(all_parameters.start_setup == full) {
                  if(Simulation.num_cell_types[cancer] < 1) {
                      std::cout << "cancer eradicated";
                      break;
                  }
                  if(Simulation.num_cell_types[normal] < 1) {
                      std::cout << "normal tissue gone\n";
                      break;
                  }
               }

              prev_timepoint = next_t;
            }

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
  freq[resistant] = cell_counts[resistant];
  float total_num_cells = cell_counts[normal] + cell_counts[cancer] + cell_counts[infected] + cell_counts[resistant];
  // we want to know:
  // A: tumor eradication, e.g. only normal cells remain
  // B: tumor victory, e.g. only tumor cells remain
  // C: co-existence of the three populations
  for(size_t i = 0; i < 4; ++i) freq[i] *= 1.0f / total_num_cells;

  if(freq[resistant] < 1e-6f) {

   /* if(freq[normal] >= (1-1e-6f) && freq[cancer] <= 1e-6f && freq[infected] <= 1e-6f) {
        return "A";
    }*/
    if(freq[cancer] <= 1e-6f) {
        return "A";
    }

    if(freq[normal] <= 1e-6f     && freq[cancer] >= (1-1e-6f) && freq[infected] <= 1e-6f) {
        return "B";
    }

    return "C";
  }

  return "D";


}
