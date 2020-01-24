#include <cstring>


#include "../Simulation/node.hpp"
#include "../Simulation/simulation.hpp"
#include "config_parser.h"

// forward declaration
void read_parameters_from_ini(Param& p, const std::string file_name);
void obtain_equilibrium(simulation& Simulation, const Param& all_parameters);
std::string get_outcome(const std::array<int, 5>& cell_counts);
std::string do_analysis(Param all_parameters);

int main(int argc, char *argv[]) {

    std::cout << "Welcome to this In Silico Simulation of oncolytic tumor virotherapy\n";
    std::cout << "Copyright D. Bhatt\n";
    std::cout << "Coding by Dr. T. Janzen\n";
    std::cout << "Supervision: Prof. Dr. F.J. Weissing\n";
    std::cout << "This is version: 0.5\n";



    std::cout << "All files are to be found in this folder: \n";
    std::cout << argv[0] << "\n";

    std::string file_name = "config.ini";

    Param all_parameters;

    read_parameters_from_ini(all_parameters, file_name);
    // Berg et al. varied two parameters:
    // lambda: speed of virus replication
    // varied in [0, 100]
    // delta: infected cell death
    // varied in [0, 15]
    // these are in the model:
    // all_parameters.birth_infected;
    // all_parameters.death_infected;

    // depending on computing power, several ways of exploring can be done:
    /*
     *  option 1:
     *  for(float lambda = 0; lambda < 100; lambda += 1) {
     *    for(float delta = 0; delta < 15;  delta += 0.5) {
     *          all_parameters.birth_infected = lambda;
     *          all_parameters.death_infected = delta;
     *          std::string outcome = do_analysis(all_parameters);
     *          std::ofstream outfile("output.txt", std::ios::app);
     *          outfile << lambda << "\t" << delta << "\t" << outcome << "\n";
     *          outfile.close();
     *     }
     *   }
     *
     * option 2:
     * // use command line arguments, then run the program many different times
     * // for instance on a computer cluster
     * std::string name_of_program = argv[0];
     * all_parameters.birth_infected = argv[1];
     * all_parameters.death_infected = argv[2];
     * std::string outcome = do_analysis(all_parameters);
     * std::ofstream outfile("output.txt", std::ios::app);
     * outfile << all_parameters.birth_infected << "\t"
     *         << all_parameters.death_infected << "\t"
     *         << outcome << "\n";
     * outfile.close();
     *
     * option 3:
     * // write a script that generates separate config files for each
     * // parameter setting
    */

   // std::string outcome = do_analysis(all_parameters);
   // std::cout << all_parameters.birth_infected << "\t"
   //           << all_parameters.death_infected << "\t"
   //           << outcome << "\n";


   std::string outcome = do_analysis(all_parameters);
   std::ofstream outfile("output.txt", std::ios::app);
   outfile << all_parameters.birth_infected << "\t"
           << all_parameters.death_infected << "\t"
           << outcome << "\n";
   outfile.close();


    return 0;
}

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
  std::array<int, 5> cell_counts;

  if(all_parameters.start_setup == converge) {
      std::cout << "simulating until having reached equilibrium with normal cells\n";
      obtain_equilibrium(Simulation, all_parameters); // this obtains a fully grown grid, with normal cells
  }

  Simulation.t = 0.f;
  bool cancer_added = false;
  bool virus_added = false;
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
             virus_added == false) {
              std::cout << "adding virus!\n";
             Simulation.add_cells(infected);
             Simulation.t = 0.f;
             virus_added = true;
          }

          if(static_cast<int>(Simulation.t) - static_cast<int>(prev_t) == 1) {
            std::cout << Simulation.t << "\t" << cell_counts[normal] << "\t" << cell_counts[cancer] << "\t" << cell_counts[infected] << "\n";
            cell_counts = Simulation.count_cell_types();
            if(cell_counts[cancer] < 1 && cancer_added == true && cell_counts[infected] < 1) {
                std::cout << "Simulation stopped because the cancer is gone\n";
                break; // stop if cancer is extinct
            }
            if(cell_counts[normal] < 1 && cell_counts[infected] < 1 && virus_added == true) {
                std::cout << "Simulation stopped because cancer has taken over\n";
                break; // stop if normal and virus are extinct
            }
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

void obtain_equilibrium(simulation& Simulation, const Param& all_parameters) {

  std::vector< float > densities(10, 0);
  float prev_t = Simulation.t;
  std::array<int, 5> cell_counts = Simulation.count_cell_types();
  int count = 0;
  while(Simulation.t < all_parameters.maximum_time) {
      Simulation.update_one_step();

      if(static_cast<int>(Simulation.t) - static_cast<int>(prev_t) == 10) {
           cell_counts = Simulation.count_cell_types();
           densities[count % 10] = cell_counts[normal];
           count++;
           if(count / 10 > 1) {
               float sum_first_half = 0.f;
               float sum_second_half = 0.f;
               for(size_t i = 0; i < 5; ++i) {
                   sum_first_half  += densities[i ];
                   sum_second_half += densities[i + 5];
               }
               std::cout << Simulation.t << "\t" << sum_first_half * 0.2f
                         << "\t" << sum_second_half * 0.2f << "\n";
               if(sum_first_half >= sum_second_half) {
                   break;
               }
           }
           prev_t = Simulation.t;
      }
  }
  return;
}





std::string get_outcome(const std::array<int, 5>& cell_counts) {
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

  if(freq[normal] <= 1e-6f && freq[cancer] >= (1-1e-6f) && freq[infected] <= 1e-6f) {
      return "B";
  }

  return "C";
}




void read_parameters_from_ini(Param& p, const std::string file_name) {

  ConfigFile from_config(file_name);


  p.maximum_time = from_config.getValueOfKey<int>("maximum_time");
  p.time_adding_cancer = from_config.getValueOfKey<int>("time_adding_cancer");
  p.time_adding_virus = from_config.getValueOfKey<int>("time_adding_virus");

  p.initial_number_cancer_cells = from_config.getValueOfKey<int>("initial_number_cancer_cells");
  p.initial_number_normal_cells = from_config.getValueOfKey<int>("initial_number_normal_cells");

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
