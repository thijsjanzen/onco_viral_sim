#include <cstring>
#include <chrono>

#include "../Simulation/node.hpp"
#include "../Simulation/simulation.hpp"
#include "config_parser.h"

// forward declaration
void read_parameters_from_ini(Param& p, const std::string file_name);
void obtain_equilibrium(simulation& Simulation, const Param& all_parameters);
std::string get_outcome(const std::array<size_t, 5>& cell_counts);
std::string do_analysis(Param all_parameters);

int main(int argc, char *argv[]) {

    std::cout << "Welcome to this In Silico Simulation of oncolytic tumor virotherapy\n";
    std::cout << "Copyright 2019 - 2020, D. Bhatt, T. Janzen & F.J. Weissing\n";
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
          << outcome                       << "\n";
  outfile.close();
  return 0;
}
