#include <cstring>


#include "../Simulation/node.hpp"
#include "../Simulation/simulation.hpp"
#include "config_parser.h"


void read_parameters_from_ini(Param& p, const std::string file_name);

int main(int argc, char *argv[])
{
    std::string file_name = "config.ini";


    Param all_parameters;

    read_parameters_from_ini(all_parameters, file_name);


    return 0;
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

 p.resolution = from_config.getValueOfKey<int>("resolution");

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

  return;
}
