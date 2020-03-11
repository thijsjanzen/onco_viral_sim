#ifndef ANALYSIS_H
#define ANALYSIS_H

#include "parameters.hpp"
#include "config_parser.h"
#include "simulation.hpp"

#include <array>
#include <string>


std::string do_analysis(Param all_parameters);
std::string get_outcome(const std::array<size_t, 5>& cell_counts);
void read_parameters_from_ini(Param& p, const std::string file_name);



#endif // ANALYSIS_H
