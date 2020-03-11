#include "../Simulation/node.hpp"
#include "../Simulation/simulation.hpp"
#include "../Simulation/analysis.hpp"

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <cassert>

  // Test the analysis module



BOOST_AUTO_TEST_CASE(virus_no_birth)
{
  std::clog << "Testing zero variance partitioning...\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.birth_infected = 0.f;
  all_parameters.birth_cancer = 0.f;

  std::string outcome = do_analysis(all_parameters);
  BOOST_CHECK_EQUAL(outcome, "A");

  //BOOST_CHECK_EQUAL(Xst({ 0.0, 0.0, 0.0 }, { 10u, 10u, 20u }), 0.0);
}

