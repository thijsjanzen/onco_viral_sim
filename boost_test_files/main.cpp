#define BOOST_TEST_MODULE Verify_code
#include <boost/test/included/unit_test.hpp>

#include "../Simulation/node.hpp"
#include "../Simulation/simulation.hpp"
#include "../Simulation/analysis.hpp"

BOOST_AUTO_TEST_CASE( check_full_2 )
{
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.birth_infected = 0.2f;
  all_parameters.death_infected = 2.0f;

  // GROW SETUP
  all_parameters.start_setup = full;

  std::string outcome = do_analysis(all_parameters);
  BOOST_CHECK_EQUAL(outcome, "B");

  all_parameters.use_voronoi_grid = true;
  outcome = do_analysis(all_parameters);
  BOOST_CHECK_EQUAL(outcome, "B");
}

BOOST_AUTO_TEST_CASE( check_grow_2 )
{
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.birth_infected = 0.2f;
  all_parameters.death_infected = 2.0f;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 10000;
  all_parameters.time_adding_cancer = 1000;
  all_parameters.time_adding_virus = 2000;

  std::string outcome = do_analysis(all_parameters);
  BOOST_CHECK_EQUAL(outcome, "B");

  all_parameters.use_voronoi_grid = true;
  outcome = do_analysis(all_parameters);
  BOOST_CHECK_EQUAL(outcome, "B");
}

BOOST_AUTO_TEST_CASE( check_basic_functionality )
{
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;

  // set birth infected and birth cancer to 0
  all_parameters.birth_infected = 0.f;
  all_parameters.birth_cancer = 0.f;

  // FULL SETUP
  all_parameters.start_setup = full;

  std::string outcome = do_analysis(all_parameters);
  BOOST_CHECK_EQUAL(outcome, "A");

  all_parameters.use_voronoi_grid = true;
  outcome = do_analysis(all_parameters);
  BOOST_CHECK_EQUAL(outcome, "A");
}

BOOST_AUTO_TEST_CASE( check_grow )
{
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.birth_infected = 0.f;
  all_parameters.birth_cancer = 0.f;

  // GROW SETUP
  all_parameters.start_setup = grow;

  std::string outcome = do_analysis(all_parameters);
  BOOST_CHECK_EQUAL(outcome, "A");

  all_parameters.use_voronoi_grid = true;
  outcome = do_analysis(all_parameters);
  BOOST_CHECK_EQUAL(outcome, "A");
}








