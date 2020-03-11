#define BOOST_TEST_MODULE Verify_code
#include <boost/test/included/unit_test.hpp>

#include "../Simulation/node.hpp"
#include "../Simulation/simulation.hpp"
#include "../Simulation/analysis.hpp"





BOOST_AUTO_TEST_CASE( Verification_of_function )
{
  //BOOST_TEST( true /* test assertion */ );
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.birth_infected = 0.f;
  all_parameters.birth_cancer = 0.f;



  std::string outcome = do_analysis(all_parameters);
  BOOST_CHECK_EQUAL(outcome, "A");
}
