#define BOOST_TEST_MODULE Verify_code
#include <boost/test/included/unit_test.hpp>

#include "../Simulation/node.hpp"
#include "../Simulation/simulation.hpp"
#include "../Simulation/analysis.hpp"



BOOST_AUTO_TEST_CASE( birth_death )
{
  std::cout << "testing birth and death\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);
  std::cout << "starting simulation\n";

  Simulation.t = 0.f;


  // test update_rates

  Simulation.test_update_rates();
  for(auto c : {normal, cancer, infected, resistant}) {
   BOOST_CHECK_EQUAL(Simulation.death_prob[c].get_total_sum(), 0.0);
   BOOST_CHECK_EQUAL(Simulation.growth_prob[c].get_total_sum(), 0.0);
   BOOST_CHECK_EQUAL(Simulation.get_rates(c), 0.0);
  }

  Simulation.test_change_cell_type(0, normal);
  Simulation.test_update_rates();
  BOOST_CHECK_EQUAL(Simulation.get_rates(1), all_parameters.death_normal);

  Simulation.test_change_cell_type(3, cancer);
  Simulation.test_update_rates();
  BOOST_CHECK_EQUAL(Simulation.get_rates(3), all_parameters.death_cancer);

  Simulation.test_change_cell_type(0, empty);
  Simulation.test_change_cell_type(3, empty);
  Simulation.test_update_rates();




  // test change_chell_type
  // add a normal cell, then kill it

  for(auto c : {normal, cancer, infected, resistant}) {
    std::cout << c << "\n";
    Simulation.test_change_cell_type(5000, c);

    BOOST_CHECK_EQUAL(Simulation.world[5000].get_cell_type(),
                      c);

    BOOST_CHECK_EQUAL(Simulation.death_prob[c].get_value(5000),
                      1.0);

    Simulation.test_change_cell_type(5000, empty);

    BOOST_CHECK_EQUAL(Simulation.world[5000].get_cell_type(),
                      empty);

    BOOST_CHECK_EQUAL(Simulation.death_prob[c].get_value(5000),
                      0.0);
   }

  // test do_event

  std::vector< cell_type > v = {normal, empty, cancer, empty,
                                infected, empty, resistant, empty};

   for(size_t i = 0; i < 7; ++i) {
    Simulation.test_event(i);

    BOOST_CHECK_EQUAL(Simulation.world[0].get_cell_type(),
                      v[i]);

    if (v[i] != empty) {
        BOOST_CHECK_EQUAL(Simulation.death_prob[v[i]].get_value(0),
                          1.0);
    }
   }


}






/*
BOOST_AUTO_TEST_CASE( t_cells )
{
  std::cout << "simulation with growth\n";
  std::cout << "testing t cell effect\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = center_infection;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 500;
  all_parameters.time_adding_virus = 100;

  all_parameters.t_cell_increase = 1.f;
  all_parameters.prob_infection_upon_death = 0.1f;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "C");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "C");
}

BOOST_AUTO_TEST_CASE( infect_second_time )
{
  std::cout << "simulation allowing growth from the start\n";
  std::cout << "using random infection\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = random_infection;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 500;
  all_parameters.time_adding_virus = 100;
  all_parameters.time_adding_virus_2 = 100;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "C");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  std::string outcome2 = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, outcome2);
}










BOOST_AUTO_TEST_CASE( check_grow_2 )
{
  std::cout << "simulation allowing growth from the start\n";
  std::cout << "Using weakened virus, expected tumor win\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.birth_infected = 0.2f;
  all_parameters.death_infected = 2.0f;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 100;
  all_parameters.time_adding_virus = 200;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "B");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "B");
}




BOOST_AUTO_TEST_CASE( set_infection )
{
  Param all_parameters;
  simulation Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);

  float target = 0.1f;
  Simulation.set_percent_infected(target);
  BOOST_CHECK_EQUAL(Simulation.get_parameters().percent_infected,
                    target);

  infection_routine test = random_infection;
  Simulation.set_infection_type(test);
  BOOST_CHECK_EQUAL(Simulation.get_parameters().infection_type,
                    test);
}


BOOST_AUTO_TEST_CASE( infect_random )
{
  std::cout << "simulation allowing growth from the start\n";
  std::cout << "using random infection\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = random_infection;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 500;
  all_parameters.time_adding_virus = 100;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "C");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  std::string outcome2 = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, outcome2);
}

BOOST_AUTO_TEST_CASE( infect_periphery )
{
  std::cout << "simulation allowing growth from the start\n";
  std::cout << "using random infection\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = periphery_infection;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 500;
  all_parameters.time_adding_virus = 100;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "C");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  std::string outcome2 = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, outcome2);
}



BOOST_AUTO_TEST_CASE( infect_all )
{
  std::cout << "simulation with growth\n";
  std::cout << "testing INFECT ALL\n";
  std::cout << "using random infection\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = random_infection;
  all_parameters.percent_infected = 1.f;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 500;
  all_parameters.time_adding_virus = 100;

  all_parameters.initial_number_normal_cells = 5000;
  all_parameters.initial_number_cancer_cells = 100;
  all_parameters.birth_cancer = 0.6f;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "A");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "A");
}

BOOST_AUTO_TEST_CASE( infect_center )
{
  std::cout << "simulation with growth\n";
  std::cout << "testing INFECT CENTER\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.infection_type = center_infection;

  // GROW SETUP
  all_parameters.start_setup = grow;
  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 500;
  all_parameters.time_adding_virus = 2;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "C");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "C");
}



BOOST_AUTO_TEST_CASE( test_converge )
{
  std::cout << "simulation with convergence\n";
  std::cout << "Using weakened virus, expected tumor win\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.birth_infected = 0.2f;
  all_parameters.death_infected = 2.0f;

  // GROW SETUP
  all_parameters.start_setup = converge;

  all_parameters.maximum_time = 1000;
  all_parameters.time_adding_cancer = 1000;
  all_parameters.time_adding_virus = 500;


  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "B");

  all_parameters.use_voronoi_grid = false;
  result = do_analysis(all_parameters);
  outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "B");
}

// TODO: add tests with long distance infection
BOOST_AUTO_TEST_CASE( long_distance_infection )
{

  std::cout << "simulation starting full\n";
  std::cout << "Using strong virus, expected virus win\n";
  std::cout << "Using long distance infection as well\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.prob_infection_upon_death = 0.5f;
  all_parameters.distance_infection_upon_death = 1.0f;

  // GROW SETUP
  all_parameters.start_setup = full;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "A");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  outcome = get_outcome(result);
   BOOST_CHECK_EQUAL(outcome, "A");
}

// TODO: add tests with resistant cells
BOOST_AUTO_TEST_CASE( resistance )
{
  std::cout << "Simulation with resistant cells\n";
  std::cout << "aiming for resistant persistence\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.prob_normal_infection = 0.01f;
  all_parameters.freq_resistant = 0.01f;

  all_parameters.start_setup = full;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
   BOOST_CHECK_EQUAL(outcome, "D");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "D");
}


BOOST_AUTO_TEST_CASE( check_full_2 )
{

  std::cout << "simulation starting full\n";
  std::cout << "Using weakened virus, expected tumor win\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.birth_infected = 0.2f;
  all_parameters.death_infected = 2.0f;

  // GROW SETUP
  all_parameters.start_setup = full;

  std::array<size_t, 5> result = do_analysis(all_parameters);
  std::string outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "B");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "B");
}
*/


// TODO:







