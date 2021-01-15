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

  Simulation.t = 0.f;

  // count cell types
  std::array<size_t, 5> cells = Simulation.count_cell_types();
  for(size_t i = 0; i < 4; ++i) {
     BOOST_CHECK_EQUAL(cells[i], 0);
  }
  BOOST_CHECK_EQUAL(cells[4], all_parameters.sq_num_cells *
                              all_parameters.sq_num_cells *
                              all_parameters.sq_num_cells);


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
  //  std::cout << c << "\n";
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

   for(size_t i = 0; i < 8; ++i) {
    Simulation.test_event(i);

    BOOST_CHECK_EQUAL(Simulation.world[0].get_cell_type(),
                      v[i]);

    if (v[i] != empty) {
        BOOST_CHECK_EQUAL(Simulation.death_prob[v[i]].get_value(0),
                          1.0);
    }
   }

   // test pick event
   std::array<float, 8 > vx = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
   BOOST_CHECK_EQUAL(Simulation.test_pick_event(vx, 0.0), 0);
   for(size_t i = 0; i < 8; ++i) {
    vx[i] = 1.0f;
    BOOST_CHECK_EQUAL(Simulation.test_pick_event(vx, 1.0), i);
    vx[i] = 0.0f;
   }
}

BOOST_AUTO_TEST_CASE( check_outcome )
{
  std::cout << "testing check_outcome\n";
  std::array<size_t, 5> cell_counts_normal = {100000, 0, 0, 0, 0};

  BOOST_CHECK_EQUAL(get_outcome(cell_counts_normal), "A");

  std::array<size_t, 5> cell_counts_cancer = {0, 1000000, 0, 0, 0};

  BOOST_CHECK_EQUAL(get_outcome(cell_counts_cancer), "B");

  std::array<size_t, 5> cell_counts_resistant = {0, 0, 0, 1000000, 0};

  BOOST_CHECK_EQUAL(get_outcome(cell_counts_resistant), "C");
}


BOOST_AUTO_TEST_CASE( find_central_cell )
{
  std::cout << "testing finding center cells\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);

  // first we try to find the central empty cell (50, 50), coordinate: 5000
  auto index = Simulation.find_central_cell(empty);
  auto x = Simulation.world[index].x_;
  auto y = Simulation.world[index].y_;
  BOOST_CHECK_EQUAL(x, 49);
  BOOST_CHECK_EQUAL(y, 49);

  size_t row_size = all_parameters.sq_num_cells;
  // now we add a square of cancer cells, and calculate the center
  std::vector< size_t > positions;
  for(size_t x = 10; x < 21; ++x) {
      for(size_t y = 10; y < 21; ++y) {
          size_t pos = y * row_size + x;
          Simulation.test_change_cell_type(pos, cancer);
          positions.push_back(pos);
        }
  }
  auto index2 = Simulation.find_central_cell(cancer);
  auto x2 = Simulation.world[index2].x_;
  auto y2 = Simulation.world[index2].y_;
  BOOST_CHECK_EQUAL(x2, 15);
  BOOST_CHECK_EQUAL(y2, 15);

  auto index3 = Simulation.find_central_cell(positions);
  auto x3 = Simulation.world[positions[index3]].x_;
  auto y3 = Simulation.world[positions[index3]].y_;
  BOOST_CHECK_EQUAL(x2, x3);
  BOOST_CHECK_EQUAL(y2, y3);
}

BOOST_AUTO_TEST_CASE( add_cells )
{
  std::cout << "testing adding cells\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;
  all_parameters.initial_number_normal_cells = 1000;
  all_parameters.initial_number_cancer_cells = 100;

  simulation Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);

  std::array<size_t, 5> cell_cnt = Simulation.count_cell_types();
  BOOST_CHECK_EQUAL(cell_cnt[normal], 0);

  Simulation.add_cells(normal);
  cell_cnt = Simulation.count_cell_types();
  BOOST_CHECK_EQUAL(cell_cnt[normal],
                    all_parameters.initial_number_normal_cells);

  Simulation.add_cells(cancer);
  cell_cnt = Simulation.count_cell_types();
  BOOST_CHECK_EQUAL(cell_cnt[cancer],
                    all_parameters.initial_number_cancer_cells);
}

BOOST_AUTO_TEST_CASE( setup_types)
{
  std::cout << "testing setup types\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = full;

  simulation Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);
  std::array<size_t, 5> cell_cnt = Simulation.count_cell_types();

  size_t total_num_cells = all_parameters.sq_num_cells *
                           all_parameters.sq_num_cells *
                           all_parameters.sq_num_cells;

  BOOST_CHECK_EQUAL(cell_cnt[normal],
                    total_num_cells * 0.9);
  BOOST_CHECK_EQUAL(cell_cnt[cancer],
                    total_num_cells * 0.09);
  BOOST_CHECK_EQUAL(cell_cnt[infected],
                    total_num_cells * 0.01);
}

BOOST_AUTO_TEST_CASE( ask_infect_neighbours)
{
  std::cout << "testing ask infect neighbours\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = full;
  all_parameters.distance_infection_upon_death = 1.0;
  all_parameters.prob_infection_upon_death = 100.f;

  simulation Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);

  size_t row_size = all_parameters.sq_num_cells;
  size_t x = 50;
  for(size_t y = 40; y < 61; ++y) {
      size_t pos = y * row_size + x;
      Simulation.test_change_cell_type(pos, cancer);
  }
  size_t y2 = 50;
  for(size_t x2 = 40; x2 < 61; ++x2) {
      size_t pos = y2 * row_size + x2;
      Simulation.test_change_cell_type(pos, cancer);
  }

  // now we have a cross of individuals
  size_t initial_pos = 50 * 100 + 50;
  Simulation.test_change_cell_type(initial_pos, infected);
  size_t dist = 3;
  Simulation.test_ask_infect_neighbours(dist, initial_pos);
  for(size_t x = 40; x < 61; ++x) {
    size_t pos = 50 * row_size + x;
    int dist_x = static_cast<int>(x) - static_cast<int>(50);
    if (dist_x < 0) dist_x *= -1;

    auto ct = Simulation.world[pos].get_cell_type();
    if (dist_x > static_cast<int>(all_parameters.distance_infection_upon_death)) {
      BOOST_CHECK_EQUAL(ct, cancer);
    } else {
      BOOST_CHECK_EQUAL(ct, infected);
    }
  }
}

BOOST_AUTO_TEST_CASE( random_stuff )
{
  std::cout << "testing randomizer\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 10;

  simulation Simulation(all_parameters);

  // test random numbers
  Simulation.growth_prob[0].update_entry(10, 1.0f);
  auto x = Simulation.growth_prob[0].draw_explicit(Simulation.rndgen);
  BOOST_CHECK_EQUAL(x, 10);

  Simulation.growth_prob[0].update_entry(10, 0.0f);
  Simulation.growth_prob[0].update_entry(50, 0.01f);
  auto x2 = Simulation.growth_prob[0].draw_explicit(Simulation.rndgen);
  BOOST_CHECK_EQUAL(x2, 50);

  // now we populate with equal numbers
  int total_num_cells = all_parameters.sq_num_cells * all_parameters.sq_num_cells;
  for(int i = 0; i < total_num_cells; ++i) {
    Simulation.growth_prob[0].update_entry(i, 1.0f);
  }
  // and now we draw many numbers:
  int num_zero = 0;
  int num_ten  = 0;
  int num_repl = 100000;
  for(int r = 0; r < num_repl; ++r) {
      auto xx = Simulation.growth_prob[0].draw_explicit(Simulation.rndgen);
      if(xx == 0) num_zero++;
      if(xx == 10) num_ten++;
  }

  float expected_freq = 1.f / total_num_cells;
  float freq_zero     = 1.f * num_zero / num_repl;
  float freq_ten      = 1.f * num_ten  / num_repl;


  BOOST_CHECK_CLOSE(freq_zero,
                    expected_freq, 30.f);
  BOOST_CHECK_CLOSE(freq_ten,
                    expected_freq, 30.f);
  BOOST_CHECK_CLOSE(freq_zero,
                    freq_ten, 30.f);




}


BOOST_AUTO_TEST_CASE( diffuse )
{
  std::cout << "testing diffusion\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;
  all_parameters.t_cell_increase = 10.f;
  all_parameters.evaporation = 0.0f;
  all_parameters.diffusion = 0.1f;

  simulation Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);

  size_t x = 50;
  size_t y = 50;
  size_t z = 50;
  size_t pos = x + 100 * (y + 100 * z);

  Simulation.test_increase_t_cell_concentration(pos);
  Simulation.test_diffuse();

  // now all direct neighbors of (50, 50) should have 1.

  static int relative_points[6][3] = { {-1,  0,  0},
                                       { 1,  0,  0},
                                       { 0,  1,  0},
                                       { 0, -1,  0},
                                       { 0,  0,  1},
                                       { 0,  0, -1} };


  for(int i = 0; i < 6; ++i) {
      int other_x = static_cast<int>(x) + relative_points[i][0];
      int other_y = static_cast<int>(y) + relative_points[i][1];
      int other_z = static_cast<int>(y) + relative_points[i][2];

       int other_pos = other_x + static_cast<int>(100) * (other_y  +
                                                          static_cast<int>(100) * other_z);

        if(other_pos != pos)   {

          float local_conc = Simulation.world[other_pos].t_cell_concentration ;

          BOOST_CHECK_EQUAL( static_cast<int>(100 * local_conc) ,
                             static_cast<int>(100 * 0.166));
        }
   }

  Simulation.test_diffuse();
    for(int i = 0; i < 6; ++i) {
        int other_x = static_cast<int>(x) + relative_points[i][0];
        int other_y = static_cast<int>(y) + relative_points[i][1];
        int other_z = static_cast<int>(y) + relative_points[i][2];

         int other_pos = other_x + static_cast<int>(100) * (other_y  +
                                                            static_cast<int>(100) * other_z);

          if(other_pos != pos)   {

            float local_conc = Simulation.world[other_pos].t_cell_concentration ;

            BOOST_CHECK_EQUAL( static_cast<int>(ceil(100.0 * local_conc)) ,
                               static_cast<int>(100 * 0.3f) );  // ERROR
                                           //  0.2 + ((9 - 0.2) * 0.1) / 6 -
                                           // 5 * (0.2 * 0.1) / 6
          }
     }

}

BOOST_AUTO_TEST_CASE( infect_periphery )
{
  std::cout << "testing infect periphery\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          for (size_t z = 40; z < 60; ++z) {
            size_t pos = x + 100 * (y + z * 100);
            Simulation.test_change_cell_type(pos, cancer);
          }
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  BOOST_CHECK_EQUAL(cell_counts_before[cancer],
                    400);


  Simulation.test_infect_periphery(1.0);
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  BOOST_CHECK_LT(cell_counts_after[cancer],
                 cell_counts_before[cancer]); // after < before


  BOOST_CHECK_GT(cell_counts_after[infected],  // after > before
                 cell_counts_before[infected]);

  BOOST_CHECK_EQUAL(cell_counts_after[infected],
                    2168);  // 18 * 18 * 6 + 4 * 18 + 2 * 76
}

BOOST_AUTO_TEST_CASE( infect_random)
{
  std::cout << "test infect random\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  BOOST_CHECK_EQUAL(cell_counts_before[cancer],
                    400);

  // now we randomly infect some cells
  Simulation.test_infect_random(0.1);
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  BOOST_CHECK_LT(cell_counts_after[cancer],
                 cell_counts_before[cancer]); // after < before


  BOOST_CHECK_GT(cell_counts_after[infected],  // after > before
                 cell_counts_before[infected]);

  BOOST_CHECK_EQUAL(cell_counts_after[infected],
                    40);
}

BOOST_AUTO_TEST_CASE( infect_center )
{
  std::cout << "test infect center\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  BOOST_CHECK_EQUAL(cell_counts_before[cancer],
                    400);

  // now we randomly infect some cells
  Simulation.test_infect_center(0.1f);
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  BOOST_CHECK_LT(cell_counts_after[cancer],
                 cell_counts_before[cancer]); // after < before


  BOOST_CHECK_GT(cell_counts_after[infected],  // after > before
                 cell_counts_before[infected]);

  BOOST_CHECK_EQUAL(cell_counts_after[infected],
                    40);

  size_t central_pos = Simulation.find_central_cell(cancer);
  auto ctype = Simulation.world[central_pos].get_cell_type();
  BOOST_CHECK_EQUAL(ctype, cancer);
}

BOOST_AUTO_TEST_CASE( infect_center_largest)
{
  std::cout << "test infect center\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);

  // create two squares
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  for (size_t x = 20; x < 30; ++x) {
      for (size_t y = 20; y < 30; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  BOOST_CHECK_EQUAL(cell_counts_before[cancer],
                    500);

  Simulation.test_infect_center_largest(0.1f);
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  BOOST_CHECK_LT(cell_counts_after[cancer],
                 cell_counts_before[cancer]); // after < before


  BOOST_CHECK_GT(cell_counts_after[infected],  // after > before
                 cell_counts_before[infected]);

  BOOST_CHECK_EQUAL(cell_counts_after[infected],
                    40);

  auto index = Simulation.find_central_cell(infected);
  float x1 = Simulation.world[index].x_;
  float y1 = Simulation.world[index].y_;
  BOOST_CHECK_EQUAL(x1, 50);
  BOOST_CHECK_EQUAL(y1, 50);
}


BOOST_AUTO_TEST_CASE( voronoi )
{
  std::cout << "testing voronoi\n";
  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = true;
  all_parameters.start_setup = empty_grid;

  simulation Simulation(all_parameters);

  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);

  BOOST_CHECK_EQUAL(filler.size(), all_parameters.sq_num_cells *
                                   all_parameters.sq_num_cells *
                                   all_parameters.sq_num_cells);

  for(const auto& i : Simulation.world) {
        BOOST_CHECK_GT(i.neighbors.size(), 0);
  }
}

BOOST_AUTO_TEST_CASE( infect_all_cancer )
{
  // infect_all_cancer
  std::cout << "test infect all cancer\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  BOOST_CHECK_EQUAL(cell_counts_before[cancer],
                    400);

  // now we randomly infect some cells
  Simulation.test_infect_all_cancer();
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  BOOST_CHECK_EQUAL(cell_counts_after[cancer],
                    0); // after < before


  BOOST_CHECK_GT(cell_counts_after[infected],  // after > before
                 cell_counts_before[infected]);

  BOOST_CHECK_EQUAL(cell_counts_after[infected],
                    400);
}


// add_infected
BOOST_AUTO_TEST_CASE( add_infected )
{
  std::cout << "test add infected\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  BOOST_CHECK_EQUAL(cell_counts_before[cancer],
                    400);

  // now we randomly infect some cells
  Simulation.add_infected(random_infection, 0.1f);
  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();

  BOOST_CHECK_LT(cell_counts_after[cancer],
                    cell_counts_before[cancer]); // after < before


  BOOST_CHECK_GT(cell_counts_after[infected],  // after > before
                 cell_counts_before[infected]);

  BOOST_CHECK_EQUAL(cell_counts_after[infected],
                    40);
}


// update_one_step
BOOST_AUTO_TEST_CASE( update_one_step)
{
  std::cout << "test update one step\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.death_cancer = 0.0;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;

  simulation Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  std::array<size_t, 5> cell_counts_before = Simulation.count_cell_types();
  BOOST_CHECK_EQUAL(cell_counts_before[cancer],
                    400);

  Simulation.update_one_step();

  std::array<size_t, 5> cell_counts_after = Simulation.count_cell_types();
  BOOST_CHECK_EQUAL(cell_counts_after[cancer],
                    401);
}



// infect_long_distance
BOOST_AUTO_TEST_CASE( infect_long_distance ) {
  std::cout << "test infect long distance\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = empty_grid;
  all_parameters.distance_infection_upon_death = 1.0;
  all_parameters.prob_infection_upon_death = 100.f;

  simulation Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);

  // create a square
  for (size_t x = 40; x < 60; ++x) {
      for (size_t y = 40; y < 60; ++y) {
          size_t pos = x + y * 100;
          Simulation.test_change_cell_type(pos, cancer);
      }
  }

  size_t x = 50;
  size_t y = 50;
  size_t central_pos = y * 100 + x;
  //Simulation.test_change_cell_type(central_pos, infected);
  Simulation.test_infect_long_distance(central_pos);

  static int relative_points[4][2] = { {-1, 0},
                                      {1, 0},
                                      {0, 1},
                                      {0, -1} };

  for (size_t i = 0; i < 4; ++i) {
     int x2 = x + relative_points[i][0];
     int y2 = y + relative_points[i][1];
     size_t pos = y2 * 100 + x2;
     BOOST_CHECK_EQUAL(Simulation.world[pos].get_cell_type(),
                       infected);
  }
}

BOOST_AUTO_TEST_CASE( obtain_equilibrium )
{
  std::cout << "test obtain equilibrium\n";

  Param all_parameters;
  all_parameters.sq_num_cells = 100;
  all_parameters.use_voronoi_grid = false;
  all_parameters.start_setup = grow;
  all_parameters.initial_number_normal_cells = 5000;

  simulation Simulation(all_parameters);
  std::vector< std::vector< voronoi_point > > filler;

  Simulation.initialize_network(filler);

  Simulation.obtain_equilibrium();

  std::array<size_t, 5> ctypes = Simulation.count_cell_types();

  BOOST_CHECK_GT(ctypes[normal], 0.0);
  BOOST_CHECK_GT(ctypes[normal], 4000);
  BOOST_CHECK_LT(ctypes[normal], 10000);
  BOOST_CHECK_LT(ctypes[normal], 6000);
}


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
  BOOST_CHECK_EQUAL(outcome, "A");

  all_parameters.use_voronoi_grid = true;
  result = do_analysis(all_parameters);
  outcome = get_outcome(result);
  BOOST_CHECK_EQUAL(outcome, "A");
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


BOOST_AUTO_TEST_CASE( infect_random2 )
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

BOOST_AUTO_TEST_CASE( infect_periphery2 )
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

BOOST_AUTO_TEST_CASE( infect_center2 )
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
