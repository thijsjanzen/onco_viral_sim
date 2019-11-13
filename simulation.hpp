//
//  simulation.hpp
//  Cancer_v1
//
//  Created by Thijs Janzen on 13/11/2019.
//  Copyright © 2019 Thijs Janzen. All rights reserved.
//

#ifndef simulation_hpp
#define simulation_hpp

#include <stdio.h>
#include <vector>
#include "parameters.hpp"
#include "random_thijs.hpp"
#include <array>

#include "node.hpp"

namespace model {
  class Simulation {
  public:
    Simulation(const Param& param);
    void run();
    void initialize_network();

  private:
    Param parameters;
    rnd_t rndgen;

    std::vector< std::vector< node > > world;
    std::vector< int > num_cell_types;
    void update_cell_types(); // this is very heavy, should not be ran often
    void update_rates(std::array< float, 6>& rates);
    int pick_event(const std::array< float, 6>& rates, float sum);
    void update_cell_counts(int event);
    void do_event(int event);
  };
}



#endif /* simulation_hpp */
