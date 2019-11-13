#ifndef RANDOM_THIJS_H
#define RANDOM_THIJS_H

// Copyright 2019 Thijs Janzen

#include <random>
#include <map>
#include <vector>
#include <array>
#include <algorithm>

struct rnd_t {
    std::mt19937 rndgen;

    rnd_t() {
        std::random_device rd;
        std::mt19937 rndgen_t(rd());
        rndgen = rndgen_t;
    }

    std::cauchy_distribution<float> cauch_dist =
                                std::cauchy_distribution<float>(0.f, 0.01f);
    std::bernoulli_distribution bern_dist =
                                std::bernoulli_distribution(0.01);

    std::bernoulli_distribution bern_dist_zero =
                                std::bernoulli_distribution(0.001);

    std::uniform_real_distribution<float> unif_dist =
                                std::uniform_real_distribution<float>(0, 1.0);

    int random_number(int n)    {
        if(n <= 1) return 0;
        return std::uniform_int_distribution<> (0, n - 1)(rndgen);
    }

    float uniform()    {
        return unif_dist(rndgen);
    }

    float Expon(float lambda) {
        if (lambda == 0.0) return 1e20f;
        return std::exponential_distribution<float>(lambda)(rndgen);
    }

    float normal(float m, float s) {
        return std::normal_distribution<float>(m, s)(rndgen);
    }

    void set_seed(unsigned seed)    {
        std::mt19937 new_randomizer(seed);
        rndgen = new_randomizer;
    }

    bool bernouilli() {
        return bern_dist(rndgen);
    }

    bool bernouilli_zero() {
      return bern_dist_zero(rndgen);
    }

    float cauchy() {
        return cauch_dist(rndgen);
    }

    void set_cauchy(float m, float s) {
        cauch_dist = std::cauchy_distribution<float>(m, s);
    }

    void set_bernoulli(float p) {
        bern_dist = std::bernoulli_distribution(p);
    }

    void set_bernoulli_zero(float p) {
      bern_dist_zero = std::bernoulli_distribution(p);
    }

  std::array< std::normal_distribution<float>, 5> norm_dists;

  void set_up_norm_dists(const std::array< float, 5>& avg_waiting_times) {
    for(int i = 0; i < 5; ++i) {
      float mu = avg_waiting_times[i];
      float s = avg_waiting_times[i] / 10;
      norm_dists[i] = std::normal_distribution<float>(mu, s);
    }
  }

  float normal_t(size_t task) {
    float output = norm_dists[task](rndgen);
    while(output <= 0.f) output = norm_dists[task](rndgen);
    return output;
  }





  int draw_from_dist(const std::vector< float >& v, int max) {
      int max_index = static_cast<int>(v.size());
      while (true) {
          int index = random_number(max_index);
          if (max == 0.f) return index;

          if (uniform() < (1.f * v[index] / max)) {
              return index;
          }
      }
  }

  int draw_from_dist(const std::vector< float >& v, float min, float max) {
      int max_index = static_cast<int>(v.size());

      if (max == min) return random_number(max_index);
      float max_diff = max - min;

      while (true) {
          int index = random_number(max_index);

          if (uniform() < (1.f * (v[index] - min) / max_diff)) {
              return index;
          }
      }
  }

  int draw_from_dist_max(const std::vector< float >& v) {
    auto which_is_max = std::max_element(v.begin(), v.end());
    return static_cast<int>(std::distance(v.begin(), which_is_max));
  }

};

#endif  // RANDOM_THIJS_H
