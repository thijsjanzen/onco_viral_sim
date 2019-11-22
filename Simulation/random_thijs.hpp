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

    std::uniform_real_distribution<float> unif_dist =
                                std::uniform_real_distribution<float>(0, 1.0);

    int random_number(size_t n)    {
        if(n <= 1) return 0;
        return std::uniform_int_distribution<> (0, static_cast<int>(n - 1))(rndgen);
    }

    float uniform()    {
        return unif_dist(rndgen);
    }

    float Expon(float lambda) {
        if (lambda == 0.f) return 1e20f;
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

    float cauchy() {
        return cauch_dist(rndgen);
    }

    void set_cauchy(float m, float s) {
        cauch_dist = std::cauchy_distribution<float>(m, s);
    }

    void set_bernoulli(double p) {
        bern_dist = std::bernoulli_distribution(p);
    }


  template< typename it>
  size_t draw_from_dist(it begin, it end, float max) {
    size_t max_index = end - begin;
    float mult = 1.f / max;
    while (true) {
      size_t index = static_cast<size_t>(random_number(max_index));
      if (max == 0.f) return index;
      float prob = *(begin + index) * mult;

      if (uniform() < prob) {
          return index;
      }
    }
  }
};

#endif  // RANDOM_THIJS_H