#ifndef RANDOM_THIJS_H
#define RANDOM_THIJS_H

// Copyright 2019 Thijs Janzen

#include <random>
#include <map>
#include <vector>
#include <array>
#include <algorithm>
#include "rndutils.hpp"

struct rnd_t {
    rndutils::xorshift128 rndgen_;

    rnd_t() {
        std::random_device rd;
        rndgen_ = rndutils::make_random_engine(rd());
    }

    std::cauchy_distribution<float> cauch_dist =
                                std::cauchy_distribution<float>(0.f, 0.01f);
    std::bernoulli_distribution bern_dist =
                                std::bernoulli_distribution(0.01);

    rndutils::uniform01_distribution<float> rndutil_norm;


    int random_number(size_t n)    {
        if(n <= 1) return 0;
        return std::uniform_int_distribution<> (0, static_cast<int>(n - 1))(rndgen_);
    }

    float uniform()    {
        return rndutil_norm(rndgen_);
    }

    float Expon(float lambda) {
        if (lambda == 0.f) return 1e20f;
        return std::exponential_distribution<float>(lambda)(rndgen_);
    }

    float normal(float m, float s) {
        return std::normal_distribution<float>(m, s)(rndgen_);
    }

    int poisson(double lambda) {
        return std::poisson_distribution<int>(lambda)(rndgen_);
    }


    bool bernouilli() {
        return bern_dist(rndgen_);
    }

    float cauchy() {
        return cauch_dist(rndgen_);
    }

    void set_cauchy(float m, float s) {
        cauch_dist = std::cauchy_distribution<float>(m, s);
    }

    void set_bernoulli(double p) {
        bern_dist = std::bernoulli_distribution(p);
    }

    int binomial(int n, double p) {
        return std::binomial_distribution<>(n, p)(rndgen_);
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



template <int num_bins>
struct binned_distribution {
public:
    binned_distribution() {
        row_sum = {0.f};
        total_sum = 0.f;
        rows = rndutils::mutable_discrete_distribution< int, rndutils::all_zero_policy_nothing >(row_sum.begin(), row_sum.end());
    }

   template< typename It>
   binned_distribution(It first, It last) {
     int N = std::distance(first, last);
     bin_size = N / num_bins;
     for(int i = 0; i < num_bins; ++i) {
       auto start_it = first + i * bin_size;
       auto end_it = start_it + bin_size - 1;
       dist[i] = rndutils::mutable_discrete_distribution< int, rndutils::all_zero_policy_nothing >(start_it, end_it);
       row_sum[i] = dist[i].cdf().back();
     }
     rows = rndutils::mutable_discrete_distribution< int, rndutils::all_zero_policy_nothing >(row_sum.begin(), row_sum.end());
     total_sum = rows.cdf().back();
   }

   template <typename It>
   int draw_from_dist(It first, It last, float max_val, rnd_t& r) {
     int max_index = std::distance(first, last);
     while (true) {
       int index = r.random_number(static_cast<size_t>(max_index));
       float val = *(first + index);
       if (r.uniform() < (1.f * val / max_val)) {
         return index;
       }
     }
   }

   template <typename Reng>
   int operator()(Reng& reng) const
   {
    // draw a row
    size_t row = rows(reng);
    size_t col = dist[row](reng);
    return row * bin_size + col;
   }

   template <typename It>
   int draw_explicit(It first,
                     It last,
                     rnd_t& r)
   {
       rows.mutate(row_sum.begin(), row_sum.end());
       size_t row = rows(r.rndgen_);
       size_t col;
       float frac = row_sum[row] * 1.f / bin_size;
       auto start = first + row * bin_size;
       auto end = start + bin_size - 1;
       if(frac < 0.1f) {
        dist[row].mutate(start, end);
        col = dist[row](r.rndgen_);
       } else { 
        col = draw_from_dist(start, end, 1.f, r);
       }
       int result = row * bin_size + col;
    //   assert(result < std::distance(first, last));
   //    assert(result != 0);
       return result;
   }

   template< typename It>
   void mutate(It first, It last, size_t pos) {
       auto N = std::distance(first, last);
       bin_size = N / num_bins;
       size_t row = static_cast<size_t>(pos / bin_size);

       auto start_it = first + row * bin_size;
       auto end_it = start_it + bin_size - 1;

       dist[row].mutate(start_it, end_it);
       row_sum[row] = dist[row].cdf().back();
       rows.mutate(row_sum.begin(), row_sum.end());
   }

   template< typename It>
   void mutate_all(It first, It last) {
       auto N = std::distance(first, last);
       bin_size = N / num_bins;
       for(int i = 0; i < num_bins; ++i) {
         auto start_it = first + i * bin_size;
         auto end_it = start_it + bin_size - 1;
         dist[i].mutate(start_it, end_it);
         row_sum[i] = dist[i].cdf().back();
       }
       rows.mutate(row_sum.begin(), row_sum.end());
   }

   float get_total_sum() {
       return total_sum;
   }



   template< typename It>
   void update_row_sum(It first, int pos, float new_val) {
       float old_val = *(first + pos);
       if(old_val == new_val) return;
       size_t row = static_cast<size_t>(pos / bin_size);
       float adjust = new_val - old_val;
       row_sum[row] += adjust;
       if(row_sum[row] < 0) row_sum[row] = 0;
       total_sum += adjust;
 //      assert(row_sum[row] >= 0);
   }

   template< typename It>
   void update_entry(It first, It last, int pos, float new_val) {
       int N = std::distance(first, last);
       bin_size = N / num_bins;
       update_row_sum(first, pos, new_val);
   }

   template< typename It>
   void update_row_cdf(It first) {
       for(int i = 0; i < num_bins; ++i) {
           It start = first + i * bin_size;
           It end = start + bin_size - 1;
           row_sum[i] = std::accumulate(start, end, 0.0);
       }

       rows.mutate(row_sum.begin(), row_sum.end());
       total_sum = rows.cdf().back();
   }

 private:
    int bin_size;
    float total_sum;
    std::array<float, num_bins> row_sum;
    std::array< rndutils::mutable_discrete_distribution< int, rndutils::all_zero_policy_nothing >, num_bins > dist;
    rndutils::mutable_discrete_distribution< int, rndutils::all_zero_policy_nothing > rows;
};



#endif  // RANDOM_THIJS_H
