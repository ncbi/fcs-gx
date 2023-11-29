/*
-----------------------------------------------------------------------------
                             PUBLIC DOMAIN NOTICE
                 National Center for Biotechnology Information

  This software is a "United States Government Work" under the terms of the
  United States Copyright Act.  It was written as part of the author's official
  duties as a United States Government employees and thus cannot be copyrighted.
  This software is freely available to the public for use. The National Library
  of Medicine and the U.S. Government have not placed any restriction on its use
  or reproduction.

  Although all reasonable efforts have been taken to ensure the accuracy and
  reliability of this software, the NLM and the U.S. Government do not and
  cannot warrant the performance or results that may be obtained by using this
  software. The NLM and the U.S. Government disclaim all warranties, expressed
  or implied, including warranties of performance, merchantability or fitness
  for any particular purpose.

  Please cite NCBI in any work or product based on this material.

-----------------------------------------------------------------------------
*/
#pragma once
#include <vector>
#include <utility>
#include <cassert>

namespace gx
{

// https://en.wikipedia.org/wiki/Kahan_summation_algorithm
struct kahan_accumulator_t
{
    double m_curr_sum = 0;
    double m_err_term = 0;

    kahan_accumulator_t(double x)
    {
        m_curr_sum = x;
        m_err_term = 0;
    }

    kahan_accumulator_t& operator+=(double x)
    {
        x += m_err_term;
        const double orig_sum = m_curr_sum;

        m_curr_sum += x;
        const double added_x = m_curr_sum - orig_sum;

        m_err_term = x - added_x;
        return *this;
    }

    operator double() const
    {
        return m_curr_sum;
    }
};

// Implementation of Kadane and Ruzzo-Tompa algos.
template<typename Score = float, typename Pos = uint32_t>
struct max_subarray_sum
{
    struct node_t
    {
        // [start, end), corresponding to values in array I in reference implementation.
        std::pair<Pos, Pos> ivl = {};

        // Cumulative scores, corresponding to values in L and R arrays in reference implementation.
        // ::first is up to, not including ivl.first and second is up to, not including ::second.
        // The sum of values in the interval is therefore ::second - ::first.
        std::pair<Score, Score> score = {}; 

        Score get_score() const
        {
            return score.second - score.first;
        }

        bool operator==(const node_t& other) const { return ivl == other.ivl && score == other.score; };

    }; //__attribute__((packed));
    static_assert(sizeof(node_t) == (sizeof(Pos) + sizeof(Score)) * 2, "");
    using nodes_t = std::vector<node_t>;


    // Kadane's algorithm https://en.wikipedia.org/wiki/Maximum_subarray_problem
    template<typename F>
    static auto kadane(size_t n, F score_fn/*defined on [0, n)*/) -> node_t
    {
        static_assert(std::is_same_v<decltype(score_fn(0)), Score>, "");
        static_assert(std::is_signed_v<Score>, "");
        static_assert(sizeof(Score) >= 4, ""); // sanity-check.

        node_t curr = {};
        node_t best = {};

        using accumulator_t = std::conditional_t<std::is_floating_point_v<Score>, kahan_accumulator_t, Score>;
        accumulator_t total = 0;

        for (Pos i = 0; i < n; ++i) {
            total += score_fn(i);

            curr.score.second = (Score)total;
            curr.ivl.second = i + 1;

            if (best.get_score() < curr.get_score()) {
                best = curr;
            }
 
            if (curr.get_score() < 0) {
                curr = node_t{ {i + 1, i + 1}, {total, total} };
            }
        }

        return best;
    }


    // https://en.wikipedia.org/wiki/Ruzzo%E2%80%93Tompa_algorithm
    // https://homes.cs.washington.edu/~ruzzo/papers/maxseq.pdf
    // https://ieeexplore.ieee.org/document/6182645
    template<typename F>
    static auto ruzzo_tompa(size_t n, F score_fn/*defined on [0, n)*/, nodes_t arr = {}) -> nodes_t
    {
        arr.clear();

        static_assert(std::is_same_v<decltype(score_fn(0)), Score>, "");
        static_assert(std::is_signed_v<Score>, "");
        static_assert(sizeof(Score) >= 4, ""); // sanity-check.

        using accumulator_t = std::conditional_t<std::is_floating_point_v<Score>, kahan_accumulator_t, Score>;
        accumulator_t total = 0;

        size_t k = 0; // index of the last ivl in arr.

        const auto get_maxj = [&]
        {
            for (size_t j = k - 1; j < k; --j) {
                if (arr[j].score.first < arr[k].score.first) {
                    return j;
                }
            }
            return size_t(-1);
        };

        const auto update_loop_invariant = [&](Pos i)
        {
            for (auto maxj = get_maxj();
                 maxj != size_t(-1) && arr[maxj].score.second < arr[k].score.second;
                 maxj = get_maxj())
            {
                arr[maxj].ivl.second = i + 1;
                arr[maxj].score.second = (Score)total;
                k = maxj;
            }
            ++k;
        };

        for (Pos i = 0; i < n; ++i) {
            const Score s = score_fn(i);
            total += s;

            if (s > 0) {
                arr.resize(k + 1);
                arr[k] = node_t{{ i, i + 1 }, { total - s, total }};
                update_loop_invariant(i);
            }
        }

        arr.resize(k);
        return arr;
    }

    static void test()
    {
        const auto vec = std::vector<int>{{5, 15, -30, 10, -50, 40, -5, 10, -1000, 1}};
        const auto nodes = ruzzo_tompa(vec.size(), [&](size_t i){ return vec[i]; });
        const auto best_node =  kadane(vec.size(), [&](size_t i){ return vec[i]; });

        assert((nodes == nodes_t{ {{0, 2},  {0    , 20    }}, 
                                  {{3, 4},  {-10  , 0     }}, 
                                  {{5, 8},  {-50  , -5    }}, 
                                  {{9, 10}, {-1005, -1004 }} }));
        assert((best_node == nodes[2]));
    }
};

}
