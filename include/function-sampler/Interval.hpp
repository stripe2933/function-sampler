#pragma once

#include <cassert>

namespace function_sampler{
    struct Interval{
        double low;
        double high;

        constexpr double length() const noexcept { 
            assert(high > low);
            return high - low; 
        }
        constexpr bool empty() const noexcept { return low >= high; }
        constexpr bool contains(double x) const noexcept { return x >= low && x <= high; }
    };
};