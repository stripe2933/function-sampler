#pragma once

#include <cassert>
#include <algorithm>

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

        constexpr Interval &operator+=(double delta) noexcept { low += delta; high += delta; return *this; }
        constexpr Interval &operator-=(double delta) noexcept { low -= delta; high -= delta; return *this; }
        constexpr Interval operator+(double delta) const noexcept { return { low + delta, high + delta }; }
        constexpr Interval operator-(double delta) const noexcept { return { low - delta, high - delta }; }

        constexpr static Interval intersection(const Interval &interval1, const Interval &interval2) noexcept{
            return { std::max(interval1.low, interval2.low), std::min(interval1.high, interval2.high) };
        }
    };

    constexpr static Interval UNIT_INTERVAL { 0, 1 };
    constexpr static Interval REAL_INTERVAL { -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity() }; // TODO: can be simplified
};