#pragma once

#include "Point.hpp"
#include <cmath>

namespace function_sampler{
    struct Segment{
        Point point1;
        Point point2;

        constexpr double length() const noexcept { return std::hypot(point1.x - point2.x, point1.y - point2.y); }
    };
};