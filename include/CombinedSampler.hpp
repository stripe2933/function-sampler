#pragma once

#include <random>
#include <stack>
#include <optional>
#include <numbers>
#include <concepts>

#include "concept/Point.hpp"
#include "concept/Interval.hpp"
#include "Point.hpp"
#include "Interval.hpp"

namespace function_sampler{
    template <std::floating_point NumberType = double, typename PointType = Point, typename IntervalType = Interval> requires
        concepts::Point<PointType, NumberType> &&
        concepts::Interval<IntervalType, NumberType>
    class CombinedSampler{
    private:
        // type definitions.
        using number_type = NumberType;
        using point_type = PointType;
        using interval_type = IntervalType;
        using function_type = std::function<number_type(number_type)>;

        // constants.
        constexpr static double PI = std::numbers::pi_v<double>;

        /**
         * @brief Check if a function is singular at a given point.
         * 
         * @param function function to check
         * @param value value to check singularity
         * @return true the function is singular at the given point.
         * @return false the function is nonsingular at the given point.
         */
        bool is_singular(const function_type& function, number_type x) const noexcept{
            const auto h = numerical_threshold / 2;
            const auto xl2 = x - 2 * h, xl1 = x - h, xr1 = x + h, xr2 = x + 2 * h;

            for (const auto xi : { xl2, xl1, x, xr1, xr2 }){
                const auto y = function(xi);

                if (std::abs(y) > infinity_threshold || std::isnan(y)){
                    return true;
                }
            }

            // LR discontinuity check: if LR > 0.8 -> discontinue
            const auto fr = 3 * function(x) - 4 * function(xr1) + function(xr2),
                       fl = 3 * function(x) - 4 * function(xl1) + function(xl2);
            const auto lr = std::abs(fr * fr - fl * fl) / (fr * fr + fl * fl);

            return lr > 0.8;
        }

        /**
         * @brief Calculate a complementary of an acute angle (i.e. pi - (acute angle)) consisted by the three points in order, \p p1 - \p p2 - \p p3 .
         * 
         * @param p1 point 1
         * @param p2 point 2
         * @param p3 point 3
         * @return number_type the complementary of the acute angle in radian unit
         */
        static number_type refinement(const point_type& p1, const point_type& p2, const point_type& p3) noexcept{
            const auto x1 = p1.x - p2.x, y1 = p1.y - p2.y,
                       x2 = p3.x - p2.x, y2 = p3.y - p2.y;
            const auto dot = x1 * x2 + y1 * y2,
                       det = x1 * y2 - y1 * x2;
            
            const auto angle = std::abs(std::atan2(det, dot));
            return PI - angle;
            
        }

        /**
         * @return std::optional<number_type> If the discontinuity detected, return the x value of discontinuity. Otherwise, return std::nullopt.
         */
        template <typename OutputIt>
        std::optional<number_type> init_sampling(const function_type& function, const interval_type& domain, OutputIt output) const noexcept{
            for (const auto x : { domain.low, domain.high }){
                if (is_singular(function, x)){
                    return x;
                }
            }

            const point_type pa { domain.low, function(domain.low) },
                             pb { domain.high, function(domain.high) };
            *output = pa;
            auto result = recursive_sampling(function, pa, pb, 1, output);
            if (result){
                return result;
            }
            *output = pb;
            return std::nullopt;
        }

        /**
         * @return std::optional<number_type> If the discontinuity detected, return the x value of discontinuity. Otherwise, return std::nullopt.
         */
        template <typename OutputIt>
        std::optional<number_type> recursive_sampling(const function_type& function, const point_type& pa, const point_type& pb, uint depth, OutputIt output) const noexcept{
            const auto x_diff = pb.x - pa.x;
            if (depth > max_depth || x_diff < numerical_threshold){
                return std::nullopt;
            }

            // random number generation
            static std::random_device rd;
            static std::mt19937 gen { rd() };
            static std::uniform_real_distribution<double> dis { 0.45, 0.55 };

            const auto r1 = dis(gen), r2 = dis(gen), r3 = dis(gen);
            const auto x1 = pa.x + 0.5 * r1 * x_diff,
                       x2 = pa.x + r2 * x_diff,
                       x3 = pa.x + 1.5 * r3 * x_diff; // four quarters of [pa.x, pb.x] : [pa.x, x1], [x1, x2], [x2, x3], [x3, pb.x]

            for (const auto xi : { x1, x2, x3 }){
                if (is_singular(function, xi)){
                    return xi;
                }
            }

            const auto y1 = function(x1), y2 = function(x2), y3 = function(x3);
            const point_type p1 { x1, y1 }, p2 { x2, y2 }, p3 { x3, y3 };
            const auto a1 = refinement(pa, p1, p2),
                       a2 = refinement(p1, p2, p3),
                       a3 = refinement(p2, p3, pb);

            if (a1 > refinement_threshold || depth <= min_depth){
                auto result = recursive_sampling(function, pa, p1, depth + 1, output);
                if (result){
                    return result;
                }
            }
            *output = p1;
            if (a1 > refinement_threshold || a2 > refinement_threshold || depth <= min_depth){
                auto result = recursive_sampling(function, p1, p2, depth + 1, output);
                if (result){
                    return result;
                }
            }
            *output = p2;
            if (a2 > refinement_threshold || a3 > refinement_threshold || depth <= min_depth){
                auto result = recursive_sampling(function, p2, p3, depth + 1, output);
                if (result){
                    return result;
                }
            }
            *output = p3;
            if (a3 > refinement_threshold || depth <= min_depth){
                auto result = recursive_sampling(function, p3, pb, depth + 1, output);
                if (result){
                    return result;
                }
            }

            return std::nullopt;
        }

    public:
        number_type numerical_threshold = 1e-4; // the lower bound of a sampling interval; If the interval's length is less than this, stop sampling.
        uint min_depth = 3; // the lower bound of a sampling depth; if depth is deeper than this, stop sampling.
        uint max_depth = 8; // the upper bound of a sampling depth; if depth is shallower than this, continue sampling even if the other conditions. (ex. refinement, ...)
        number_type refinement_threshold = 0.05; // the upper bound of the adjusting points' refinement; if the refinement is greater than this, continue sampling.
        number_type infinity_threshold = 1e8; // the criterion for determining infinity
        uint max_singularity_count = 500; // the max value of singularity count in a sampling; If the singularity count is more than this, stop sampling and throw error (exceptions::too_many_singularity).

        /**
         * @brief Sample a function points from a domain.
         * 
         * @param function the function to sample
         * @param domain the domain of the function
         * @return a vector that consists of continuous points
         */
        template <typename FunctionPP> requires
            std::is_convertible_v<FunctionPP, function_type>
        std::vector<std::vector<point_type>> sample_points(FunctionPP&& function, const interval_type& domain) const{
            std::stack<interval_type> interval_stack;
            interval_stack.push(domain);

            std::vector<std::vector<point_type>> result;
            while (interval_stack.empty() == false){
                const auto xrange = interval_stack.top();
                interval_stack.pop();

                std::vector<point_type> points;
                auto sampling_result = init_sampling(function, xrange, std::back_inserter(points));
                if (sampling_result){
                    if (result.size() > max_singularity_count){
                        result.clear();
                        throw std::runtime_error { "Too many singularity detected" };
                    }
                    else{
                        auto discontinuity_x = sampling_result.value();
                        /*
                            divide the range with discontinuity_x

                            if discontinuity_x is near to the xrange.low -> process the range as [xrange.low, discontinuity_x) (same for high)
                            otherwise                                    -> [xrange.low, discontinuity_x), (discontinuity_x, xrange.high]
                        */
                        if (xrange.low <= xrange.high && xrange.length() > numerical_threshold){
                            if (xrange.low <= discontinuity_x && std::abs(discontinuity_x - xrange.low) <= numerical_threshold){
                                interval_stack.emplace(interval_type { xrange.low + numerical_threshold, xrange.high });
                            }
                            else if (xrange.high >= discontinuity_x && std::abs(discontinuity_x - xrange.high) <= numerical_threshold){
                                interval_stack.emplace(interval_type { xrange.low, xrange.high - numerical_threshold });
                            }
                            else if (xrange.low < discontinuity_x && discontinuity_x < xrange.high && std::abs(discontinuity_x - xrange.low) > numerical_threshold && std::abs(discontinuity_x - xrange.high) > numerical_threshold){
                                interval_stack.emplace(interval_type { xrange.low, discontinuity_x - numerical_threshold });
                                interval_stack.emplace(interval_type { discontinuity_x + numerical_threshold, xrange.high });
                            }
                        }
                    }
                }
                else{
                    result.emplace_back(std::move(points));
                }
            }

            return result;
        }

        /**
         * @brief Generate line segments from sampled points.
         * 
         * @param sampled_points sampled points
         * @param output the output iterator to store generated segments
         * @param adapter the transformation function (point_type -> VertexType)
         */
        template <typename VertexType, typename OutputIt, typename AdapterPP> requires
            std::is_convertible_v<AdapterPP, std::function<VertexType(point_type)>> // for perfect forwarding
        static void generate_segments(const std::vector<std::vector<point_type>>& sampled_points, OutputIt output, AdapterPP&& adapter) noexcept{
            for (const auto& continuous_points : sampled_points){
                *output = std::move(adapter(continuous_points.front()));

                for (auto it = continuous_points.cbegin() + 1; it != continuous_points.cend() - 1; ++it){
                    const auto& vertex = adapter(*it);
                    *output = vertex; // if using std::move, the original object must be invalidated
                    *output = vertex;
                }

                *output = std::move(adapter(continuous_points.back()));
            }
        }

        /**
         * @brief Generate line segments from sampled points.
         * 
         * @param sampled_points sampled points
         * @param output the output iterator to store generated segments
         * 
         * @note This function calls an overloaded generated_segments function with adapter = std::identity { }.
         */
        template <typename OutputIt>
        static void generate_segments(const std::vector<std::vector<point_type>>& sampled_points, OutputIt output) noexcept{
            generate_segments(sampled_points, output, std::identity { });
        }

        /**
         * @brief Generate line segments from sampled points.
         * 
         * @param sampled_points sampled points
         * @param vertices the output vector to store generated segments
         * @param adapter the transformation function (point_type -> VertexType)
         * 
         * @note This function is faster than generate_segments(const vector<vector<point_type>>, OutputIt) because it reserves the output vertice count.
         */
        template <typename VertexType, typename AdapterPP> requires
            std::is_convertible_v<AdapterPP, std::function<VertexType(point_type)>> // for perfect forwarding
        static void generate_segments(const std::vector<std::vector<point_type>>& sampled_points, std::vector<VertexType>& vertices, AdapterPP&& adapter) noexcept{
            // vertices are two endpoints of the line segment
            // each vertice count of continuous_points is 2 * continuous_points.size() - 2
            std::size_t vertex_count = std::reduce(sampled_points.cbegin(), sampled_points.cend(), 0, [](std::size_t counter, const auto& vec) { return counter + 2 * vec.size() - 2; });
            vertices.reserve(vertices.size() + vertex_count);

            generate_segments<VertexType>(sampled_points, std::back_inserter(vertices), std::forward<std::function<VertexType(point_type)>>(adapter));
        }

        /**
         * @brief Generate line segments from sampled points.
         * 
         * @param sampled_points sampled points
         * @param vertices the output vector to store generated segments
         * @param adapter the transformation function (point_type -> VertexType)
         * 
         * @note This function is faster than generate_segments(const vector<vector<point_type>>, OutputIt) because it reserves the output vertice count. This function calls an overloaded generated_segments function with adapter = std::identity { }.
         */
        template <typename VertexType = point_type>
        static void generate_segments(const std::vector<std::vector<point_type>>& sampled_points, std::vector<VertexType>& vertices) noexcept{
            generate_segments(sampled_points, vertices, std::identity { });
        }
    };
};