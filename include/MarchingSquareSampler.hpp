#pragma once

#include <bitset>
#include <concepts>

#include "concept/Point.hpp"
#include "concept/Interval.hpp"
#include "Point.hpp"
#include "Interval.hpp"

namespace function_sampler{
    template <std::floating_point NumberType = double, typename PointType = Point, typename IntervalType = Interval> requires
        concepts::Interval<IntervalType, NumberType>
    class MarchingSquareSampler{
    public:
        // type definitions.
        using number_type = NumberType;
        using point_type = PointType;
        using interval_type = IntervalType;
        using function_type = std::function<number_type(number_type, number_type)>;

        struct SampleResult{
            interval_type xbound;
            interval_type ybound;
            unsigned int x_sampling_count;
            unsigned int y_sampling_count;
            std::vector<number_type> z_values;
        };

        // fields.
        unsigned int horizontal_sampling_count;
        unsigned int vertical_sampling_count;

    private:
        /**
         * @brief A number divides a line segment of a given interval in the ratio ( \p index : \p subdivision_count - \p index ).
         * 
         * @param interval the interval which consists the line segment
         * @param subdivision_count the number of interval division (including start point and end point)
         * @param index the index of internally dividing point, this can be 0... \p subdivision_count
         * @return number_type internally dividing point
         */
        static constexpr number_type divide_interval(const interval_type& interval, uint subdivision_count, uint index) noexcept{
            assert(interval.length() > 0 && index <= subdivision_count);
            return (interval.low * (subdivision_count - index) + interval.high * index) / subdivision_count;
        }

        /**
         * @brief Get the corner bits of the given function values
         * 
         * @param tl function value of the top left point of the grid cell
         * @param tr function value of the top right point of the grid cell
         * @param br function value of the bottom right point of the grid cell
         * @param bl function value of the bottom left point of the grid cell
         * @param threshold the cantour value (k for f(x,y) = k)
         * @return std::bitset<4> The Marching Square corner bits
         * 
         * @note Be aware that the order of the parameter is tl-tr-br-bl (clockwise).
         */
        static constexpr std::bitset<4> get_corner_bits(number_type tl, number_type tr, number_type br, number_type bl, number_type threshold) noexcept{
            return { static_cast<unsigned long>(((tl > threshold) << 3) + ((tr > threshold) << 2) + ((br > threshold) << 1) + (bl > threshold)) };
        }

        /**
         * @brief Get the corner bits of the given grid cell
         * 
         * @param sample_result a sample result fetched by \p sample_points function
         * @param row row index of the grid cell (0... \p sampled_points.y_xampling_count - 1)
         * @param column column index of the grid cell (0... \p sampled_points.x_xampling_count - 1)
         * @param threshold the cantour value (k for f(x,y) = k)
         * @return std::bitset<4> The Marching Square corner bits
         * 
         * @note This function calls \p get_corner_bits(number_type,number_type,number_type,number_type,number_type) internally.
         */
        static constexpr std::bitset<4> get_corner_bits(const SampleResult& sample_result, unsigned int row, unsigned int column, number_type threshold) noexcept{
            assert(row < sample_result.x_sampling_count);
            assert(column < sample_result.y_sampling_count);
            return get_corner_bits(
                sample_result.z_values[sample_result.x_sampling_count * row       + column    ], 
                sample_result.z_values[sample_result.x_sampling_count * row       + column + 1], 
                sample_result.z_values[sample_result.x_sampling_count * (row + 1) + column + 1],
                sample_result.z_values[sample_result.x_sampling_count * (row + 1) + column    ], 
                threshold
            );
        }

        template <typename OutputIt>
        static void fetch_isolines(std::bitset<4> corner_bits, const point_type& left_mid, const point_type& top_mid, const point_type& right_mid, const point_type& bottom_mid, OutputIt output) noexcept{
            switch (corner_bits.to_ulong()){
                case 0b1000: case 0b0111:
                    *output = left_mid; *output = top_mid; break;
                case 0b0100: case 0b1011:
                    *output = top_mid; *output = right_mid; break;
                case 0b0010: case 0b1101:
                    *output = right_mid; *output = bottom_mid; break;
                case 0b0001: case 0b1110:
                    *output = left_mid; *output = bottom_mid; break;
                case 0b1100: case 0b0011:
                    *output = left_mid; *output = right_mid; break;
                case 0b0110: case 0b1001:
                    *output = top_mid; *output = bottom_mid; break;
                case 0b0101:
                    *output = top_mid; *output = right_mid; *output = left_mid; *output = bottom_mid; break;
                case 0b1010:
                    *output = top_mid; *output = right_mid; *output = left_mid; *output = bottom_mid; break;
                default: // no isoline for 0b0000, 0b1111
                    break;
            }
        }

    public:
        MarchingSquareSampler(unsigned int horizontal_sampling_count, unsigned int vertical_sampling_count) noexcept : horizontal_sampling_count { horizontal_sampling_count }, vertical_sampling_count { vertical_sampling_count } { }

        SampleResult sample_points(const function_type& function, const interval_type& xdomain, const interval_type& ydomain) const noexcept{
            std::vector<number_type> z_values;
            z_values.reserve(horizontal_sampling_count * vertical_sampling_count);

            #pragma omp declare reduction (merge : std::vector<number_type> : omp_out.insert(omp_out.end(), omp_in.cbegin(), omp_in.cend()))
            #pragma omp parallel for reduction(merge: z_values)
            for (int i = 0; i < vertical_sampling_count; ++i){
                const auto y = divide_interval(ydomain, vertical_sampling_count - 1, i);
                for (int j = 0; j < horizontal_sampling_count; ++j){
                    const auto x = divide_interval(xdomain, horizontal_sampling_count - 1, j);

                    z_values.emplace_back(function(x, y));
                }
            }

            return SampleResult { xdomain, ydomain, horizontal_sampling_count, vertical_sampling_count, std::move(z_values) };
        }

        template <typename VertexType, typename OutputIt, typename AdapterPP> requires
            std::is_convertible_v<AdapterPP, std::function<VertexType(point_type)>> // for perfect forwarding
        static void generate_segments(const SampleResult& sample_result, number_type threshold, OutputIt output, AdapterPP&& adapter) noexcept {
            for (int i = 0; i < sample_result.y_sampling_count - 1; ++i){
                const auto top    = divide_interval(sample_result.ybound, sample_result.y_sampling_count - 1, i    ),
                           bottom = divide_interval(sample_result.ybound, sample_result.y_sampling_count - 1, i + 1),
                           v_mid  = (bottom + top) / 2;

                for (int j = 0; j < sample_result.x_sampling_count - 1; ++j){
                    const auto left  = divide_interval(sample_result.xbound, sample_result.x_sampling_count - 1, j    ),
                               right = divide_interval(sample_result.xbound, sample_result.x_sampling_count - 1, j + 1),
                               h_mid = (left + right) / 2;

                    const auto corner_bits = get_corner_bits(sample_result, i, j, threshold);

                    const point_type top_mid    { h_mid, top    },
                                     bottom_mid { h_mid, bottom },
                                     left_mid   { left , v_mid  },
                                     right_mid  { right, v_mid  };

                    fetch_isolines(corner_bits, left_mid, top_mid, right_mid, bottom_mid, output);
                }
            }
        }

        template <typename OutputIt>
        static void generate_segments(const SampleResult& sample_result, number_type threshold, OutputIt output) noexcept {
            generate_segments<point_type>(sample_result, threshold, output, std::identity { });
        }
    };
};