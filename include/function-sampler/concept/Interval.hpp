#pragma once

namespace function_sampler::concepts{
    template <typename T, typename NumberType>
    concept Interval = requires(T t){
        { t.low } -> std::convertible_to<NumberType>;
        { t.high } -> std::convertible_to<NumberType>;
        { t.length() } -> std::convertible_to<NumberType>;
    };
};