#pragma once

namespace function_sampler::concepts{
    template <typename T, typename NumberType>
    concept Point = requires(T t){
        { t.x } -> std::convertible_to<NumberType>;
        { t.y } -> std::convertible_to<NumberType>;
    };
};