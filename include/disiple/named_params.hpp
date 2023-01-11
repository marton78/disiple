#pragma once

namespace disiple {

    // See test/named_params.cpp for usage
    template <typename ListArgs, typename... Params>
    struct Parameters;

    template <typename... Args>
    struct List;

    template <typename T, template <T> class Tag>
    struct Required;

    template <typename T, template <T> class Tag, T Default>
    struct Optional;

    
    // Errors

    template <typename T, template <T> class Tag>
    struct ParameterNotFound;

    template <typename... Args>
    struct UnusedArguments;


    // Implementation of ExtractRequired and ExtractOptional

    // Helper templates to extract parameters from arguments. The two lists are:
    //    * ListArgs: all elements yet to be processed
    //    * ListDone: all elements except for the one that has been found (by default empty)
    template <typename ListArgs, typename T, template <T> class Tag, typename ListDone = List<>>
    struct ExtractRequired;
    
    template <typename ListArgs, typename T, template <T> class Tag, T Default, typename ListDone = List<>>
    struct ExtractOptional;

    // If the current element in ListArgs is the Tag we're looking for, we're done
    template <typename T, template <T> class Tag, T Value, typename... Tail, typename... Done>
    struct ExtractRequired<List<Tag<Value>, Tail...>, T, Tag, List<Done...>>
    {
        using This = Tag<Value>;
        using Rest = List<Done..., Tail...>;
    };

    template <typename T, template <T> class Tag, T Value, T Default, typename... Tail, typename... Done>
    struct ExtractOptional<List<Tag<Value>, Tail...>, T, Tag, Default, List<Done...>>
    {
        using This = Tag<Value>;
        using Rest = List<Done..., Tail...>;
    };

    // Otherwise, the current element in ListArgs is not what we're looking for.
    // Remove it from ListArgs, append it to ListDone, and continue with the next element
    template <typename Head, typename... Tail, typename T, template <T> class Tag, typename... Done>
    struct ExtractRequired<List<Head, Tail...>, T, Tag, List<Done...>>
    : ExtractRequired<List<Tail...>, T, Tag, List<Done..., Head>>
    {};

    template <typename Head, typename... Tail, typename T, template <T> class Tag, T Default, typename... Done>
    struct ExtractOptional<List<Head, Tail...>, T, Tag, Default, List<Done...>>
    : ExtractOptional<List<Tail...>, T, Tag, Default, List<Done..., Head>>
    {};

    // If we have processed all elements in ListArgs, we're done
    template <typename T, template <T> class Tag, typename... Done>
    struct ExtractRequired<List<>, T, Tag, List<Done...>>
    : ParameterNotFound<T, Tag>
    {};

    template <typename T, template <T> class Tag, T Default, typename... Done>
    struct ExtractOptional<List<>, T, Tag, Default, List<Done...>>
    {
        using This = Tag<Default>;
        using Rest = List<Done...>;
    };


    // Implementation of Parameters

    // If the current parameter is Optional, extract it and continue processing
    template <typename... Args, typename T, template <T> class Tag, T Default, typename... Tail>
    struct Parameters<List<Args...>, Optional<T, Tag, Default>, Tail...>
    : ExtractOptional<List<Args...>, T, Tag, Default>::This
    , Parameters<typename ExtractOptional<List<Args...>, T, Tag, Default>::Rest, Tail...>
    {};

    // If the current parameter is Required, extract it and continue processing
    template <typename... Args, typename T, template <T> class Tag, typename... Tail>
    struct Parameters<List<Args...>, Required<T, Tag>, Tail...>
    : ExtractRequired<List<Args...>, T, Tag>::This
    , Parameters<typename ExtractRequired<List<Args...>, T, Tag>::Rest, Tail...>
    {};

    // If there are still Args left, but all Params have been processed, throw an error
    template <typename... Args>
    struct Parameters<List<Args...>>
    : UnusedArguments<Args...>
    {};

    // If there are no Args left and all Params have been processed, all good
    template <>
    struct Parameters<List<>>
    {};

}
