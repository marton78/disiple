#pragma once

namespace disiple {

    // some common named parameters
    template <int N> struct Channels { enum { channels = N }; };
    template <int N> struct Stages   { enum { stages = N }; };
    template <int N> struct Length   { enum { length = N }; };

    // See test/named_params.cpp for usage
    template <typename ListArgs, typename... Params>
    struct Parameters;

    template <typename... Args>
    struct List;

    template <typename T, template <T> class Tag>
    struct RequiredValue;

    template <typename T, template <T> class Tag, T Default>
    struct OptionalValue;

    template <template <typename> class Tag>
    struct RequiredType;

    template <template <typename> class Tag, typename Default>
    struct OptionalType;

    
    // Errors

    template <typename T, template <T> class Tag>
    struct RequiredValueNotFound;

    template <template <typename> class Tag>
    struct RequiredTypeNotFound;

    template <typename... Args>
    struct UnusedArguments;


    // Implementation of Extract

    // Helper template to extract parameters from arguments. The two lists are:
    //    * ListArgs: all elements yet to be processed
    //    * ListDone: all elements except for the one that has been found (by default empty)
    template <typename ListArgs, typename Parameter, typename ListDone = List<>>
    struct Extract;
    
    // If the current element in ListArgs is the Tag we're looking for, we're done
    template <typename T, template <T> class Tag, T Value, typename... Tail, typename... Done>
    struct Extract<List<Tag<Value>, Tail...>, RequiredValue<T, Tag>, List<Done...>>
    {
        using This = Tag<Value>;
        using Rest = List<Done..., Tail...>;
    };

    template <typename T, template <T> class Tag, T Value, T Default, typename... Tail, typename... Done>
    struct Extract<List<Tag<Value>, Tail...>, OptionalValue<T, Tag, Default>, List<Done...>>
    {
        using This = Tag<Value>;
        using Rest = List<Done..., Tail...>;
    };

    template <template <typename> class Tag, typename Type, typename... Tail, typename... Done>
    struct Extract<List<Tag<Type>, Tail...>, RequiredType<Tag>, List<Done...>>
    {
        using This = Tag<Type>;
        using Rest = List<Done..., Tail...>;
    };

    template <template <typename> class Tag, typename Type, typename Default, typename... Tail, typename... Done>
    struct Extract<List<Tag<Type>, Tail...>, OptionalType<Tag, Default>, List<Done...>>
    {
        using This = Tag<Type>;
        using Rest = List<Done..., Tail...>;
    };

    // Otherwise, the current element in ListArgs is not what we're looking for.
    // Remove it from ListArgs, append it to ListDone, and continue with the next element
    template <typename Head, typename... Tail, typename Parameter, typename... Done>
    struct Extract<List<Head, Tail...>, Parameter, List<Done...>>
    : Extract<List<Tail...>, Parameter, List<Done..., Head>>
    {};

    // If we have processed all elements in ListArgs, we're done
    template <typename T, template <T> class Tag, typename... Done>
    struct Extract<List<>, RequiredValue<T, Tag>, List<Done...>>
    : RequiredValueNotFound<T, Tag>
    {};

    template <typename T, template <T> class Tag, T Default, typename... Done>
    struct Extract<List<>, OptionalValue<T, Tag, Default>, List<Done...>>
    {
        using This = Tag<Default>;
        using Rest = List<Done...>;
    };

    template <template <typename> class Tag, typename... Done>
    struct Extract<List<>, RequiredType<Tag>, List<Done...>>
    : RequiredTypeNotFound<Tag>
    {};

    template <template <typename> class Tag, typename Default, typename... Done>
    struct Extract<List<>, OptionalType<Tag, Default>, List<Done...>>
    {
        using This = Tag<Default>;
        using Rest = List<Done...>;
    };



    // Implementation of Parameters


    // Extract the current parameter and continue processing
    template <typename... Args, typename Parameter, typename... Tail>
    struct Parameters<List<Args...>, Parameter, Tail...>
    : Extract<List<Args...>, Parameter>::This
    , Parameters<typename Extract<List<Args...>, Parameter>::Rest, Tail...>
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
