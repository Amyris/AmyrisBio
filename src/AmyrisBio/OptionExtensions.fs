/// Extension functions to the Option module to give the Option type the composable functions it
/// sorely needs, specifically around defaults when None.
namespace Amyris.OptionExtensions

module Option =
    /// Flatten a nested option by mapping
    /// Some(Some(a)) -> Some(a)
    /// Some(None) -> None
    /// None -> None
    let flatten =
        function
        | Some a -> a
        | None -> None

    /// Return the wrapped value if it is Some, otherwise return d.
    let getOrElse d = 
        function
        | Some a -> a
        | None -> d

    /// Return a if a.isSome; otherwise, return b,
    let orElse b a =
        match a with
        | Some _ -> a
        | None -> b