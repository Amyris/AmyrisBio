namespace System
open System.Reflection

[<assembly: AssemblyTitleAttribute("AmyrisBio")>]
[<assembly: AssemblyProductAttribute("AmyrisBio")>]
[<assembly: AssemblyDescriptionAttribute("Amyris FSharp BioLib")>]
[<assembly: AssemblyVersionAttribute("2.0.4")>]
[<assembly: AssemblyFileVersionAttribute("2.0.4")>]
do ()

module internal AssemblyVersionInformation =
    let [<Literal>] Version = "2.0.4"
    let [<Literal>] InformationalVersion = "2.0.4"
