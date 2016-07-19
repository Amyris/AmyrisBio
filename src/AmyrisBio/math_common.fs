namespace Amyris.Bio 
/// math functions that are ubiquitous throughout our other math libraries
module math_common = 

    open Microsoft.FSharp.Math

    // some constants out to 200+ units variable precision arithmetic 
    let E = 2.71828182845904523536028747135266249775724709369995 // System.Math.E higher vpa
    let PI = 3.1415926535897932384626433832795028841971693993751058 // System.Math.PI higher vpa
    let Epsilon = System.Double.Epsilon // generic epsilon that comes with fsharp

    /// mean value for an array of arbitrary numerical type
    let inline mean (x:'T[]) = (x |> Array.reduce (+)) |> float |> (fun sum -> sum/float(x.Length))

    /// sum squared of an array of arbitrary numerical type
    let inline ssq x = Array.map (fun x->x*x) x |> mean

    /// variance of an array of arbitrary numerical type
    let inline var x = (mean x) |> (fun y -> y*y) |> (-) (ssq x)

    /// standard deviation of an array of arbitrary numerical type
    let inline sd x = sqrt(var x)

    /// median of an array of arbitrary (numerical) type
    let inline median (x:'T[]) : float =
        let xl = Array.length x
        let eo = xl%2 // even/odd
        let x2 = Array.copy x |> Array.sortWith (fun a b -> if a<b then -1 elif a > b then 1 else 0)
        if eo=1 then x2.[(xl-1)/2] |> float // middle number if odd length
        else ( float(x2.[xl/2]) + float(x2.[(xl/2)+1]) )/2.0 // average of middle two if even length

    /// factoring integer k; k should be an integer for this function to make sense - else use gamma function
    /// generic form is provided because when fact (k:float) will produce pretty good results for a larger range of k than (k:int)
    let inline fact k =
        let one = LanguagePrimitives.GenericOne<_>
        let rec _fact (k:'T) (t:'T) = 
            if k <= one then t else _fact (k-one) (t*k)
        _fact k one

    /// highest common factor of two integers or floats a and b
    let inline hcf a b = 
        let rec hcff a b = 
            if a=LanguagePrimitives.GenericZero<_> then b 
            elif a<b then hcff a (b-a)
            else hcff (a-b) b 
        hcff a b