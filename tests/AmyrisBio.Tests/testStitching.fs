module testStitching
open System
open NUnit.Framework
open Amyris.Bio


[<Test>]
let testSeqSlice() =
    let x = [| 1; 2; 3; 4 |]

    let checkPair s e =
        Assert.AreEqual(x.[s..e], (utils.arraySliceSeq x s e) |> Array.ofSeq)
        ()

    let pairsToCheck = [(0, 0); (0, 1); (0, 3); (1, 3); (3, 3); (3,1)]

    for (s, e) in pairsToCheck do
        checkPair s e