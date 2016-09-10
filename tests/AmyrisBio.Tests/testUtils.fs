module testUtils
open System
open NUnit.Framework
open Amyris.Bio.utils

[<Test>]
let testPaddedZip() =
    let smaller = [1;2;];
    let larger = [1;2;3;4;];

    let correctZip =
        [Some(1), Some(1); Some(2), Some(2); None, Some(3); None, Some(4)]

    let check x y correct =
        Assert.AreEqual( (zipWithPad x y) |> List.ofSeq, correct)

    check smaller larger correctZip

    check [] [] []
    check [1] [] [Some(1), None]
    check [] [1] [None, Some(1)]