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

[<Test>]
let testFormat60() =
    Assert.AreEqual("", format60 [||])
    Assert.AreEqual("A", format60 [|'A'|])
    Assert.AreEqual("AT", format60 [|'A'; 'T'|])

    let A59 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

    let A60 = A59 + "A"
    let A60TrailingNewline = A60 + "\n"
    let A60Arr = [|for c in A60 -> c|]

    let A61Arr = Array.append A60Arr [|'A'|]
    let A61Result = A60 + "\nA"

    Assert.AreEqual(A60TrailingNewline, format60 A60Arr)
    Assert.AreEqual(A61Result, format60 A61Arr)