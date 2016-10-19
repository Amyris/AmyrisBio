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
    let check a b = Assert.AreEqual(a, format60 b)
    check "" [||]
    check "A" [|'A'|]
    check "AT" [|'A'; 'T'|]

    let aArray n = Array.create n 'A'

    let A59Arr = aArray 59
    let A59 = arr2seq A59Arr

    let A60Arr = aArray 60
    let A60 = arr2seq A60Arr

    let A61Arr = aArray 61
    let A61Result = A60 + "\nA"
    let A120Result = A60 + "\n" + A60
    let A121Result = A120Result + "\nA"

    check A59 A59Arr
    check A60 A60Arr
    check A61Result A61Arr
    check A120Result (aArray 120)
    check A121Result (aArray 121)