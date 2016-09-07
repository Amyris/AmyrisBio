module testStitching
open System
open NUnit.Framework
open Amyris.Bio
open Amyris.Bio.Stitching

// 100 random DNA base pairs for testing
let randomDna100BPLower =
    "gacgcacactcccttccttgaaaacgcacaatcatacaactgggcacataatgcgtacgcccatctaatacatccaactctctaggtcctgttcaagagc"
let randomDna100BPUpper = randomDna100BPLower.ToUpper()

[<Test>]
let testSeqSlice() =
    let x = [| 1; 2; 3; 4 |]

    let checkPair s e =
        Assert.AreEqual(x.[s..e], (utils.arraySliceSeq x s e) |> Array.ofSeq)
        ()

    let pairsToCheck = [(0, 0); (0, 1); (0, 3); (1, 3); (3, 3); (3,1)]

    for (s, e) in pairsToCheck do
        checkPair s e

[<Test>]
let testDnaType() =
    let d = Dna(randomDna100BPLower)
    // Check to make sure DNA was uppercased
    Assert.AreEqual(randomDna100BPUpper, d.str)

    Assert.Throws<System.Exception>(fun () -> Dna("BADCHARSYO") |> ignore) |> ignore

    // Make sure revcomp is properly cyclically linked
    let rc = d.RevComp()

    Assert.AreSame(d.RevComp(), rc)
    Assert.AreSame(rc.RevComp(), d)

let rcStr (s: string) = s.ToCharArray() |> biolib.revComp |> String
   
// more random DNA snippets
// 30 bp linkers
let linker0 = "TGGAAGAGCACCCTCCACTTGGTCAAGTGA"
let linker1 = "TATCCTCGTAAGGCAAGCTCGTACCGTCAT"
// 75 bp overlap region
let overlapRegion = "TCATGCGGAAGGGGTAAGACCATTAGAAGTAGGGATAGTCCCGAACCTCACTTACCACTCCCAATAAGGGATCCC"
let overlapRegionRC = rcStr overlapRegion
// 80 bp tails
let tail0 = "AATAAATCTGTCGTAGTAACCGGCTTCAACGACCCGTACAGGTGGCACTTCAGGAGGGGCCCGCAGGGAGGAAGTTTTCT"
let tail1 = "GCTATTCGTGGCCGTTCGTGGTAACTAGTTGCGTTCCTAGCCACTACAATTGTTTCTAAGCCGTGTAATGAGAACAACCA"

[<Test>]
let testOverlapStitching() =
    // Test that a reasonable things stitch together
    let seq0 = Dna(tail0 + overlapRegion + linker0)
    let seq1 = Dna(tail1 + overlapRegionRC + linker1)
    let linker0, linker1 = Dna(linker0), Dna(linker1)


    let correctSequence = tail0 + overlapRegion + (rcStr tail1)

    let stitchedSeq = overlapStitchWithLinkers seq0 seq1 linker0 linker1 defaultSearchParameters

    Assert.AreEqual(stitchedSeq.str, correctSequence)

    // Test that stitching still works with no linkers
    let seq0NoLinker = Dna(tail0 + overlapRegion)
    let seq1NoLinker = Dna(tail1 + overlapRegionRC)

    let emptyDna = Dna("")

    let stitchedSeqNoLinkers =
        overlapStitchWithLinkers
            seq0NoLinker seq1NoLinker emptyDna emptyDna defaultSearchParameters
    
    Assert.AreEqual(stitchedSeqNoLinkers.str, correctSequence)

    // Test that a sequence doesn't stitch to itself
    let err = Assert.Throws<System.Exception>(fun () ->
        overlapStitchWithLinkers seq0 seq0 linker0 linker0 defaultSearchParameters
        |> ignore)

    Assert.That(err.Message.Contains("No exact overlap found for input sequences."))