module testStitching
open System
open NUnit.Framework
open Amyris.Dna
open Amyris.Bio
open Amyris.Bio.Stitching
open Chessie.ErrorHandling

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


let assertStitchedTo s correct =
    match s with
    | Ok(Stitchable(stitchedSeq), msgs) ->
        Assert.AreEqual(stitchedSeq.str, correct)
        Assert.IsEmpty(msgs)
    | x -> Assert.Fail(sprintf "%A" x)

let assertUnstitchable s =
    match s with
    | Ok(Unstitchable, msgs) ->
        Assert.IsEmpty(msgs)
    | x -> Assert.Fail(sprintf "%A" x)

[<Test>]
let testOverlapStitching() =
    // Test that a reasonable things stitch together
    let req0 =
       {seq0 = Dna(tail0 + overlapRegion + linker0);
        margin0 = Dna(linker0);
        seq1 = Dna(tail1 + overlapRegionRC + linker1);
        margin1 = Dna(linker1);
        searchParams = defaultSearchParameters}

    let correctSequence = tail0 + overlapRegion + (rcStr tail1)

    assertStitchedTo (overlapStitchWithMargins req0) correctSequence

    // Test that stitching still works with no linkers
    let emptyDna = Dna("")

    let req1 =
       {seq0 = Dna(tail0 + overlapRegion);
        margin0 = emptyDna;
        seq1 = Dna(tail1 + overlapRegionRC);
        margin1 = emptyDna;
        searchParams = defaultSearchParameters}

    assertStitchedTo (overlapStitchWithMargins req1) correctSequence
   
    // Test that a sequence doesn't stitch to itself
    let reqSelf = {req0 with seq1 = req0.seq0; margin1 = req0.margin0}

    assertUnstitchable (overlapStitchWithMargins reqSelf)