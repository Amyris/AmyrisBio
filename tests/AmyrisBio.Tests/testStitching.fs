module testStitching
open System
open NUnit.Framework
open Amyris.Dna
open Amyris.Bio
open Amyris.Bio.Stitching
open Amyris.ErrorHandling

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


let assertStitchedTo (s:Result<Dna option,_>) correct =
    match s with
    | Ok(Some(stitchedSeq), msgs) ->
        Assert.AreEqual(stitchedSeq.str, correct)
        Assert.IsEmpty(msgs)
    | x -> Assert.Fail(sprintf "%A" x)

let assertUnstitchable s =
    match s with
    | Ok(None, msgs) ->
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


[<Test>]
let testMismatchedTails() =

    let check a b correctA correctB =
        let aTail, bTail = mismatchedTails (Dna(a)) (Dna(b))
        Assert.AreEqual(Dna(correctA), aTail)
        Assert.AreEqual(Dna(correctB), bTail)

    check "atcggctatgcagtg" "atcggctCatgcagtg" "atgcagtg" "Catgcagtg"
    check "" "a" "" "a"
    check "aa" "a" "a" ""
    check "" "" "" ""
    check "atcg" "atcg" "" ""

[<Test>]
let testEqualWithIndel() =
    let check a b correct =
        // check both argument orderings for paranoia
        Assert.AreEqual(correct, equalWithIndel (Dna(a)) (Dna(b)))
        Assert.AreEqual(correct, equalWithIndel (Dna(b)) (Dna(a)))

    check "atcgcgatcgtacgacaagta" "atcgcgatcgtacgacaagta" Equal // identical seqs
    check "atcgcgatcAgtacgacaagta" "atcgcgatcgtacgacaagta" EqualWithIndel // easy case, one indel in the middle
    check "atcg" "atAcg" EqualWithIndel
    check "atcg" "atcgA" NotEqual // reject single appended BP
    check "Catcg" "atcg" NotEqual // reject single prepended BP
    check "AAA" "AAA" Equal
    check "A" "A" Equal
    check "ATC" "AGC" NotEqual // reject single BP mutations
    check "atcggcatcagc" "aCcggcaTtcagc" NotEqual // reject if mutation present with indel

let defaultReq s = {s = s; searchParams = defaultLoopoutSearchParams}

let assertSingleLoopout s res =
    match computeLoopoutScar (defaultReq s) with
    | Ok(Some({scar0 = ls; scar1=None}), msgs) ->
        Assert.AreEqual(res, ls)
        Assert.IsEmpty(msgs)
    | x ->
        Assert.Fail(sprintf "Wrong result: %A" x)

let assertMultiLoopout s (res: Dna * Dna) =
    match computeLoopoutScar (defaultReq s) with
    | Ok(Some {scar0 = ls0; scar1 = Some(ls1)}, msgs) ->
        Assert.AreEqual(res, (ls0, ls1))
        Assert.IsEmpty(msgs)
    | x ->
        Assert.Fail(sprintf "Wrong result: %A" x)

let assertNoLoopout s =
    match computeLoopoutScar (defaultReq s) with
    | Ok(None, msgs) ->
        Assert.IsEmpty(msgs)
    | x ->
        Assert.Fail(sprintf "Wrong result: %A" x)

let assertLoopoutError s =
    match computeLoopoutScar (defaultReq s) with
    | Bad([msg]) ->
        Assert.That(msg.Contains("Multiple possible loopout matches found:"))
    | x ->
        Assert.Fail(sprintf "Wrong result: %A" x)

[<Test>]
let testLoopoutSimple() =

    let repeat60BP = "AGCACCCTCCACAAGGTCAAGTGGTATCCTGGTAAGGTAAGCTCGTACCGTGATTCATGC"
    let loopedOutRegion = "GACAGGGGTAAGACCATCAGTAGTAGGGATAGTGCCAAACCTCACTCACCACTGCCAATAAGGGGTCCTTACCTGAAGAATAAGTGTCAGCCAGTGTAAC"

    let s = Dna(repeat60BP + loopedOutRegion + repeat60BP)

    assertSingleLoopout s (Dna(repeat60BP))

    let repeat60BPAddIndel = "AGCACCCTCCACAAGGTCAAGTGGTATCaCTGGTAAGGTAAGCTCGTACCGTGATTCATGC"

    let sWithIndel = Dna(repeat60BP + loopedOutRegion + repeat60BPAddIndel)

    assertMultiLoopout sWithIndel (Dna(repeat60BP), Dna(repeat60BPAddIndel))

    let randomDna200BP = "CCGATGAGGAACCCAAAAGGCGAACCGGGCCAGACAACCCGGCGGTATCGCACTCAAAGCCGGGACACGACGCGTCACAGCCGGTAAGAGTAACCCCGGAGTGAAGACCTATGGGGCTGGATAAAACTGCCGTGGTAACCGCCTTCAACAACCCGAATACGTGGCACTTCAGGAGGCGCCCGGAGGGGGGATGTTTTCTA"

    assertNoLoopout (Dna(randomDna200BP))