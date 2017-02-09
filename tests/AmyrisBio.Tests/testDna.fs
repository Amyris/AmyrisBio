module testDna
open System
open NUnit.Framework
open Amyris.Dna

// 100 random DNA base pairs for testing
let randomDna100BPLower =
    "gacgcacactcccttccttgaaaacgcacaatcatacaactgggcacataatgcgtacgcccatctaatacatccaactctctaggtcctgttcaagagc"
let randomDna100BPUpper = randomDna100BPLower.ToUpper()

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

    // Check indexing
    Assert.AreEqual(d.[0], 'G')
    Assert.AreEqual(d.[1], 'A')
    Assert.AreEqual(d.[99], 'C')

    // Check slicing
    Assert.AreEqual(d.[0..2], Dna("GAC"))
    Assert.AreEqual(d.[..2], Dna("GAC"))
    Assert.AreEqual(d.[98..], Dna("GC"))

    // Check splitting
    Assert.AreEqual(
        (Dna("gacgcacactc"), Dna("ccttccttgaaaacgcacaatcatacaactgggcacataatgcgtacgcccatctaatacatccaactctctaggtcctgttcaagagc")),
        d.Split(10))

    // Check that hash codes two two identical DNA objects are equivalent.
    let d1 = Dna(randomDna100BPLower)
    let d2 = Dna(randomDna100BPLower)
    Assert.AreEqual(d1, d2)
    let dnaSet = Set.ofList [d1; d2]
    Assert.AreEqual(1, dnaSet.Count)
    Assert.AreEqual(Operators.hash d1, Operators.hash d2)
    Assert.AreEqual(d1.GetHashCode(), d2.GetHashCode())

[<Test>]
let testConcat() =
    let pieces = [Dna("AAAGGTCCAT"); Dna("AGCACGTACA"); Dna("TCGCAACCTG")]
    let d = DnaOps.concat pieces
    Assert.AreEqual(d, Dna("AAAGGTCCATAGCACGTACATCGCAACCTG"))