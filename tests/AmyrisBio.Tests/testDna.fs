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