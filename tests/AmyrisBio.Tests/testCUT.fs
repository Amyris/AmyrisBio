module testCUT

open NUnit.Framework

open Amyris.Bio.IO

let sampleCodonTable = "UUU 22.8( 29999)  UCU 20.7( 27211)  UAU 12.3( 16236)  UGU  5.3(  7030)
UUC 19.3( 25484)  UCC 14.8( 19542)  UAC 21.5( 28314)  UGC  6.9(  9073)
UUA  5.5(  7251)  UCA  9.7( 12776)  UAA  0.0(     0)  UGA  0.0(     0)
UUG 19.4( 25535)  UCG 12.8( 16886)  UAG  0.0(     0)  UGG 12.9( 16999)

CUU 15.9( 20941)  CCU 18.2( 24001)  CAU 10.0( 13195)  CGU  8.9( 11715)
CUC 15.4( 20318)  CCC 12.4( 16275)  CAC 11.8( 15584)  CGC  8.5( 11159)
CUA  8.7( 11435)  CCA 10.3( 13565)  CAA 11.9( 15684)  CGA 19.3( 25402)
CUG 30.3( 39887)  CCG  4.0(  5324)  CAG 27.4( 36122)  CGG  5.9(  7802)

AUU 34.7( 45685)  ACU 21.2( 27947)  AAU 15.6( 20529)  AGU  7.3(  9606)
AUC 19.0( 25040)  ACC 18.5( 24424)  AAC 24.6( 32356)  AGC  8.7( 11446)
AUA  2.8(  3751)  ACA  8.2( 10811)  AAA 12.0( 15869)  AGA  7.6( 10073)
AUG 21.7( 28569)  ACG  5.7(  7576)  AAG 47.2( 62206)  AGG  3.5(  4577)

GUU 17.3( 22761)  GCU 27.0( 35509)  GAU 26.3( 34683)  GGU 17.6( 23162)
GUC 16.4( 21577)  GCC 24.8( 32729)  GAC 33.8( 44500)  GGC 14.8( 19436)
GUA 12.7( 16681)  GCA 16.9( 22206)  GAA 25.9( 34101)  GGA 26.9( 35503)
GUG 25.7( 33834)  GCG  8.2( 10850)  GAG 39.4( 51875)  GGG  5.3(  6953)
"

(* Leucine encoding codons
CTG	30.3	0.318277311
TTG	19.4	0.203781513
CTT	15.9	0.167016807
CTC	15.4	0.161764706
CTA	8.7	    0.091386555
TTA	5.5	    0.057773109


Total: 95.2

*)


type TestCUT() = class
    do
        ()

    let load() = Amyris.Bio.IO.CodonUsage.loadCodonTableFromString sampleCodonTable

    let AssertSameEpsilon (a:double) (b:double) =
        if abs (a-b) > 0.00001 then 
            Assert.Fail(sprintf "comparison of %f and %f differs by %f" a b (abs (a-b)))

    [<Test>]
    member x.LoadFromString() =
        let p = load()
        Assert.AreEqual(p.Count,64)
        Assert.AreEqual(p.["GGT"],17.6)
    
    [<Test>]
    member x.PrepCutSimpleFreqLoaded() =
        let p = load()
        let p2 = Amyris.Bio.IO.CodonUsage.prepCUT 0.0 100 p // No filtering
        Assert.AreEqual(p2.Codon("GGT").Value.freq,17.6)

    [<Test>]
    member x.ThreeStopCodons() =
        let p = load()
        let p2 = Amyris.Bio.IO.CodonUsage.prepCUT 0.0 100 p// No filtering
        Assert.AreEqual(p2.AA('*').Value.Length,3) // Should be three stop codons

    [<Test>]
    member x.TestRel1() =
        let p = load()
        let p2 = Amyris.Bio.IO.CodonUsage.prepCUT 0.0 100 p// No filtering
        AssertSameEpsilon (p2.Codon("CTG").Value.relFreq1) 0.31827731 // CTG should be approx 31% of codons encoding leucine

    [<Test>]
    member x.TestRel2() =
        let p = load()
        let p2 = Amyris.Bio.IO.CodonUsage.prepCUT 0.2 100 p // Only codons >=20%
        Assert.AreEqual(p2.AA('L').Value.Length,2) // Only two codons >=20% of usage
        AssertSameEpsilon (p2.Codon("CTG").Value.relFreq1) 0.31827731 // CTG should still be approx 31% of codons encoding leucine
        AssertSameEpsilon (p2.Codon("CTG").Value.relFreq2) 0.609657948 // but CTG should be approx 61% of codons remaining encoding leucine

    [<Test>]
    member x.TestRank() =
        let p = load()
        let p2 = Amyris.Bio.IO.CodonUsage.prepCUT 0.2 100 p // Only codons >=20%
        Assert.AreEqual(p2.Codon("TTG").Value.rank,2) // TTG is the second ranked 
        Assert.AreEqual(p2.Codon("CTG").Value.rank,1) // TTG is the second ranked 
end

