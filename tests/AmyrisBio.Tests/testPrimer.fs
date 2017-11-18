module TestPrimer

open Amyris.Bio.primercore
open NUnit.Framework

[<TestFixture>]
type TestPrimer() = class     
    let common = [
         "AAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "AAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "AAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "AAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "AAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "AAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "AAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "AAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "AAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "AAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "AATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "ATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "TCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "CTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "TTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "TAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "AATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "ATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "TAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "AGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATACATATATAAAG" ;
         "GATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATACATATATAAAG" ;
         "ATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATACATATATAAAG" ;
         "TTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATACATATATAAAG" ;
         "TAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATACATATATAAAG" ;
         "AATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATACATATATAAAG" ;
         "ATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATACATATATAAAG" ;
         "TTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATACATATATAAAG" ;
         "TTAAACAGTATATGTACAGTTTTATATATATATATATATATATATACATATATAAAG" ;
         "TAAACAGTATATGTACAGTTTTATATATATATATATATATATATACATATATAAAG" ;
         "AAACAGTATATGTACAGTTTTATATATATATATATATATATATACATATATAAAG" ;
         "AAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATACATATATAAAG" ;
        ]

    let templates1 = 
        [ "GAAATCTGTACCAACCGTATAGGTGAAAGAGACCCTTTATGG" ;
          // "ATGCCATTCGTTGTTCCTAGAAGAAACCGTTCTTTGT" ;  // interesting case, but dimer check overrides GC clamp so will fail
          "ATGGTCATTGCTGAAGTTCCTAAATTAGCCTCTGCC" ;
          "AAAAAAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATAATATATATATATTATATATATATA" ;
          "AAAAAAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATAATATATATATATTATATATATATACATATATAAAG" ;
        ]@common

    let templates2 = 
        [ 
          "AAAAAAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATAATATATATATATTATATATATATA" ;
          "AAAAAAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATAATATATATATATTATATATATATACATATATAAAG" ;
          "AAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATTATATATATATATATATATAATATATATATATATATATATATATATACATATATAAAG" ;
          "AAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATTATATATATATATATAATATATATATATATATATATATATATATACATATATAAAG" ;
          "AAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATACATATATAAAG" ;
          "AAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATACATATATAAAG" ;
          "AAAAAATCTTAATAGATTAATTTAAACAGTATATGTACAGTTTTATATATATATATATATATATATATATATATATATATATATATATATATATATATACATATATAAAG" ;
          "CGGTTGGGCTTAACTTTAAAGAAAAAAGTTGAGATTAGATTTATTGTGTT"  
        ]@common

    let makeTask (template:string) = 
        {tag = "PR";
         temp = (template.ToCharArray())
         align = LEFT;
         strand = TOP;
         offset = 0;
         targetTemp = 60.0<C>;
         sequencePenalties = None}

    let testPenalty = {
        tmPenalty = 1.0; 
        tmMaxDifference = 5.0<C> ;
        positionPenalty = 5.0 ; 
        lengthPenalty = 3.0 ; 
        polyLengthThreshold = 4; 
        polyPenalty = 10. ; 
        threePrimeUnstablePenalty = 5.0 ; 
        ATPenalty=3.0<C> ; 
        targetLength=20 ;
        maxLength = 80 ;
        minLength = 20 ; 
        monovalentConc = mM2M 50.0<mM>;
        primerConc = uM2M 0.25<uM> ; 
        divalentConc = mM2M 1.5<mM> ; 
        templateConc = uM2M 0.01<uM> ; 
        dNTPConc = uM2M 0.0<uM> ;}

    do
        ()
    [<Test>]
    /// Test that with just 20 bp of template, the design can walk right up to the edge
    member __.Test20Min() =
        let pcrParams = {defaultParams with minLength = 20 }
        let taskFwd : OligoTask = { tag = "PF"  ; temp = "CCAGTTACCCCGACCCCGGC".ToCharArray() ; 
                                    align = ANCHOR.LEFT ; strand = TOP;offset = 0 ; targetTemp = 70.0<C>;
                                    sequencePenalties = None}
        let design = oligoDesign false pcrParams taskFwd
        Assert.IsTrue(design.IsSome)
        Assert.IsTrue(design.Value.oligo.Length=20)

    [<Test>]
    /// Test the DNA sequence context around melting temperature region doesn't affect Tm prediction
    member __.TestTempContext() =
        let x = "CGCTCGTCCAACGCCGGCGGACCT".ToCharArray()
        let y = "CGCTCGTCCAACGCCGGCGGACCTA".ToCharArray()
        let z = "CGCTCGTCCAACGCCGGCGGACCTC".ToCharArray()
        let zz = "CGCTCGTCCAACGCCGGCGGACCTG".ToCharArray()
        
        let a = temp Amyris.Bio.primercore.defaultParams x 24
        for s in [ y;z;zz] do
            let b = temp Amyris.Bio.primercore.defaultParams s 24
            if abs(b-a) > 1.0E-5<C> then
                Assert.Fail(sprintf "%s[.23] and \n%s[..23] give Tms of %f and %f respectively" 
                                            (Amyris.Bio.utils.arr2seq x) 
                                            (Amyris.Bio.utils.arr2seq s) (a/1.0<C>) (b/1.0<C>))


    [<Test>]
    /// Test GC rich template can make minimum length primer
    member __.TestGCRich() =
        let template="CTGCCGGCGACGTGGAGCGTCCGATTGTGACGCGCCTGAGCAACCCGGGCACGGTGCTGCGCGAGTCGTGCGACGCCTCACTGCTGGTGCAGGCCATCATCGACGCCATCGTCGACCTGGCCGTGCCCCTGACGGCCGCGTACAACGACGT"
        let p = oligoDesignWithCompromise false testPenalty (makeTask template)
        Assert.IsTrue(p.IsSome)

    [<Test>]
    /// Test that GC clamp finding is behaving sensibly by tempting designer with long AT rich regions
    member __.TestGCClamp() =
        for template in templates1 do
            // set false->true for diagnostics debugging test cases
            match oligoDesign false testPenalty (makeTask template) with
            | None -> 
                Assert.Fail (sprintf "TestGClamp failed to make a primer for %s" template)
            | Some design ->
                match design.oligo.[design.oligo.Length-1] with
                | 'G' | 'C' -> () // good
                | x -> Assert.Fail( // bad
                        sprintf "TestGCClamp yielded oligo ending in '%c' not G or C for template %s oligo=%s length=%d temp=%A" 
                            x 
                            template 
                            (Amyris.Bio.utils.arr2seq design.oligo)
                            design.oligo.Length
                            design.temp
                       )
    [<Test>]
    /// Test that we don't form a primer dimer under tempting circumstances
    member __.TestPrimerDimer() =
        let verbose = false
        for template in templates2 do
            if verbose then printfn "Testing %s for primer dimer avoidance" template
            if verbose then printfn "designing oligo for %s " template
            match oligoDesign false testPenalty (makeTask template) with
            | None -> 
                if verbose then printfn "failed design for %s " template
                Assert.Fail "TestPrimerDimer failed to make a primer"
            | Some design ->
                if verbose then printfn "succeeded design %A %s for %s " design.temp (Amyris.Bio.utils.arr2seq design.oligo) template
                let tail = Amyris.Bio.primercore.longestTailTailOverlap design.oligo design.oligo
                if verbose then printfn "tail is %d for %s" tail template
                Assert.IsTrue(tail<4)

end
