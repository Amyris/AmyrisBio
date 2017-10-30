module TestPrimer

open Amyris.Bio.primercore
open NUnit.Framework

[<TestFixture>]
type TestPrimer() = class     
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
        let template="CTGCCGGCGACGTGGAGCGTCCGATTGTGACGCGCCTGAGCAACCCGGGCACGGTGCTGCGCGAGTCGTGCGACGCCTCACTGCTGGTGCAGGCCATCATCGACGCCATCGTCGACCTGGCCGTGCCCCTGACGGCCGCGTACAACGACGT".ToCharArray()
        let pen={lengthPenalty = 3.0;
             tmPenalty = 1.0;
             tmMaxDifference = 5.0<C>;
             positionPenalty = 5.0;
             polyLengthThreshold = 4;
             polyPenalty = 10.0;
             threePrimeUnstablePenalty = 5.0;
             ATPenalty = 3.0<C>;
             maxLength = 60;
             minLength = 20;
             targetLength = 20;
             monovalentConc = 0.05<M>;
             divalentConc = 0.0015<M>;
             primerConc = 2.5e-07<M>;
             templateConc = 1e-08<M>;
             dNTPConc = 0.0<M>}

        let task = { tag = "PR";
                     temp = template
                     align = CENTERLEFT;
                     strand = TOP;
                     offset = 0;
                     targetTemp = 60.0<C>;
                     sequencePenalties = None}

        let p = oligoDesignWithCompromise false pen task
        Assert.IsTrue(p.IsSome)

    [<Test>]
    /// Test that we don't form a primer dimer under tempting circumstances
    member __.TestPrimerDimer() =
        let template ="CGGTTGGGCTTAACTTTAAAGAAAAAAGTTGAGATTAGATTTATTGTGTT"
        let pen= {
            tmPenalty = 1.0; 
            tmMaxDifference = 5.0<C> ;
            positionPenalty = 5.0 ; 
            lengthPenalty = 3.0 ; 
            polyLengthThreshold = 4; 
            polyPenalty = 10. ; 
            threePrimeUnstablePenalty = 5.0 ; 
            ATPenalty=3.0<C> ; 
            targetLength=20 ;
            maxLength = 60 ;
            minLength = 20 ; 
            monovalentConc = mM2M 50.0<mM>;
            primerConc = uM2M 0.25<uM> ; 
            divalentConc = mM2M 1.5<mM> ; 
            templateConc = uM2M 0.01<uM> ; 
            dNTPConc = uM2M 0.0<uM> ;}

        let task = { tag = "PR";
                     temp = (template.ToCharArray())
                     align = LEFT;
                     strand = TOP;
                     offset = 0;
                     targetTemp = 60.0<C>;
                     sequencePenalties = None}

        match oligoDesign false pen task with
        | None -> Assert.Fail "TestPrimerDimer failed to make a primer"
        | Some design ->
            let tail = Amyris.Bio.primercore.longestTailTailOverlap design.oligo design.oligo
            Assert.IsTrue(tail<4)

end
