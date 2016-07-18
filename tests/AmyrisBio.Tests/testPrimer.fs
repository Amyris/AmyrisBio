module testPrimer

open Amyris.primercore
open NUnit.Framework


[<TestFixture>]
type TestPrimer() = class     
    do
        ()
    [<Test>]
    /// Test that with just 20 bp of template, the design can walk right up to the edge
    member x.Test20Min() =
        let pcrParams = {defaultParams with minLength = 20 }
        let taskFwd : OligoTask = { tag = "PF"  ; temp = "CCAGTTACCCCGACCCCGGC".ToCharArray() ; 
                                    align = ANCHOR.LEFT ; strand = TOP;offset = 0 ; targetTemp = 70.0<C>;
                                    sequencePenalties = None}
        let design = oligoDesign false pcrParams taskFwd
        Assert.IsTrue(design.IsSome)
        Assert.IsTrue(design.Value.oligo.Length=20)

    [<Test>]
    /// Test the DNA sequence context around melting temperature region doesn't affect Tm prediction
    member x.TestTempContext() =
        let x = "CGCTCGTCCAACGCCGGCGGACCT".ToCharArray()
        let y = "CGCTCGTCCAACGCCGGCGGACCTA".ToCharArray()
        let z = "CGCTCGTCCAACGCCGGCGGACCTC".ToCharArray()
        let zz = "CGCTCGTCCAACGCCGGCGGACCTG".ToCharArray()
        
        let a = temp Amyris.primercore.defaultParams x 24
        for s in [ y;z;zz] do
            let b = temp Amyris.primercore.defaultParams s 24
            if abs(b-a) > 1.0E-5<C> then
                Assert.Fail(sprintf "%s[.23] and \n%s[..23] give Tms of %f and %f respectively" 
                                            (Amyris.utils.arr2seq x) 
                                            (Amyris.utils.arr2seq s) (a/1.0<C>) (b/1.0<C>))
end
