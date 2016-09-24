module testSmithWaterman
open System
open NUnit.Framework

let run (a:string) (b:string) =
    let a' = a.Replace("-","")
    let b' = b.Replace("-","")
    let ok,alignA,alignB,_ = Amyris.Bio.smithWaterman._smithWaterman false None 10 20 (a'.ToCharArray()) (b'.ToCharArray())
    if alignA <> a || alignB <> b then
        Assert.Fail(sprintf "mismatch checking alignment:\n%s\n%s\n%s\n%s" a b alignA alignB)
    ()

//    let a = "TGATCTACTGACGTATCTGATGCATCACTGATCTA"    
//    let b = "TGATCTACTGACGTATCTGATGCATCACTGATCTA"
[<Test>]
let T1_simpleSW() =
    let a = "TGATCTACTGACGTAT-TGATGCATCACTGATCTA"    
    let b = "TGATCTACT-ACGTATCTGATGCATCACTGATCTA"
    run a b

[<Test>]
let T2_leadingGaps1() =
    let a = "----------ACGTATCTGATGCATCACTGATCTA"    
    let b = "TGATCTACTGACGTATCTGATGCATCACTGATCTA"
    run a b
let T3_leadingGaps2() =
    let b = "----------ACGTATCTGATGCATCACTGATCTA"    
    let a = "TGATCTACTGACGTATCTGATGCATCACTGATCTA"
    run a b


[<Test>]
let T4_trailingGaps1() =
    let a = "TGATCTACTGACGTATCTGATGC------------"    
    let b = "TGATCTACTGACGTATCTGATGCATCACTGATCTA"
    run a b

[<Test>]
let T5_trailingGaps2() =
    let b = "TGATCTACTGACGTATCTGATGC------------"    
    let a = "TGATCTACTGACGTATCTGATGCATCACTGATCTA"
    run a b

[<Test>]
let T6_interiorGaps1() =
    let a = "TGATCTACTG------------------TGATCTA"    
    let b = "TGATCTACTGACGTATCTGATGCATCACTGATCTA"
    run a b

[<Test>]
let T7_interiorGaps2() =
    let b = "TGATCTACTG------------------TGATCTA"    
    let a = "TGATCTACTGACGTATCTGATGCATCACTGATCTA"
    run a b


(*
testSmithWaterman.test:
mismatch checking alignment:
TGATCTACTG------------------TGATCTA  7 * 4 - 18 = 10
TGATCTACTGACGTATCTGATGCATCACTGATCTA


TGATCTACTGT-G-ATCT-A---------------  6 * 4 - 18 - 1= 5
TGATCTACTGACGTATCTGATGCATCACTGATCTA
*)

[<Test>]
let test() =
    let a = "TGATCTACTG------------------TGATCTA"    
    let b = "TGATCTACTGACGTATCTGATGCATCACTGATCTA"
    let a' = a.Replace("-","")
    let b' = b.Replace("-","")
    let ok,alignA,alignB,_ = Amyris.Bio.smithWaterman._smithWaterman true None 10 40 (a'.ToCharArray()) (b'.ToCharArray())
    if alignA <> a || alignB <> b then
        Assert.Fail(sprintf "mismatch checking alignment:\n%s\n%s\n%s\n%s" a b alignA alignB)
    ()
