module testHashing

open Amyris.Bio.MerHash
open NUnit.Framework

let rc (s:string) = Amyris.Bio.biolib.revComp (s.ToCharArray()) |> Amyris.Bio.utils.arr2seq
let mer1 = "ATGTCATGCTACCCGTGTAATCTGATCATGTA"
let mer1RC = rc mer1

let template2 = "ATGTCATGCTACCCGTGTAATCTGATCATGTAC"
let query2 =     "TGTCATGCTACCCGTGTAATCTGATCATGTAC"
let query2RC = rc query2


let template3 = "ATGTCATGCTACCCGTGTAATCTGATCATGTACTTGACTAGCTAGCTGACGTATGCTACTGACTGAGCACGTATCTACGTACGGCGCGATCTGATCATCTGACACACACACTAGCTACTATCCCAGTTTTTTG"
let query3 =     "GTATCTACGTACGGCGCGATCTGATCATCTGA"
let query3RC = rc query3


        
[<Test>]
let T1_hashOne() =
    let h = hashTemplate32 32 mer1
    Assert.AreEqual((ONCE,0) ,search32 h mer1)

[<Test>]
let T1_hashOneRC() =
    let h = hashTemplate32 32 mer1
    Assert.AreEqual((ONCE,0),search32 h mer1RC)


[<Test>]
let T2_hashOne() =
    let h = hashTemplate32 32 template2
    Assert.AreEqual((ONCE,0) ,search32 h query2)

[<Test>]
let T2_hashOneRC() =
    let h = hashTemplate32 200 template2
    Assert.AreEqual((ONCE,0),search32 h query2)

[<Test>]
let T3_hashOne() =
    let h = hashTemplate32 200 template3
    Assert.AreEqual((ONCE,0) ,search32 h query3)

[<Test>]
let T3_hashOneRC() =
    let h = hashTemplate32 200 template3
    Assert.AreEqual((ONCE,0),search32 h query3)