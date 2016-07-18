module testSuffixTree
open NUnit.Framework

open System.IO
open Amyris.SuffixTree
open Amyris.utils
open System.Diagnostics
let rng = new System.Random()

[<TestFixture>]
type TestNCk() = class 
    
    do
        ()
    [<Test>]
    member x.Test() =
        let sw = Stopwatch()
        sw.Start()

        let randBase () = match rng.Next() % 4 with
                            | 0 -> 'G'
                            | 1 -> 'A'
                            | 2 -> 'C'
                            | 3 -> 'T'
                            | _ -> failwith "inconceivable"

        let master = Array.init 10000 (fun _ -> randBase())

        let oneSample() =
                let f = rng.Next() % (master.Length-1)
                let t = f+rng.Next() % (master.Length-f-1)
                master.[f..t]


        let reads = seq {
                        for i in {0..200} do
                            yield oneSample()
                        } |> Array.ofSeq

        // Make a tree of all reads separated by | symbols and all reads reversed for good measure
        let allReads = reads |> Array.fold (fun (accum:ResizeArray<char>) b -> 
                                                                            accum.Add('|') 
                                                                            accum.AddRange(b)
                                                                            accum
                                                                            ) (new ResizeArray<char>()) |> Array.ofSeq
                                                                        
        printfn "Building STree with %d chars" allReads.Length
        let tree = SuffixTree(arr2seq allReads)
        tree.buildFwdChain() // Enable counting of hits
        // Assert.IsTrue( 
        printfn "Done!  %f\n" sw.Elapsed.TotalMilliseconds

end

