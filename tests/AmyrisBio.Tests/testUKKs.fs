module testUKKs

open Amyris.Bio
open Amyris.Bio.math_common
open Amyris.Bio.math_bhg
open Amyris.Bio.math_spf
open Amyris.Bio.math_stat
open Amyris.Bio.SuffixTree
open NUnit.Framework
open System.IO   
open System

let testString2 = "Marry, well said; very well said. Look you, sir,
Inquire me first what Danskers are in Paris;
And how, and who, what means, and where they keep,
What company, at what expense; and finding
By this encompassment and drift of question
That they do know my son, come you more nearer
Than your particular demands will touch it:
Take you, as 'twere, some distant knowledge of him;
As thus, 'I know his father and his friends,
And in part him: ' do you mark this, Reynaldo?"

let testString1 = "BANANAS"

[<TestFixture>]
type TestUKKBasic() = class     
    do
        ()
    [<Test>]
    member x.Test() =
        let tree = SuffixTree(testString1)
    
        printfn "Dump tree"
        tree.dump()
        printfn "done.. testing"
    
        Assert.IsTrue(tree.Contains("NAN"))
        Assert.IsTrue(tree.Contains("S"))
        Assert.IsTrue(not (tree.Contains("Q")))
        Assert.IsTrue(not (tree.Contains("BANS")))
        Assert.IsTrue(tree.Contains(""))
        Assert.IsTrue(tree.Contains("BANANAS"))
end

[<TestFixture>]
type TestUKKLonger() = class
    [<Test>]
    member x.Test() =
        let test = testString2
        printfn "Build long string"
        let tree = SuffixTree(test)
        printfn "Test all possible substrings present"
        for i in {0..test.Length-1} do
            for j in {i..test.Length-1} do
                let s = test.Substring(i,j-i+1)
                Assert.IsTrue(tree.Contains(s))
            
        printfn "Make some strings that don't appear probably and ensure they aren't found" 
        let r = new System.Random()
    
        let mutable count = 0
        for i in {0..test.Length-1} do
            for j in {i..test.Length-1} do
                let s = test.Substring(i,j-i+1).ToCharArray()
                let a = r.Next(s.Length)
                let b = r.Next(s.Length-1)
                let c = (a + 1 + b) % (s.Length)
                let tmp = s.[a]
                s.[a] <- s.[c]
                s.[c] <- tmp
                let s' = new string(s)
                //printf "%s %s" s' (if testString.Contains s' then "yes" else "no")
                count <- count + 1
                Assert.IsTrue(tree.Contains(s') = test.Contains(s'))
        printfn "done longer string test (%d substrings searched)" count
end

[<TestFixtureAttribute>]
type UkkLocationTest() = class
    let locationTest test queries =
        printfn "Build test string"
    
        let tree = SuffixTree(test)
        printfn "build location indexes"
        tree.buildFwdChain()
    
        // Traverse tree from root
        printfn "Testing traversal from root"
        let a1 = tree.findAllPositions 0 |> Seq.sort |> Array.ofSeq
        let a2 = {1..a1.Length} |> Array.ofSeq
        let z = a1=a2
        Assert.IsTrue(z)
    
    [<Test>]
    member x.test1() =
        locationTest testString1 ["N" ; "NA"; "S"; "BANANAS" ; ""]

    [<Test>]
    member x.test2() =
        locationTest testString2 ["who" ; "what" ; "where" ; "a"]

end

[<TestFixture>]
type UkkQuery() = class
    let testQuery test qy =
        let tree = SuffixTree(test)
        printfn "build location indexes"
        tree.buildFwdChain()
        printf "testing all loc %s\n" qy
        let res2 = tree.FindAll qy
        let rec slowFind (i:int) (s:string) (q:string) res = 
                    if i >= s.Length then res else
                        let w = s.IndexOf(q,i)
                        if w = -1 then res else
                            slowFind (w+1) s q (w::res)
                        
        
        printf "%A\n" res2
        printfn ""
        let res3 =  (slowFind 0 test qy [] |> List.rev) |> Array.ofList
        printf "%A\n" res3
        printfn ""
        let z = res2=res3
        Assert.IsTrue(z)
    let locationTest test queries =
        for query in queries do
                testQuery test query

    [<Test>]
    member x.Test1() =
        locationTest testString1 ["N" ; "NA"; "S"; "BANANAS" ; ""]

    [<Test>]
    member x.Test2() =
        locationTest testString2 ["who" ; "what" ; "where" ; "a"]
        
end

[<TestFixture>]
type  UkkSaveLoad() = class

    [<Test>]
    member x.saveLoadTest() =
        printf "build tree once\n"
        let t = new SuffixTree(testString2)
        t.buildFwdChain()
        let a1 = t.FindAll("you")
        if File.Exists("testTreeSave") then
            File.Delete("testTreeSave")
        printf "Saving tree\n"
        t.save("testTreeSave")

        let t2 = new SuffixTree("")
        printf "Loading tree\n"

        t2.load("testTreeSave")
        printf "Build fwd chain\n"
        t2.buildFwdChain()
        printf "Find again\n"
        let a2 = t2.FindAll("you")
        
        printf "%A\n" a1
        printf "%A\n" a2
        let z = a1=a2
        Assert.IsTrue(z)
        if File.Exists("testTreeSave") then
            File.Delete("testTreeSave")
end


[<TestFixture>]
type  UkkCountTest() = class

    [<Test>]
    member x.countCheck() =
        printf "build tree once\n"
        let t = new SuffixTree(testString2)
        t.buildLeafCounts()

        let words = testString2.Split([| ' ' ; '\n' ; '\r' ; ',' |],StringSplitOptions.RemoveEmptyEntries)
        let slowCount (w:string) =
            let rec count (i:int) n =
                match testString2.IndexOf(w,i) with
                    | -1 -> n
                    | _ as v-> 
                        count (v+1) (n+1)
            count 0 0

        for w in words do
            let slow = slowCount w
            let fast = t.CountAll w
            if slow <> fast then
                Assert.Fail(sprintf "Slow=%d Fast=%d for word %s" slow fast w)

    [<Test>]
    member x.slowFastComparison() =
        printf "build tree once\n"
        let t = new SuffixTree(testString2)
        t.buildLeafCounts()

        let alphabet = testString2.ToCharArray() |> Set.ofSeq |> Array.ofSeq 

        let words = testString2.Split([| ' ' ; '\n' ; '\r' ; ',' |],StringSplitOptions.RemoveEmptyEntries)

        for w in words do
            let slow = t.CountAllSlow(w,alphabet)
            let fast = t.CountAll w
            if slow <> fast then
                Assert.Fail(sprintf "Slow=%d Fast=%d for word %s" slow fast w)
        
end

[<TestFixture>]
type  UkkDiskTree() = class

    let makeAndSave input =
        let t = new SuffixTree(input)
        t.save("testDiskSearch")
        t
        
    let cleanUp() =
        if File.Exists("testDiskSave") then
            File.Delete("testDiskSave")    


    [<Test>]
    member x.ContainsCheckSimple1() =
        try
            let t = makeAndSave testString1 

            use tree= new SuffixTreeDisk("testDiskSearch")
            Assert.IsTrue(tree.Contains("NAN"))
        finally
            cleanUp()
    [<Test>]
    member x.ContainsCheckSimple2() =
        try
            let t = makeAndSave testString1 
            use tree= new SuffixTreeDisk("testDiskSearch")
            Assert.IsTrue(not (tree.Contains("Q")))
        finally
            cleanUp()

    [<Test>]
    member x.ContainsCheck() =
        try
            let t = makeAndSave testString1 

            use tree= new SuffixTreeDisk("testDiskSearch")
            Assert.IsTrue(tree.Contains("NAN"))
            Assert.IsTrue(tree.Contains("S"))
            Assert.IsTrue(not (tree.Contains("Q")))
            Assert.IsTrue(not (tree.Contains("BANS")))
            Assert.IsTrue(tree.Contains(""))
            
        finally
            cleanUp()
                
    [<Test>]
    member x.diskVMemCheck() =
        try
            let t = makeAndSave testString2
        
            let alphabet = testString2.ToCharArray() |> Set.ofSeq |> Array.ofSeq 

            use d = new SuffixTreeDisk("testDiskSearch")

            let words = testString2.Split([| ' ' ; '\n' ; '\r' ; ',' |],StringSplitOptions.RemoveEmptyEntries)
            for w in words do
                let disk = d.CountAllSlow(w,alphabet)
                let fast = t.CountAllSlow(w,alphabet)
                if disk <> fast then
                    Assert.Fail(sprintf "Disk=%d Fast=%d for word %s" disk fast w)
        finally
            cleanUp()


end
