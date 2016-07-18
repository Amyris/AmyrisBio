module testPQ

open NUnit.Framework

[<TestFixture>]
type TestOPQ() = class 
    let items = [ 10 ; 2 ; 4; 7 ; 1 ;8 ; 1 ;3 ; 9 ;3;7;9 ]
    let itemsSorted = items |> List.sortWith (compare) |> List.rev
    
    let itemsF = [ 10.2 ; 2.2 ; 4.1; 7.1 ; 1.1 ;8.7 ; 1.1 ;3.9 ; 9.9 ;3.89999;7.2;9.9 ]
    let itemsSortedF = itemsF |> List.sortWith (compare) |> List.rev
    
    let setup() =
         let pq = Amyris.pqueue.PQ<int,int>()
         items |> List.iteri (fun i v -> pq.Insert(v,i))
         pq
    let setupF() =
         let pq = Amyris.pqueue.PQ<float,int>()
         itemsF |> List.iteri (fun i v -> pq.Insert(v,i))
         pq
        
    do
        ()

    [<Test>]
    member x.PushPop() =
        let pq = setup()
        let out = seq {   while not (pq.Empty) do
                                yield pq.Pop()

                        } |> List.ofSeq
        let outKeys = out |> List.map (fst)


        Assert.AreEqual(outKeys,itemsSorted)

    [<Test>]
    /// Grab more items than are available
    member x.PeekTooMuch() =
        let pq = setup()
        let x = pq.PeekN (items.Length+10) |> List.map (fst)
        Assert.AreEqual(x,itemsSorted)


    [<Test>]
    member x.Peek1() =
        let pq = setup()
        let x = pq.PeekN 1 |> List.map (fst)
        Assert.AreEqual(x,[10])

    [<Test>]
    member x.Peek3() =
        let pq = setup()
        let x = pq.PeekN 3 |> List.map (fst)
        Assert.AreEqual(x,[10;9;9])

    [<Test>]
    member x.Peek3F() =
        let pq = setupF()
        let x = pq.PeekN 3 |> List.map (fst)
        Assert.AreEqual(x,[10.2;9.9;9.9])

    [<Test>]
    member x.Peek5() =
        let pq = setup()
        let x = pq.PeekN 5 |> List.map (fst)
        Assert.AreEqual(x,[10;9;9;8;7])

    [<Test>]
    member x.PeekAll() =
        let pq = setup()
        let x = pq.PeekN items.Length |> List.map (fst)
        Assert.AreEqual(x,itemsSorted)


    [<Test>]
    member x.PeekEmpty() =
        let p = Amyris.pqueue.PQ<int,int>()
        let x = p.PeekN 10
        Assert.AreEqual(x,[])
        let x = p.PeekN 1
        Assert.AreEqual(x,[])

    [<Test>]
    member x.Peek100of10000() =
        let rng = new System.Random()
        
        let pq = Amyris.pqueue.PQ<double,int>()
        let data = [| for i in {0..9999} -> rng.NextDouble() |]
        for d in data do
            pq.Insert(d,0)
        let top100 = pq.PeekN 100 |> List.map (fst)
        Array.sortInPlace data
        let dataRev = data |> Array.rev
        Assert.AreEqual(dataRev.[..99],Array.ofList top100)

    [<Test>]
    member x.PeekAllF() =
        let pq = setupF()
        let x = pq.PeekN items.Length |> List.map (fst)
        Assert.AreEqual(x,itemsSortedF)
    [<Test>]
    member x.Simple10() =
        let rng = new System.Random()
        // Simple insert 10, remove 10
        let pq = Amyris.pqueue.PQ<int,int>()
        for i in {0..10} do
            pq.Insert(i,i*i)
        
        pq.Print()
        for i in {10..-1..0} do
            let k,v = pq.Pop()
            Assert.AreEqual(k,i)  


    [<Test>]
    member x.Simple10Min() =
        let rng = new System.Random()
        // Simple insert 10, remove 10
        let pq = Amyris.pqueue.PQMin<int,int>()
        for i in {0..10} do
            pq.Insert(i,i*i)
        
        for i in {0..10} do
            let k,v = pq.Pop()
            Assert.AreEqual(k,i)  

    [<Test>]
    member x.Random100() =
        let rng = new System.Random()
        let pq = Amyris.pqueue.PQ<int,int>()
        let arr = {0..99} |> Array.ofSeq
        for i in {0..98} do
            let r = (rng.Next() % (99-i)) + i + 1
            let a,b = arr.[r],arr.[i]
            arr.[i] <- a
            arr.[r] <- b
         
        for a in arr do
            pq.Insert(a,a*a)
      
        for i in {99..-1..0} do
            let k,v = pq.Pop()
            Assert.AreEqual(k,i)
            Assert.AreEqual(v,i*i)
    
    [<Test>]
    member x.Random100Min() =
        let rng = new System.Random()
        let pq = Amyris.pqueue.PQMin<int,int>()
        let arr = {0..99} |> Array.ofSeq
        for i in {0..98} do
            let r = (rng.Next() % (99-i)) + i + 1
            let a,b = arr.[r],arr.[i]
            arr.[i] <- a
            arr.[r] <- b
         
        for a in arr do
            pq.Insert(a,a*a)
      
        for i in {0..99} do
            let k,v = pq.Pop()
            Assert.AreEqual(k,i)
            Assert.AreEqual(v,i*i)
    [<Test>]
    member x.Mixed100000() =
        let rng = new System.Random()
        let pq = Amyris.pqueue.PQ<float,int>()
        let mutable myCount = 0
        for i in {0..100000} do
            if rng.NextDouble() > 0.25 || pq.Empty then
                pq.Insert(rng.NextDouble(),0)
                myCount <- myCount+1
                Assert.AreEqual(myCount,pq.Count)
            else
                myCount <- myCount-1
                pq.Pop() |> ignore
                Assert.AreEqual(myCount,pq.Count)
        // Drain the tree now
        let rec drainTree n last =
            if pq.Empty then n
            else
                let v,_ = pq.Pop()
                match last with
                    | None -> drainTree (n+1) (Some(v))
                    | Some(v') ->
                        Assert.GreaterOrEqual(v',v)
                        drainTree (n+1) (Some(v))
        
        drainTree 0 None |> ignore
        Assert.AreEqual(0,pq.Count)

end
