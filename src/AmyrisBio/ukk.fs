namespace Amyris.Bio
open System.Collections.Generic
/// Ukkonens suffix tree  algorithm
module SuffixTree =

    open System.IO

    [<Struct>]
    type Edge(_firstChar:int,_lastChar:int,_startNode:int,_endNode:int) = 
        member x.firstChar = _firstChar
        member x.lastChar = _lastChar
        member x.startNode = _startNode
        member x.endNode = _endNode 

    [<Struct>]
    type Node(_sl:int) = 
        member x.suffixLink = _sl 

    type Suffix = { origin : int ; firstChar : int ; lastChar : int }


    [<Struct>]
    type NodeTips(_offset:int, _count :sbyte) =
        member x.offset = _offset
        member x.count = _count


    /// A suffix is canonical if there are no more chars to add beyond its origin,  which is represented with lastChar < firstChar
    let suffixIsExplicit (s:Suffix) = s.firstChar > s.lastChar
    let suffixIsImplicit (s:Suffix) = s.firstChar <= s.lastChar


    type SearchResult = { found : bool ; edge : Edge ; qMatched : int ; eMatched : int}
    //type EdgeKeyOld = { n : int ; c : char}


    type private CountFlag =  RED | YELLOW | BLACK

    [<Struct>] 
    type EdgeKey(_n:int,_c:char) =  
        member v.n = _n
        member v.c = _c
        member x.Compare(e2:EdgeKey) = compare _n e2.n 

    /// Create a suffix tree for searching strings rapidly
    type SuffixTree(text:string) = class
        let T : char [] ref = ref (text.ToCharArray())
        let nodes = new ResizeArray<Node>([ Node(-1)])
        let edges = new  Dictionary<EdgeKey,Edge>()
        let mutable leafCounts : int array option = None

        // These are optional data structures for forward traversing the tree
        // without a specific query in mind.
        let mutable fwdChainBuilt = false
        let nodeTips : NodeTips [] ref = ref [||]
        let charTips : char [] ref = ref [||]

        let ensureFwdChain() =
            if not fwdChainBuilt then 
                failwith "ERROR: buildFwdChain must be run first before using this function"
        
        // Make a canonical form of a suffix, so that
        // the origin be the closest parent node to the end of the
        // tree    
    
        let findEdge startNode firstChar = (edges).[EdgeKey(startNode,firstChar)]
      
        let Canonize (s:Suffix) =
                if suffixIsExplicit s then s 
                else
                    let rec _canonize (edge:Edge) (s:Suffix) = 
                        let edgeSpan = edge.lastChar - edge.firstChar // # chars in edge
                        if edgeSpan > s.lastChar - s.firstChar then  // If suffix ends in this span, then done
                            s
                        else
                            // Suffix spans this edge, move suffix past edge, move node to end node
                            let s' = {firstChar = s.firstChar + edgeSpan + 1 ; origin = edge.endNode ; lastChar = s.lastChar}
                            // If still characters left, then dive into next edge
                            if s'.firstChar <= s'.lastChar then
                                _canonize (findEdge edge.endNode (!T).[s'.firstChar]) s'
                            else
                                // Otherwise go one more round, but I think we're done right?
                                _canonize edge s'
                        
                    _canonize (findEdge s.origin  ((!T).[s.firstChar]))  s
        let addEdge (e:Edge) = 
            (edges).Add( EdgeKey(e.startNode,(!T).[e.firstChar]),e) 
    
        let splitEdge (e0:Edge) (s:Suffix) =
            // Existing e0
            // O --> ABCD --> O
            //
            //   e1      newNode    e2
            // O --> AB --> O --> CD --> O

            // Where does this suffix want to split the edge?
            let splitPoint = s.lastChar - s.firstChar + 1
        
            // Remove old edge 1
            // edges <- edges.Remove({n=e0.startNode;c=(!T).[e0.firstChar]})
            edges.Remove(EdgeKey(e0.startNode,(!T).[e0.firstChar])) |> ignore
        
            // Add the new node we need
            (nodes).Add( Node(s.origin)) // Add new node
            let newNode = (nodes).Count-1
            // Build new edge
            //let e1 = { firstChar = e0.firstChar ; lastChar = e0.firstChar + s.lastChar - s.firstChar ; startNode = s.origin ; endNode = newNode}
            //edges <- addEdge e1
            addEdge (Edge(e0.firstChar,e0.firstChar + s.lastChar - s.firstChar,s.origin,newNode))

        
            // Build other new edge (replacement for old e0)
            //let e2 = { firstChar = e0.firstChar+ splitPoint; lastChar = e0.lastChar ; startNode = newNode ; endNode = e0.endNode}
            //edges <- addEdge e2
            addEdge (Edge(e0.firstChar+ splitPoint,e0.lastChar,newNode,e0.endNode))

            //e1.endNode // Return the split point
            newNode
        
        // Build the tree
        let addPrefix (active : Suffix) (lastChar:int) =  
            let testChar = (!T).[lastChar]
            //printfn "%d %d" lastChar edges.contents.Count
            let rec _update (active : Suffix) (lastParent:int) =
                /// Cleanup done only after recursion is finished
                let rec finalize (active:Suffix) parent lastParent =
                    if lastParent > 0 then
                        nodes.[lastParent] <- Node(parent)
                
                    Canonize {origin = active.origin ; firstChar =active.firstChar ; lastChar  = active.lastChar+1 }
            
                // Start update    
                let parent = active.origin
                if suffixIsExplicit active then
                    // Does this node have an edge starting with test char
                    if edges.ContainsKey(EdgeKey(active.origin,testChar)) then finalize active parent lastParent 
                    else
                        // Edge not found:  Need to add an edge here with this new character
                        //createEdgeThenUpdate active parent lastParent    
                    
                        nodes.Add(Node(if lastParent > 0 then parent else -1))
                        addEdge (Edge(lastChar,(!T).Length-1,parent,nodes.Count-1))

                        if lastParent > 0 then
                            nodes.[lastParent] <- Node(parent)
                        
                        let lastParent' = parent
                        if active.origin = 0 then  
                            let sC = Canonize {origin = 0 ; firstChar = active.firstChar+1 ; lastChar = active.lastChar} 
                            _update sC lastParent'
                        else
                            let sC =
                                {origin = (nodes).[active.origin].suffixLink;
                                 firstChar = active.firstChar;
                                 lastChar = active.lastChar}
                                |> Canonize
                            _update sC lastParent'        
                    
                else // implicit node
                    // Find the active edge
                    let edge = edges.[EdgeKey(active.origin,(!T).[active.firstChar])]
                
                    let span = active.lastChar - active.firstChar
                    // Is implicit node's next char same as test char
                    if (!T).[edge.firstChar + span + 1] = testChar then finalize active parent lastParent 
                    else
                        // oops, need to split this puppy, then pass in the new midpoint node to update
                        let parent = splitEdge edge active
                        //createEdgeThenUpdate active parent lastParent
                    
                        nodes.Add( Node(if lastParent > 0 then parent else -1))                    
                        //let newEdge = { firstChar = lastChar ; lastChar = (!T).Length-1 ; startNode = parent ; endNode = (nodes).Count-1 }
                        //edges <- addEdge newEdge
                        addEdge (Edge(lastChar,(!T).Length-1,parent,nodes.Count-1))

                        if lastParent > 0 then
                            nodes.[lastParent] <- Node(parent)
                        
                        let lastParent' = parent
                        if active.origin = 0 then  
                            let sC = Canonize {origin = 0; firstChar = active.firstChar+1; lastChar = active.lastChar} 
                            _update sC lastParent'
                        else
                            let sC =
                                {origin = (nodes).[active.origin].suffixLink;
                                 firstChar = active.firstChar;
                                 lastChar = active.lastChar}
                                |> Canonize
                            _update sC lastParent'
                    
            _update active -1
        do 
            // Build tree, feeding in successive prefixes of the target text to addPrefix, passing active point from call to call    
            let activeStart = {origin = 0 ; firstChar = 0; lastChar = -1 }  // The initial active prefix
            Seq.fold (fun (ap,i) _(*c*) -> ( (addPrefix ap i),i+1) ) (activeStart,0) (!T) |> ignore

        /// Build links to enable arbitrary tree traversal without
        /// a particular search sequence
        member this.buildFwdChain() =
            // We want a list of the leading character for edges sorted by node
            // and accompanying node info pointing into the table

            let counts : sbyte array = Array.zeroCreate nodes.Count
            // Calculate how many edges come off each node
            for e in edges do
                let n = e.Key.n
                assert(n>=0 && n <= counts.Length-1)
                counts.[n] <- counts.[n] + 1y

            let starts = Array.init nodes.Count (fun _ -> 0)
            for i in {1..nodes.Count-1} do
                starts.[i] <- starts.[i-1] + int(counts.[i-1])
                                                    
            let seen : int array = Array.zeroCreate nodes.Count
            let chars :char array = Array.zeroCreate nodes.Count

            for e in edges do
                let n = e.Key.n
                chars.[starts.[n] + seen.[n]] <- e.Key.c
                seen.[n] <- seen.[n] + 1
                assert(seen.[n] <= int(counts.[n]))

            
            nodeTips := Array.zip starts counts |> Array.map (fun (n,c) -> NodeTips(n,c))  // starts _nodeTips
            charTips :=  chars // _charTips
            fwdChainBuilt<-true
            ()

        /// Label every node with a count of the number of leaves in the subtree (self included)
        member this.buildLeafCounts() =
            // Iteratively propagate a wave of leaf counts from root upwards till every tree is labeled


            /// Number of outstanding edges for each node
            let edgeCount = Array.init nodes.Count (fun _ -> 0uy)
            let leafCountLocal = Array.init nodes.Count (fun _ -> 0)
            let color = Array.init nodes.Count (fun _ ->  BLACK) // Nodes start out black (inactive)

            for edge in edges do
                let node = edge.Key.n
                if edgeCount.[node] = 254uy then
                    failwithf "ERROR: unsupported - more than 255 character alphabets"
                else
                    edgeCount.[node]<-edgeCount.[node]+1uy

            /// Starting leaves are nodes with no edges leading out of them
            for i in {0..edgeCount.Length-1} do
                if edgeCount.[i] = 0uy then
                    leafCountLocal.[i] <- 1 // Self leaf
                    color.[i]<-RED // Flag - leaves are active
                    
            let rec label () =
                if color |> Array.exists (fun v -> v=RED) then
                    // One more round
                    for edge in edges do
                        let child = edge.Value.endNode
                        if color.[child] = RED then
                            let parent = edge.Value.startNode
                            leafCountLocal.[parent] <- leafCountLocal.[parent]+leafCountLocal.[child] // Add child's leaf count to parent
                            edgeCount.[parent] <- edgeCount.[parent]-1uy // One less outstanding child to propagate now
                            if edgeCount.[parent] = 0uy then color.[parent]<-YELLOW // Ready to activate
                    color
                    |> Array.iteri (fun i v ->
                        if v = RED then color.[i] <- BLACK elif v=YELLOW then color.[i] <- RED)
                    label() // Recurse
            label()

            leafCounts<-Some(leafCountLocal)
                    
        /// Given a starting node, calculate the lengths of all paths to tips    
        /// The length is effectively the position of the start of the sequence
        /// measured from the end of the string.
        /// Requires that  buildFwdChain be run once first
        member this.findAllPositions node =
            ensureFwdChain()
            let rec walk todo (results : ResizeArray<int>) =
                match todo with
                | [] -> results //no more work to do
                | (node,depth)::tl ->
                    // How many children does this node have?
                    let tips = (!nodeTips).[node]
                    if tips.count = 0y then
                        // I'm a leaf!  record myself
                        results.Add(depth)
                        walk tl results
                    else // I'm interior - fan out guys..
                        let newTasks =
                            seq {for i in tips.offset..tips.offset+int(tips.count)-1 do
                                    let c = (!charTips).[i]
                                    let e = (edges).[EdgeKey(node,c)]
                                    yield (e.endNode, (depth + e.lastChar - e.firstChar + 1))}
                            |> List.ofSeq
                        walk (newTasks @ tl) results
                                
            let res = new ResizeArray<int>()
            walk [(node,0)] res
        
        /// Simple tree size/string size     
        member this.Count = (nodes).Count
        member this.EdgeCount = (edges).Count // should be node count - 1
        member this.Chars = !T
        member this.LeafCount(n:int) =
            match leafCounts with
            | None -> failwith "ERROR: run buildLeafCount first"
            | Some(lc) -> lc.[n]
    
        /// Find offsets of all instances of str in original string from start of string
        /// Required buildFwdChain be run once first for the tree
        member this.FindAll (str:string) =
            let where = this._search str
            if where.found = false then [||]
            else
                let depths = this.findAllPositions where.edge.endNode
                let stemLen = where.edge.lastChar - where.edge.firstChar + 1 - where.eMatched
                // Convert depths of tree to a position in the array
                depths.ToArray()
                |> Array.map (fun d -> (!T).Length - (d + stemLen + where.qMatched))
                |> Array.sort
    
        member private this._search (str:string) =
            let q = str.ToCharArray()
            let rec _find queryIdx node (lastEdge:Edge) =
                if queryIdx = q.Length then
                    {found = true;
                     edge = lastEdge; 
                     qMatched = q.Length; 
                     eMatched = lastEdge.lastChar - lastEdge.firstChar + 1}
                else
                    // Is there a suitable edge leading in the right direction?
                    match edges.TryGetValue(EdgeKey(node,q.[queryIdx])) with
                    | false,_ -> { found = false ; edge = lastEdge ; qMatched = queryIdx ; 
                                    eMatched = lastEdge.lastChar-lastEdge.firstChar+1} // nope
                    | true,edge ->
                        // Yes - does query terminate in this edge or not?  and if
                        // so does it match
                        let eLen = edge.lastChar - edge.firstChar+1
                        let qLen = q.Length - queryIdx
                        let matchLen = min eLen qLen
                    
                        let s1 = (!T).[edge.firstChar..edge.firstChar + matchLen-1] 
                        let s2 = q.[queryIdx..queryIdx+matchLen-1]

                        let rec matchCount (a:char []) (b:char []) n =
                            if n = a.Length then n
                            else if a.[n] = b.[n] then
                                matchCount a b (n+1)
                            else n

                        let m = matchCount s1 s2 0
                        if m < matchLen then 
                            { found = false ; edge = edge ; qMatched = queryIdx+m ; 
                                    eMatched = m}(*no!*)  
                        else
                            // Did we exhaust query?
                            if matchLen = qLen then
                               {found = true; 
                                edge = edge ; 
                                qMatched = queryIdx+m ; 
                                eMatched = m} (*yes*)
                            else 
                                // No - need to follow edge
                                _find (queryIdx+eLen) edge.endNode edge
            _find 0 0 (Edge(-1,-1,0,0)) //  {firstChar = -1  ; lastChar = -1 ; startNode =0 ; endNode =0 }
        
        member this.CountAll(str:string) = 
            let res = this._search str
            match leafCounts with
                | None -> failwithf "ERROR: need to run buildLeafCount first before using this function"
                | Some(lc) ->
                    if res.found then
                        lc.[res.edge.endNode]
                    else
                        0
        /// Slower implementation of counting# instances that doesn't require pre-running buildLeafCount 
        member this.CountAllSlow(str:string,alphabet:char array) =
            let rec count (n:int) =
                alphabet
                |> Array.map (fun c ->
                    match edges.TryGetValue(EdgeKey(n,c)) with
                    | true,value -> count value.endNode
                    | false,_ -> 0)
                |> Array.sum
                |> max 1 (* if there are no leaf nodes below, we are a leaf which raises count to 1 *)

            let res = this._search str
            if res.found then // Count all leaves below the last node visited
                count res.edge.endNode
            else 0 // No leaves, nothing found

        member this.Contains (str:string) : bool = 
            let res = (this._search str)
            res.found
    
        member this.LongestMatch (str:string) =
            let res = this._search str
            res.qMatched
                   
        member this.dump () =
            for c in (!T) do
                printf "%c" c
            printfn ""
            edges |> Seq.iter (fun kv ->
                let ek = kv.Key
                let edge = kv.Value
                printfn "edge  %d %c -> %d/%d fr=%d -> to=%d "
                    ek.n ek.c edge.firstChar edge.lastChar edge.startNode edge.endNode) 
            nodes |> Seq.iteri (fun i n ->
                printfn "node %d   sLink=%d" i n.suffixLink)
            ()
        /// Write suffix tree in binary format
        member this.save(file:string) =
            use dump = new FileStream(file,FileMode.Create, FileAccess.Write, FileShare.None)  
            use bw = new BinaryWriter(dump)
        
            // Write header
            bw.Write([| 'S' ; 'T' ; 'R' ; 'E' |])

            // Write version
            bw.Write([| '2' ; '0' ; '0' ; '0' |])

            // Write size of reference string
            bw.Write((!T).Length)

            // Write # nodes
            bw.Write(nodes.Count)

            // Write # edges
            bw.Write(edges.Count)

            // Write reference string
            bw.Write(!T)

            
            (* Technically don't need node data later for things like searches but save it anyway *)
            //for n in nodes do
            //    bw.Write(n.suffixLink)

            
            // This is slow ..  Could do better with radix style sort
            let countEdges = Array.init nodes.Count (fun _ -> 0)
            let seenEdges = Array.init nodes.Count (fun _ -> 0)
            let edgeOffset = Array.init nodes.Count (fun _ -> 0)

            // Count # of edges leaving each node
            edges |> Seq.iter (fun e -> countEdges.[e.Value.startNode] <- countEdges.[e.Value.startNode]+1)
            // Work out offset for each node in edge array
            let i,j = countEdges |> Array.fold (fun (i,j) n -> edgeOffset.[i]<- j ; (i+1,j+n) ) (0,0)
            assert(i=nodes.Count)
            assert(j=edges.Count)
            let edgesToSave = Array.init edges.Count (fun _ -> Edge(0,0,0,0) )

            for e in edges do
                edgesToSave.[edgeOffset.[e.Value.startNode]+seenEdges.[e.Value.startNode]] <- e.Value
                seenEdges.[e.Value.startNode] <- seenEdges.[e.Value.startNode] + 1
               
            // Write edges
            for e in edgesToSave do
                bw.Write(e.startNode)
                bw.Write(e.endNode)
                bw.Write(e.firstChar)
                bw.Write(e.lastChar)
            
            // Write edge index
            for eo in edgeOffset do
                bw.Write(eo)

            bw.Close()
            dump.Close()
        /// Load suffix tree from filename
        member private this.loadV1(file:string) =
            //
            use f = new FileStream(file,FileMode.Open, FileAccess.Read, FileShare.None)  
            use br = new BinaryReader(f)
            let header = br.ReadChars(4)
            if header <> [| 'S' ; 'T' ; 'R'; 'E' |] then failwith "file not in suffixtree format"
            // Read text that was indexed
            let numT = br.ReadInt32()
            T := (Array.create numT '?')
        
            assert(br.Read(!T,0,numT) = numT)
        
            // Read number of nodes
            let N = br.ReadInt32()
        
            nodes.Clear()
            // Read actual nodes
            for i in {1..N} do
                nodes.Add( Node(br.ReadInt32()))
            
            //edges <- Map.empty
            edges.Clear()
            let E = br.ReadInt32()
        
            for i in {1..E} do
                let c = br.ReadChar()
                let n = br.ReadInt32()
                let start = br.ReadInt32()
                let endN = br.ReadInt32()
                let first = br.ReadInt32()
                let last = br.ReadInt32()
                edges.Add( EdgeKey(n,c),Edge(first,last,start,endN))
            ()

        /// Load suffix tree from filename
        member private this.loadV2(file:string) =
            //
            use f = new FileStream(file,FileMode.Open, FileAccess.Read, FileShare.None)  
            use br = new BinaryReader(f)
            let header = br.ReadChars(4)
            if header <> [| 'S' ; 'T' ; 'R'; 'E' |] then failwith "file not in suffixtree format"
            let version = br.ReadChars(4)
            assert(version = [| '2' ; '0' ; '0' ; '0' |])

            // Read text size that was indexed
            let numT = br.ReadInt32()
            T := (Array.create numT '?')
            
            // Read number of nodes
            let N = br.ReadInt32()
            nodes.Clear()
            for i in 0..N-1 do
                nodes.Add(Node(0)) // Dummy data - sadly we need to do this because code relies on Nodes.Count to determine nodeCount
            // Read number of edges
            let E = br.ReadInt32()
            edges.Clear()
            printfn "3"
        
            // -----------------------------------------
            // Read actual stuff now
            // ------------------------------------------
            // Read text
            let t = br.Read(!T,0,numT)
            assert( t = numT)
        
            // Read edges

            for i in {1..E} do
                let start = br.ReadInt32()
                let endN = br.ReadInt32()
                let first = br.ReadInt32()
                let last = br.ReadInt32()
                edges.Add( EdgeKey(start,(!T).[first]),Edge(first,last,start,endN))
            // For now ignore the edge index - next E int32s
            
            ()

        member this.load(path:string) =
            use inF = File.OpenRead(path)
            let buffer = Array.init 8 (fun _ -> 0uy)
            inF.Read(buffer,0,8) |> ignore
            let bufferT = buffer |> Array.map (char) |> utils.arr2seq
            if bufferT.StartsWith("STRE") then
                inF.Close()
                match bufferT with
                | "STRE2000" -> this.loadV2(path)
                | _ -> this.loadV1(path)
            else
                failwith "ERROR: bad suffix tree binary file"
    end

    /// Disk based tree access scheme that doesn't have a high startup cost loading tree into memory
    type SuffixTreeDisk(path:string) = class
        let size = FileInfo(path).Length
//        let mapName =
//            sprintf "gslc_st_%d_%s_%s"
//                (System.Diagnostics.Process.GetCurrentProcess().Id)
//                (System.DateTime.Now.ToLongDateString())
//                (System.DateTime.Now.ToLongTimeString())
        
        let openMM path = 
            // TODO - this is not going to share nicely between processes.  We should ideally share the
            // actual memory map between concurrent compiler sessions.  At the moment, the second user
            // will get a locking error ;(
            (*
            Example with shared access
                 return MemoryMappedFile.CreateFromFile(
               //include a readonly shared stream
               File.Open(path, FileMode.Open, FileAccess.Read, FileShare.Read),
               //not mapping to a name
               null,
               //use the file's actual size
               0L, 
               //read only access
               MemoryMappedFileAccess.Read, 
               //not configuring security
               null,
               //adjust as needed
               HandleInheritability.None,
               //close the previously passed in stream when done
               false);
            *)
            //MemoryMappedFiles.MemoryMappedFile.CreateFromFile(path,FileMode.Open,mapName,size)
            let fileStream = File.Open(path,FileMode.Open,FileAccess.Read,FileShare.Read)
            let mapName = null // not mapping to a name
            let capacity = 0L // use the file's actual size
            let accessType = MemoryMappedFiles.MemoryMappedFileAccess.Read // read only access
            let inheritability = HandleInheritability.None // not configuring security
            let leaveOpen = false
            MemoryMappedFiles.MemoryMappedFile.CreateFromFile(
                                                    fileStream,
                                                    mapName,
                                                    capacity,
                                                    accessType,
                                                    inheritability,
                                                    leaveOpen)

        let mm = openMM path

        let va = mm.CreateViewAccessor(0L,0L,MemoryMappedFiles.MemoryMappedFileAccess.Read)
        let headerBytes : byte array = Array.zeroCreate 8
        
        do
            va.ReadArray<byte>(0L,headerBytes,0,8) |> ignore
            if headerBytes <>( [| 'S' ; 'T'; 'R' ; 'E' ; '2' ; '0' ; '0' ; '0' |] |> Array.map (byte) ) then 
                failwithf "ERROR: incorrect header type %A" headerBytes

        let textLen : int32 = va.Read(8L)
        let nodeCount : int32 = va.Read(12L)
        let edgeCount :int32 = va.Read(16L)
        let textStart = 12L (* counts *)  + 8L (* header *)
        let edgeStart = textStart + int64 textLen
        let edgeOffStart = edgeStart + (int64 edgeCount) * 16L
        do
            if edgeOffStart + (int64 nodeCount)*4L <> size then
                failwithf "ERROR: suffixdisktree %s fails size check\nedgeCount=%d\nnodeCount=%d\ntextLen=%d\nedgeStart=%d\nedgeOffStart=%d\nsize=%d" 
                    path edgeCount nodeCount textLen
                    edgeStart edgeOffStart size


        let readChar(i:int) =
            let offset = (int64 i) + textStart
            let c = va.ReadByte(offset) |> char
            c

        let readChars(i:int) (j:int) =
            //let offset = (int64 i) + textStart
            let buffer : byte array = Array.zeroCreate (j-i+1)
            let read = va.ReadArray(textStart+(int64 i),buffer,0,(j-i+1)) 
            assert(read=(j-i+1))
            buffer |> Array.map (char)

        let readEdge(i:int) =
            let offset = (int64 edgeStart) + (int64 i)* 16L
            if offset >= size then
                failwithf "ERROR: readEdge offset for i=%d is %d which is past end of file size=%d" i offset size
            let startNode = va.ReadInt32(offset)
            let endNode = va.ReadInt32(offset+4L)
            let firstChar = va.ReadInt32(offset+8L)
            let endChar = va.ReadInt32(offset+12L)
            Edge(firstChar,endChar,startNode,endNode)

        let readEdgeOffset(i:int) =
            let offset = (int64 edgeOffStart) + (int64 i)*4L
            if offset >= size then
                failwithf "ERROR: readEdge offset for node=%d is %d which is past end of file size=%d" i offset size
            va.ReadInt32(offset)

        let tryFindEdge (n:int) (c:char) =
            let startOffset = readEdgeOffset n
            
            let rec find (i:int) =
                if i = edgeCount then None
                else
                    let e = readEdge i
                    if e.startNode <> n then None
                    else
                        if readChar (e.firstChar ) = c then Some(e) else find (i+1)
            find startOffset

        let search (str : string) =
            let q = str.ToCharArray()
            let rec _find queryIdx node (lastEdge:Edge) =
                if queryIdx = q.Length then
                    {found = true; 
                     edge = lastEdge; 
                     qMatched = q.Length ; 
                     eMatched = lastEdge.lastChar - lastEdge.firstChar + 1} 
                else
                    // Is there a suitable edge leading in the right direction?
                    match tryFindEdge node (q.[queryIdx]) with
                        | None ->
                            {found = false ; 
                             edge = lastEdge ; 
                             qMatched = queryIdx ; 
                             eMatched = lastEdge.lastChar-lastEdge.firstChar+1} // nope
                        | Some(edge) ->
                            // Yes - does query terminate in this edge or not?  and if
                            // so does it match
                            let eLen = edge.lastChar - edge.firstChar+1
                            let qLen = q.Length - queryIdx
                            let matchLen = min eLen qLen
                        
                            let s1 = readChars edge.firstChar (edge.firstChar + matchLen-1)
                            let s2 = q.[queryIdx..queryIdx+matchLen-1]

                            let rec matchCount (a:char []) (b:char []) n =
                                if n = a.Length then n
                                else if a.[n] = b.[n] then
                                    matchCount a b (n+1)
                                else n

                            let m = matchCount s1 s2 0
                            if m < matchLen then 
                                { found = false ; edge = edge ; qMatched = queryIdx+m ; 
                                        eMatched = m}(*no!*)  
                            else
                                // Did we exhaust query?
                                if matchLen = qLen then
                                    {found = true ;
                                     edge = edge ; 
                                     qMatched = queryIdx+m ; 
                                     eMatched = m} (*yes*)
                                else 
                                    // No - need to follow edge
                                    _find (queryIdx+eLen) edge.endNode edge
            _find 0 0 (Edge(-1,-1,0,0)) //  {firstChar = -1  ; lastChar = -1 ; startNode =0 ; endNode =0 }
        
        member this.Contains (str:string) : bool = 
                let res = (search str)
                res.found

        /// Slower implementation of counting# instances that doesn't require pre-running buildLeafCount 
        member this.CountAllSlow(str:string,alphabet:char array) =
            let rec count (n:int) =
                alphabet
                |> Array.map (fun c ->
                    match tryFindEdge n c with
                    | Some(value) -> count value.endNode
                    | None -> 0)
                |> Array.sum
                |> max 1 (* if there are no leaf nodes below, we are a leaf which raises count to 1 *)

            let res = search str
            if res.found then // Count all leaves below the last node visited
                count res.edge.endNode
            else 0 // No leaves, nothing found

        interface System.IDisposable with
            member x.Dispose() =
                va.Dispose()
                mm.Dispose()
    end
