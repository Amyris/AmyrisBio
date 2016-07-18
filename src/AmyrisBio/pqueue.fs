/// Priority queue implementation

namespace Amyris

/// PriorityQueue implementation useful for event driven simulation
module pqueue =
    open System

    /// Node for a PQ priority queue
    type PQNode<'T,'T2> =
    | Leaf of 'T*'T2
    | Node of PQNode<'T,'T2> * 'T * 'T2 * PQNode<'T,'T2>
    | Empty

    /// Priority queue with Key/Value types.
    type PQGeneric<'T,'T2 when 'T : comparison (* and 'T2 : equality  *) > (lt_op,gt_op,gteq_op) = class
        let mutable tree = Empty
        let mutable count = 0
        let peek tree =
            match tree with
            | Empty -> failwith "Trying to peek empty node"
            | Node(_,k,_,_) -> k
            | Leaf(k,_) -> k

        member x.Count with get() = count

        /// Insert Key/Value pair into priority queue
        member x.Insert(k:'T,v:'T2) =
            let rec insert  t k v = 
                match t with
                | Empty -> Leaf(k,v)
                | Leaf(k1,v1) as l ->
                    if lt_op k1 k then Node(l,k,v,Empty)
                    else Node(Empty,k1,v1,Leaf(k,v))
                | Node(left,kMid,vMid,Empty) ->
                    if gteq_op kMid k then Node(left,kMid,vMid,Leaf(k,v))
                    else Node(left,k,v,Leaf(kMid,vMid))
                | Node(Empty,kMid,vMid,right) ->
                    if gteq_op kMid k then Node(Leaf(k,v),kMid,vMid,right)
                    else Node(Leaf(kMid,vMid),k,v,right)
                | Node(left,kMid,vMid,right) when lt_op kMid k -> // New value higher than this node
                    if lt_op (peek right) (peek left) then Node(left,k,v,insert right kMid vMid)
                    else Node(insert left kMid vMid,k,v,right)
                | Node(left,kMid,vMid,right) -> // New value lower than this node, insert into one subtree
                    if lt_op (peek right) (peek left) then Node(left,kMid,vMid,insert right k v)
                    else Node(insert left k v,kMid,vMid,right)

            tree <- insert tree k v
            count <- count + 1
        /// Ascii representation of tree contents
        member x.Print() =
            let pad n = "                                                                                                               ".Substring(0,n)
            let rec print d tree =
                match tree with
                | Empty -> printf "%sEmpty\n" (pad d)
                | Leaf(k,_) -> printf "%sLeaf(%A)\n" (pad d) k
                | Node(l,k,_,r) -> 
                    printf "%sNode(\n" (pad d)
                    print (d+4)  l
                    printf "%s,%A,\n" (pad d) k
                    print (d+4) r
                    printf "%s)\n" (pad d)
            print 0 tree
            printf "\n"

        /// Pull highest Key/Value pair from queue
        member x.Pop() : ('T*'T2) =
            let rec promote tree =
                match tree with 
                | Empty -> failwith "ERROR: pulling from empty PQ"
                | Leaf(k,v) -> k,v,Empty
                | Node(Empty,k,v,r) -> k,v,r
                | Node(l,k,v,Empty) -> k,v,l
                | Node(l,k,v,r) ->
                    if gt_op (peek r)  (peek l) then 
                        let k',v',r' = promote r
                        k,v,Node(l,k',v',r')
                    else
                        let k',v',l' = promote l
                        k,v,Node(l',k',v',r)
                        
            let k,v,newTree = promote tree
            tree <- newTree
            count <- count - 1
            k,v
        /// Is priority queue empty?
        member x.Empty with get() = match tree with | Empty -> true | _ -> false 

        /// Look at highest value if it exists. Returns Some/None
        member x.Peek with get() =
                            match tree with
                            | Empty -> None
                            | Leaf(k,v) -> Some(k,v)
                            | Node(_,k,v,_) -> Some(k,v)

        /// Look at the top N entries
        member x.PeekN(n) =
            let rec expand seekingN front result =
                let proc x = x |> List.rev |> Seq.take (min n x.Length) |> List.ofSeq
                    
                if seekingN <=0 then proc result
                else
                    match front with
                    | [] -> proc result
                    | _ ->
                        // Find largest value across current front
                        let maxValue =
                            front
                            |> List.fold
                                (fun mx v -> 
                                    match mx,v with
                                    | None,Empty -> None
                                    | None,Leaf(k,_) -> Some(k)
                                    | None,Node(_,k,_,_) -> Some(k)
                                    | Some(mv),Empty -> Some(mv)
                                    | Some(mxv),Leaf(k,_) -> if gt_op k mxv then Some(k) else Some(mxv)
                                    | Some(mxv),Node(_,k,_,_) -> if gt_op k mxv then Some(k) else Some(mxv))
                                None

                        match maxValue with
                        | None -> result // No actual values in remaining front
                        | Some(mv) ->
                            let values,newFront =
                                front
                                |> List.fold
                                    (fun (values,front) x -> 
                                        match x with
                                        | Empty -> (values,front)
                                        | Leaf(k,v) when gteq_op k mv -> (k,v)::values,front
                                        | Leaf(_,_) -> values,x::front
                                        | Node(l,k,v,r) when gteq_op k mv -> (k,v)::values,l::r::front
                                        | Node(_,_,_,_) -> values,x::front)
                                    ([],[])

                            expand (seekingN-values.Length) newFront (values@result)      
            expand n [tree] []                                   
    end

    type PQ<'T,'T2 when 'T : comparison (* and 'T2 : equality  *) >() = class
        inherit PQGeneric<'T,'T2>( (<),(>),(>=) )
    end

    type PQMin<'T,'T2 when 'T : comparison (* and 'T2 : equality  *) >() = class
        inherit PQGeneric<'T,'T2>( (>),(<),(<=) )
    end