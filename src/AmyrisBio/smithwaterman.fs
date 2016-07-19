namespace Amyris.Bio
/// optimized banded smithwaterman alignment
module smithWaterman =
    open System
    open Amyris.Bio.MerHash
    open System.IO
    open utils
    type SWDirection = UP | LEFT | DIAG |NOWHERE

    [<Struct>]
    type SWCell(_s:int,_direc:SWDirection) = // { s : int ; direc : SWDirection }
        member x.s = _s
        member x.direc = _direc

    type SwScores = 
        { matchScore:int ; deletionScore:int ; mismatchScore:int}    
    let defaultSwScores = {matchScore= 4; deletionScore= -1; mismatchScore= -1} 

    let _smithWaterman debug (swScores:SwScores option) N W (seq1 : char []) (seq2 : char [])  =
        
        let costs = match swScores with 
                        | None -> defaultSwScores
                        | Some(c) -> c
        
        // make lookup table of mers from sequence1 with offsets for each mer
        
        //let time0 = DateTime.Now
        let count = ref 0
    
        let rec enc count i t =
            if count = 0 then t
            else
                enc (count-1) (i+1) (t*4uL + encDNA seq1.[i])
                
        let merTable =
            seq { for i in {0..seq1.Length-N-1} do
                    //let mer = seq1.[i..i+N-1]
                    let n = enc N i 0uL // encBases mer
                    yield n,i}
            |> Map.ofSeq // Possible issue here with duplicate mers which will be preferentially mapped towards end

        let rec enc2 count i t =
            if count = 0 then t
            else
                enc2 (count-1) (i+1) (t*4uL + encDNA seq2.[i])
        // Find mers in seq2 and generate relative offsets
        let relOffset =
            seq2
            |> Array.mapi (fun i _ -> 
                if i < seq2.Length-N then
                    let merI = enc2 N i 0uL //  mer
                    match merTable.TryFind(merI) with
                    | Some(x) -> Some(i-x)
                    | None -> None
                else
                    None)

        // Work out median offset of seq2 from seq1 based on hits
        let hits =
            relOffset
            |> Array.choose id
            |> Array.sort
                        
        if hits.Length = 0 then
            //printfn "No alignment"
            false, (seq1 |> arr2seq), (seq2 |> arr2seq), 0
        else
            let medianOffset = hits.[hits.Length/2]
            // Pad out the sequences so they roughly align
            let s1,s2 = 
                if medianOffset > 0 then
                    let s1 = Array.concat [ Array.create medianOffset ' '  ; seq1 ]
                    (s1 , seq2)
                else
                    let s2 =  Array.concat [ Array.create (-medianOffset) ' ' ; seq2 ]
                    seq1, s2
            let maxLen = max s1.Length s2.Length
        
        
            //let time1 = DateTime.Now    
            // Width of band on matrix

            let s1b = Array.concat [ Array.create (W/2) ' ' ; s1 ; Array.create (maxLen - s1.Length) ' ']
            let s2b = Array.concat [ s2 ; Array.create (maxLen - s2.Length) ' ' ; Array.create (W/2) ' ']
        
        
            let L = s1b.Length

            let cell = Array2D.create L W (SWCell(-999999,NOWHERE))
        
            let del = costs.deletionScore // -1
            let matchCost = costs.matchScore // 4
            let mismatchCost = costs.mismatchScore //-1
        
            // Fill out matrix
            cell
            |> Array2D.iteri
                (fun y x _ ->
                    // Calc 3 scores for 3 possible directions
                    let left,upLeft,up = 
                        match x,y with
                        | 0,0 -> 0,0,0  // top left
                        | 0,_ when x < W-1 -> -999999,(cell.[y,x].s + del) ,(cell.[y-1,x+1].s + del) // Up and left but not last cell on row
                        | _,_ when x = W-1 -> (cell.[y,x-1].s + del),-999999,-999999 // Last cell in row, just left
                        | _,0 -> (cell.[y,x-1].s + del),0,del // Top row
                        | _,_ -> (cell.[y,x-1].s+del), (cell.[y-1,x].s) ,(cell.[y-1,x+1].s+del) // general case
                
                
                    // Final determination of the diagonal cost depends on the two sequences
                    let finalUpleft = 
                        // Assuming we're still in the sequence we're matching
                        if x+y < s1b.Length then
                            match s1b.[x + y],s2b.[y] with
                            | _,' ' -> upLeft
                            | ' ', _ -> upLeft
                            | a,b when a=b -> upLeft + matchCost
                            | a,b when a<>b -> upLeft + mismatchCost
                            | _ -> upLeft - matchCost // This shouldn't happen really but penalize..
                        else
                            upLeft
                
                                    
                    // Update cell value and note direction
                    cell.[y,x] <-
                        if finalUpleft > max left up then SWCell(finalUpleft,DIAG)
                        elif left > up then SWCell(left,LEFT)
                        else SWCell(up,UP)
                    count := !count + 1)
        
            //Dump matrix
            if debug then
                use outF = new StreamWriter("swdump.txt")
                outF.Write("\t")
                Array.iter (fun c -> outF.Write(sprintf "%c\t" c)) s1b
                outF.WriteLine()

                cell
                |> Array2D.iteri (fun y x (v:SWCell) -> 
                    if x = 0 then
                        outF.Write(sprintf "%c" s2b.[y])
                        for j in { 0 .. y} do
                            outF.Write("\t")
                    outF.Write(sprintf "%d\t" v.s)
                    if x = W-1 then
                        outF.Write("\n"))

            //let time2 = DateTime.Now    
            // Generate alignment, starting with highest value on last row of alignment
            let mutable best = -999
            let mutable bestI = -1
            let mutable bestJ = -1
            for i in {0..W-1} do
                if cell.[L-1,i].s > best then
                    best <-  cell.[L-1,i].s
                    bestI <- i
                    bestJ <- L-1
            // Walk up that column
            for j in {L-1-W .. L-1} do
                // How far to index into row as we move up column
                let thisI = bestI+(L-1)-j
                //printfn "Walk.. %d %d %d" thisI j cell.[j,thisI].s
                if thisI < W && cell.[j,thisI].s > best then
                    best <- cell.[j,thisI].s
                    bestJ <- j
                    bestI <- thisI
                
            if debug then printfn "Starting at %d,%d -> %d" bestI bestJ best
        
            let rec traceback x y (res1:ResizeArray<char>,res2:ResizeArray<char>) =
                // Explore the best alignment following the highest score at each step
                if y < 0 then
                    for b in Array.rev (s1b.[..x-1]) do
                        if b <> ' ' then
                            res1.Add(b)
                            res2.Add('-')
                    (res1,res2) 
                else
                    let c1 = if x+y < s1b.Length then s1b.[x+y] else ' '
                
                    match cell.[y,x].direc with
                    | NOWHERE -> 
                        // Hit top row, but still want to go back to origin to
                        // generate full alignment.  Pad out remaining bases
                        for b in Array.rev (s1b.[..x-1]) do
                            if b <> ' ' then
                                res1.Add(b)
                                res2.Add('-')
                        (res1,res2 ) // Hit top row
                    | UP -> 
                        res1.Add('-')
                        res2.Add(s2b.[y])
                        traceback (x+1) (y-1) (res1,res2) // Y insert
                    | LEFT -> 
                        res1.Add(c1)
                        res2.Add('-')
                        traceback (x-1) y  (res1,res2)
                    | DIAG -> 
                        res1.Add(c1)
                        res2.Add(s2b.[y])
                        traceback x (y-1) (res1,res2)  // match
                   
        
            let align1 : ResizeArray<char>  = new ResizeArray<char>()
            let align2 : ResizeArray<char>  = new ResizeArray<char>()
            let _,_ = traceback bestI bestJ (align1,align2)  
            let rev (s:string) =
                s.ToCharArray()
                |> Array.rev
                |> Array.map (fun x -> if x = ' ' then '-' else x)
                |> arr2seq
        
            //printfn "%s\n%s" align1 align2                  
            //let time3 = DateTime.Now    
        
            let ret = 
                (true,
                 (align1 |> Array.ofSeq |> arr2seq |> rev),
                 (align2 |> Array.ofSeq |> arr2seq |> rev),
                 best)
            //let time4 = DateTime.Now    
//            let t1 = time1 - time0
//            let t2 = time2 - time1
//            let t3 = time3 - time2
//            let t4 = time4-time3
            //printfn "SW: %d %f %f %f %f total=%f cells=%d" (align1.Count) (t1.TotalSeconds) (t2.TotalSeconds) (t3.TotalSeconds) 
            //                    (t4.TotalSeconds) ((time4-time0).TotalSeconds) !count
            ret
    
    /// Count mismatches between 2 aligned sequences. 
    /// Sequences align1 and align2 must be of same length
    let countMismatches (align1:string) (align2:string) = 
        if align1.Length <> align2.Length then
            failwithf "Both sequence inputs must be of same length."
        let l = Seq.zip align1 align2 |> Seq.findIndex (fun (_,b) -> b <> '-')
        let rr = Seq.zip align1 align2 |> List.ofSeq |> List.rev |>  List.findIndex (fun (_,b) -> b <> '-')
        let r = align2.Length - rr-1
        let aa = (align1.Substring(l,r-l+1))
        let bb = (align2.Substring(l,r-l+1))
        Seq.zip aa bb  |> Seq.fold (fun c (a,b) -> if a <> b then c+1 else c) 0
    
    /// Retrieve the position of the mismatches between 2 aligned sequences. 
    /// Sequences align1 and align2 must be of same length
    let getMismatchPositions (align1:string) (align2:string) =
        Array.zip (align1.ToCharArray()) (align2.ToCharArray())
        |> Array.mapi (fun i (a1, a2) -> i,a1=a2)
        |> Array.filter (fun (_,b) -> not b)
    
    /// Get coordinate translation from the alignment to the original sequence
    /// Returns an array of the same length as alignment, with the position of the base
    /// in the original sequence. If not present in the original sequence (insertion, out 
    /// of alignment), None is returned.
    let convertCoordAln2Seq (sequenceStartOffset:int) (alignment:string)  = 
        let alnArr = alignment.ToCharArray()
        let coordMap:(int option []) = Array.init (alnArr.Length) (fun _ -> None)
        let rec genCoordMapping alnCoord prevCoord = 
            if alnCoord < alnArr.Length then
                let newSeqCoord = 
                    if alnArr.[alnCoord] <> '-' then
                        coordMap.[alnCoord] <- Some (prevCoord + 1)
                        prevCoord + 1
                    else
                        prevCoord
                genCoordMapping (alnCoord + 1) newSeqCoord
        if alnArr.Length >= 1 then
           genCoordMapping 0 (sequenceStartOffset-1)
        coordMap

    /// Get coordinate translation from the original sequence to the alignment
    /// Returns an array of the same length as the original sequence, with the position of the base
    /// in the alignment. If not present in the alignment, None is returned.
    let convertCoordSeq2Aln (sequenceStartOffset:int) (sequenceLength:int) (alignment:string)  = 
        let alnArr = alignment.ToCharArray()
        let coordMap:(int option []) = Array.init (sequenceLength) (fun _ -> None)
        let rec genCoordMapping alnCoord prevCoord = 
            if alnCoord < alnArr.Length then
                let newSeqCoord = 
                    if alnArr.[alnCoord] <> '-' then
                        coordMap.[prevCoord + 1] <- Some alnCoord
                        prevCoord + 1
                    else
                        prevCoord
                genCoordMapping (alnCoord + 1) newSeqCoord
        if alnArr.Length >= 1 then
           genCoordMapping 0 (sequenceStartOffset-1)
        coordMap

    //[<Obsolete>]
    /// Banded Smith Waterman align two sequences, seeding alignment with shared N mers
    /// suitable for similar sequences with relatively straight alignment. 
    /// N - size of mer
    /// seq1 - char array sequence 1
    /// seq2 - char array sequence 2    
    /// (Window set to 24)
    let smithWaterman N (seq1 : char []) (seq2 : char []) =
        let W = 24
        let (aligned, a1, a2, _) = _smithWaterman false None N W seq1 seq2
        (aligned, a1, a2)

    //[<Obsolete>]
    /// Banded Smith Waterman align two sequences, seeding alignment with shared N mers
    /// suitable for similar sequences with relatively straight alignment. 
    /// N - size of mer (decrease for more dissimilar sequences)
    /// W - width of alignment window (increase if alignment wanders)
    /// seq1 - char array sequence 1
    /// seq2 - char array sequence 2    
    let smithWatermanW N W (seq1 : char []) (seq2 : char []) =
        let (aligned, a1, a2, _) = _smithWaterman false None N W seq1 seq2
        (aligned, a1, a2)
    
    /// Banded Smith Waterman align two sequences, seeding alignment with shared N mers
    /// suitable for similar sequences with relatively straight alignment. 
    /// N - size of mer (decrease for more dissimilar sequences)
    /// W - width of alignment window (typically 24, increase if alignment wanders)
    /// seq1 - char array sequence 1
    /// seq2 - char array sequence 2    
    let smithWatermanWscore N W (seq1 : char []) (seq2 : char []) =
        _smithWaterman false None N W seq1 seq2


    /// Find the alignment boundaries of two aligned sequences. Both sequences should have the same length. 
    /// This method is typically applied on the output of smithWaterman
    let findAlignmentBoundaries (a1:string) (a2:string) =
        if a1.Length <> a2.Length then failwith "ERROR: Aligned sequences must have the same length"
        let rec leadingPadCount (a:string) i = 
            if i = a.Length then i else
                if a.[i] = '-' then leadingPadCount a (i+1) else i
        let rec trailingPadCount (a:string) i = 
            if i = 0 then 0 else
                if a.[i] = '-' then trailingPadCount a (i-1) else i
        let left = max (leadingPadCount a1 0) (leadingPadCount a2 0) 
        let right = min (trailingPadCount a1 (a1.Length-1)) (trailingPadCount a2 (a1.Length-1))
        (left, right)

    /// Find the coordinates of the aligned sequence in the original query sequence. 
    let findAlignmentCoordInOriginalSequence (alignment:string) (sequence:char []) = 
        let seqString = arr2seq sequence
        let alignmentNoGap = alignment.Replace("-","")
        let left = seqString.IndexOf alignmentNoGap
        if left = -1 then
            failwith "ERROR: alignment sequence is not included in the query sequence. Check that the alignment really corresponds to this query sequence."
        let right = left + alignmentNoGap.Length - 1 
        (left, right)

    