namespace Amyris.Bio

/// Optimized lookup for searches of DNA sequencew against whole chromosomes
module SequenceSearch =

    open Amyris.Bio.biolib
    open Amyris.Bio.utils
    open Amyris.Bio.SuffixTree
    open Amyris.Bio.smithWaterman

    // A - First method- use suffix trees to find exact matches of substrings
    // 1- Find exact matches of substrings

    type seqMatch = {queryL:int; queryR:int; targetL:int; targetR:int; targetFwd:bool; targetChr:string}
    /// Find all substring exact matches of one sequence on a chromosome, on fwd and reverse direction, down to a length of matchLengthCutoff (usually 50) 
    ///     findSubstringMatchesOnChromosome (matchLengthCutoff:int) (targetSeq:char []) (querySeq:char [])
    ///  targetSeq must be the chromosome. 
    let findSubstringMatchesOnChromosome (matchLengthCutoff:int) (targetSeq:char []) (querySeq:char []) =
        if querySeq.Length < matchLengthCutoff then 
            sprintf "WARNING: query sequence shorter than the match length cutoff %d. No match returned"
                matchLengthCutoff |> stderr.Write
            []
        else
            // Create a string with forward$revcomp  ie AATCG$CGATT 
            let chromString = Array.append (Array.append targetSeq [|'$'|]) (revComp targetSeq) |> arr2seq
            let chromST = SuffixTree(chromString)
            chromST.buildFwdChain()
            let sequenceString = querySeq |> arr2seq

            /// Function to translate coordinates on chromString to real coordinates on the chromosome
            let getChromCoordinateFromChromString (l:int) (r:int) =
                if l > r then failwith "ERROR: left should be < right"
                if l < targetSeq.Length && r > targetSeq.Length then
                    failwith "ERROR: left coordinate on forward and right on reverse part of chromosome string. Both should be on the same part."
                if l = targetSeq.Length || r = targetSeq.Length then
                    failwith "ERROR: coordinate cannot include the separator character of the chromosome string"
                if r < targetSeq.Length then 
                    (l,r,true)
                elif l > targetSeq.Length then 
                    (chromString.Length+1-r, chromString.Length+1-l, false)
                else
                    failwith "ERROR: issue with coordinates. This should have been captures previously in the function"

            /// recurse on sequence substrings, from the left. For each substring, find the longest match and store it if it is longer than the threshold.
            /// If a match is a longest match is reported, mask the reported sequence and search again for a potential shorter match. 
            let rec findExactMatches
                    (matches:seqMatch list)
                    (targetSeqString:string)
                    (targetSeqST:SuffixTree)
                    (sequenceLeft:int) =
                if (sequenceString.Length - sequenceLeft) < matchLengthCutoff then 
                    // we are at the end of the sequence, stop here
                    matches
                else
                    // look for longest match
                    let longestMatchLength = targetSeqST.LongestMatch(sequenceString.[sequenceLeft..])
                    if longestMatchLength < matchLengthCutoff then 
                        // longest match is shorter than cutoff, go to next substring
                        findExactMatches matches targetSeqString targetSeqST (sequenceLeft+1)
                    else
                        // retrieve all matches 
                        let allMatches =
                            targetSeqST.FindAll(sequenceString.[sequenceLeft..(sequenceLeft+longestMatchLength-1)])
                        let allMatchesCoord =
                            allMatches
                            |> Seq.map (fun l -> getChromCoordinateFromChromString l (l+longestMatchLength-1) )
                            |> Seq.map (fun (tl, tr, tf) ->
                                {queryL=sequenceLeft;
                                 queryR=(sequenceLeft+longestMatchLength-1);
                                 targetL=tl; 
                                 targetR=tr; 
                                 targetFwd=tf; 
                                 targetChr=""})
                            |> List.ofSeq
                        // mask their sequence on the target string, and build a new suffix tree
                        let maskedTargetSeqString =
                            allMatches
                            |> Seq.fold
                                (fun (prevSeq:string) l ->
                                    if l-1 < 0 then
                                        ((Array.create longestMatchLength '#' |> arr2seq)
                                        + prevSeq.[(l+longestMatchLength)..])
                                    else
                                        (prevSeq.[..(l-1)]
                                        + (Array.create longestMatchLength '#' |> arr2seq)
                                        + prevSeq.[(l+longestMatchLength)..])
                                )
                                targetSeqString

                        let maskedTargetSeqST = SuffixTree(maskedTargetSeqString)
                        maskedTargetSeqST.buildFwdChain()
                        // search for the same substring again
                        findExactMatches (allMatchesCoord @ matches) maskedTargetSeqString maskedTargetSeqST sequenceLeft
            findExactMatches [] chromString chromST 0

    /// Find all substring exact matches of the chromosomes of the genome, on fwd and reverse directions, down to a length of matchLengthCutoff (usually 50)
    let findSubstringMatchesOnGenome
            (matchLengthCutoff:int)
            (genome:Map<string, char []>)
            (sequence:char []) =
        genome
        |> Seq.map (fun pk ->
            (findSubstringMatchesOnChromosome matchLengthCutoff pk.Value sequence
            |> Seq.map (fun m -> {m with targetChr=pk.Key})))
        |> Seq.concat |> List.ofSeq

    // 2- Group substring matches and infer full region of similarity

    /// Type holding filtered matches. queryCoord2targetCoord is a function giving an approximate mapping between query coordinates and target coordinates, up to the given tolerance.
    type FilteredMatch = 
        {tolerance:int;
         queryCoord2targetCoord:((int*int)->(int*int)); 
         chr:string; 
         fwd:bool; 
         inferredL:int; 
         inferredR:int; 
         matches:(seqMatch list); 
         SwScoreAndLength:(int*int*string*string) option}

    /// translate a query sequence coordinate to a target sequence coordinate, given reference query and target coordinates, and target direction
    let query2targetCoord queryRefL (targetRefL,targetRefR) fwd (queryL,queryR) =
        if fwd then
            (targetRefL + (queryL - queryRefL), targetRefL + (queryR - queryRefL))
        else 
            (targetRefR - (queryR - queryRefL), targetRefR - (queryL - queryRefL))

    /// group substring matches together, given some tolerance on how much submatches can be offset.
    /// groupMatches (tolerance:int) (queryLength:int) (matches:seqMatch list)
    let groupMatches (tolerance:int) (queryLength:int) (matches:seqMatch list) =

        let isInFilteredMatch (myMatch:seqMatch) (fM:FilteredMatch) = 
            let (predL, predR) = fM.queryCoord2targetCoord (myMatch.queryL, myMatch.queryR)
            fM.chr=myMatch.targetChr &&
            fM.fwd=myMatch.targetFwd &&
            (abs(predL - myMatch.targetL) < fM.tolerance) &&
            (abs(predR - myMatch.targetR) < fM.tolerance)

        let rec _groupMatches (filteredMatches:FilteredMatch list) (matches:seqMatch list) = 
            match matches with 
            | [] -> filteredMatches
            | currMatch::remMatches ->
                if (filteredMatches.Length > 0 &&
                        filteredMatches
                        |> Seq.map (isInFilteredMatch currMatch)
                        |> Seq.reduce (||))
                then 
                    // current matches fall in an already filtered match
                    let newFilteredMatches =
                        filteredMatches
                        |> List.map (fun fM ->
                            if isInFilteredMatch currMatch fM then
                                {fM with matches=currMatch::fM.matches}
                            else fM)
                    _groupMatches newFilteredMatches remMatches
                else 
                    let q2tCoord = query2targetCoord currMatch.queryL (currMatch.targetL, currMatch.targetR) currMatch.targetFwd
                    let (infL,infR)=q2tCoord (0,queryLength)
                    let newFilteredMatches =
                        {tolerance=tolerance; 
                         queryCoord2targetCoord=q2tCoord; 
                         chr=currMatch.targetChr; 
                         fwd=currMatch.targetFwd; 
                         inferredL=infL; 
                         inferredR=infR; 
                         matches=[currMatch]; 
                         SwScoreAndLength=None}
                        ::filteredMatches

                    _groupMatches newFilteredMatches remMatches
        _groupMatches [] matches

    // 3- confirm similarity by local SW alignment

    ///// Find the alignment boundaries of two aligned sequences. Both sequences should have the same length. 
    ///// This method is typically applied on the output of smithWaterman
    //let findAlignmentBoundaries (a1:string) (a2:string) =
    //    if a1.Length <> a2.Length then failwith "ERROR: Aligned sequences must have the same length"
    //    let rec leadingPadCount (a:string) i = 
    //        if i = a.Length then i else
    //            if a.[i] = '-' then leadingPadCount a (i+1) else i
    //    let rec trailingPadCount (a:string) i = 
    //        if i = 0 then 0 else
    //            if a.[i] = '-' then trailingPadCount a (i-1) else i
    //    let left = max (leadingPadCount a1 0) (leadingPadCount a2 0) 
    //    let right = min (trailingPadCount a1 (a1.Length-1)) (trailingPadCount a2 (a1.Length-1))
    //    (left, right)

    let getSwAlignmentScoreAndLength (query:char []) (chromosome:char []) (filteredMatch:FilteredMatch) =
        let target =
            let start = (max 0 (filteredMatch.inferredL-filteredMatch.tolerance))
            let stop = (min (chromosome.Length-1) (filteredMatch.inferredR+filteredMatch.tolerance))
            (if filteredMatch.fwd then id else revComp) (chromosome.[start..stop])

        let ((*aligned*)_,a1,a2, score) = smithWatermanWscore 20 24 query target
        let (left, right) = findAlignmentBoundaries a1 a2
        (score, right-left+1, a1, a2)
    
    let evaluateFilteredMatch
            (query:char [])
            (genome:Map<string, char []>)
            (filteredMatches:FilteredMatch list) = 
        filteredMatches
        |> List.map (fun fM ->
            let (score, length, queryAln, genomeAln) =
                getSwAlignmentScoreAndLength query genome.[fM.chr] fM
            {fM with SwScoreAndLength=Some (score,length, queryAln, genomeAln)})


