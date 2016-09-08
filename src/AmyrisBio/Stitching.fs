/// Algorithms for stitching genetic constructs
module Amyris.Bio.Stitching

/// Domain type for immutable DNA sequences.
/// Interally backed by a char array.
/// Ensures that contents have been uniformly converted to uppercase if validation is performed.
/// Lazy memoized access to string representation as well as reverse compliment version of itself.
/// This data structure is absolutely not thread safe.
type Dna private (asArray: char [], rc: Dna option) =

    // Private mutable values to permit memoization.
    let mutable asString: string option = None
    let mutable revComp: Dna option = rc

    new(sequence, ?validate) =
        let doValidate = match validate with | Some(v) -> v | None -> true
        if doValidate then
            let badChars =
                sequence
                |> Seq.filter (fun c -> not (biolib.isDnaBaseStrict c))
                |> Array.ofSeq
            if badChars.Length > 0 then
                failwithf "Found bad chars in DNA string: %A" badChars
        
        let seqArr =
            sequence
            |> (if doValidate then Seq.map (fun (c: char) -> System.Char.ToUpper(c)) else id)
            |> Array.ofSeq
        Dna(seqArr, None)

    with
    /// Return the char array representation of this DNA payload.
    /// Client code shouldn't need to call this.
    member x.arr = asArray
    /// Return the string array representation of this DNA payload.
    /// Client code shouldn't need to call this.
    member x.str =
        match asString with
        | Some(s) -> s
        | None ->
            let s = new System.String(asArray)
            asString <- Some(s)
            s

    /// Return a sequence viewing a slice of this DNA sequence.
    /// Both start and end indices are includsive to correspond to the semantics of F# slices.
    /// Unlike an array or string slice, this is represented as a sequence and has minimal
    /// memory footprint and does not copy the slice to a new data structure.
    member x.Subseq(?startInd, ?endInd) =
        let starti = match startInd with | Some(s) -> s | None -> 0
        let endi = match endInd with | Some(s) -> s | None -> asArray.Length - 1
        utils.arraySliceSeq asArray starti endi

    /// Return True if this DNA sequence ends with the provided sequence.
    member x.EndsWith(other: Dna) = x.str.EndsWith(other.str)

    /// Return True if this DNA sequence starts with the provided sequence.
    member x.StartsWith(other: Dna) = x.str.StartsWith(other.str)

    member x.Length = asArray.Length

    /// Get the reverse compliment of this sequence.
    member x.RevComp() =
        match revComp with
        | Some(rc) -> rc
        | None ->
            let rcArr = biolib.revComp asArray
            let rcDna = new Dna(rcArr, Some(x))
            revComp <- Some(rcDna)
            rcDna

type OverlapSearchParameters =
   {minOverlap: uint64 option;
    maxOverlap: uint64 option}

let defaultMinOverlap = 50UL;

let defaultSearchParameters =
   {minOverlap = Some(defaultMinOverlap);
    maxOverlap = None}

/// Compute the DNA sequence expected to result from homologous recombination.
/// This function assumes both sequences are provided in the sense direction,
/// and will internally take the reverse compliment of seq1 to find the overlap.
/// The returned sequence will read in the sense direction of seq0.
/// Any linker tails must be explicitly provided or stitching will fail.
/// This function will check for exact overlap within the provided window,
/// starting in the suggested hint direction to check.
let overlapStitchWithLinkers
        (seq0: Dna) (seq1: Dna) (seq0Linker: Dna) (seq1Linker: Dna) searchParams =
    let minOverlap = int (match searchParams.minOverlap with | Some(m) -> m | None -> defaultMinOverlap)
        
    /// Check to make sure that:
    /// - the specified linkers are actually the tails of the sequence.
    /// - the sequences minus linkers satisfy a minimum length requirement.
    let checkSeq (toCheck: Dna) (linker: Dna) which =
        if toCheck.Length - linker.Length < minOverlap then
            failwithf
                "Sequence %i minus linker is not long enough to satify a minimum overlap of %i"
                which
                minOverlap
        if not(toCheck.EndsWith(linker)) then
            failwithf
                "Sequence %i does not end with the specified linker."
                which

    checkSeq seq0 seq0Linker 0
    checkSeq seq1 seq1Linker 1

    let seq0ActiveLength = seq0.Length - seq0Linker.Length
    let seq1ActiveLength = seq1.Length - seq1Linker.Length

    /// Length of the smaller active non-linker region out of the two sequences
    let activeLength = min seq0ActiveLength seq1ActiveLength

    /// Last index of seq0 that isn't part of the linker
    let seq0End = seq0ActiveLength - 1

    let seq1RC = seq1.RevComp()

    /// First index of seq1RC that isn't part of the linker
    let seq1RCStart = seq1Linker.Length

    /// If max specified, constrain to maximum length we have available.
    /// Otherwise, use full available overlap region.
    let maxOverlap =
        match searchParams.maxOverlap with
        | Some(m) -> min (int m) activeLength
        | None -> activeLength

    ///
    let rec _tryAlign len =
        let seq0Slice = seq0.Subseq(seq0End - len + 1, seq0End)
        let seq1Slice = seq1RC.Subseq(seq1RCStart, seq1RCStart + len - 1)
        // If slices are equal, we have aligned the two stitches.
        if utils.compareSliceSeqs seq0Slice seq1Slice then
            Some(len)
        elif len >= maxOverlap then
            None
        else
            _tryAlign (len+1)
    match _tryAlign minOverlap with
    | None -> failwith "No exact overlap found for input sequences."
    | Some(l) ->
        // Return a new Dna sequence from the computed overlap
        // Get the tail of seq1 minus linker and overlap region
        let stitchedSeq = Seq.append (seq0.Subseq(0, seq0End)) (seq1RC.Subseq(seq1RCStart + l))
        // skip validation as incoming DNA should already be validated
        Dna(stitchedSeq, false)

type LoopoutSearchParams =
   {prefixSearchLength: int;
    minOverlap: uint64;
    maxOverlap: uint64;
    allowAnomaly: bool}

let defaultLoopoutSearchParams =
   {prefixSearchLength = 8;
    minOverlap = 50UL;
    maxOverlap = 500UL;
    allowAnomaly = false}

/// Find a pair of candiates indices that bookend a potential direct repeat
/// Allow matches to be within +-1 bp to allow for up to one SNP.
let findBookendCandidates (s: Dna) searchParams start =
    let snippetLen = searchParams.prefixSearchLength
    let tail = s.str.Substring(s.Length - snippetLen)
    let head = s.str.Substring(0, snippetLen)

    let minIndex = (int searchParams.minOverlap) - snippetLen
    let maxIndex = (int searchParams.maxOverlap) - snippetLen

    /// Given a candidate match of tail to head, see if head is present at the mirror position.
    let checkMirrorMatch index =
        ()

    let rec findCloseMatch (start: int) =
        // Starting at start, find the next occurrence of the tail snippet
        match s.str.IndexOf(tail, start) with
        | -1 -> None // No matches left, we're done
        | i when i < minIndex ->
            // Found a candidate but it implies insufficient overlap.  Continue.
            findCloseMatch (i+1)
        | i when i > maxIndex ->
            // Found a candidate but it implies too large of an overlap.  Abort.
            None
        | i ->
            // Potentially good candiate, see if its a loose match
            None
    ()
         
/// Given a sequence that may loop itself out, compute the looped-out sequence.
/// The loopout repeat sections must be at the very beginning and very end of the sequence.
/// This algorithm optionally tolerates up to one SNP or mutation in the two repeats.
let computeLoopoutSequence (sequence: Dna) =

    ()