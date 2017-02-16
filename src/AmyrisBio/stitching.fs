/// Algorithms for stitching genetic constructs simulating molecular biology operations
module Amyris.Bio.Stitching
open biolib
open Amyris.Dna
open Amyris.ErrorHandling
open utils


type OverlapSearchParameters =
   {minOverlap: uint64 option;
    maxOverlap: uint64 option}

let defaultMinOverlap = 50UL;

let defaultSearchParameters =
   {minOverlap = Some(defaultMinOverlap);
    maxOverlap = None}

type OverlapStitchRequest =
   {seq0: Dna; seq1: Dna; margin0: Dna; margin1: Dna; searchParams: OverlapSearchParameters}

type OverlapStitchResult = Option<Dna>

/// Compute the DNA sequence expected to result from homologous recombination.
/// This function assumes both sequences are antiparallel and provided in the sense direction,
/// seq0 ------------------------------------LLLLLLL>
///                       <LLLLL------------------------------------ seq1
///                             overlapoverla
/// and will internally take the reverse complement of seq1 to find the overlap.
/// The returned sequence will read in the sense direction of seq0.
/// Any margin tails (for example, linkers) must be explicitly provided or stitching will fail.
/// This function will check for exact overlap within the provided window.
let overlapStitchWithMargins req : Result<OverlapStitchResult, string> =
    let minOverlap = int (defaultArg req.searchParams.minOverlap defaultMinOverlap)
        
    /// Check to make sure that:
    /// - the specified margins are actually the tails of the sequence.
    /// - the sequences minus margins satisfy a minimum length requirement.
    let checkSeq (toCheck: Dna) (margin: Dna) which =
        let errs = ResizeArray()
        if toCheck.Length - margin.Length < minOverlap then
            errs.Add(sprintf
                "Sequence %i minus margin is not long enough to satify a minimum overlap of %i"
                which
                minOverlap)
        if not(toCheck.EndsWith(margin)) then
            errs.Add(sprintf
                "Sequence %i does not end with the specified margin."
                which)
        errs

    /// Check both sequences, agglomerating errors from both.
    let checkSeqs req =
        let errs = checkSeq req.seq0 req.margin0 0
        errs.AddRange(checkSeq req.seq1 req.margin1 1)
        if errs.Count <> 0 then
            Bad(errs |> List.ofSeq)
        else
            ok(req)


    let tryStitch req =

        let seq0, seq1, margin0, margin1 = req.seq0, req.seq1, req.margin0, req.margin1

        let seq0ActiveLength = seq0.Length - margin0.Length
        let seq1ActiveLength = seq1.Length - margin1.Length

        /// Length of the smaller active non-margin region out of the two sequences
        let activeLength = min seq0ActiveLength seq1ActiveLength

        /// Last index of seq0 that isn't part of the margin
        let seq0End = seq0ActiveLength - 1

        let seq1RC = seq1.RevComp()

        /// First index of seq1RC that isn't part of the margin
        let seq1RCStart = margin1.Length

        /// If max specified, constrain to maximum length we have available.
        /// Otherwise, use full available overlap region.
        let maxOverlap =
            match req.searchParams.maxOverlap with
            | Some(m) -> min (int m) activeLength
            | None -> activeLength

        /// Try to align subsections of each sequence recursing over length.
        let rec tryAlign len =
            let seq0Slice = seq0.Subseq(seq0End - len + 1, seq0End)
            let seq1Slice = seq1RC.Subseq(seq1RCStart, seq1RCStart + len - 1)
            // If slices are equal, we have aligned the two stitches.
            if compareSlices seq0Slice seq1Slice then
                Some(len)
            elif len >= maxOverlap then
                None
            else
                tryAlign (len+1)

        match tryAlign minOverlap with
        | None -> ok None
        | Some(l) ->
            // Return a new Dna sequence from the computed overlap
            // Get the tail of seq1 minus margin and overlap region
            let stitchedSeq = Seq.append (seq0.Subseq(0, seq0End)) (seq1RC.Subseq(seq1RCStart + l))
            // skip validation as incoming DNA should already be validated
            ok (Some(Dna(stitchedSeq, false)))

    req |> (checkSeqs >=> tryStitch)

type LoopoutSearchParams =
   {prefixSearchLength: int;
    minOverlap: uint64;
    maxOverlap: uint64}

let defaultLoopoutSearchParams =
   {prefixSearchLength = 8;
    minOverlap = 50UL;
    maxOverlap = 500UL}

type LoopoutScars = {scar0: Dna; scar1: Dna option}

type LoopoutResult = Option<LoopoutScars>

/// Given two sequences, return the tails of the two sequences following
/// any mismatch in character.
let mismatchedTails (a: Dna) (b: Dna) =
    let mutable foundMismatch = false
    let aTail, bTail =
        seq {
            for x, y in zipWithPad a b do
                if foundMismatch then yield (x, y)
                elif x <> y then
                    foundMismatch <- true
                    yield (x, y)}
        |> Array.ofSeq
        |> Array.unzip

    let dnaify x = Dna(Seq.choose id x, false)

    (dnaify aTail, dnaify bTail)

type SeqCompareResultWithIndel =
    | NotEqual
    | Equal
    | EqualWithIndel

/// Determine whether two Dna segments are the same allowing for a single indel.
/// The indel cannot be the first or last BP in the sequence.  Note that this restriction
/// implies that some problematic edge cases will be declared false when they might be
/// true.
let equalWithIndel (a: Dna) (b: Dna) =
    /// diffLen = 0 implies no indel, +1 means missing BP in b, -1 means missing BP in a
    let diffLen = a.Length - b.Length

    if (abs diffLen) > 1 then NotEqual
    elif diffLen = 0 then (if a = b then Equal else NotEqual)
    else
        let aTail, bTail = mismatchedTails a b

        let longer, shorter =
            if aTail.Length > bTail.Length then aTail, bTail else bTail, aTail

        // Reject some edge cases
        // Longer sequence has extra BP appended or doesn't match
        if longer.Length = 1 then NotEqual

        // Longer sequence has extra BP prepended or doesn't match
        elif longer.Length = (if diffLen = 1 then a else b).Length then NotEqual
        elif compareSlices (longer.Subseq(1)) shorter then EqualWithIndel
        else NotEqual


type LoopoutRequest = {s: Dna; searchParams: LoopoutSearchParams}

/// Validate and massage loopout search parameters.
let validateLoopoutParameters (req: LoopoutRequest) =
    // Sequence needs to be at least 2x the size of the min loopout construct
    if req.s.Length < 2 * int(req.searchParams.minOverlap) then
        fail (sprintf
            "Sequence is too short to loop out; min repeat size is %i, sequence is %i bp."
            (req.searchParams.minOverlap)
            (req.s.Length))
    else
        // Restrict the max overlap search size to seq length / 2
        // Might want a far more conservative limit here, like some kind of min looped out seq length
        ok {req with searchParams = {req.searchParams with maxOverlap = uint64(req.s.Length / 2)}}


/// Compute the loopout scar(s) for a sequence, if it satisfies these requirements:
/// The direct repeat segments must be at the very beginning and very end of the sequence.
/// The direct repeats may be identical or differ by a single indel.
/// Returns one or two sequences (two are returned if the repeat contained a indel).
let computeLoopoutScarPostValidation (req: LoopoutRequest) =

    let s, searchParams = req.s, req.searchParams

    let snippetLen = searchParams.prefixSearchLength
    let tail = s.str.Substring(s.Length - snippetLen)
    let head = s.str.Substring(0, snippetLen)

    let minIndex = (int searchParams.minOverlap) - snippetLen
    let maxIndex = (int searchParams.maxOverlap) - snippetLen

    /// get the sequence start indices at mirror position that match the snippet sequence
    let getMirrorMatches ind =
        let mirrorIndexStart = s.Length - ind - tail.Length - 1
        // prefer no offset to +- 1
        let offsets = [0; 1; -1]

        let checkMatch offset =
            let start = mirrorIndexStart + offset
            let finish = start + head.Length - 1
            compareSlices head (s.Subseq(start, finish))

        offsets 
        |> List.filter checkMatch
        |> List.map (fun offset -> offset + mirrorIndexStart)

    /// Given a candidate mirror match, check if the sequences are equivalent.
    let checkMatch leftStart rightStart =
        // leftStart is where the match was found, full seq includes the tail segment
        let leftSeq = s.[..leftStart + tail.Length - 1]
        let rightSeq = s.[rightStart..]
        match equalWithIndel leftSeq rightSeq with
        | NotEqual -> None
        | Equal -> Some {scar0 = leftSeq; scar1 = None}
        | EqualWithIndel -> Some {scar0 = leftSeq; scar1 = Some(rightSeq)}

    /// Search for candidate locations where this diagram holds true:
    /// HEADSEQ...n bp...TAILSEQ... ...HEADSEQ...n+-1 bp...TAILSEQ
    /// If candidate matches are found, check them and return the loopout scar(s) if they match.
    let rec search (start: int) =
        // Starting at start, find the next occurrence of the tail snippet
        match s.str.IndexOf(tail, start) with
        | -1 -> ok None // No matches left, we're done
        | i when i < minIndex ->
            // Found a candidate but it implies insufficient overlap.  Continue.
            search (i+1)
        | i when i > maxIndex ->
            // Found a candidate but it implies too large of an overlap.  Abort.
            ok None
        | i ->
            // Potentially good candiate, see if its a loose match
            // Check for head snippet at the mirror position +- 1 bp
            match getMirrorMatches i with
            | [] -> search (i+1) // no mirror match, keep looking
            | x -> // found some matches
                match x |> List.choose (checkMatch i) with
                | [] -> search (i+1) // no good matches
                | [scars] -> ok (Some scars) // one perfect match
                | x -> fail (sprintf "Multiple possible loopout matches found: %A" x)
    
    search 0

/// Compute the loopout scar(s) for a sequence, if it satisfies these requirements:
/// The direct repeat segments must be at the very beginning and very end of the sequence.
/// The direct repeats may be identical or differ by a single indel.
/// Returns one or two sequences (two are returned if the repeat contained a indel).
let computeLoopoutScar (req: LoopoutRequest) =
    req |> (validateLoopoutParameters >> bind computeLoopoutScarPostValidation)