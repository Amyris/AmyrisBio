/// Algorithms for stitching genetic constructs
module Amyris.Bio.Stitching
open biolib
open Amyris.Dna
open Chessie.ErrorHandling


type OverlapSearchParameters =
   {minOverlap: uint64 option;
    maxOverlap: uint64 option}

let defaultMinOverlap = 50UL;

let defaultSearchParameters =
   {minOverlap = Some(defaultMinOverlap);
    maxOverlap = None}

type OverlapStitchRequest =
   {seq0: Dna; seq1: Dna; margin0: Dna; margin1: Dna; searchParams: OverlapSearchParameters}

type OverlapStitchResult =
    | Stitchable of Dna
    | Unstitchable

/// Compute the DNA sequence expected to result from homologous recombination.
/// This function assumes both sequences are provided in the sense direction,
/// and will internally take the reverse compliment of seq1 to find the overlap.
/// The returned sequence will read in the sense direction of seq0.
/// Any margin tails (for example, linkers) must be explicitly provided or stitching will fail.
/// This function will check for exact overlap within the provided window,
/// starting in the suggested hint direction to check.
let overlapStitchWithMargins req =
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
            fail(errs)
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
            if utils.compareSliceSeqs seq0Slice seq1Slice then
                Some(len)
            elif len >= maxOverlap then
                None
            else
                tryAlign (len+1)

        match tryAlign minOverlap with
        | None -> ok Unstitchable
        | Some(l) ->
            // Return a new Dna sequence from the computed overlap
            // Get the tail of seq1 minus margin and overlap region
            let stitchedSeq = Seq.append (seq0.Subseq(0, seq0End)) (seq1RC.Subseq(seq1RCStart + l))
            // skip validation as incoming DNA should already be validated
            ok (Stitchable(Dna(stitchedSeq, false)))

    req |> (checkSeqs >> bind tryStitch)

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