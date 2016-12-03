namespace Amyris.Dna
open Amyris.Bio
open System
open System.Collections
open System.Collections.Generic
open biolib
open utils

/// Is a given DNA sequence allowed to have ambiguous bases or not?
type SequenceSemantics =
    | Explicit
    | AllowAmbiguousBases

 // NOTE: this domain type was added quite some time after much of this package was developed.
// As such, most functions are written in terms of either char [] or string, while only a few
// newer ones are written in terms of Dna.  Ideally, new development should leverage this type
// to enforce some assumptions such as valid DNA content and to provide efficiency through
// immutable seq views and memoization.
/// Domain type for immutable DNA sequences.
/// Interally backed by a char array.
/// Ensures that contents have been uniformly converted to uppercase if validation is performed.
/// Lazy memoized access to string representation as well as reverse compliment version of itself.
/// This data structure is absolutely not thread safe due to lazy memoization.  A function
/// that fills in all the memoizations would make the data structure thread-safe after it
/// returns.
type Dna private (asArray:char [], rc: Dna option, mode: SequenceSemantics) =
    // TODO: refactor this class to use ImmutableArray rather than Array.  At the moment
    // clients could shoot themselves in the foot by mutating the array contents, or feeding
    // it to a function that mutates it.

    // Private mutable values to permit memoization.
    let mutable asString: string option = None
    let mutable revCompPartner: Dna option = rc

    new(sequence, ?validate, ?mode) =
        let doValidate = defaultArg validate true
        let mode = defaultArg mode Explicit
        let baseCheck = match mode with | AllowAmbiguousBases -> isDnaBase | Explicit -> isDnaBaseStrict
        if doValidate then
            let badChars =
                sequence
                |> Seq.filter (fun c -> not (baseCheck c))
                |> Array.ofSeq
            if badChars.Length > 0 then
                failwithf "Found bad chars in DNA string: %A" badChars
        
        let seqArr =
            sequence
            |> (if doValidate then Seq.map (fun (c: char) -> System.Char.ToUpper(c)) else id)
            |> Array.ofSeq
        Dna(seqArr, None, mode)


    with
    member x.SemanticMode = mode
    /// Return the char array representation of this DNA payload.
    /// Client code shouldn't need to call this except for interop with older functions.
    member x.arr = asArray
    /// Return the string array representation of this DNA payload.
    /// Client code shouldn't need to call this except for interop with older functions.
    member x.str =
        match asString with
        | Some(s) -> s
        | None ->
            let s = new System.String(asArray)
            asString <- Some(s)
            s

    // Implement convenience interfaces for various standard things

    override x.ToString() = x.str

    // Implement F# indexing to return a single char.
    member x.Item(i) = x.arr.[i]

    // Implement F# slicing syntax to return another Dna type via copy.
    member x.GetSlice(start: int option, finish: int option) =
        let start = defaultArg start 0
        let finish = defaultArg finish (asArray.Length-1)
        Dna(asArray.[start..finish], None, mode)

    interface IEnumerable<char> with
        member x.GetEnumerator() = (Seq.cast<char> asArray).GetEnumerator()
    interface IEnumerable with
        member x.GetEnumerator() = asArray.GetEnumerator()

    interface IEquatable<Dna> with
        member x.Equals(other) = x.arr = other.arr

    // Implement equality between similiar types
    override x.Equals(other) =
        match other with
        | :? Dna as o -> x.arr = o.arr
        | _ -> false

    interface IComparable with
        member x.CompareTo(other) =
            match other with
            | null -> 1
            | :? Dna as o -> compare x.arr o.arr
            | _ -> invalidArg "other" "not an instance of Dna"
    interface IComparable<Dna> with
        member x.CompareTo(other) = compare x.arr other.arr

    // Implement hash as hash of the underlying sequence, ignoring memoization.
    override x.GetHashCode() = asArray.GetHashCode()

    /// Return a view of a slice of this DNA sequence.
    /// Start and finish are both inclusive indices.
    member x.Subseq(start, ?finish: int) =
        let finish = defaultArg finish (asArray.Length-1)
        ArraySegment(asArray, start, finish-start+1)

    /// Split this DNA about a particular base pair.
    /// The index provided will be the last base pair of the first half of the split.
    member x.Split(endOfFirstPart) =
        x.[..endOfFirstPart], x.[endOfFirstPart+1..]

    /// Return True if this DNA sequence ends with the provided sequence.
    member x.EndsWith(other: Dna) = x.str.EndsWith(other.str)

    /// Return True if this DNA sequence starts with the provided sequence.
    member x.StartsWith(other: Dna) = x.str.StartsWith(other.str)

    /// Return True if this DNA sequence contains the provided sequence.
    member x.Contains(other: Dna) = x.str.Contains(other.str)

    /// Return the index of the first appearence of this substring.
    member x.IndexOf(substr: string) = x.str.IndexOf(substr)

    /// Return the index of the first appearence of this character.
    member x.IndexOf(c: char) = x.str.IndexOf(c)

    member x.Length = asArray.Length

    /// Get the reverse compliment of this sequence.
    member x.RevComp() =
        match revCompPartner with
        | Some(rc) -> rc
        | None ->
            let rcArr = revComp asArray
            let rcDna = new Dna(rcArr, Some(x), mode)
            revCompPartner <- Some(rcDna)
            rcDna

module DnaOps =
    /// Determine what combined mode a set of Dna sequences should have.
    let private combinedMode (seqs: seq<Dna>) =
        seqs
        |> Seq.fold
            (fun accumMode s ->
                match accumMode, s.SemanticMode with
                | AllowAmbiguousBases, _ -> AllowAmbiguousBases
                | _, AllowAmbiguousBases -> AllowAmbiguousBases
                | Explicit, Explicit -> Explicit)
            Explicit
    /// Concatentate a sequence of Dna types into a single, new Dna type.
    /// If any of the sequences allow ambiguous bases, the resulting Dna type will as well.
    let concat (seqs: seq<Dna>) =
        seqs
        |> Seq.map (fun s -> s.arr)
        |> Array.concat
        |> fun d -> Dna(d, false, combinedMode seqs)

    /// Append a Dna type to a second Dna type.
    let append (a: Dna) (b: Dna) =
        Array.append a.arr b.arr
        |> fun d -> Dna(d, false, combinedMode [a; b])
    
    /// Pipelineable call to reverse complement a piece of DNA.
    let revComp (a: Dna) = a.RevComp()

    /// Helper function to return reverse complement depending on a conditional.
    let revCompIf condition (a: Dna) = if condition then a.RevComp() else a

    // These functions are here to enable clients to use functions exclusively in the Dna type domain
    // without needed to manually reach into the DNA type.  Adding aliases here rather than calling
    // the original versions of these functions should ease refactoring of the base functions later.

    let codon2aa (codon: Dna) = codon2aa codon.arr

    /// Translate DNA to AA 3 bases at a time.
    let translate (s: Dna) = translate s.arr

module DnaConstants =
    let codons =
        biolib.codons
        |> Array.map (fun codon -> Dna(codon))