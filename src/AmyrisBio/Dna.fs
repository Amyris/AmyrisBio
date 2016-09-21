namespace Amyris.Dna
open Amyris.Bio
open System
open System.Collections
open System.Collections.Generic
open biolib
open utils

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
type Dna private (asArray:char [], rc: Dna option) =
    // TODO: refactor this class to use ImmutableArray rather than Array.  At the moment
    // clients could shoot themselves in the foot by mutating the array contents, or feeding
    // it to a function that mutates it.

    // Private mutable values to permit memoization.
    let mutable asString: string option = None
    let mutable revCompPartner: Dna option = rc

    new(sequence, ?validate) =
        let doValidate = match validate with | Some(v) -> v | None -> true
        if doValidate then
            let badChars =
                sequence
                |> Seq.filter (fun c -> not (isDnaBaseStrict c))
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
        Dna(asArray.[start..finish], None)

    interface IEnumerable<char> with
        member x.GetEnumerator() = (Seq.cast<char> asArray).GetEnumerator()
    interface IEnumerable with
        member x.GetEnumerator() = asArray.GetEnumerator()

    // Implement equality between similiar types
    override x.Equals(other) =
        match other with
        | :? Dna as o -> x.arr = o.arr
        | _ -> false

    // TODO: might want to implement hash code by hashing the sequence

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

    member x.Length = asArray.Length

    /// Get the reverse compliment of this sequence.
    member x.RevComp() =
        match revCompPartner with
        | Some(rc) -> rc
        | None ->
            let rcArr = revComp asArray
            let rcDna = new Dna(rcArr, Some(x))
            revCompPartner <- Some(rcDna)
            rcDna

[<AutoOpen>]
module DnaOps =
    let concat (seqs: seq<Dna>) =
        seqs
        |> Seq.map (fun s -> s.arr)
        |> Array.concat
        |> fun d -> Dna(d, false)