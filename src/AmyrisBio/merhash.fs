namespace Amyris.Bio
/// Rapid searches for matching or repetitive sequences using Mer hashing
module MerHash =
    open Amyris.Bio.biolib
    open System
    open System.IO
    let mutable randBase = 0

    (*
    // Example
    
    // First, choose a sufficiently large hash size to make sure that enough slots will be available. If there are not enough slots available, 
    // the merhash.insert function will loop infinitely, trying to find an empty slot that does not exist...
    // Suggestion for setting hashsize: first prime number after sequence size. For yeast genome: 14999981 works fine . 
    let hash = MerHash(hSize)

    // Pattern for processing a DNA sequence dna (char array)
    // NB - assumes 32 bit mer size, otherwise a bit mask is needed after bit shift

    let init = seq { for i in 0..mer-1 -> uint64(encDNA (dna.[i])) <<< ((mer-i-1)*2) } |> Seq.sum
    let initRC = seq { for i in 0..mer-1 -> uint64(encCompDNA (dna.[i])) <<< ((mer-i-1)*2) } |> Seq.sum

    let totalMerCount = dna.Length-mer+1

    let rec procDNA count i (v:uint64) (v2:uint64) =
        hash.insert(min v v2)
        if count = 0 then ()
        else
            procDNA (count-1) (i+1) 
                ((v <<< 2) + uint64(match dna.[i] with |'G' -> 0uy | 'A' -> 1uy | 'T' -> 2uy | 'C' -> 3uy | _ -> failwithf "Bad base '%c'" (dna.[i]) ))
                ((v2 <<<2) + uint64(match dna.[i] with |'C' -> 0uy | 'T' -> 1uy | 'A' -> 2uy | 'G' -> 3uy | _ -> failwithf "Bad base '%c'" (dna.[i]) ))

    procDNA (totalMerCount-1) (mer) init initRC
    *)


    let encDNA (b : char) =
        match b with
            | 'G' | 'g' -> 0uL
            | 'A' | 'a' -> 1uL
            | 'T' | 't' -> 2uL
            | 'C' | 'c' -> 3uL
            | _ ->
                randBase <- randBase + 1
                (randBase % 4) |> uint64

    let encCompDNA (b : char) =
        match b with
            | 'G' | 'g' -> 3uL
            | 'A' | 'a' -> 2uL
            | 'T' | 't' -> 1uL
            | 'C' | 'c' -> 0uL
            | _ ->
                randBase <- randBase + 1
                (randBase % 4)|>uint64
    /// Encode a whole DNA string to an int
    let encBases (seq : char []) = Array.fold (fun a b -> a*4uL + (encDNA b) ) 0uL seq




    type HState = EMPTY | USED | USEDTWICE | SKIPPED | SKIPPEDTWICE
    type HSearchResult = ONCE | TWICE | NONE

    type MerHash(hSize:int) = class 
        let hash1 = Array.create hSize 0uL
        let hash2 = Array.create hSize EMPTY
        let hFunc (v:uint64) = (v % (uint64 hSize)) |> int32

        do
            //printf "hSize=%d" (uint64 hSize)

            ()

        member x.insert(v:uint64) =
            let rec place i =
                if i = hSize then
                    place 0
                else
                    match hash2.[i] with
                        | EMPTY ->
                            hash1.[i] <-v 
                            hash2.[i] <- USED

                        | USED when hash1.[i] = v ->
                            hash2.[i] <- USEDTWICE
                        | USED when hash1.[i] <> v ->
                            // skip past this bucket
                            hash2.[i] <- SKIPPED 

                            place (i+1)
                        | USEDTWICE when hash1.[i] = v ->
                            ()
                        | USEDTWICE when hash1.[i] <> v ->
                            // skip past this bucket
                            hash2.[i] <- SKIPPEDTWICE
                            place (i+1)
                        | SKIPPED when hash1.[i] = v ->
                            hash2.[i] <- SKIPPEDTWICE

                
                        | SKIPPED when hash1.[i] <> v ->
                            // skip past this bucket
                            place (i+1)
                        | SKIPPEDTWICE when hash1.[i] = v ->
                            // Stop counting at two
                            ()
                        | SKIPPEDTWICE when hash1.[i] <> v ->
                            // skip past this bucket
                            place (i+1)
                        | _ -> failwithf "ERROR: unexpected hash state %A %A\n" (hash2.[i]) (hash1.[i] = v)
            place (hFunc v)


        member x.check(v) =
            let rec search i =
                if hash1.[i] = v then
                    match hash2.[i] with
                        | EMPTY -> NONE
                        | USED -> ONCE
                        | USEDTWICE -> TWICE
                        | SKIPPED -> ONCE
                        | SKIPPEDTWICE -> TWICE
                else match hash2.[i] with
                        | SKIPPED -> search (i+1)
                        | SKIPPEDTWICE -> search (i+1)
                        | _ -> NONE


            search (hFunc v)
    
        member x.save(file:string) =
            use dump = new FileStream(file,FileMode.Create, FileAccess.Write, FileShare.None)  
            use bw = new BinaryWriter(dump)

            bw.Write([| 'M' ; 'H' ; 'S'; 'H' |]) // header
            bw.Write(hSize)
            hash1 |> Array.iter (fun v -> bw.Write(v))
            hash2 |> Array.iter (fun v -> bw.Write(match v with | EMPTY -> 0uy | USED ->1uy | USEDTWICE ->2uy| SKIPPED ->3uy
                                                                | SKIPPEDTWICE -> 4uy))
            ()
        member x.load(file:string) =
            use loadF = new FileStream(file,FileMode.Open, FileAccess.Read, FileShare.None)  
            use br = new BinaryReader(loadF)

            assert(br.ReadChars(4) = [| 'M' ; 'H' ; 'S'; 'H' |])
            let hSize' = br.ReadInt32()
            assert(hSize' = hSize)

            for i in {0..hSize-1} do
                hash1.[i] <- br.ReadUInt64()
            
            for i in {0..hSize-1} do
                hash2.[i] <- match br.ReadByte() with
                                    | 0uy -> EMPTY
                                    | 1uy -> USED
                                    | 2uy -> USEDTWICE
                                    | 3uy -> SKIPPED
                                    | 4uy -> SKIPPEDTWICE
                                    | _ as x -> failwithf "Invalid hash value %d\n" x


    end
