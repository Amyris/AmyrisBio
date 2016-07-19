namespace Amyris.Bio

///
/// Simulation of DNA / epcr operations in silico
///
module dnaops =
    open biolib
    open utils
    open System

    /// Contains the sequence around and inside an ePCR result
    /// fwdStart = primerOffset into product
    /// template = product and surrounding DNA
    /// Product gives just the dna product
    /// chrFrom and chrTo are the first and last bases of the product relative to original
    /// primers, i.e chrFrom > chrTo if result was flipped around
    type EPCRResult =
        {template : string ; 
         chr: string ; 
         chrFrom : int ; 
         chrTo : int ; 
         chrFwd : bool ; 
         fwdStart : int ; 
         revEnd : int } with
        /// Actual PCR product
        member x.Product = x.template.[x.fwdStart..x.revEnd] 

    /// Locate a primer in the genome and get a window around it.
    /// Only one match per chromosome is reported
    [<Obsolete>]
    let findPrimer (genome:Map<string,char array>) (p:char array) upMargin downMargin=
        let ps = (arr2seq p).ToUpper()
        let psr = (revComp p |> arr2seq).ToUpper()
        seq {for pk in genome do
                let chr = pk.Value |> arr2seq
                match chr.IndexOf(ps)  with
                | -1 ->
                    match chr.IndexOf(psr) with
                    | -1 -> ()
                    | i -> 
                        let substringStart = max (i-downMargin-psr.Length) 0
                        yield
                           (pk.Key,
                            (i+psr.Length-1),
                            false,
                            chr.Substring(
                                    substringStart,
                                    min (downMargin+upMargin+psr.Length) (chr.Length-substringStart))
                                .ToCharArray()
                            |> revComp |> arr2seq)
                | i ->
                    let substringStart = max (i-upMargin) 0
                    yield
                       (pk.Key,
                        i,
                        true,
                        chr.Substring(
                            substringStart,
                            min (upMargin+downMargin) (chr.Length-substringStart)).ToUpper())
        }

    /// Locate a primer in the genome and get a window around it.
    /// Allows multiple matches on the same chromosome to be reported/
    let findPrimersAllowingMultipleMatches
            (genome:Map<string,char array>)
            (p:char array)
            upMargin
            downMargin =
        let ps = (arr2seq p).ToUpper()
        let psr = (revComp p |> arr2seq).ToUpper()

        /// Returns: chrName, chrPos, chrFwd, DNA region in same orientation as primer with
        ///margin on each side, primer start position on DNA region
        let rec generatePrimerMatches (chrList:(string*string)list) offset isForward = seq {
            match chrList with 
            | (chrName, chrSeq)::tail -> 
                if isForward then 
                    match chrSeq.[offset..].IndexOf(ps) with 
                    | -1 -> yield! generatePrimerMatches chrList 0 false
                    | i ->  let substringStart = max (offset+i-upMargin) 0
                            let primerStartInSubstring = min (offset+i) upMargin
                            yield
                                chrName,
                                offset+i,
                                true,
                                chrSeq.Substring(
                                        substringStart,
                                        min (primerStartInSubstring+ps.Length+downMargin) (chrSeq.Length-substringStart))
                                    .ToUpper(),
                                primerStartInSubstring
                            if offset+i+1+ps.Length < chrSeq.Length then
                                yield! generatePrimerMatches chrList (offset+i+1) isForward
                else 
                    match chrSeq.[offset..].IndexOf(psr) with 
                    | -1 -> yield! generatePrimerMatches tail 0 true
                    | i ->  let substringStart = max (offset+i-downMargin) 0
                            let primerStartInSubstringRev = min (offset+i) downMargin
                            let primerStartInSubstring = min (chrSeq.Length-1-(offset+i+psr.Length-1)) upMargin
                            yield
                                chrName,
                                (offset+i+psr.Length-1),
                                false,
                                (chrSeq.Substring(
                                        substringStart,
                                        min (primerStartInSubstringRev+psr.Length+upMargin) (chrSeq.Length-substringStart))
                                    .ToCharArray() |> revComp |> arr2seq),
                                    primerStartInSubstring
                            if offset+i+1+psr.Length < chrSeq.Length then
                                yield! generatePrimerMatches chrList (offset+i+1) isForward
            | [] -> ()
            }
        generatePrimerMatches
            (genome |> Seq.map (fun pk -> (pk.Key, pk.Value |> arr2seq)) |> List.ofSeq)
            0
            true
  
    /// Given a genome map and fwd/rev primers, return all possible PCR products. 
    let ePCRAllowMultipleProducts (genome:Map<string,char array>) (fwd:string) (rev:string) (maxPcrLength:int) =
        let revRev = rev.ToCharArray() |> revComp |> arr2seq

        let hits =
            findPrimersAllowingMultipleMatches genome (fwd.ToCharArray()) 1000 (1000+maxPcrLength)
            |> Seq.choose (fun (chr,cOff,cFwd,s,primerOffset) -> 
                // commenting out this following line. Coordinate of fwd primer is now 
                // directly grabbed from the findPrimersAllowingMultipleMatches function. 
                // In some rare cases, it could happen that the primer had a 2nd match in 
                // the upstream region, which resulted in getting the wrong primer Offset
                //let primerOffset = s.ToUpper().IndexOf(fwd.ToUpper())
                assert(s.ToUpper().[primerOffset..primerOffset+fwd.Length-1] = fwd.ToUpper())
                match s.IndexOf(revRev.ToUpper()) with
                | -1 -> None
                | i when i > (100+primerOffset) -> 
                    Some({template = s ;
                          chr = chr; 
                          chrFrom = cOff ; 
                          chrTo =
                              (if cFwd then cOff+i-primerOffset+rev.Length-1
                               else cOff-(i-primerOffset+rev.Length-1));
                          chrFwd = cFwd;
                          fwdStart = primerOffset;
                          revEnd = i+rev.Length-1})
                | _ -> None)
            |> List.ofSeq

        // handle the special case where fwd = rev. Remove duplicate PCR products from hits.
        if fwd <> rev then
            hits
        else
            let rec dedup toProcess processed = 
                match toProcess with 
                | [] -> processed
                | pcr::tail -> 
                    if processed |> List.exists (fun p ->
                        p.chr=pcr.chr && 
                        (p.chrFrom=pcr.chrTo && p.chrTo=pcr.chrFrom || 
                         p.chrFrom=pcr.chrFrom && p.chrTo=pcr.chrTo ))
                    then
                        dedup tail processed
                    else
                        dedup tail (pcr::processed)
            dedup hits []

        
    /// Given a genome map and fwd/rev primers, see if there is
    /// a single unique product. Same as ePCR except the maximal length can be tuned
    let ePCR2 (genome:Map<string,char array>) (fwd:string) (rev:string) (maxPcrLength:int) =
        let hits = ePCRAllowMultipleProducts genome fwd rev maxPcrLength
        match hits with
        | [] -> None
        | [ single ] -> Some(single)
        | _ -> None


    /// Given a genome map and fwd/rev primers, see if there is
    /// a single unique product.  
    /// PCR length limit at 10kb 
    let ePCR (genome:Map<string,char array>) (fwd:string) (rev:string) =
        ePCR2 genome fwd rev 10000

