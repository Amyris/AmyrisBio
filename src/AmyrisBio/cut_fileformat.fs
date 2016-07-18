namespace Amyris.IO

/// Codon Usage table file formats
module CodonUsage =
    open System
    open System.Text.RegularExpressions
    open Amyris
    open Amyris.utils
    open Amyris.sgd
    open Amyris.biolib
    open System.Collections.Generic
    open System.IO
    (*

    UUU 26.1(170666)  UCU 23.5(153557)  UAU 18.8(122728)  UGU  8.1( 52903)
    UUC 18.4(120510)  UCC 14.2( 92923)  UAC 14.8( 96596)  UGC  4.8( 31095)
    UUA 26.2(170884)  UCA 18.7(122028)  UAA  1.1(  6913)  UGA  0.7(  4447)
    UUG 27.2(177573)  UCG  8.6( 55951)  UAG  0.5(  3312)  UGG 10.4( 67789)

    CUU 12.3( 80076)  CCU 13.5( 88263)  CAU 13.6( 89007)  CGU  6.4( 41791)
    etc

    *)
    let private loadCodonCommon (strSeq:seq<string>) =
        strSeq
        |> Seq.choose (fun line ->
            if line.Trim().Length = 0 then None
            else
                Some([line.[0..15] ; line.[18..18+15] ; line.[36..36+15] ; line.[54..54+15] ]))
        |> List.concat
        |> List.map (fun col ->
            let m = Regex.Match(col,@"([AUCG][AUCG][AUCG]) *(\d*.\d*)\( *\d*\)")
            if not m.Success then 
                failwithf "ERROR: parsing codon table entry '%s'" col
            else
                m.Groups.[1].Value.Replace("U","T"),float(m.Groups.[2].Value))
        |> Map.ofList

    let loadCodonTableFromString (str:string) =
        str.Split([| '\n' ; '\r' |])
        |> Seq.filter (fun line -> line.Trim() <> "")
        |> loadCodonCommon
                
    let loadCodonTable (path:string) =
        loadCodonCommon (eachLineIn path)
                          
    /// Given an amino acid sequence and codon frequency, pick a dumb optimal codon sequence                          
    let codonOpt (codonTable:Map<string,float>) (prot:string) =
        prot
        |> Seq.map (fun (aa:char) -> // find best codon
            codonTable
            |> Seq.choose (fun kv ->
                if biolib.codon2aa (kv.Key.ToCharArray()) = aa then Some(kv.Value, kv.Key)
                else None)
            |> List.ofSeq
            |> List.sortWith (fun (v1,_) (v2,_) -> compare v2 v1)
            |> List.head
            |> snd)
        |> Seq.concat
        |> Array.ofSeq
        |> arr2seq

    /// computes codon usage for a genome from its REF folder
    /// WARNING: For each gene, computation of codon frequency is stopped as soon as a stop codon is found. 
    /// This may be wrong for mitochondrial genes, which often use TGA for tryptophane..
    /// A feature filter should be provided, for instance (fun f -> true) to use all features or 
    /// (fun f -> f.chr <> 17 && f.featType="ORF" && (not (f.description.ToLower().Contains("dubious")))) 
    /// to only use only non mitochondrial ORFs that are not dubious.
    let computeCodonFrequencies (featureFilter:Feature->bool) (refFolder:String) = 
        let orgName = baseName refFolder
        let (genome, annotation) =
            (Amyris.biolib.readReference (opj refFolder (sprintf "%s.fsa" orgName)),
             loadFeatures (opj refFolder (sprintf "%s_features.tab" orgName)))
        
        // stop codons
        let stops = Set([|"TGA"; "TAA"; "TAG"|])
        // initialize codon frequencies table
        let codonFreqTable = new Dictionary<string,int>()
        for c1 in [|'T';'C';'A';'G'|] do
            for c2 in [|'T';'C';'A';'G'|] do
                for c3 in [|'T';'C';'A';'G'|] do
                    codonFreqTable.[arr2seq [|c1;c2;c3|]] <- 0
        
        // recurse through the sequence
        let rec countCodons(seq:char[]) = 
            if seq.Length < 3 then 
                if seq.Length <> 0 then 
                    printfn "Warning, gene length is not a multiple of 3"
            else
                let codon = seq.[0..2] |> arr2seq
                codonFreqTable.[codon] <- codonFreqTable.[codon] + 1
                if seq.Length > 3 then 
                    if stops.Contains(codon) then
                        printfn "Warning, stop codon found before the end of the sequence"
                    else
                        countCodons seq.[3..]

        // iterate over all genes
        annotation 
        |> Array.filter featureFilter
        |> Array.iter (fun f ->
            printfn "processing %s %s" f.id f.gene
            // grab the underlying coding sequences
            buildGeneStructure annotation f.id
            |> Array.map (fun cds -> // get each CDS sequence
                (if cds.feature.fwd then id else revComp)
                    (genome.[string cds.feature.chr].[cds.feature.l..cds.feature.r]))
            |> Array.concat  //combine all CDS
            |> countCodons) // count the codons

        codonFreqTable

    /// Output the computed codon usage in the proper text format (http://www.kazusa.or.jp/codon/)
    let writeCodonTable (filename:string) (codonFreqTable:Dictionary<string,int>) = 
        /// total number of codons
        let totalCount = float ( codonFreqTable |> Seq.sumBy (fun pk -> pk.Value) )
        /// output one result in the format:  UUU 26.1(170666) 
        let outputOne(codon:string) =
            let count =
                if codonFreqTable.ContainsKey(codon) then codonFreqTable.[codon]
                else 0
            sprintf "%s%5.1f(%6d)"
                (codon.ToUpper().Replace("T","U")) (((float count)*1000.)/totalCount) count
        use outF = new StreamWriter(filename) 
        for c1 in [|'T';'C';'A';'G'|] do
            for c3 in [|'T';'C';'A';'G'|] do
                let line =
                    String.Join("  ", [|'T';'C';'A';'G'|]
                    |> Array.map (fun c2 -> outputOne (arr2seq [|c1;c2;c3|])))
                outF.WriteLine line
            outF.WriteLine ""
        outF.Close()


    /// Output the computed codon usage in the proper text format (http://www.bioinformatics.org/sms2/rev_trans.html)
    /// If provided as None, genetic code will the one in Amyris.biolib.codon2aa (standard code)
    let writeCodonTableToGCG
            (filename:string)
            (geneticCode:(char[]->char) option)
            (codonFreqTable:Dictionary<string,int>) = 
        let gencode = match geneticCode with | Some gc -> gc | None -> codon2aa
        /// total number of codons
        let totalCount = float ( codonFreqTable |> Seq.sumBy (fun pk -> pk.Value) )
        /// total count per amino acid
        let totalCountAA = new Dictionary<string,int>()
        for codon in codonFreqTable do 
            let aa = codon.Key.ToCharArray() |> gencode |> aaLetterToTrigram
            if totalCountAA.ContainsKey(aa) then
                totalCountAA.[aa] <- totalCountAA.[aa] + codon.Value
            else
                totalCountAA.[aa] <- codon.Value
        /// output one result in the format:  UUU 26.1(170666) 
        let outputOne(codon:char []) =
            let codonStr = arr2seq codon
            let count =
                if codonFreqTable.ContainsKey(codonStr) then codonFreqTable.[codonStr]
                else 0
            let aa = codon |> gencode |> aaLetterToTrigram
            let codonAaFrac =
                if totalCountAA.ContainsKey(aa) then (float count)/(float (totalCountAA.[aa]))
                else 1.0

            sprintf "%s     %s%13.2f%10.2f%10.2f"
                aa (codonStr.ToUpper()) (float count) (((float count)*1000.)/totalCount) codonAaFrac

        use outF = new StreamWriter(filename) 
        outF.WriteLine("AmAcid  Codon      Number    /1000     Fraction   ..")
        for c1 in [|'G';'A';'T';'C'|] do
            for c2 in [|'G';'A';'T';'C'|] do
                for c3 in [|'G';'A';'T';'C'|] do
                    outF.WriteLine (outputOne ([|c1;c2;c3|]))
                outF.WriteLine ""
        outF.Close()
    
    /// Details of an individual codon showing rank, relFreq1 (proportion of codons using this codon for this aa)
    /// relFreq2 - same as relFreq1 but for non filtered AAs
    /// freq is original freq per 10,000 codon score
    type CodonLU =
        {codon : char array; 
         aa : char ; 
         freq : float ; 
         relFreq1 : float ; 
         relFreq2 : float ; 
         rank : int}

    /// Codon table optimized for lookup and selection
    type CodonLookup = { byAA: Map<char,CodonLU list>; byCodon : Map<string,CodonLU> } with
        member x.Choose(rInit:float,c:char) = 
            let rec pick (r:float) (l:CodonLU list) =
                match l with
                | [] -> failwithf "ERROR: ran out of codons in CodonLookup::Choose"
                | h::_ when h.relFreq2 > r -> h.codon
                | [h] -> h.codon
                | h::tl -> pick (r-h.relFreq2) tl
            pick
                rInit
                (match x.byAA.TryFind(c) with
                 | Some(x) -> x
                 | None -> failwithf "ERROR: in CodonLookup, unable to find '%c' in byAA" c)

        member x.Codon(s:string) = x.byCodon.TryFind(s)
        member x.AA(aa:char) = x.byAA.TryFind(aa)
    
    /// Prepare a codon lookup table to make it more useful,
    /// calculating relative codon frequency and filtering codons if needed
    /// Takes a  a minFrequency and maxRank to retain and finally a string:float map for codon table,
    /// Rank is the nth codon.  e.g. maxRank of 3 means no 4th ranked codons are retained
    /// regardless of frequency
    let prepCUT minFreq maxRank (codonTable:Map<string,float>) =
        // We want an amino acid -> relative frequency map.
        // Start by getting aa -> codon list

        let epsilon = 1e-06
        // Messy but basically collects serial groups of the same amino acid
        // making a map of amino acid ->  [ codon, codon, codon]
        let aaLookup1 = 
            match 
                codonTable
                |> Seq.map (fun pk ->
                    {aa = codon2aa (pk.Key.ToCharArray());
                     codon = pk.Key.ToCharArray();
                     freq = pk.Value;
                     relFreq1 = -1.0;
                     relFreq2 = -1.0; rank = -666})
                |> List.ofSeq
                |> List.sortWith (fun a b -> compare a.aa b.aa) // collect similar amino acids
                |> List.fold
                    (fun (last,aaLookup:Map<char,CodonLU list>,current) v -> 
                        match last with
                        | None -> Some(v),aaLookup,[v] // start current run of similar amino acids
                        | Some(l) when l.aa = v.aa -> Some(v),aaLookup,v::current // continue current run of similar amino acids
                        | Some(l) -> Some(v),aaLookup.Add(l.aa,current),[v]) // end one run of amino acids, start a new one
                    (None,Map.empty,[])
                with
            | None,lookup,_ -> lookup // Might never have started a group
            | Some(last),lookup,current -> lookup.Add(last.aa,current) // finish up last group

        // Sort the codon choices for each amino acid from least used to most used to make selection later
        // easier.  Eliminate any super low frequency amino acids.  Calculate relative frequency as well 
        // leave in least to most frequent order
        let aaL =
            aaLookup1
            |> Map.map (fun _ codonList -> 
                let firstPass =
                    codonList
                    |> List.sortWith (fun c1 c2 -> compare c1.freq c2.freq) 
                    |> List.mapi (fun i c -> { c with rank = codonList.Length-i})
                    // Nothing removed this pass
                    //|> List.filter (fun c -> (* c.freq > minFreq && *) c.rank <= maxRank) // remove low frequency codons

                let totalFreq =
                    firstPass
                    |> List.map (fun c -> c.freq |> max epsilon)
                    |> List.sum

                let secondPass =
                    firstPass
                    |> List.map (fun c -> {c with relFreq1 = (c.freq |> max epsilon) / totalFreq}) // fill in relative frequency
                    |> List.filter ( fun c -> c.relFreq1 >= minFreq && c.rank <= maxRank)

                let secondTotalFreq =
                    secondPass
                    |> List.map (fun c -> c.freq |> max epsilon)
                    |> List.sum

                // Recalculate the relative frequency from the new group
                secondPass
                |> List.map (fun c -> {c with relFreq2 = (c.freq |> max epsilon) / secondTotalFreq}))
        /// Get all codons and also index by codon -> properties
        let individuals =
            aaL
            |> Seq.map (fun kv -> kv.Value)
            |> Seq.concat
            |> Seq.map (fun codon -> (arr2seq codon.codon),codon)
            |> Map.ofSeq
        {byAA = aaL; byCodon = individuals}     