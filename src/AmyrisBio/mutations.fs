// Tools to handle SNPs, Indels.

namespace Amyris
open utils
open biolib
open System.Text.RegularExpressions
open System.Collections.Generic
/// Extracting mutation data from pileup/VCF format files
module mutations =

    /// Single nucleotide polymorphism | Multiple nucleotide polymorphism | Insertion-deletions
    type MutationType = | SNP | MNP | INDEL   

    type MutationCoverage = {reffwd:int ; refrev:int ; altfwd:int; altrev:int}

    type MutationSamtoolsPval = 
        {strandBias:float ; baseQBias:float ; mapQBias:float ; tailDistBias:float }

    /// parse depth from vcf info field
    let getDepthFromVcfInfo (vcfInfo:string) = 
        let reg = new Regex("DP=([0-9]+)")
        let m = reg.Match(vcfInfo)
        if m.Success then
            Some m.Groups.[1].Value
        else 
            None
    
    /// parse dp4 from vcf info field
    let getDP4FromVcfInfo (vcfInfo:string) = 
        let reg = new Regex("DP4=(\d+),(\d+),(\d+),(\d+)")
        let m = reg.Match(vcfInfo)
        if m.Success then
            Some
                {reffwd= int m.Groups.[1].Value ;
                 refrev= int m.Groups.[2].Value ;
                 altfwd= int m.Groups.[3].Value ;
                 altrev= int m.Groups.[4].Value}
        else 
            None

    /// parse pv4 from vcf info field
    let getPV4FromVcfInfo (vcfInfo:string) = 
        // see http://msdn.microsoft.com/en-us/library/bs2twtah.aspx for help on regex construct
        let floatRegexpGroupMatch = "(\d+(?:\.\d*)?(?:[eE][+-]?\d+)?)"
        let reg =
            new Regex(
                sprintf "PV4=%s,%s,%s,%s"
                    floatRegexpGroupMatch floatRegexpGroupMatch floatRegexpGroupMatch floatRegexpGroupMatch)

        let m = reg.Match(vcfInfo)
        if m.Success then
            Some
                {strandBias= float m.Groups.[1].Value ;
                 baseQBias= float m.Groups.[2].Value ; 
                 mapQBias= float m.Groups.[3].Value ;
                 tailDistBias= float m.Groups.[4].Value}
        else 
            None

    /// parse mapqual from vcf info field
    let getMapQualFromVcfInfo (vcfInfo:string) = 
        let reg = new Regex("MQ=(\d+)")
        let m = reg.Match(vcfInfo)
        if m.Success then
            Some (int m.Groups.[1].Value)
        else 
            None

    /// Create a type for storing mutations, should originate from a VCF reader, or aligner function
    /// Should additional data on the mutation be stored here ? Eg quality measures, coverage, experiments identifying it, etc...
    type VcfMutation = 
        {chr: string ;
         refpos: int ; 
         refseq: char [] ; 
         altseq: char [] ; 
         qual: float ; 
         mapqual: int ; 
         dp4: MutationCoverage ; 
         pv4: MutationSamtoolsPval option}

        /// Get the type of mutation
        member m.MutationType = 
            if m.refseq.Length = m.altseq.Length then
                if m.refseq.Length = 1 then
                    SNP
                else 
                    MNP
            else
                INDEL
        
        /// get the number of reads calling either reference or alternate allele with sufficient base quality
        member m.CallDepth =
            m.dp4.reffwd + m.dp4.refrev + m.dp4.altfwd + m.dp4.altrev
        
        /// get the ratio:
        /// (number of reads calling the alternate allele with sufficient base quality)
        ///  / 
        /// (number of reads calling either reference or alternate allele with sufficient base quality)
        member m.Mutness = 
            float(m.dp4.altfwd + m.dp4.altrev) / (float m.CallDepth)

        /// Many functions from the aligner pipeline do not support MNPs.
        /// This functions creates individual SNPs from one MNP, with the same qual, coverage, etc. 
        member m.DivideMnpIntoSnps() =
            if m.MutationType <> MNP then
                failwith "This method only handles mutations of type MNP"
            else
                Array.zip m.refseq m.altseq 
                |> Array.mapi (fun i (ref, alt) ->
                    {m with refpos = m.refpos + i; refseq = [|ref|]; altseq = [|alt|]})
                |> List.ofArray

        /// Applies the mutation to the provided sequence, and returns the mutated sequence.
        /// The coordinate offset is used when the original sequence is only a subpart of the sequence that was used to provide the SNP coordinates.
        /// For instance, is refpos of the mutation is relative to the chromosome, and originalSeq is only a gene sequence on that chromosome, coordinateOffset
        /// should be the position of the first nucleotide of the gene on the chromosome.
        member m.ApplyToSequence(originalSeq:char [], coordinateOffset:int) = 
            // first, check that the mutation falls within the original sequence, coordinate-wise
            if (coordinateOffset > m.refpos) ||
               ((originalSeq.Length + coordinateOffset) < (m.refpos + m.refseq.Length))
            then
                failwith "Coordinates of the mutation fall outside the provided sequence"

            if not (originalSeq.[(m.refpos-coordinateOffset)..(m.refpos + m.refseq.Length-1-coordinateOffset)] = m.refseq)
            then 
                failwith "Mutation reference sequence does not match the provided reference sequence: %s vs %s"
                    originalSeq.[(m.refpos-coordinateOffset)..(m.refpos + m.refseq.Length-1-coordinateOffset)] m.refseq
            let (before, after) =
                match (m.refpos - coordinateOffset) with 
                | 0 -> ([||], originalSeq.[(m.refseq.Length)..])
                | x when (x = (originalSeq.Length-m.refseq.Length-1)) ->
                    (originalSeq.[..(x-1)], [||])
                | x -> (originalSeq.[..(x-1)], originalSeq.[(x+m.refseq.Length)..])

            after |> Array.append m.altseq |> Array.append before

        override m.ToString() =
            sprintf "%s\t%d\t.\t%s\t%s\t%f\t.\tMQ=%d;DP4=%d,%d,%d,%d%s" 
                m.chr m.refpos (arr2seq m.refseq) (arr2seq m.altseq) m.qual m.mapqual m.dp4.reffwd m.dp4.refrev m.dp4.altfwd m.dp4.altrev 
                (match m.pv4 with Some pv4 -> sprintf ";PV4=%f,%f,%f,%f" pv4.strandBias pv4.baseQBias pv4.mapQBias pv4.tailDistBias| None -> "")
    

    /// Simple reader of VCF file. In case of multiple alternate sequences, only the first one is used.
    /// Warning: we use 0-based coordinates, VCF uses 1-based coordinates
    let loadMutationsFromVCFfile(vcfFileName:string) = 
        eachLineIn vcfFileName 
        |> Seq.filter (fun line -> (line.Length>0) && (line.[0]<>'#')) 
        |> Seq.map (fun line ->
            let fields = line.Split('\t')

            let refseq = fields.[3].ToCharArray()
            refseq
            |> Seq.iter (fun c ->
                if not (isDnaBase c) then failwithf "Inproper VCF entry \n%s" line)

            let altseqvariants = fields.[4].Split(',')

            if altseqvariants.Length > 1 then
                printfn "WARNING: several variants provided in VCF entry \n%s\n. Only the first one is kept."
                    line

            let altseq = altseqvariants.[0].ToCharArray()
            altseq
            |> Seq.iter (fun c ->
                if not (isDnaBase c) then failwithf "Inproper VCF entry %s" line)

            let qual = 
                match fields.[5] with
                | "." | "inf" -> 0.
                | _ -> float fields.[5]

            let dp4 =
                match getDP4FromVcfInfo fields.[7] with 
                | Some dp4 -> dp4
                | None -> {reffwd= 0 ; refrev= 0 ; altfwd= 0 ; altrev= 0}

            let pv4 = getPV4FromVcfInfo fields.[7] 
            let mq =
                match getMapQualFromVcfInfo fields.[7] with 
                | Some mq -> mq
                | None -> 0
            {chr=fields.[0];
             refpos= (int(fields.[1])-1); 
             refseq= refseq; 
             altseq= altseq; 
             qual= qual ; 
             mapqual= mq; 
             dp4= dp4 ; 
             pv4= pv4} ) 

    /// sort a list of VcfMutations 
    let sortVcfMutations (mutations:VcfMutation list) =
        mutations |> List.sortBy (fun vcf -> int(vcf.chr)*10000000 + vcf.refpos)


    /// Some tests on mutations
    let vcfTest() = 
        let originalSeq = "ATCGGTCATGTCGTGCCGGTAATGATGCCTCCTCAACTACGCTACGACT".ToCharArray()
        {chr="1" ; refpos=6 ; refseq= [|'C'|] ; altseq= [|'A'|] ; qual= 100. ; pv4=None; dp4= {reffwd= 0 ; refrev= 0 ; altfwd= 0 ; altrev= 0}; mapqual=20} |> fun m -> printfn "Mutation %A applied to \n%s gives \n%s" m (arr2seq originalSeq) (arr2seq (m.ApplyToSequence(originalSeq, 0)))
        {chr="1" ; refpos=0 ; refseq= [|'A'|] ; altseq= [|'G'|] ; qual= 100. ; pv4=None; dp4= {reffwd= 0 ; refrev= 0 ; altfwd= 0 ; altrev= 0}; mapqual=20} |> fun m -> printfn "Mutation %A applied to \n%s gives \n%s" m (arr2seq originalSeq) (arr2seq (m.ApplyToSequence(originalSeq, 0)))
        {chr="1" ; refpos=47 ; refseq= [|'C';'T'|] ; altseq= [|'A';'G';'G'|] ; qual= 100. ; pv4=None; dp4= {reffwd= 0 ; refrev= 0 ; altfwd= 0 ; altrev= 0}; mapqual=20} |> fun m -> printfn "Mutation %A applied to \n%s gives \n%s" m (arr2seq originalSeq) (arr2seq (m.ApplyToSequence(originalSeq, 0)))
        {chr="1" ; refpos=25047 ; refseq= [|'C';'T'|] ; altseq= [|'A';'G';'G'|] ; qual= 100. ; pv4=None; dp4= {reffwd= 0 ; refrev= 0 ; altfwd= 0 ; altrev= 0}; mapqual=20} |> fun m -> printfn "Mutation %A applied to \n%s gives \n%s" m (arr2seq originalSeq) (arr2seq (m.ApplyToSequence(originalSeq, 25000)))
            
        let mutations = loadMutationsFromVCFfile @"C:\seq\YeastDeltaTest\PROCESSED\test2_new\test2.var.raw.vcf"
        printfn "Loaded mutations %A" mutations

    
    /// Parameters for filtering mutations of a VCF file. Same parameters as the vcfutils.pl program provided with samtools. Entries
    ///  *Q: INT    minimum RMS mapping quality for SNPs
    ///  *d: INT    minimum read depth
    ///  *D: INT    maximum read depth
    ///  *a: INT    minimum number of alternate bases
    ///  *w: INT    SNP within INT bp around a gap to be filtered
    ///  *W: INT    window size for filtering adjacent gaps
    ///  *p1: FLOAT  min P-value for strand bias (given PV4)
    ///  *p2: FLOAT  min P-value for baseQ bias
    ///  *p3: FLOAT  min P-value for mapQ bias
    ///  *p4: FLOAT  min P-value for end distance bias
    ///  *e: FLOAT  min P-value for HWE (plus F<0)
    type vcfFilterParams = 
      { Q:int;
        d:int;
        D:int;
        a:int;
        w:int;
        W:int;
        p1:float;
        p2:float;
        p3:float;
        p4:float;
        e:float }
        
        /// default settings
        /// Q=10 ; d=2 ; D=10000000 ; a=2 ; w=3 ; W=10 ; p1=1e-4 ; p2=1e-100 ; p3=0. ; p4=1e-4 ; e=1e-4
        static member Default =
            {Q=10 ; d=2 ; D=10000000 ; a=2 ; w=3 ; W=10 ; p1=1e-4 ; p2=1e-100 ; p3=0. ; p4=1e-4 ; e=1e-4}


    type VcfFilterStatus = 
        | DEPTH_TOO_LOW
        | DEPTH_TOO_HIGH
        | ALT_DEPTH_TOO_LOW
        | MAPPING_QUAL_TOO_LOW
        | FAILED_PV
        | FAILED_HWE
        | BETTER_NEIGHBORING_MUTATION
        | NEIGHBORING_INDEL
        | PASS
        override st.ToString() =
            match st with 
            | DEPTH_TOO_LOW -> "coverage_too_low" 
            | DEPTH_TOO_HIGH -> "coverage_too_high" 
            | ALT_DEPTH_TOO_LOW -> "alternate_base_coverage_too_low"
            | MAPPING_QUAL_TOO_LOW -> "mapping_quality_too_low" 
            | FAILED_PV -> "failed_Pvalues_thresholds"
            | FAILED_HWE -> "failed_HWE"
            | BETTER_NEIGHBORING_MUTATION -> "better_mutation_in_neighborhood"
            | NEIGHBORING_INDEL -> "indel_in_neighborhood"
            | PASS -> "pass_filter"

    type stagingEntry = 
        { score:int;
          vtype:MutationType;
          mutable fltTag:VcfFilterStatus;
          indelSpan:int;
          vcf:VcfMutation }
 
    /// flush out of the staging lists the elements that are outside the window.
    /// If filt_tag is = 0, than element is output to outputVCF
    let flushStagingL(stagingL:stagingEntry list, currChr:string, currPos:int, param:vcfFilterParams) = 
        let (outputVCF, newStagingL) =
            stagingL
            |> List.partition (fun s ->
                ((s.vcf.chr <> currChr) || (s.vcf.refpos + s.indelSpan + (max param.W param.w) < currPos))) 
            |> fun (a, b) ->
                let a' = a |> List.choose (fun s ->
                    if (s.fltTag = PASS) then Some(s.vcf)
                    else
                        printfn "%s\t%s" (s.fltTag.ToString()) (s.vcf.ToString())
                        None)
                a', b
        (outputVCF, newStagingL)

    let processVcfEntry(stagingL:stagingEntry list, vcf:VcfMutation, param:vcfFilterParams) =
        // Find type of variant
        let vtype = 
            if (vcf.refseq.Length = vcf.altseq.Length ) then
                if vcf.refseq.Length = 1 then SNP // ref and alt have same length = 1 -> SNP
                else MNP // ref and alt have same length > 1 -> MNP
            else 
                INDEL // ref and alt have different length -> indel
    
        // Get total depth and depth of variant from INFO slot
        let dpTot, dpAlt = vcf.CallDepth, vcf.dp4.altfwd + vcf.dp4.altrev
    
        let arePVsWrong() =
            match vcf.pv4 with
            | Some pv4 ->
                pv4.strandBias < param.p1
                || pv4.baseQBias < param.p2
                || pv4.mapQBias < param.p3
                || pv4.tailDistBias < param.p4
            | None -> false
    
        // not implemented
        let isHWEWrong() = false

        // score is defined so that quality has priority over depth of variant (see vcfutils.pl)
        let score = int(vcf.qual) * 1000 + dpAlt  

        let hasStagingListBetterEntry() =
            match vtype with
            | INDEL ->
                // identify SNP and MNP in the staging list that are in the window of influence, and discard them 
                stagingL
                |> List.filter (fun se ->
                    se.vtype <> INDEL
                    && se.fltTag = PASS
                    && (se.vcf.refpos + se.indelSpan + param.w >= vcf.refpos))
                |> Seq.iter (fun se -> (se.fltTag <- NEIGHBORING_INDEL))

                // identify Indels in the staging list that are in the window of influence, and partition them in list of higher and lower score
                let (betterScore, lowerScore) = 
                    stagingL
                    |> List.filter (fun se ->
                        se.vtype = INDEL
                        && se.fltTag = PASS
                        && (se.vcf.refpos + se.indelSpan + param.W >= vcf.refpos))
                    |> List.partition (fun se -> se.score > score)

                // discard those with lower score 
                lowerScore |> Seq.iter (fun se -> se.fltTag <- BETTER_NEIGHBORING_MUTATION)
                // if an entry with better score is found, than our staging list has a better entry
                betterScore.Length > 0
            | _ ->
                // Check if there is a bigger MNP calling the exact same variant, and choose the one with the higher score
                let (betterScore, lowerScore) = 
                    stagingL
                    |> List.filter (fun se ->
                        se.vtype <> INDEL
                        && se.fltTag = PASS
                        && (se.vcf.refpos + se.indelSpan >= vcf.refpos))
                    |> List.partition (fun se -> se.score > score)

                lowerScore |> Seq.iter (fun se -> (se.fltTag <- BETTER_NEIGHBORING_MUTATION))

                // check if there are indels in the staging list that are not flagged within the window of the SNP/MNP
                let indelInWindow =
                    stagingL
                    |> List.filter (fun se ->
                        se.vtype = INDEL
                        && se.fltTag = PASS
                        && (se.vcf.refpos + se.indelSpan + param.w >= vcf.refpos ))

                betterScore.Length > 0 || indelInWindow.Length > 0


        let flt = 
            if (dpTot >= 0 && dpTot < param.d) then DEPTH_TOO_LOW
            elif (dpTot >= 0 && dpTot > param.D) then DEPTH_TOO_HIGH
            elif (dpAlt >= 0 && dpAlt < param.a) then ALT_DEPTH_TOO_LOW
            elif (vcf.mapqual >= 0 && vcf.mapqual < param.Q) then MAPPING_QUAL_TOO_LOW
            elif arePVsWrong() then FAILED_PV
            elif isHWEWrong() then FAILED_HWE
            elif hasStagingListBetterEntry() then
                match vtype with 
                | INDEL -> BETTER_NEIGHBORING_MUTATION
                | SNP|MNP -> NEIGHBORING_INDEL
            else PASS
    
        stagingL@[{score=score; vtype=vtype; fltTag = flt; indelSpan = (vcf.refseq.Length - 1); vcf=vcf}] 
        

    /// filters variants using the same criteria as vcfutils.pl, provided by samtools
    let filterVcf
            (stagingL:stagingEntry list, outputL:VcfMutation list, param:vcfFilterParams)
            (vcf:VcfMutation)
        : (stagingEntry list*VcfMutation list*vcfFilterParams) =

        if (vcf.refseq = [|'N'|] || vcf.altseq = [|'.'|]) then
            (stagingL, outputL, param)
        else 
            // clear out of range elements in staging list
            let (newOutputL, flushedStagingL) = flushStagingL(stagingL, vcf.chr, vcf.refpos , param)
            let processedStaging = processVcfEntry(flushedStagingL, vcf, param)
            (processedStaging, List.append outputL newOutputL, param)

    /// Filters the mutations stored in the given VCF file, using parameters specified with a vcfFilterParams object.
    let filterVcfList  (param:vcfFilterParams) (mutations:VcfMutation list) = 
        let (stagingFinal, outputTmp, _) = 
            mutations
                |> sortVcfMutations                                
                |> List.fold filterVcf ([], [], param) 
        outputTmp@(fst (flushStagingL(stagingFinal, "", 1000000000, param)) )
