namespace Amyris
/// routines for using some Saccharomyces Genome Database formats,
module sgd =

    open System.Collections.Generic
    open System.IO
    open System
    open System.Text
    open utils
    open biolib

    let lookup =
       [("ref|NC_001133|",1); ("ref|NC_001134|",2) ; ("ref|NC_001135|",3); ("ref|NC_001136|", 4);
        ("ref|NC_001137|",5); ("ref|NC_001138|",6); ("ref|NC_001139|",7); ("ref|NC_001140|",8); 
        ("ref|NC_001141|",9); ("ref|NC_001142|",10); ("ref|NC_001143|",11); ("ref|NC_001144|",12);
        ("ref|NC_001145|",13); ("ref|NC_001146|",14); ("ref|NC_001147|",15);("ref|NC_001148|",16);
        ("ref|NC_001224|",17) ]

    let chrLookup = Map.ofSeq lookup 

    /// Should be replaced by Amyris.biolib.readReference
    [<Obsolete>]
    let readReference (path:string) =
        use f = new StreamReader(path)
        //let chrs = File.ReadAllText(path).Split([|'>'|]) |> Seq.filter (fun seq -> seq.Length > 0 )
        let chrs = (f.ReadToEnd()).Split([|'>'|]) |> Seq.filter (fun seq -> seq.Length > 0 )
        let splitEntry (s:string) =
            let firstEOL = s.IndexOf('\n')
            let header = s.Substring(0,firstEOL)
            let firstName = header.Split([|' '|]).[0]
            let seq = s.Substring(firstEOL).ToCharArray() // |> Array.filter (fun x -> x <> '\n')
            let rec nlCount i c =
                if i = seq.Length then c 
                elif seq.[i] = '\n' then nlCount (i+1) (c+1) else nlCount (i+1) c
            let nl = nlCount 0 0
            let newSeq = Array.create (seq.Length-nl) ' ' 
            let rec transfer i j =
                if i = seq.Length then ()
                else
                    if seq.[i] = '\n' then transfer (i+1) j
                    else
                        newSeq.[j] <- seq.[i]
                        transfer (i+1) (j+1)
            transfer 0 0
            (firstName,newSeq)
        
        f.Close()
        //let chrSeq = new HashMultiMap<string,char [] >()  
        
        let chrSeq = Seq.map (splitEntry) chrs  |> Map.ofSeq //Seq.iter (fun (name,seq) -> chrSeq.Add(name,seq))
        chrSeq
    /// Map integer chromosome ids to names 
    let chromNames2Ints =
       [("ref|NC_001133|",1); ("ref|NC_001134|",2); ("ref|NC_001135|",3); ("ref|NC_001136|", 4);
        ("ref|NC_001137|",5); ("ref|NC_001138|",6); ("ref|NC_001139|",7); ("ref|NC_001140|",8); 
        ("ref|NC_001141|",9); ("ref|NC_001142|",10); ("ref|NC_001143|",11); ("ref|NC_001144|",12);
        ("ref|NC_001145|",13); ("ref|NC_001146|",14); ("ref|NC_001147|",15);("ref|NC_001148|",16);
        ("ref|NC_001224|",17)]

    // Get reference sequences
    // let chrSeq = readReference yeastRef
    //  let chrSeq = readReference yeastRef

    /// chromosome number to name lookup
    let chrLookupInt2Str = chromNames2Ints |> List.map (fun (name,id) -> (id,name)) |> Map.ofSeq // new HashMultiMap<int,string>()
    //List.iter (fun (name,id) -> chrLookupInt2Str.Add( id,name )) chromNames2Ints 

    /// Take a chromosome name that is either an int or a genbank ref and return the canonical int
    let canonicalChromName (x:string) = if x.StartsWith("ref|NC") then chrLookup.[x] else int(x)
    
    /// Make a new hash multi map keyed on integers, renaming
    /// the chromosomes if we have legacy genbank NA identifiers we'd
    /// like to get rid of
    let renameChroms (chroms : Map<string,char []>) =
        //let n = new HashMultiMap<int,int>()
        //let x = n.
        chroms
        |> Seq.map (fun p ->
            if p.Key.StartsWith("ref|NC") then
                (chrLookup.[p.Key], p.Value)
            else
                (int(p.Key), p.Value))
        |> Map.ofSeq
    //--------------------------------------
 
    /// One yeast genome feature
    type Feature =
        {chr : int ; 
         l : int ;
         r : int ;
         description : string ;
         featType : string ;
         fwd : bool ;
         gene : string ;
         sysName : string ;
         id : string}
       
    let cmp a b =
    
        let al = a.chr * 10000000 + a.l
        let ar = a.chr * 10000000 + a.r
        let bl = b.chr * 10000000 + b.l
        let br = b.chr * 10000000 + b.r
        if ar < bl then -1
        else if al > br then 1
        else 0               
    
    /// Parses one line of SGD tab-delimited feature files
    let parseSgdFeatureLine (line:string) =
        let cols = line.Split([|'\t'|],StringSplitOptions.None )
        try
            if cols.Length >= 11 && cols.[8] <> "" then
                let chr = int ( cols.[8])
                let l = int (cols.[9]) // -1
                let r = int (cols.[10]) // -1  NEED -1 for 1 based coordinate system
                let sysName = if cols.[3] <> "" then cols.[3] else cols.[6] // systematic name
                // CDFs and ORFs have the gene name in different columns, and we really want to group this together
                // so select it accordingly
                // S000000069|ORF|Verified|YAR014C|BUD14||chromosome 1||1|168866|166743|C||2004-07-20|2004-01-27|1996-07-31|Protein involved in bud-site selection, Bud14p-Glc7p complex functions as a cortical regulator of dynein; diploid mutants display a random budding pattern instead of the wild-type bipolar pattern
                // S000031457|CDS|        |       |     ||YAR014C     ||1|168866|166743|C||2004-07
                let gene =
                    if cols.[6].StartsWith("chromosome") then 
                        if cols.[4] <> "" then cols.[4]
                        else cols.[3]
                    else cols.[6] // changed 3->4

                let d = gene + "/" + cols.[4] + "/" + (if cols.Length >= 16 then cols.[15] else "") + "/"
                    
                let fwd =
                    match cols.[11] with 
                    | "W" -> true
                    | "C" -> false
                    | _ ->
                        if cols.[1] <> "ARS" then
                            printfn "WARNING: Feature strand not understood for feature %s: '%s'. Using start and end coordinates to determine strand."
                                sysName cols.[11]
                        l <= r

                if d.Trim() <> "" then
                    let f =
                       {chr = chr ;
                        l = min l r; 
                        r = max l r ; 
                        description = d.Trim() ; 
                        featType = cols.[1] ; 
                        gene = gene.ToUpper() ; 
                        fwd = fwd ; 
                        sysName = sysName.ToUpper() ; 
                        id = cols.[0]}
                    //printfn "Adding %s : %s" gene f.description
                    //printfn "Appending %s" line
                    Some(f)
                else    
                    //nqprintfn "Discarding %s gene=%s cols6=%s cols3=%s\n" line gene cols.[6] cols.[3]
                    None
            else
                None
        with x -> 
            //if cols.Length >8 then 
                //printfn "Skipping %s"  line
            None
                
    /// Load up chromosome features for yeast 
    let loadFeatures path =
        (*
        1.   Primary SGDID (mandatory)
        2.   Feature type (mandatory)
        3.   Feature qualifier (optional)
        4.   Feature name (optional)
        5.   Standard gene name (optional)
        6.   Alias (optional, multiples separated by |)
        7.   Parent feature name (optional)
        8.   Secondary SGDID (optional, multiples separated by |)
        9.   Chromosome (optional)
        10.  Start_coordinate (optional)
        11.  Stop_coordinate (optional)
        12.  Strand (optional)
        13.  Genetic position (optional)
        14.  Coordinate version (optional)
        15.  Sequence version (optional)
        16.  Description (optional)
        *)
    
        //let final = eachLineIn path (*SGDFeaturePath*) |> Seq.fold (addLine) featStore  |> Array.ofSeq
        let final = File.ReadAllLines( path) (*SGDFeaturePath*) |> Array.choose (parseSgdFeatureLine) 
    
        // let final = eachLineIn SGDFeaturePath |> Seq.fold (addLine) featStore 
        Array.Sort(final,cmp)
        final 


    /// Convert a Feature to a string, for saving in a tab delimited text file.
    let feature2fileline (f:Feature) = 
        (*
            1.   Primary SGDID (mandatory)
            2.   Feature type (mandatory)
            3.   Feature qualifier (optional)
            4.   Feature name (optional)
            5.   Standard gene name (optional)
            6.   Alias (optional, multiples separated by |)
            7.   Parent feature name (optional)
            8.   Secondary SGDID (optional, multiples separated by |)
            9.   Chromosome (optional)
            10.  Start_coordinate (optional)
            11.  Stop_coordinate (optional)
            12.  Strand (optional)
            13.  Genetic position (optional)
            14.  Coordinate version (optional)
            15.  Sequence version (optional)
            16.  Description (optional)
        *)
        let cols =
          [|f.id ;
            f.featType; 
            "" ; 
            f.sysName ; 
            f.gene ; 
            "" ; 
            sprintf "chromosome %d" f.chr ;
            "" ; 
            sprintf "%d" f.chr ; 
            sprintf "%d" (if f.fwd then f.l else f.r) ; 
            sprintf "%d" ( if f.fwd then f.r else f.l) ;  
            (if f.fwd then "W" else "C");
            "" ; 
            "2/17/1970" ; 
            "2/17/1970" ; 
            f.description|]
        String.Join("\t", cols)

    /// Write a sequence of features to a file 
    let dumpFeatures (features:seq<Feature>) (fileName:string) =
        use outF = new StreamWriter(fileName)
        features |> Seq.iter (fun f -> outF.Write(sprintf "%s\n" (feature2fileline f)))

    /// Write a sequence of features to a string, with tab delimited fields 
    let dumpFeaturesToString (features:Feature seq) = 
        let sb = StringBuilder()
        features |> Seq.iter (fun f -> (sb.Append(feature2fileline f)|> ignore); (sb.Append("\n") |> ignore)) 
        sb.ToString()


    let rec chop (arr : Feature []) l r f =

        if r-l <2 then
            l
        else
            let x = (l+r) / 2
            let mid = if x > l then x else (l+1)
            //printfn "Chop l=%d mid=%d r=%d  (%d->%d)" l mid r arr.[mid].l arr.[mid].r
            match (cmp arr.[mid] f) with
                | -1 -> chop arr mid r f
                | 1 -> chop arr l mid f
                | 0 -> mid
                | _ -> failwith "Illegal comparison result"
            
    let findFeats (featStore : Feature []) (feat:Feature) =
        let start = chop featStore 0 (featStore.Length-1) feat
        let overlaps =
            [for i in (max 0 (start-100)) .. (min (featStore.Length-1) (start+100)) do 
                if featStore.[i].l <= feat.l && featStore.[i].r >= feat.l  && featStore.[i].chr = feat.chr then
                    yield featStore.[i]] 

        let featArray = Array.ofList overlaps                        
        if featArray.Length > 0 then featArray
        else
            // If nobody overlapped us strictly, still grab the nearest location and then
            // any adjacent features with the same coordinates
            // printf "start=%d\n" start
            let a,b =
                match cmp feat featStore.[start] with
                | -1 -> // start is to right of sought coordinate
                    if start > 0 then
                        featStore.[start-1],featStore.[start]
                    else
                        featStore.[start],featStore.[start]
                        

                | 1 -> // start is to left of sought coordinate
                    if start < featStore.Length-1 then
                        featStore.[start],featStore.[start+1]
                    else
                        featStore.[start],featStore.[start]
                | _ -> failwith "impossible case in findFeats"
            
            if b.chr <> feat.chr && a.chr<> feat.chr then
                failwithf "ERROR: in findFeats, a.chr=%d  b.chr=%d and feat.chr=%d"
                    a.chr b.chr feat.chr
            //assert(b.chr = feat.chr || a.chr = feat.chr) // One must be on correct chromosome
        
            // Choose final feature based on proximity (and ensure on correct chromosome to handle end cases)
            let f =
                if a = b then a else
                if a.chr <> feat.chr then b
                else if b.chr <> feat.chr then a else
                    let aDist = abs (feat.l - a.r)
                    let bDist = abs (b.l - feat.r)
                    //assert(aDist >=0)
                    //assert(bDist >=0)
                    if aDist < bDist then a else b

            featStore |> Array.filter (fun f2 -> f2.chr = f.chr && f.l = f2.l && f.r = f2.r)
   
   
    /// Find gene nearest to a chromosome/position and report the gene and position of base rel to gene
    let nearestGene features c p =
        let f =
           {chr = c ; 
            l = p; 
            r = p ; 
            description = "" ; 
            featType = "" ; 
            gene = "" ; 
            fwd = true ; 
            sysName = "" ; 
            id = ""}
        // Find all matching features and take the one with the longest description
        let res,_ =
            findFeats features f
            |> Array.fold 
                (fun (best,length) f -> 
                    if  (best.featType = "ORF"  && f.featType = "ORF"  && f.description.Length > length) ||
                        (best.featType <> "ORF"  && f.featType <> "ORF"  && f.description.Length > length) ||
                        (best.featType <> "ORF"  && f.featType = "ORF"  )
                    then 
                        (f,f.description.Length) 
                    else 
                        (best,length))
                (f,0)

        // Locate mutation relative to ORF.. 100,000 + downstream, -ve upstream                        
        let featRelPos = 
            if res.fwd then
                if p < res.l then p-res.l
                else if p > res.r then
                    p-res.r + 100000
                else p-res.l
            else
                if p > res.r then res.r-p
                else if p < res.l then
                    res.l-p + 100000
                else
                    res.r-p     
        res,featRelPos


    /// Loads chromosome sequences and features from a reference folder 
    /// Warning, not all information is loaded for features..
    let loadGenomeFolder(strainFolder:string) = 
        if not (Directory.Exists strainFolder) then
            failwithf "ERROR: strain genome folder %s does not exist" strainFolder 
    
        let strainId = baseName strainFolder
        let (fastaFile, featureFile) = ((opj strainFolder strainId+".fsa"), (opj strainFolder strainId+"_features.tab"))
        if not (File.Exists fastaFile) then
            failwithf "ERROR: strain fasta file %s does not exist" fastaFile 
        if not (File.Exists featureFile) then
            failwithf "ERROR: strain feature %s does not exist" featureFile 

        (Amyris.biolib.readReference fastaFile,  loadFeatures featureFile)


    /// phase uses the GFF3 convention, ie nb of bases to remove before reaching the first codon of the CDS
    type CdsInGene = 
        {phase:int ; feature:Feature}

    /// Retrieves the CDSs corresponding to a gene identified by its ORF id, and returns them in their 
    /// respective position in the gene (from gene start to gene end) with information on their phase. 
    /// CDS from genes are found using the feature.gene slot, which should be equal to the gene feature.sysName slot
    /// This means there should not be any duplication of systematic names of ORFs.
    let buildGeneStructure (features:Feature []) (geneId:string) = 
        let geneL = features |> Seq.filter (fun f -> f.id=geneId && f.featType<>"CDS")
        if Seq.length geneL > 1 then failwithf "ERROR: Several features with same ID %s" geneId
        let geneFeat = Seq.head geneL
        let cdsL =
            features
            |> Seq.filter (fun f -> f.featType="CDS" && f.gene=geneFeat.sysName )
            |> List.ofSeq
            |> List.sortBy (fun f -> if geneFeat.fwd then f.l else -1 * f.l )

        // various checks on CDS
        let rec _buildGeneStructure (geneStructure:CdsInGene []) cdsList prevCoord phaseForNextCds = 
            match cdsList with 
            | [] -> geneStructure
            | f::_ ->   // first perform some checks on the CDS (included in the ORF, non overlapping)
                if (f.l < geneFeat.l ||
                    f.r > geneFeat.r ||
                    f.chr <> geneFeat.chr ||
                    (f.l < prevCoord && f.fwd) ||
                    (f.r > prevCoord && (not f.fwd)) ||
                    f.fwd <> geneFeat.fwd)
                then 
                    failwithf "ERROR: improper CDS structure for gene %A. CDS triggering the error: %A" geneFeat f
                else
                    let nextPhase =
                        match ((f.r - f.l + 1 - phaseForNextCds) % 3) with 
                        | 1 -> 2
                        | 2 -> 1
                        | 0 -> 0
                        | _ -> failwith "ERROR: cannot occur"
                    _buildGeneStructure
                        (Array.append geneStructure [|{phase=phaseForNextCds; feature=f}|])
                        (List.tail cdsList)
                        (if geneFeat.fwd then f.r else f.l)
                        nextPhase

        _buildGeneStructure [||] cdsL (if geneFeat.fwd then geneFeat.l else geneFeat.r) 0
    
    /// More elaborated store of feature, with several indexes to retrieve Features by type, names, etc...
    /// Use the GetIndicesBy... to retrieve the indices of features. 
    /// Indices are returned as int list (previous idea: returned as Set so that they can be easily combined, but set operations to combine criteria are 
    /// not as efficient as doing a first filter using the index for the most stringent filter (eg search for a feature by name), then filter down for 
    /// other criteria using a simple Array.filter (to keep only CDS for instance ) )
    /// Designed to be used for QUERYING features, not modifying them
    type FeatureStore(feats:Feature [])= 
        let geneIndex = new Dictionary<string,int list>()
        //let featTypeIndex = new Dictionary<string,Set<int>>()
        //let sysNameIndex = new Dictionary<string,Set<int>>()
        let idIndex = new Dictionary<string,int list>()
        //let allIndex = Set([0..feats.Length-1])
        do
            // create indexes
            let addElementToIndex (index:Dictionary<'T,int list>) el i =
                if index.ContainsKey(el) then
                    index.[el] <- i::index.[el]
                else
                    index.[el] <- [i]
            feats
            |> Array.iteri (fun i f ->
                addElementToIndex geneIndex f.gene i
                //addElementToIndex featTypeIndex f.featType i
                //addElementToIndex sysNameIndex f.sysName i
                addElementToIndex idIndex f.id i ) 

        member x.Features = feats
        //member x.AllIndex = allIndex
        member x.GetIndicesByGene(gene:string) = geneIndex.[gene] 
        //member x.GetIndicesByFeatType(feattype:string) = featTypeIndex.[feattype] 
        //member x.GetIndicesBySysName(sysname:string) = sysNameIndex.[sysname] 
        member x.GetIndicesById(id:string) = idIndex.[id] 
        member x.GetFeatureSubset(indexes:seq<int>) = 
            [|for i in indexes -> x.Features.[i] |]

    /// Retrieves the CDSs corresponding to a gene identified by its ORF id, and returns them in their 
    /// respective position in the gene (from gene start to gene end) with information on their phase. 
    /// CDS from genes are found using the feature.gene slot, which should be equal to the gene feature.sysName slot
    /// This means there should not be any duplication of systematic names of ORFs.
    /// This version uses FeatureStore structure, and is much more efficient then the regular buildGeneStructure when number of features > several hundreds (rough assessment)
    let buildGeneStructureWithFeatStore (features:FeatureStore) (geneId:string) = 
        let geneL =
            features.GetFeatureSubset(features.GetIndicesById(geneId))
            |> Array.filter (fun f -> f.featType <> "CDS")

        if Array.length geneL > 1 then failwithf "ERROR: Several features with same ID %s" geneId

        let geneFeat = geneL.[0]
        let cdsL = 
            features.GetFeatureSubset(features.GetIndicesByGene(geneFeat.sysName)) 
            |> List.ofArray 
            |> List.filter (fun f -> f.featType = "CDS")
            |> List.sortBy (fun  f -> if geneFeat.fwd then f.l else -1 * f.l ) 

        // various checks on CDS
        let rec _buildGeneStructure (geneStructure:CdsInGene []) cdsList prevCoord phaseForNextCds = 
            match cdsList with 
            | [] -> geneStructure
            | f::_ -> // first perform some checks on the CDS (included in the ORF, non overlapping)
                if (f.l < geneFeat.l ||
                    f.r > geneFeat.r ||
                    f.chr <> geneFeat.chr || 
                    (f.l < prevCoord && f.fwd) || 
                    (f.r > prevCoord && (not f.fwd)) || 
                    f.fwd <> geneFeat.fwd)
                then 
                    failwithf "ERROR: improper CDS structure for gene %A. CDS triggering the error: %A" geneFeat f
                else
                    let nextPhase =
                        match ((f.r - f.l + 1 - phaseForNextCds) % 3) with 
                        | 1 -> 2
                        | 2 -> 1
                        | 0 -> 0
                        | _ -> failwith "ERROR: cannot occur"
                    _buildGeneStructure
                        (Array.append geneStructure [|{phase=phaseForNextCds; feature=f}|])
                        (List.tail cdsList)
                        (if geneFeat.fwd then f.r else f.l)
                        nextPhase

        _buildGeneStructure [||] cdsL (if geneFeat.fwd then geneFeat.l else geneFeat.r) 0



      //////////////////////////
     // -- export in GFF3 -- //
    //////////////////////////

    /// Output features in gff3, and handle properly gene structure
    /// GFF3 format: http://www.sequenceontology.org/gff3.shtml
    let features2gff3stream (featuresArr:Feature []) (gff3stream:StreamWriter) =
        // ,=; characters must be URL-escaped in GFF3
        let urlEscape (inputString:string) =
            inputString.Replace("%", "%25").Replace(",", "%2C").Replace("=", "%3D").Replace(";", "%3B").Replace("\t","%09")
        
        let features = List.ofArray featuresArr
        let featureStore = FeatureStore featuresArr

        gff3stream.Write "##gff-version 3\n"
        // 1- retrieve genes and write them to gff3
        let genes, notgenes =
            features
            |> List.partition (fun f -> f.featType = "ORF" || f.featType = "transposable_element_gene")

        for f in genes do    
            let gffFields =
                [|string f.chr; 
                  "."; 
                  (if f.featType="ORF"then "gene" else f.featType); 
                  string (f.l+1); 
                  string (f.r+1); 
                  "." ; 
                  (if f.fwd then "+" else "-") ; 
                  "." ; 
                  (sprintf "ID=%s;Name=%s;Note=%s"
                    (urlEscape f.sysName) (urlEscape f.gene) (urlEscape f.description))
                |]
            gff3stream.Write (sprintf "%s\n" (String.Join("\t", gffFields)) )
        
        // 2- retrieve CDS underlying the ORFs and write them to gff3
        for f in genes do
            let geneStruct = buildGeneStructureWithFeatStore featureStore f.id
            geneStruct 
            |> Array.iteri (fun i cds ->    
                let gffFields =
                    [|string cds.feature.chr; 
                      "."; 
                      "CDS"; 
                      string (cds.feature.l+1) ; 
                      string (cds.feature.r+1); 
                      "." ; 
                      (if cds.feature.fwd then "+" else "-") ; 
                      string cds.phase ; 
                      (sprintf "ID=%s_CDS%s;Parent=%s;Note=%s"
                        (urlEscape cds.feature.sysName) (string i) (urlEscape f.sysName) (urlEscape cds.feature.description))
                    |]

                gff3stream.Write (sprintf "%s\n" (String.Join("\t", gffFields))) 
                )

        // 3- retrieve CDS that are not parented to any ORF
        notgenes
        |> List.filter (fun f ->
            f.featType = "CDS" &&
            not (genes |> List.exists (fun orf -> orf.sysName=f.gene)))
        |> List.iter (fun f -> 
            let gffFields2 =
                [|string f.chr; 
                  "."; 
                  "CDS"; 
                  string (f.l+1) ; 
                  string (f.r+1); "." ; 
                  (if f.fwd then "+" else "-") ; 
                  "0" ; 
                  (sprintf "ID=%s_CDS1;Parent=%s;Note=%s"
                    (urlEscape f.sysName) (urlEscape f.gene) (urlEscape f.description))
                |]

            gff3stream.Write (sprintf "%s\n" (String.Join("\t", gffFields2)))
            printfn "Retrieving CDS with no parent ORF: %s" (String.Join("\t", gffFields2)))
        
        // 4- dump all other types of features (except introns)
        notgenes
        |> List.filter (fun f ->
            f.featType <> "ORF" && f.featType <> "CDS" && f.featType <> "intron") 
        |> List.iter (fun f ->
            let gffFields =
                [|string f.chr; 
                  "."; 
                  f.featType; 
                  string (f.l+1) ; 
                  string (f.r+1); "." ; 
                  (if f.fwd then "+" else "-") ; 
                  "." ; 
                  (sprintf "ID=%s_%s;Name=%s;Note=%s"
                    (urlEscape f.sysName) (f.featType.Replace(" ","")) (urlEscape f.gene) (urlEscape f.description))
                |]

            gff3stream.Write (sprintf "%s\n" (String.Join("\t", gffFields) )))

    /// Output features in gff3, and handle properly gene structure
    /// GFF3 format: http://www.sequenceontology.org/gff3.shtml
    let features2gff3 (features:Feature []) (gff3File:string) =
        use outF = new StreamWriter(gff3File)
        features2gff3stream features outF
        outF.Close()

      /////////////////////////////
     // -- Validate features -- //
    /////////////////////////////

    /// possible issues found in gene annotation
    type GeneIssue = 
        PREMATURE_STOP | LENGTH_NOT_MULT_OF_3 | STRANGE_START_CODON | MISSING_STOP_CODON
        override x.ToString() =
            match x with
            | PREMATURE_STOP -> "PREMATURE_STOP"
            | LENGTH_NOT_MULT_OF_3 -> "LENGTH_NOT_MULT_OF_3"
            | STRANGE_START_CODON -> "STRANGE_START_CODON"
            | MISSING_STOP_CODON -> "MISSING_STOP_CODON"
    
    /// Check different things on a gene sequence: start codon, proper stop, premature stops, length multiple of 3.
    /// geneticCode is usually the codon2aa function from Amyris.biolib. gene input is usually provided by builGeneStructure
    let checkGene
            (geneticCode:char []->char)
            (genome:Dictionary<string,char[]>)
            (gene:CdsInGene [])
        : (GeneIssue list) = 

        let transcript = 
            gene 
            |> Array.map (fun cds ->
                if cds.feature.fwd then genome.[string cds.feature.chr].[cds.feature.l..cds.feature.r]
                else revComp genome.[string cds.feature.chr].[cds.feature.l..cds.feature.r]) 
            |> Array.concat
        
        
        let rec checkGene (transcript:char[]) (issues:GeneIssue list) = 
            match transcript with 
            | [|_;_;_|] ->
                if geneticCode transcript <> '*' then MISSING_STOP_CODON::issues
                else issues
            | [|_;_|] | [|_|] -> LENGTH_NOT_MULT_OF_3::issues
            | _ ->
                if geneticCode transcript.[0..2] = '*' then
                    checkGene transcript.[3..] (PREMATURE_STOP::issues)
                else checkGene transcript.[3..] issues

        if transcript.Length < 3 then 
            [LENGTH_NOT_MULT_OF_3]
        else
            let initIssues =
                if geneticCode transcript.[0..2] <> 'M' then [STRANGE_START_CODON]
                else []
            checkGene transcript initIssues