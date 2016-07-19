namespace Amyris.Bio
/// Legacy binary format (Amyris) for storing genomes, coverage and SNPs. 
module Goo =
    open System.IO

    let gooVerbose = false

    // Genome open organization
    type ChromLoc = { id : int ; off : int ; len : int}
    type SNP = { chr : byte ; off : int; (*frBase: byte ; toBase:char ;*) _frTo : byte ; fwd:uint16 ; rev : uint16}
        with
        member x.frBase() =
            match (x._frTo &&& 3uy) with
            | 0uy -> 'A'
            | 1uy -> 'C'
            | 2uy -> 'G'
            | _ -> 'T'

        member x.toBase() = 
            match ((x._frTo >>>2)&&& 3uy) with
            | 0uy -> 'A'
            | 1uy -> 'C'
            | 2uy -> 'G'
            | _ -> 'T'
                        
    type GOO = { variants : SNP [] ; depth : uint16 [] ; offSets : ChromLoc[]}

    // let mapping2Base = [| 'A' ; 'C' ; 'G' ; 'T' |]

    let decodeVariant8 (i:int) (off:int) =
        let chr = i &&& 255
        let b2 = (i >>> 8) &&& 255
        //let wt = mapping2base.[(b2 &&& 3)]
        //let alt = mapping2base.[(b2 >>> 2)&&&3]
        let frTo = byte(b2 &&& 15) // bottom 4 bits
        let fwd = (i>>> 16) &&& 255
        let rev = (i>>> 24) &&& 255
    
        { chr = byte(chr); off = off ; _frTo = frTo; fwd = uint16 fwd ; rev = uint16 rev}

    let decodeVariant16 (i:uint32) (i2:uint16) (off:int) =
        let chr = i &&& (255ul)
    
        let eFwd = (i2 &&& 255us) <<< 8
        let eRev = i2 &&& (255us <<< 8)
        //let b2 = (i >>> 8) &&& (255ul)
        //let wt = mapping2base.[(b2 &&& (3ul))|> int32]
        //let alt = mapping2base.[(b2 >>> 2)&&&3ul |> int32]
        let frTo = (i>>>8) &&& (15ul) |> byte
        let fwd = (i>>> 16) &&& 255ul + uint32(eFwd)
        let rev = (i>>> 24) &&& 255ul + uint32(eRev)
        { chr = byte(chr); off = off ; _frTo = frTo (* frBase = wt ; toBase = alt *) ; fwd = uint16(fwd) ; rev = uint16(rev)}
    
    type Goo(file:string) = class   
        let _loadFromBR (br:BinaryReader) (prog:string -> unit) =
            let header = br.ReadChars(4)
            prog "Reading header"
            if header <> [| 'G' ; 'O' ; 'O'; 'F' |] then
                (sprintf "file %s not in goo format header=%s" file (header.ToString()) |> failwith )
        
            let version = br.ReadChars(4) 
            let major,minor =
                match version with 
                | [| '0' ;'1'; '0'; '0' |] -> 1,0
                | [| '0' ;'1'; '1'; '0' |] -> 1,10
                | _ -> 0,0 
                            
            if major = 0 then failwith "incompatible version"     
        
            let numBases = br.ReadInt32()
            let numVariants = br.ReadInt32()
            let numChroms = br.ReadInt32()
            sprintf "Loading %d chroms %d bases %d variants" numChroms numBases numVariants |> prog
        
            // Optional read in resolution for depth
            let resolution = if major > 1 || minor >= 10 then
                                  br.ReadByte() |> int else 8
                        
            if gooVerbose then printfn "Reading chrom info"
            let offsets =
                seq { for _ in {0..numChroms-1}
                        -> { id = br.ReadInt32() ; off = br.ReadInt32() ; len = br.ReadInt32() }}
                |> Array.ofSeq

            if gooVerbose then 
                for c in offsets do
                    printf "id=%d off=%d len=%d\n" c.id c.off c.len
                
            prog "reading depth info"
            
            let depth =
                match resolution with
                | 8 -> br.ReadBytes(numBases) |> Array.map (uint16)
                | 16 -> Array.init numBases (fun _ -> br.ReadUInt16())
                | _ -> failwith "bad res"
        
            prog "reading variant data"
        
            let variants =
                match resolution with
                | 8 ->
                    Array.init numVariants (fun _ -> decodeVariant8 (br.ReadInt32()) (br.ReadInt32())) 
                    |> Array.filter (fun s -> s.fwd + s.rev >= 4us)
                | 16 ->
                    Array.init numVariants (fun _ -> decodeVariant16 (br.ReadUInt32()) (br.ReadUInt16()) (br.ReadInt32()))
                    |> Array.filter (fun s -> s.fwd + s.rev >= 4us)
                | _ -> failwith (sprintf "Illegal resolution %d" resolution)
        
            prog "sorting variants"                
            let sortedVariants =
                variants
                |> Array.sortWith (fun v1 v2 ->
                    match compare v1.chr v2.chr with 
                    | -1 -> -1 | 1 -> 1 | 0 -> compare v1.off v2.off
                    | _ -> failwith "bad sort out" )

            br.Close()
            prog "processed coverage"
            sortedVariants,depth,offsets
    
        let _load (file:string) =
            let load = new FileStream(file,FileMode.Open, FileAccess.Read, FileShare.Read)  
            let br = new BinaryReader(load)
            let x =_loadFromBR br (fun _ -> ())
            load.Close()
            x
        
        let mutable variants: SNP array = [||]
        let mutable depth : uint16 array = [||]
        let mutable offsets : ChromLoc array =[||]
        let mutable avgDepth : Map<int,float> = Map.empty
    
        /// Calc average depths post load
        let calcDepths() =
            avgDepth <-
                seq {
                    for cl in offsets do
                        //let avgDepth = depth.[cl.off .. cl.off+cl.len-1] |> Array.fold (fun total v -> total+float(v)) 0.0 |> fun t -> t/float(cl.len)
                        let avgDepth =
                            depth.[cl.off .. cl.off+cl.len-1]
                            |> Array.fold
                                (fun total v -> total+int64(v))
                                0L
                            |> fun t -> double(t)/double(cl.len)
                            
                        yield cl.id,avgDepth
                } |> Map.ofSeq

        do
            if file <> "" then
                let a,b,c = _load(file) 
                variants <- a
                depth <- b
                offsets <- c
                calcDepths()

        member this.GetAvgDepth(chrom:int) = avgDepth.[chrom]
    
        member this.LoadFromBinaryReader(br:BinaryReader) (prog : string -> unit) =
            let a,b,c = _loadFromBR br prog
            variants <- a
            depth <- b
            offsets <- c
            br.Close()
            calcDepths()
           
        member this.LoadFromBytes(b:byte array) (prog: string -> unit) =
            let ms = new MemoryStream(b)
            let br = new BinaryReader(ms)
            this.LoadFromBinaryReader br prog
        
            br.Close()
                        
        member this.getDepthRange(chr,l,r) =
            let x = offsets |> Array.find (fun o -> o.id = chr)
            if r >= x.len then failwithf "Index of %d exceeds chrom %d length of %d" r chr x.len
            let off = x.off
            depth.[off+l..off+r]
        
        /// Find all variation records between a range    
        member this.getVariants(chrI,l,r) =
            let chr = byte(chrI)
            let rec check i j res =
                //printf "Checking i=%d j=%d %d/%d  %d/%d\n" i j (variants.[i].chr) (variants.[i].off) (variants.[j].chr) (variants.[j].off) 
                if (variants.[i].chr > chr) || (variants.[i].chr = chr && variants.[i].off > r) then res else
                if (variants.[j].chr < chr) || (variants.[j].chr = chr && variants.[j].off < l) then res else
                    if j-i <= 10 then
                        seq {for v in variants.[i..j] do
                                if v.chr = chr && v.off >=l && v.off <=r then yield v}
                        |> List.ofSeq
                    else
                        let k = (i+j)/2
                        (check i k res)@(check (k+1) j res)
            check 0 (variants.Length-1) []    
        member this.getChromLen id =
            match (offsets |> Array.tryFind (fun cl -> cl.id = id)) with
            | None -> failwith "goo: request for len for nonexistent id=%d"
            | Some(cl) -> cl.len

    end