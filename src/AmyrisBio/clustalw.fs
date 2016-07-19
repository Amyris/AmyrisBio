namespace Amyris.Bio.IO

/// parser for CLUSTAL file format
module clustalw =
    type AlignedSeq = { name : string ; alignment : string}
    type ClustalAlign = { seq : AlignedSeq array ; conservation : string }

    /// Parse a string containing a clustalw formatted output
    let parse (s:string) =
        let lines = s.Replace("\r","").Split([| '\n' |]) |> List.ofArray
        let header = "CLUSTAL"
        if lines.Length = 0 || not (lines.Head.StartsWith(header)) then
            failwithf "ERROR: no header line found"
        else
            let names =
                lines.Tail.Tail.Tail
                    |> Seq.takeWhile (fun line -> not (line.StartsWith(" "))) 
                    |> Seq.map (fun line -> line.Split([| ' '|]).[0])
                    |> List.ofSeq

            let rec splitBlock
                    (lines:string list)
                    (margin:int option)
                    (thisBlock:string list)
                    (blocks:(string list*string) list) =
                match lines with
                | [] when thisBlock = [] -> blocks
                | [] -> failwithf "Unexpected EOF in parse clustalw"
                | hd::tl when hd.Length > 0 && hd.[0]<> ' ' -> 
                    let spaceOff = hd.IndexOf(' ')
                    let rec findNextNonSpace i =
                        if i = hd.Length then None
                        elif hd.[i] = ' ' then findNextNonSpace (i+1)
                        else Some(i)
                    let margin'= match findNextNonSpace spaceOff with | None -> margin | Some(v) -> Some(v)
                    splitBlock tl margin' (hd::thisBlock) blocks
                | a::tl when a.StartsWith(" ") -> 
                    // End of block a is the conservation line, b is the gap between blocks
                    let blocks' =
                       (thisBlock
                        |> List.rev
                        |> List.map (fun line -> line.[margin.Value..]),a.[margin.Value..])
                       ::blocks
                    match tl with
                    | [] -> blocks'
                    | hd::tl2 ->
                        assert(hd.Trim() = "")
                        splitBlock tl2 None [] blocks'
                | _ -> failwithf "Unexpected pattern in parse clustalw"

            let alignLines,consLine = 
                splitBlock (lines.Tail.Tail.Tail) None [] []
                |> List.rev 
                |> List.fold
                    (fun (blocksJoined:string list, consLineJoined:string) 
                         (block:string list,consLine:string) ->
                        List.zip blocksJoined block
                        |> List.map (fun (a,b) ->a+b),consLineJoined+consLine)
                    ([for _ in names -> ""],"")
            {seq =
                List.zip names alignLines
                |> List.map (fun (a,b) -> {name = a ; alignment = b})
                |> Array.ofList;
             conservation=consLine}

