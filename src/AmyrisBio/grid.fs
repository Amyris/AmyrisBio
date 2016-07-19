namespace Amyris.Bio

/// A tool for importing tab delimited data with a header row allowing name based
/// lookup of columns
module grid =
    open System.IO
    open System

    type TabRow = { row: int ; gridRef : TabGrid ref} with
        static member (?) (tabRow:TabRow, name : string) = (!tabRow.gridRef).Cell(tabRow.row,name)
        member x.Col (name : string) = (!x.gridRef).Cell(x.row,name)

    and 
     TabGrid(lines:string [] []) = class
        
        do if lines.Length = 0 then failwith "ERROR: expected at least one header line"

        let headers = lines.[0] |> Array.mapi (fun i v -> (v,i)) |> Map.ofSeq
        do
            ()

        member x.Headers = lines.[0] // |> Seq.map (fun pk -> pk.Key) |> Array.ofSeq
        member x.Row(j:int) = lines.[j+1]
        member x.Cell(i:int,j:int) = lines.[j].[i]
        member x.Item(j:int) = lines.[j+1]
        member x.Col(name:string) = 
            let i = headers.[name]
            lines.[1..] |> Array.map (fun row -> row.[i])

        member x.Cell(row:int,name:string) = 
            let i = headers.[name]
            lines.[row+1].[i]

        member x.Rows = seq { for i in 0..lines.Length-2 -> {row = i ; gridRef = ref x } }

        /// Parse a tab delimited file to a TabGrid
        new(path:string) = 
            let lines =
                File.ReadAllLines(path)
                |> Array.map (fun row -> row.Split([| '\t' |],StringSplitOptions.None))
            TabGrid(lines)
    end



    //let (?) (tabRow:TabRow) (name : string) = (!tabRow.gridRef).Cell(tabRow.row,name)


(*

    //let x = t.Row(5)
    //let z = t.[10]
    // let col = t.Col("total")
    // let cell = t.Cell(10,5) // row 10, col 5 zero based
    // let cell = t.Cell(10,"total") // row 10, col 5 zero based
    // for row in t.Rows do
    //     printf "total = %s" row?total

    // Example
    let tg = TabGrid(@"c:\data\GI\deltalog.txt")  // tab separated file with text header row at top

    tg.Rows   ->  sequence generator that can be iterated or filtered
    tg.Row -> get one row
    tg.Col   -> get one row

    Row objects can be accessed with the ? operator to look up columns based on the header name


    for row in tg.Rows |> Seq.filter (fun row -> row?situation = "root" && row?unclean <> "n") do
                let comment = row?comment


    One thing I couldn't get working is overloading the ? operator to work on both the whole table (to extract a column) and a row (as it currently does),
    So I've used it to pull out elements of an extracted row.   Everything is handled as a string, so do any type conversion you need yourself,


*)