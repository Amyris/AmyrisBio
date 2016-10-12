namespace Amyris.Bio.IO
module pileup =

    // http://en.wikipedia.org/wiki/Pileup_format

    // Lots of features in the dots field we need to watch out for
    // Indels.   e.g. +1t  +2AC -1t -2?
    // 62331	2198	A	793	,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t.+1T.+1T.+1T.+1T.+1T.+1T.+1T..+1T.+1T.+1T.+1T..+1T.+1T..+1T.+1T..+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T..+1T..+1T.+1T.+1T.+1T.+1T.+1T.+1T..+1T.+1T.+1T.....$.+1T.+1T.+1T.+1T..+1T.+1T.+1T..+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T..+1T.+1T.+1T.+1T.+1T..$.+1T.+1T..+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T...+1T.+1T.+1T..+1T.+1T.+1T.+1T..+1T.$.+1T.+1T,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,$,+1t,+1t,+1t,+1t,$,+1t,,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t.+1T.+1T.+1T.+1T.+1T.+1T,+1t,+1t.+1T.+1T.+1T.+1T.+1T.+1T.+2TT.+1T.+1T.+1T.+1T,+1t,+2tt,+1t.+1T,+1t,+1t,+1t,+1t,+1t.+1T,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t.+1T,+1t.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t.+1T.+1T,+1t.+1T.+1T.+1T,+1t,+1t.+1T.+1T,+1t.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T..+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t.+1T.+1T.+1T.+1T,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t.+1T.+1T.+1T.+1T.+1T.+1T.+1T,+1t,+1t,+1t.+1T,+1t,+1t,+1t,+1t,+1t,+1t.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t.+1T,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t.+2TT.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T,+1t,+1t,+1t,+1t,+1t,+1t.+1T.+1T.+1T.+1T,+1t,+1t.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t.+1T.+1T.+1T.+1T.+1T,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t.+1T.+1T,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T.+1T,+1t,+1t,+1t,+1t.+1T.+1T,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t.+1T.+1T,+1t,,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t,+1t	HHGHHHHHHHHHHHHHHHHH;H0HHF;FHHHHGG3HFEFHHEFHG<GFHBFEHGH1HBFF1HHHHEHBHHEEHH3CGHEHHGHG3F1FFF1EHGHBHGH3FGHGFFBHHHH1GGHHHHH1FHGHHHFHHH3HGHAHHFHHFHEFHHHH;HGHHHHHHGHHGHHFHHHHGHGCHGHH1HFHGHAHGFHHHHC3HHHFF3GGHH3HGGFHHHHHHGHGHFGHHGHHHGGGG;HHFHHHFHGHHHHHGFHHHHHHH0BHGHHAHG;GHHHH?HHHHFHHGHHHHHA0H;FHHHHHH1HHHFHHHFHCHHHHGFFHHHHHFHGFFHHGHHGHHGHFHHHHF41HHHB4HH1HHHH4HHHFCHHHHHHH4HGHHHC4HHHGFFGGHH<H;FHHHCGGH0HHH0GH0GHHHGFH;HHHHHHGGBCHHHHGGFGGHHHHGH;HGHH?HHHGHHCFHF1HHCHHHGHAHHHH1HHHBHHGGGHHGGHB@HHHFHFAHH1FHGGFHHHHHHFHGHHH0HHGH;HHHFHHHCEGH?H1HHHHHHHHAFAFHHGHHHBHFHHHHHHHHHH;HHHHHGHGGGGCGGGGFGFGGHHHHGHFGGBHHFFGFGGGGGGFGGGFGGFHHH;HHHF9HGHHHHHGGGGGFFH9GHHGFHHGHHFHEHHGHH9GHHHGHHGGGCGFGGGCCAFHHHHHH;HFHHGH<EHGGHHHFHHGHHHHH9H<BFFFFCFFFB0HHHFFGGFHHHHHHHGGHGH9FHH9FHHHCHHG9FHHHGHF9HH90HHGFFHHHHHFHHHHHHHGFHHHHHHGH
    // 
    // Simple substitutions   tttttTTTTT

    /// Things we can identify in a pileup string
    type DotFeat =
            | FwdSNP of char
            | RevSNP of char
            | Indel of int*string

    /// Active pattern to parse a series of digits
    let rec (|Digits|_|) = function 
        | c::tl when c >= '0' && c <= '9' ->
        
            let rec aux (total:int) = function
                | c::tl when c >= '0' && c <= '9' -> aux (total*10+((int c) - int '0')) tl
                | _ as x -> Some(total,x)
            aux 0 (c::tl)
        | _ -> None

    let rec take (l:char list) (res:char list) (n:int) =
        if n = 0 then (List.rev res),l
        else
            take (l.Tail) (l.Head::res) (n-1)
        
    
    /// Match a series of characters in a pileup file, parsing one character each time called
    let rec private parseOne  start stop fwd rev (out : DotFeat list) (dots:char list)=
        match dots with
        | ','::tl  -> parseOne start stop fwd (rev+1) out tl
        |'.'::tl -> parseOne start stop (fwd+1) rev out tl
        | '$'::tl -> parseOne start (stop+1) fwd rev  out tl
        | x::tl  when x = 'A' || x = 'T' || x = 'C' || x = 'G' || x = 'N' -> parseOne start stop (fwd+1) rev (FwdSNP(x)::out) tl
        | x::tl  when x = 'a' || x = 't' || x = 'c' || x = 'g' || x = 'n' -> parseOne start stop fwd (rev+1) (RevSNP(x)::out) tl
        | '*'::tl -> parseOne start stop fwd rev out tl // asterisk appears as placeholder for a multibase deletion documented on a previous line, we can ignore
        | '^'::_::tl -> parseOne (start+1) stop fwd rev out  tl // not storing mapping qual char which is in second position or start read information
        | '+'::Digits(d,rem) -> // +2AC  type insertion 
                    let bases,rem2 = take rem [] d
                    match bases.[0] with
                        | 'A' | 'T' | 'C' | 'G' -> parseOne start stop (fwd+1) rev (Indel(d,new string(bases |> Array.ofList))::out) rem2 
                        | 'a' | 't' | 'c' | 'g' -> parseOne start stop fwd (rev+1) (Indel(d,new string(bases |> Array.ofList))::out) rem2 
                        |_ -> parseOne start stop fwd rev (Indel(d,new string(bases |> Array.ofList))::out) rem2 

        | '-'::Digits(d,rem) -> // -1C   deletion
                    let bases,rem2 = take rem [] d
                    match bases.[0] with
                        | 'A' | 'T' | 'C' | 'G' ->  parseOne  start stop (fwd+1) rev (Indel(-d,new string(bases |> Array.ofList) )::out) rem2 
                        | 'a' | 't' | 'c' | 'g' -> parseOne  start stop fwd (rev+1) (Indel(-d,new string(bases |> Array.ofList) )::out) rem2 
                        | _ ->  parseOne  start stop fwd rev (Indel(-d,new string(bases |> Array.ofList) )::out) rem2 
               
        | [] -> out,start,stop,fwd,rev
        | _ -> failwithf "ERROR: parsing dot string %s" (new string(Array.ofList  dots))
 
     /// Turn a pileup character string into features e.g. 
     ///  +1C,+1c..AAaa  ->  Indel(1,'C'),SNPFwd('A')  SNPFwd('A')  SNPRev('A') etc    
    let parse (dots:string) = parseOne 0 0 0 0 [] (dots.ToCharArray() |> List.ofArray )