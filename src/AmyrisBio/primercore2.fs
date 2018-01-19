namespace Amyris.Bio
open utils
open System
/// Tools for oligonucleotide design
module primercore =
    [<Measure>] type M // Molar
    [<Measure>] type uM // micromolar
    [<Measure>] type mM // millimolar
    [<Measure>] type nM // nanomolar
    [<Measure>] type K // kelvin
    [<Measure>] type C // Celsius
    [<Measure>] type Cal // Calories
    let inline K2C (k:float<K>) = (k*1.0<C/K>) - 273.5<C>
    let inline C2K (c:float<C>) = (c*1.0<K/C>) + 273.5<K>
    let inline mM2M (c:float<mM>) = c*0.001<M/mM>
    let inline uM2M (c:float<uM>) = c/1000000.0<uM/M>
    let inline nM2M (c:float<nM>) = c/1000000000.0<nM/M>

    // Table 2: http://pubs.acs.org/doi/pdf/10.1021/bi702363u
    let private a = 3.92E-5<1/K>
    let private b = -9.11E-6<1/K>
    let private c = 6.26E-5<1/K>
    let private d = 1.42E-5<1/K>
    let private e = -4.82E-4<1/K>
    let private f = 5.25E-4<1/K>
    let private g = 8.31E-5<1/K>

    /// fGC - proportion of GC bases in full primer array
    let gc (primer : char array) =
        let n =
            primer
            |> Array.fold
                (fun gc c ->
                    match c with
                    | 'G' | 'g' | 'C' | 'c' -> gc+1
                    | _ -> gc)
                0
            |> float
        n / (float primer.Length)
    
    /// fGC - proportion of GC bases in first N bases of primer
    let gcN (primer : char array) N = 
        let rec count i gc =
            if i = N then (float gc)/(float N)
            else
                match primer.[i] with
                | 'G' | 'g' | 'C' | 'c' -> count (i+1) (gc+1)
                | _ -> count (i+1) gc

        count 0 0

    /// Take natural log of a molar constant producing a unitless value
    let inline logMolar (x:float<M>) = log ( x / 1.0<M>)

    //let targetLength = 20
    // Previous distance based penalty for not ending on a G or C
    //let ATPenalty = 5

    ///<summary>Parameters used to design primers. These include penalties used to 
    ///rank primers using the various criteria and settings for the calculation of Tm. 
    ///</summary>
    ///<typeparam name="lengthPenalty"> Penalty parameter for difference from targeted primer length. The higher the parameter the lower the penalty. </typeparam> 
    ///<typeparam name="tmPenalty"> Penalty parameter for difference from targeted Tm. The higher the parameter the lower the penalty. </typeparam> 
    ///<typeparam name="positionPenalty"> Penalty parameter for difference from targeted position. The higher the parameter the lower the penalty. </typeparam> 
    ///<typeparam name="polyLengthThreshold"> Upper admissible length for polynucleotide. </typeparam> 
    ///<typeparam name="polyPenalty"> Penalty for going beyong the polynucleotide length limit. </typeparam> 
    ///<typeparam name="threePrimeUnstablePenalty"> Penalty for having an unstable 3' end. </typeparam> 
    ///<typeparam name="ATPenalty"> Penalty for the primer ending on a A or T. </typeparam> 
    ///<typeparam name="maxLength"> Maximum primer length. </typeparam> 
    ///<typeparam name="minLength"> Minimum primer length. </typeparam> 
    ///<typeparam name="targetLength"> Target primer length. </typeparam> 
    ///<typeparam name="monovalentConc"> Monovalent cation concentration (Na+, K+) in mol/L. </typeparam> 
    ///<typeparam name="divalentConc"> Divalent cation concentration (Mg2+, Ca2+) in mol/L. </typeparam> 
    ///<typeparam name="primerConc"> Primer concentration in mol/L. </typeparam> 
    ///<typeparam name="templateConc"> Template concentration in mol/L. </typeparam> 
    ///<typeparam name="dNTPConc"> dNTP concentration in mol/L. </typeparam> 
    type PrimerParams =
       {lengthPenalty : float ; 
        tmPenalty : float ; 
        tmMaxDifference : float<C> ;
        positionPenalty : float ; 
        polyLengthThreshold : int ; 
        polyPenalty : float ; 
        threePrimeUnstablePenalty : float ; 
        ATPenalty :float<C> ;
        maxLength : int ; 
        minLength : int ; 
        targetLength : int ;
        monovalentConc : float<M> ; 
        divalentConc : float<M> ; 
        primerConc : float<M> ; 
        templateConc : float<M> ; 
        dNTPConc : float<M>}
    
    ///<summary>Parameters used to design primers. These include penalties used to 
    ///rank primers using the various criteria and settings for the calculation of Tm. 
    ///</summary>
    ///<typeparam name="lengthPenalty"> Penalty parameter for difference from targeted primer length. The higher the parameter the lower the penalty. </typeparam> 
    ///<typeparam name="tmPenalty"> Penalty parameter for difference from targeted Tm. The higher the parameter the lower the penalty. </typeparam> 
    ///<typeparam name="positionPenalty"> Penalty parameter for difference from targeted position. The higher the parameter the lower the penalty. </typeparam> 
    ///<typeparam name="polyLengthThreshold"> Upper admissible length for polynucleotide. </typeparam> 
    ///<typeparam name="polyPenalty"> Penalty for going beyong the polynucleotide length limit. </typeparam> 
    ///<typeparam name="threePrimeUnstablePenalty"> Penalty for having an unstable 3' end. </typeparam> 
    ///<typeparam name="ATPenalty"> Penalty for the primer ending on a A or T. </typeparam> 
    ///<typeparam name="maxLength"> Maximum primer length. </typeparam> 
    ///<typeparam name="minLength"> Minimum primer length. </typeparam> 
    ///<typeparam name="targetLength"> Target primer length. </typeparam> 
    ///<typeparam name="monovalentConc"> Monovalent cation concentration (Na+, K+) in mol/L. </typeparam> 
    ///<typeparam name="divalentConc"> Divalent cation concentration (Mg2+, Ca2+) in mol/L. </typeparam> 
    ///<typeparam name="primerConc"> Primer concentration in mol/L. </typeparam> 
    ///<typeparam name="templateConc"> Template concentration in mol/L. </typeparam> 
    ///<typeparam name="dNTPConc"> dNTP concentration in mol/L. </typeparam> 
    let defaultParams =
       {tmPenalty = 1.0; 
        tmMaxDifference = 5.0<C> ;
        positionPenalty = 5.0 ; 
        lengthPenalty = 3.0 ; 
        polyLengthThreshold = 4; 
        polyPenalty = 10. ; 
        threePrimeUnstablePenalty = 5.0 ; 
        ATPenalty=3.0<C> ; 
        targetLength=20 ;
        maxLength = 60 ;
        minLength = 20 ; 
        monovalentConc = mM2M 50.0<mM>;
        primerConc = uM2M 0.25<uM> ; 
        divalentConc = mM2M 1.5<mM> ; 
        templateConc = uM2M 0.01<uM> ; 
        dNTPConc = uM2M 0.0<uM> ;}

    /// Final oligo design result with the oligo, its melting temp and offset into the template
    type Oligo = { tag : string ; oligo : char array ; temp : float<C> ; offset : int }

    type GeneInfo = { systematic : string ; name : string ; sequence : char [] }

    type ANCHOR = LEFT | RIGHT | CENTERLEFT | CENTERRIGHT
    type STRAND = TOP | BOTTOM

    /// Inputs to the primer designer. 
    /// temp: template sequence
    /// align: how to align the primer with respect to the template
    /// strand: on what strand to search the primer
    /// offset: 
    /// sequencePenalties: a float [] option of the same size as temp, giving nucleotide specific penalties. If None, an array of 0. will be used.
    ///    These penalties are only used in the CENTERLEFT and CENTERRIGHT ANCHOR types.
    type OligoTask =
       {tag : string ;
        temp : char [] ; 
        align :ANCHOR ; 
        strand : STRAND ; 
        offset : int ; 
        targetTemp : float<C> ; 
        sequencePenalties: float [] option}
    
    type OligoDesign = { tag : string ; oligo : string ; temp : float<C> ; note : string ; pos : int}

    type PrimerHit = { chr : int ; sysName : string ; tag : string ; pos : int ; fwd : bool}
    type PCRHit = { p1 : PrimerHit ; p2 : PrimerHit }

    (*
    Primer: 0.25 uM
    Mon+: 50 mM (although it is actually 60 mM, because there is also 10 mM Tris)
    Mg: 1.5 mM
    *)

    let base2GC b =
        match b with
        | 'G' | 'C' | 'g' | 'c' -> 1
        | 'A' | 'a' | 'T' | 't' -> 0
        | _ -> failwithf "bad base %c in base2GC" b

    type HS = { h : float<Cal/M> ; s: float<Cal/K/M>}
    /// Santa lucia PNAS Feb 17 1998 vol 95 no 3. 1460-1465
    let HSpair (a:char) (b:char) : HS =
        match a,b with
        | 'A','A' | 'T','T'     -> {h= -7.9<Cal/M>; s= -22.2<Cal/M/K> }
        | 'A','G' | 'C','T'     -> {h= -7.8<Cal/M>; s= -21.0<Cal/M/K>}
        | 'A','T'               -> {h= -7.2<Cal/M>; s= -20.4<Cal/M/K> }
        | 'A','C'| 'G','T'      -> {h= -8.4<Cal/M>; s= -22.4<Cal/M/K>}
        | 'G','A' | 'T','C'     -> {h= -8.2<Cal/M>; s= -22.2<Cal/M/K> }
        | 'G','G' | 'C','C'     -> {h= -8.0<Cal/M>; s= -19.9<Cal/M/K> }
        | 'G','C'               -> {h= -9.8<Cal/M>; s= -24.4<Cal/M/K>}
        | 'T','A'               -> {h= -7.2<Cal/M>; s= -21.3<Cal/M/K> }
        | 'T','G' | 'C','A'     -> {h= -8.5<Cal/M>; s= -22.7<Cal/M/K>}
        | 'C','G'               -> {h= -10.6<Cal/M>; s= -27.2<Cal/M/K>}
        | 'G','N'               -> {h= -8.600<Cal/M>; s= -22.225<Cal/M/K> }
        | 'A','N'               -> {h= -7.825<Cal/M>; s= -21.500<Cal/M/K> }
        | 'T','N'               -> {h= -7.950<Cal/M>; s= -22.100<Cal/M/K> }
        | 'C','N'               -> {h= -8.725<Cal/M>; s= -22.700<Cal/M/K> }
        | 'N','G'               -> {h= -8.725<Cal/M>; s= -22.700<Cal/M/K> }
        | 'N','A'               -> {h= -7.950<Cal/M>; s= -22.100<Cal/M/K> }
        | 'N','T'               -> {h= -7.825<Cal/M>; s= -21.500<Cal/M/K> }
        | 'N','C'               -> {h= -8.600<Cal/M>; s= -22.225<Cal/M/K> }
        | 'N','N'               -> {h= -8.275000<Cal/M>; s= -22.131250<Cal/M/K> }
        | _ -> failwith "bad nucpair"

    let deltaG (a:char) (b:char) =
        match a,b with
        | 'A','A' | 'T','T'     -> -1.0
        | 'A','T'               -> -0.88
        | 'T','A'               -> -0.58
        | 'C','A' | 'T','G'     -> -1.45
        | 'G','T' | 'A','C'     -> -1.44
        | 'C','T' | 'A','G'     -> -1.28
        | 'G','A' | 'T','C'     -> -1.30
        | 'C','G'               -> -2.17
        | 'G','C'               -> -2.24
        | 'G','G' | 'C','C'     -> -1.84
        | 'G','N'               -> -1.705
        | 'N','G'               -> -1.685
        | 'A','N'               -> -1.15
        | 'N','A'               -> -1.0825
        | 'T','N'               -> -1.0825
        | 'N','T'               -> -1.15
        | 'C','N'               -> -1.685
        | 'N','C'               -> -1.705
        | 'N','N'               -> -1.405625
        | _ -> failwithf "bad nucpair /%c/%c/" a b    

    let deltaG2 (a:char) (b:char) =
        match a,b with
        | 'A','A' | 'T','T'     -> -1.9
        | 'A','T'               -> -1.5
        | 'T','A'               -> -0.9
        | 'C','A' | 'T','G'     -> -1.9
        | 'G','T' | 'A','C'     -> -1.3
        | 'C','T' | 'A','G'     -> -1.6
        | 'G','A' | 'T','C'     -> -1.6
        | 'C','G'               -> -3.6
        | 'G','C'               -> -3.1
        | 'G','G' | 'C','C'     -> -3.1
        | 'G','N'               -> -2.275
        | 'N','G'               -> -2.55
        | 'A','N'               -> -1.575
        | 'N','A'               -> -1.575
        | 'T','N'               -> -1.575
        | 'N','T'               -> -1.575
        | 'C','N'               -> -2.55
        | 'N','C'               -> -2.275
        | 'N','N'               -> -1.99375
        | _ -> failwithf "bad nucpair /%c/%c/" a b   
    let initG2 (c:char) =
        match c with
        | 'G' | 'C' -> -1.96 // 5.0
        | 'A' | 'T' -> -1.96 - 0.05 // 6.0
        | 'N' -> (-1.96 - 0.05)/2.0 // 6.0
        | _  -> failwithf "bad nuc %c in initG2" c
             
    /// 3' end of primer is special.  We would like it to be
    /// unpromiscuous, and a little harder to bind than the rest
    /// of the primer, so we don't get non specific priming            
    let threePrimeStable verbose (o : char array) =
        if o.Length < 6 then
            failwithf "oligo of length %d too short for threePrimeStability check" o.Length
        let L = o.Length
        let g3 =
            seq {for i in {L-5..L-2} do
                    if verbose then printf "i=%d dg=%f\n" i (deltaG2 (o.[i]) (o.[i+1]))
                    yield deltaG2 (o.[i]) (o.[i+1])}
            |> Seq.sum

        let g3Final = g3 + initG2 (o.[L-1])
    
        if verbose then printf "deltaG 3' = %f\n" g3Final
        
        g3Final > -9.0

    /// Test if an oligo is self complementary
    let selfComp (oligo:char array) length =
        let rec checkOne l r =
            l > r ||
            match oligo.[l],oligo.[r] with
            | 'G','C' | 'C','G' | 'A','T' | 'T','A' -> checkOne (l+1) (r-1)
            | _ -> false

        if length%2=1 then false else checkOne 0 (length-1)

    /// NN implementation with 1M Na+ concentration
    /// http://www.pnas.org/content/95/4/1460.full
    /// http://en.wikipedia.org/wiki/DNA_melting
    let wikiTemp1000 (primerA:float<M>) (primerB:float<M>) (o : char array) N =
        // Nucleic Acids Research, 1996, Vol. 24, No. 22 4501–4505
        //let endBase2 (a:char) = match a with 'A' | 'T' -> {h=0.6 ; s= -9.0}  | 'C' | 'G' -> {h=0.6;s= -9.0} | _ -> failwith (sprintf "bad start base '%c'" a)
        let endBase1 (a:char) =
            match a with
            | 'A' | 'T' -> {h=2.3<Cal/M> ; s=4.1<Cal/M/K>}
            | 'C' | 'G' -> {h=0.1<Cal/M>;s= -2.8<Cal/M/K>}
            | 'N' -> {h=1.2<Cal/M>;s= 0.65<Cal/M/K>} // Take average for unknown case
            | _ -> failwith (sprintf "bad start base '%c'" a)

        let rec score i hs =
            if i = N then 
                let hs' = endBase1 (o.[i-1])
                {h=hs.h + hs'.h ; s = hs.s + hs'.s}
            else
                let hs' = HSpair (o.[i-1]) (o.[i]) 
                score (i+1) {h= hs.h+hs'.h ; s= hs.s + hs'.s } // +(h+h') (s+s')
        let hSF = endBase1 (o.[0]) // 0.0,0.0 
        let hs = 
            match N with
                | 1 -> hSF // hS,hF
                | 2 ->
                    let hSF2 = endBase1 (o.[1])  // Start/finish
                    {h=hSF.h+hSF2.h ;  s= hSF.s+hSF2.s}
                | _ -> score 1 hSF
        let R = 1.98722<Cal/M/K> // ideal gas constant

        // Santalucia
        // Proc.Natl.Acad. Sci.USA Vol. 95, pp. 1460–1465, February 1998 Biochemistry
        // eq (3)
        // Self complementary molecules get Ct = primerConcentration
        // Non self complementary case gets primerConc/4 if strands are in equal concentration
        // or (Ca-Cb/2) if different concentrations were Ca and Cb strands respectively and Ca > Cb
        let comparable (a:float<M>) (b:float<M>) = abs(a-b)/(a+b) < 0.1
        let symmetrical = selfComp o N 
        let Ct =
            if symmetrical then primerA
            elif comparable primerA primerB then primerA / 4.0
            else ( (max primerA primerB)-(min primerA primerB)/2.0)
               
        let hs' = { h = hs.h * 1000.0 ; s = hs.s } // This magic is due to deltaH units being in 100Cal/mol and deltaS being in 0.1Cal/Mol 
        hs'.h/(
                ( hs'.s+(if symmetrical then -1.4<Cal/K/M> else 0.0<Cal/K/M>)+(R*logMolar(Ct) ) )
                     (* /1000.0  *)  (* This magic is due to deltaH units being in 100Cal/mol and deltaS being in 0.1Cal/Mol *)
             )  (*|> K2C *)

    /// Calc temp given a GC and oligo length
    let gcN2Temp gc n = 
        if n < 14 then (4*gc) + 2 * (n-gc)*2
        else (64.9+41.0*(float(gc)-16.4)/float(n) ) |> int

    let monoSaltCorrect (mon:float<M>) (Ca:float<M>) (Cb:float<M>) (o : char array) N =
        let tm1000 = wikiTemp1000 Ca Cb o N
        let fgc = gcN o N
        //http://www.owczarzy.net/Tm_for_duplexes-IDT_Tech.pdf
        let correction = (4.29*fgc-3.95)*1E-5<1/K>*logMolar(mon)+9.40E-6<1/K>*logMolar(mon)*logMolar(mon)
        let tmMonRecip = (1.0/tm1000) + correction
        1.0/tmMonRecip |> K2C

    // eq 16 table 2
    let eq16Table2 (div:float<M>) (primerA:float<M>) (primerB:float<M>) (dNTP:float<M>) (primer: char array) (N:int) =
        let mgConc = max 0.0<M> (div-dNTP) // NTPs bing divalent ions, so remove that from further consideration
        let l = logMolar mgConc
        let tm1000 = wikiTemp1000 primerA primerB primer N
    
        // fGC is the fraction of residues that are G or C
        let fGC = gcN primer N

        let mgRecip = 1.0<1> /tm1000 + a + b*l + fGC * (c+d*l) + (1.0/(2.0 * float (N-1)) * (e+f*l+g*l*l) )
        1.0/mgRecip |> K2C

    // eq 16 table 2
    let eq16Table2Mod (div:float<M>) (mon:float<M>) (Ca:float<M>) (Cb:float<M>) (dNTP:float<M>) (primer: char array) (N:int) =
        let mgConc = max 0.0<M> (div-dNTP) // NTPs bind divalent ions, so remove that from further consideration
        let l = logMolar mgConc
        let lm = logMolar mon
        let tm1000 = wikiTemp1000 Ca Cb primer N
    
        let modA  = 3.92E-5<1/K> * (0.843 - 0.352 * sqrt(mon / 1.0<M> (* UNITHACK *) ) * lm)
        let modD  = 1.42E-5<1/K> *(1.279-4.03E-3*lm-8.03E-3*lm*lm)
        let modG  = 8.31E-5<1/K> *(0.486-0.258*lm+5.25E-3*lm*lm*lm)

        // fGC is the fraction of residues that are G or C
        let fGC = gcN primer N

        let mgRecip = 1.0 / tm1000 + modA + b*l + fGC * (c+modD*l) + (1.0/(2.0 * float (N-1)) * (e+f*l+modG*l*l) )
        1.0/mgRecip |> K2C

    /// Get oligo Tm for the first N bases.
    let _temp (p:PrimerParams) (oligo: char array) (N:int) =
        // http://pubs.acs.org/doi/pdf/10.1021/bi702363u
        // figure 9
        if p.monovalentConc < 1e-8<M> then  
            // eq16 table 2
            1,999.0,eq16Table2 p.divalentConc p.primerConc p.templateConc p.dNTPConc oligo N // case 1
        else 
            let r = sqrt(p.divalentConc*1.0<M> (* UNITHACK *) ) / p.monovalentConc
            if r < 0.22 then 
                // use monovalent salt correct eq 4
                2,r,monoSaltCorrect p.monovalentConc p.primerConc p.templateConc oligo N // case 2
            else
                if r < 6.0 then 
                    // Modified eq 16 table 2
                    // with a,d,g varying according to eq 18,19,20
                    3,r,eq16Table2Mod p.divalentConc p.monovalentConc p.primerConc p.templateConc p.dNTPConc oligo N // case 3
                else
                    // eq 16 table 2
                    4,r,eq16Table2 p.divalentConc p.primerConc p.templateConc p.dNTPConc oligo N // case 4
    
    /// Get oligo Tm for the first N bases.
    let temp (p:PrimerParams) (oligo: char array) (N:int) =
        let _,_,x = _temp p oligo N
        x

    /// Cut out the region from fr -> to and extend
    // Given oligo array from to gcInit offset , return oligo, temp, offset
    let cutToGC (debug:bool) (existingTemp : float<C>) (p:PrimerParams) (s:char[]) f t _ (*startingGC*) offset =
        //let startingN = t-f + 1
        let rec findGCFwd t' (*gc' *) =
            if t' < 0 || t' >= s.Length then
                failwithf "ERROR: primercord findGC array bounds exception t'=%d\n" t'
            match s.[t'] with
            |'G' | 'g' | 'C' | 'c' when threePrimeStable false (s.[..t'])-> 
                (t', temp p (s.[f..t']) (t'-f+1)) // end on G/C
                    
            | _ when t' < (s.Length-1) && t'-f+1 < p.maxLength -> findGCFwd (t'+1)
            | _ -> (t, (temp p s.[f..t] (t-f+1))) // fall back on original oligo if we run out of S without finding G

        let rec findGCRev t' =
            if t' < 0 || t' >= s.Length then
                failwithf "ERROR: primercord findGC array bounds exception t'=%d\n" t'
            match s.[t'] with
            |'G' | 'g' | 'C' | 'c' when threePrimeStable false (s.[..t']) -> 
                (t', temp p (s.[f..t']) (t'-f+1)) // end on G/C
                    
            | _ when t' > 1  && t'-f+1 > p.minLength -> findGCRev (t'-1)  
            | _ -> (t, (temp p s.[f..t] (t-f+1))) // fall back on original oligo if we run out of S without finding G
                                           
        if p.ATPenalty < 0.0001<C> then
            // No GC optimization
            if debug then printfn "cutToGC: no atPenalty, done"
            { tag = ""; oligo = s.[f..t] ; temp = temp p (s.[f..t]) (t-f+1); offset = f+offset } 
        else
            // Is there an alternative starting point that would end on a G or C that's not too far away?                
            let altTFwd,altTempFwd = findGCFwd t
            let altTRev,altTempRev = findGCRev t

            assert(altTempFwd > 10.0<C>)
            assert(altTempRev > 10.0<C>)
            if debug then
                printfn "cutToGC: altTempFwd=%A altTFwd=%d altTempRev=%A altRFwd=%d"
                    altTempFwd altTFwd altTempRev altTRev

            // Choose the direction to a GC that was least perturbative
            let altT,altTemp =
                if abs (altTempFwd-existingTemp) < abs(altTempRev-existingTemp) then
                    altTFwd,altTempFwd
                else altTRev,altTempRev

            if abs (altTemp-existingTemp) <= p.ATPenalty then
                if debug then printfn "cutToGC: using altGC option temp=%A f=%d t=%d" altTemp f altT 

                { tag = ""; oligo = s.[f..altT] ; temp = altTemp ; offset = f+offset}
            else
                if debug then printfn "cutToGC: ignore altGC option temp=%A f=%d t=%d" temp f t 
                { tag = ""; oligo = s.[f..t] ; temp = temp p (s.[f..t]) (t-f+1); offset = f+offset } 

    let upper (s:char array) =
        [| for x in s ->
            match x with 
            | 'a' -> 'A'
            | 't' -> 'T'
            | 'g' -> 'G'
            | 'c' -> 'C'
            | _ as x -> x |]
         
    let lowerStr s =
        let l c =
            match c with 
            | 'A' -> 'a'
            | 'T' -> 't'
            | 'G' -> 'g'
            | 'C' -> 'c'
            | _ as x -> x
        Array.map (l) s

    ///
    /// oMin oMax template nextBase gcCount offset tTemp
    /// nextBase - next position to add to extending oligo 0 to start
    /// gcCount - starts at 0
    /// offset - user offset to add to final position report
    /// tTemp - target temperature    
    let rec designLeft (debug:bool) (p:PrimerParams) (sInit:char []) nextBase gcCount offset tTemp  =
        let sInitUpper = upper sInit

        let rec designLeftInternal (s:char []) nextBase gcCount lastTemp =
            let thisTemp = 
                if nextBase < p.minLength-1 then 0.0<C> 
                        else temp p s (nextBase+1)   // Calculate temp for oligo including position nextBase
            if debug then printfn " designLeft: thisTemp=%f nextBase=%d lastTemp=%f" (thisTemp/1.0<C>) nextBase (lastTemp/1.0<C>)
            if thisTemp < tTemp && nextBase < s.Length-1 then 
                // Keep going till we exceed target or hit limit
                designLeftInternal s (nextBase+1) (gcCount + base2GC s.[nextBase]) thisTemp
            else 
                if nextBase >= p.maxLength then None else 
                    if nextBase<p.minLength || abs(thisTemp-tTemp) < abs(lastTemp-tTemp) then
                        if debug then printfn "designLeft, stopping at nextBase=%d thisTemp=%f" nextBase (thisTemp/1.0<C>)
                        Some(cutToGC debug thisTemp p s 0 nextBase gcCount offset)
                    else
                        if debug then printfn "designLeft, stopping at nextBase-1=%d thisTemp=%f" (nextBase-1) (thisTemp/1.0<C>)
                        Some(cutToGC debug lastTemp p s 0 (nextBase-1) gcCount offset)

        designLeftInternal sInitUpper nextBase gcCount (0.0<C>)
    
    /// Check for the presence of mono- or dinucleotide repeats longer than N
    let hasPolyrun n (s:char[]) =
        let rec worstMonoNuclRun i lastChar thisRun worst =
            if i = s.Length then worst |> max thisRun else
            if s.[i] = lastChar then worstMonoNuclRun (i+1) lastChar (thisRun+1) worst else
            worstMonoNuclRun (i+1) (s.[i]) 1 (worst |> max thisRun)
        let rec worstDiNuclRun i lastChars thisRun worst =
            if i >= s.Length then worst |> max thisRun 
            elif i < 1 then worstDiNuclRun 1 lastChars thisRun worst 
            elif s.[i-1..i] = lastChars then
                worstDiNuclRun (i+2) lastChars (thisRun+1) worst 
            else
                worstDiNuclRun (i+2) (s.[i-1..i]) 1 (worst |> max thisRun)
        (worstMonoNuclRun 0 '?' 0 0) >=n 
            || worstDiNuclRun 1 [||] 0 0 >= n 
            || worstDiNuclRun 2 [||] 0 0 >= n

    
    [<Obsolete>]
    let hasPolyGRun n (s:char[]) =
        let rec worstRun i thisRun worst =
                    if i = s.Length then worst |> max thisRun else
                    if s.[i] = 'G' then worstRun (i+1) (thisRun+1) worst else
                    worstRun (i+1) 0 (worst |> max thisRun)
        (worstRun 0 0 0) >=n

    type GCN = {gc : int ; n : int}

    type Candidate =
       {l : int ; 
        r : int ; 
        posD : float ; 
        tempD : float ; 
        lenD : float ; 
        t: float<C> ; 
        sequenceP : float ; 
        threeP : bool ; 
        hasPoly : bool ; 
        pen : float}
    
    /// Design all compatible primers from a specified sequence.
    let designCenterMany (pen:PrimerParams) (s:char[]) (seqPen:float []) tTemp =
        let n = s.Length
        let temps = Array2D.init n n (fun _ _ -> {gc = 0 ; n = 0} )
    
        let mid = s.Length / 2
        seq {
            for l in {0..n-1} do
            for r in {0..n-1} do
                if r > l && r-l+1>=pen.minLength && r-l+1 <= pen.maxLength then
                    temps.[l,r] <-
                        if l = r then
                            {gc = base2GC s.[l] ; n = 1}
                        else 
                            {gc = temps.[l,r-1].gc + base2GC s.[r] ; n = temps.[l,r-1].n+1 } 
                                
                    let t = temp pen (s.[l..]) (r-l+1)  //  gcN2Temp (temps.[l,r].gc) (temps.[l,r].n)                
                    let tempCloseness = float (abs(t - tTemp)) / pen.tmPenalty
                                
                    // We have an ideal length in mind, how close is this solution?
                    let lenCloseness = float(abs(r-l+1-pen.targetLength)) / pen.lengthPenalty 
                                
                    // Record how close to ideal position and penalize ending on an A/T
                    let posCloseness = 
                        (float(abs(l-mid))/pen.positionPenalty) + 
                        (match s.[r] with | 'G' | 'g' | 'C' | 'c' -> 0.0 | _ -> float(pen.ATPenalty) )
                    
                    // Banned CCCC in primer
                    // http://www.pnas.org/content/93/22/12116.full.pdf
                    let polyC = (arr2seq (s.[l..r])).Contains("CCCC")

                    let hasPoly = hasPolyrun (pen.polyLengthThreshold+1) (s.[l..r])
                    let threeP = if r-l+1 > 5 then threePrimeStable false (s.[l..r]) else false
                    // get average of all nucleotide specific penalties
                    let seqP = Array.average seqPen.[l..r]
                    let p = posCloseness + 
                            tempCloseness + 
                            lenCloseness + 
                            (if threeP then 0.0 else pen.threePrimeUnstablePenalty) + 
                            (if hasPoly then pen.polyPenalty else 0.0) + 
                            seqP
                    if (abs(t - tTemp) <= pen.tmMaxDifference) && (not polyC) then
                        yield
                           {l = l ; 
                            r = r ; 
                            posD = posCloseness ; 
                            tempD = tempCloseness ; 
                            lenD = lenCloseness ; 
                            t=t ; 
                            threeP = threeP ; 
                            sequenceP = seqP ;
                            hasPoly = hasPoly ;
                            pen = p }
        }

    /// Design an oligo centered around the middle base of this sequence    
    let designCenter debug (pen:PrimerParams) (s:char[]) (seqPen:float []) offFn tTemp =
        
        let viable = designCenterMany pen s seqPen tTemp
                     
        let dummy =
           {l= -1 ; 
            r= -1; 
            posD=999999.0 ; 
            tempD=999999.0 ; 
            lenD=99999.0 ; 
            t=999999.0<C> ; 
            threeP = false ; 
            pen = 999999.0 ; 
            sequenceP= 999999.0; 
            hasPoly = true} 
        
        let best =
            viable
            |> Seq.fold
                (fun a b ->
                    if debug then
                        printf "considering l= %d r= %d len= %d temp= %A posD= %f tempD= %f len= %f 3s=%s hasPoly=%s score/pen=%f %s\n" 
                            b.l b.r (b.r-b.l+1) b.t b.posD b.tempD b.lenD (if b.threeP then "Y" else "N")
                            (if b.hasPoly then "Y" else "N")
                            b.pen
                            (arr2seq s.[b.l..b.r])
                                    
                    if a.pen <= b.pen then a else b
                
                    (*
                    if b.posD < a.posD && b.tempD < a.tempD && b.lenD < a.lenD then b // New obviously better
                    else if b.posD > a.posD && b.tempD > a.tempD && b.lenD > a.lenD then a // old obviously better
                    else 
                        // If the threeprimer end is too stable, we need to levy a penalty because we'd
                        // like it to not prime too easily
                        // threeP true means threeP end is ok
                        let threePDev (c:Candidate) = if c.threeP then 0.0 else pen.three
                        let totalDeviance = a.posD + a.tempD + a.lenD + (threePDev a)
                        let totalDeviance' = b.posD + b.tempD + b.lenD + (threePDev b)
                        if totalDeviance < totalDeviance' then a // old better in aggregate
                        else b // New better in aggregate      
                    *)
                )
                dummy       
    
        //Some(s.[best.l..best.r],best.t,offFn best.l) 
        if debug then
            printf "best l= %d r= %d len= %d temp= %A posD= %f tempD= %f lenD= %f hasPoly=%s score/Pen=%f oligo=%s\n"
                best.l best.r (best.r-best.l+1) best.t best.posD 
                best.tempD best.lenD (if best.hasPoly then "Y" else "N") best.pen
                (if best.l <> -1 then arr2seq s.[best.l..best.r] else "")
        if best.t >= 999990.0<C> then None else 
            let o = { tag = ""; oligo = s.[best.l..best.r] ; temp = best.t ; offset = offFn best.l}
            assert(o.oligo.Length <= pen.maxLength) 
            Some( o ) 


    /// Main entry point for oligo design
    /// Takes debug flag, Penalty Params and an Oligo task
    /// Returns a tag,(primer,temp,offset) option  tuple 
    let oligoDesign debug pen (o : OligoTask) =
        if pen.maxLength < 0 then
            failwithf "oligoDesign: pen.oMax must be positive, and it was %d\n" pen.maxLength
        if pen.minLength > pen.maxLength then
            failwithf "oligoDesign: pen.oMin (%d) should be less than or equal to pen.oMax (%d)\n"
                pen.minLength pen.maxLength

        // First deal with strand request.  If they want a design on
        // the bottom strand, need to flip sequence, and also reverse alignment/anchor request
        let align, s, seqPen, offFn = 
            // initialize sequence penalties, if not provided
            let seqPen = 
                match o.sequencePenalties with 
                | Some p -> p
                | None -> Array.init o.temp.Length (fun _ -> 0.)
            
            match o.strand with
                | TOP -> o.align, (o.temp) , seqPen, (fun x -> o.offset+x)
                | BOTTOM ->
                    let align =
                        match o.align with
                        | LEFT -> RIGHT
                        | RIGHT -> LEFT
                        | CENTERRIGHT -> CENTERLEFT 
                        | CENTERLEFT -> failwith "should not get centerleft"
                    align, biolib.revComp (o.temp), Array.rev seqPen, (fun x -> o.offset-x)

        let s' = upper s
        match align with
        | LEFT ->
            match designLeft debug pen s' 0 0 o.offset (float o.targetTemp * 1.0<C>) with
            | None -> // If there was no solution, go with longest possible oligo
                if debug then
                    printf "No oligo in required length, so make max allowable max=%d templateLen=%d\n" 
                        pen.maxLength s'.Length
                let oligo = s'.[..(min s'.Length pen.maxLength)-1] // longest allowed
                Some( { tag=o.tag; oligo = oligo ; temp = temp pen oligo pen.maxLength; offset = o.offset } ) // calc temperature and offset is 0 plus their supplied offset
            | x -> x // 
                    
        | RIGHT -> failwith "should not get right here"
        | CENTERLEFT -> designCenter debug pen s' seqPen offFn ( (float o.targetTemp)*1.0<C>)
    
                        (*with 
                        | Some(a,b,c,d) -> Some( { tag=o.tag; oligo = a ; temp = b ; offset = c} ) 
                        | None -> None  *)

        | CENTERRIGHT -> failwith "should not get centerright here"

    (*
       // When you can't meet the default primer design criteria (e.g. min 20, 60+/- 5C tm), what should you do

       Biologist 1: I vote for preserving the minimum annealing length (20 bp or maybe 18 bp would be okay), 
                    and I’m willing to pay the cost of higher melting temps. 
       Biologist 2: I’d try to get the right melting temperature. Oligos smaller than 20 work. Plus, they are 
                    cheaper! Just make sure to add a note that operators will have to tap those plates 
                    three times while spinning counterclockwise or they will fail.

       Biologist 3: I think its ok to overshoot the target melting temp. Very anecdotal, but I have found 
                    in my personal experience that a primer with higher mp than needed had a better chance of working.
    
    *)
    type GcContent = GCHIGH|GCLOW|GCMED|GCUNKNOWN

    /// Same function as oligoDesign but with heuristics to relax the constraints if necessary to
    /// get a design
    let oligoDesignWithCompromise debug pen (o : OligoTask) =
        let rec tryWith (gcContent:GcContent) (penCurr:PrimerParams)  =
            // First just try to design an oligo
            match oligoDesign debug penCurr o with
            | None -> 
                // that didn't go so well, have we worked out gc regime yet?
                // establish if we have high or low GC since that dictates different
                // strategies
                // establish what situation we are dealing with if we haven't already
                let gcContent' = match gcContent with
                                    | GCUNKNOWN -> 
                                        let gc = gc o.temp
                                        if gc < 0.4 then GCLOW 
                                        elif gc > 0.6 then GCHIGH
                                        else GCMED
                                    | _ -> gcContent
                let pen' = 
                    match gcContent' with
                        | GCHIGH  ->
                            // in high GC we fail to reach min primer length.  Might be able to relax that.
                            // if that gets too dire, allow temp deviation to float up
                            if penCurr.minLength > 18 then 
                                if debug then 
                                    printfn "primerDesignWithCompromise: relax minLength to %d" 
                                        (penCurr.minLength-1)
                                Some {penCurr with minLength = penCurr.minLength-1}
                            elif penCurr.tmMaxDifference < 10.0<C> then
                                if debug then 
                                    printfn "primerDesignWithCompromise: increase tmMaxDiff to %fC" 
                                        (penCurr.tmMaxDifference/1.0<C>)
                                Some {penCurr with tmMaxDifference = penCurr.tmMaxDifference+1.0<C>}
                            else
                                None // run out of ideas,  fail
                        | GCMED  ->
                            // unclear if we need to allow longer or shorter primers
                            // try relaxing TmMaxDiff first
                            if penCurr.tmMaxDifference < 10.0<C> then
                                if debug then 
                                    printfn "primerDesignWithCompromise: increase tmMaxDiff to %fC" 
                                        (penCurr.tmMaxDifference/1.0<C>)
                                Some {penCurr with tmMaxDifference = penCurr.tmMaxDifference+1.0<C>}
                            elif penCurr.minLength > 18 then 
                                if debug then 
                                    printfn "primerDesignWithCompromise: relax minLength to %d" 
                                        (penCurr.minLength-1)
                                Some {penCurr with minLength = penCurr.minLength-1}
                            elif penCurr.maxLength < 30 then 
                                if debug then 
                                    printfn "primerDesignWithCompromise: relax maxLength to %d" 
                                        (penCurr.maxLength+1)
                                Some {penCurr with maxLength = penCurr.maxLength+1}
                            else
                                None // run out of ideas,  fail
                        | GCLOW  ->
                            // Low GC usually results in long primers that fail to reach their target tm
                            // try relaxing TmMaxDiff first
                            if penCurr.tmMaxDifference < 10.0<C> then
                                if debug then 
                                    printfn "primerDesignWithCompromise: increase tmMaxDiff to %fC" 
                                        (penCurr.tmMaxDifference/1.0<C>)
                                Some {penCurr with tmMaxDifference = penCurr.tmMaxDifference+1.0<C>}
                            elif penCurr.maxLength < 30 then 
                                if debug then 
                                    printfn "primerDesignWithCompromise: relax maxLength to %d" 
                                        (penCurr.maxLength+1)
                                Some {penCurr with maxLength = penCurr.maxLength+1}
                            else
                                None // run out of ideas,  fail
                        | GCUNKNOWN -> 
                            failwithf "Impossible"

                // if there are suggested params to try, go for it, otherwise fail whole
                // operation
                match pen' with
                    | None -> None
                    | Some x ->
                        // try again with new parameters
                        tryWith gcContent' x
            | x -> x // victory, return a result

        tryWith GCUNKNOWN pen

    (*
    /// Stitch two sequences left and right together with overlapping oligos creating an optional middle sequence
    /// ensure at least minOverlap bases shared and a maximum length with a target melting temperature
    type StitchSol = { fwd : char [] ; rev : char [] ; fwdOff : int ; revOff : int ; fwdTemp : int ; revTemp : int; overTemp : int ; pen : float} 
    let stitchThree debug (left:char array) (mid:char array) (right:char array) minOverlap maxLength temp =
        //
        //           o-------------------->
        //  <------------------o
        // =======================ABCDE==========================
        //
        //                           o----------->
        //                <----------------o
        // =======================ABCDE==========================
        let merged = Array.concat [ left ; mid; right] |> upper
        let merged' = Array.concat [left ; mid; right |> lowerStr]
    
        let rOff = left.Length + mid.Length // base at which right starts
    
        let targetOverTemp = 60
        let oMax = 60
    
        // TODO.... need to finish this up for stitching across a gap with sequence in it
        // 
        let attempt = seq { for i in {rOff-minOverlap .. rOff} do
                                // Starting at position i, see how much overlap is needed to get a decent stitching overlap
                                //                                          rhsEnd
                                //                       i>>>>>>>>>>>>.......
                                //                       <<<<<<<<<<<<<
                                // ============================== ======================
                                //                                ^rhsStart
                                match designLeft oMax (merged.[i..]) 0 0 0 targetOverTemp with
                                    | None -> ()
                                    | Some(fwd,fTemp,fOff) ->
                                            // Next check how far we'd need to extend to get a decent RHS primer
                                            let rhsStart = max i rOff
                                            match designLeft oMax (merged.[rhsStart..]) 0 0 0 temp with
                                                | None -> ()
                                                | Some(rhs,rhsTemp,rhsOff) ->
                                                    let rhsEnd1 = rhsStart + rhs.Length // How far we need to go out to get rhs to prime
                                                    let rhsEnd2 = i+fwd.Length-1 // How far we need to go out to get enough overlap
                                                    let rhsEnd = max rhsEnd1 rhsEnd2
                                                
                                                    let len1Pen = (rhsEnd - i + 1-60) |> max 0  |> float |> fun x -> x * 5.0
                                                
                                                    // Next check how far we'd need to extend left to get a decent LHS primer
                                                    let lhsStart = min (rOff-1) (i+fwd.Length-1)
                                                    match designLeft oMax (merged.[..lhsStart] |> revComp) 0 0 0 temp with
                                                        | None -> ()
                                                        | Some(lhs,lhsTemp,lhsOff) ->
                                                            let lhsEnd = lhsStart-lhs.Length+1
                                                            let overlapMidpoint =  i+fwd.Length/2
                                                            let len2Pen = (i+fwd.Length-1-lhsEnd+1-60) |> max 0 |> float |> fun x -> x*5.0
                                                        
                                                            let pen = float(
                                                                        abs (rhsTemp-temp) + // deviation from ideal rhs temp
                                                                        abs (lhsTemp-temp) + // deviation from lhs temp
                                                                        abs (fTemp-targetOverTemp) // deviation from overlap temp
                                                                      ) + float(abs(overlapMidpoint-rOff))/5.0 // 1 point penalty for 5 bp slide
                                                                      + float(rhsEnd-i+1+i+fwd.Length-1-lhsEnd+1)/10.0 + len2Pen + len1Pen
                                                                  
                                                                  
                                                            let fwdPrimer = merged.[i..rhsEnd]
                                                            let revPrimer = merged.[lhsEnd..rhsEnd2] |> revComp
                                                        
                                                            if debug then
                                                                printf "len1Pen=%f\nlen2Pen=%f\nrhsEnd1=%d\nrhsEnd2=%d\np1Len=%d\np2Len=%d\n" 
                                                                        len1Pen len2Pen rhsEnd1 rhsEnd2 fwdPrimer.Length revPrimer.Length
                                                        
                                                                printf "%f %d %d %d\n%s%s\n%s%s\n%s\n\n" pen lhsTemp fTemp rhsTemp 
                                                                            (pad i) (fwdPrimer |> arr2seq) (pad lhsEnd) (revPrimer |> arr2seq) 
                                                                            (merged' |> arr2seq)
                                                            yield {fwdOff = i ; fwd = fwdPrimer
                                                                    ; revOff = lhsEnd ; rev = revPrimer
                                                                    ; fwdTemp = rhsTemp ; revTemp = lhsTemp ; overTemp = fTemp
                                                                    ; pen = pen }
                               } |> List.ofSeq
    
        //if debug then
        //    for ss in attempt do
        //        printf "%f %d %d %d\n%s%s\n%s%s\n%s\n\n" ss.pen ss.revTemp ss.overTemp ss.fwdTemp 
        //                        (pad ss.fwdOff) (ss.fwd |> arr2seq) (pad ss.revOff) (ss.rev |> arr2seq) (merged' |> arr2seq)
    
        match attempt |> List.fold (fun bestOpt x ->
                                    match bestOpt with
                                        | None -> Some(x)
                                        | Some(best) ->
                                                    if x.pen < best.pen then
                                                        Some(x)
                                                    else bestOpt
                             ) None
                                                                        
                with
                | None -> None
                | Some(best) -> 
                        if debug then
                            printf "Best: pen=%f %d/%d/%d l1=%d l2=%d\n" best.pen best.revTemp best.overTemp best.fwdTemp best.fwd.Length best.rev.Length
                        Some(best.fwd,best.rev)
      *)