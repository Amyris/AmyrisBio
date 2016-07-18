
namespace Amyris
/// efficient binomial, hypergeometric and prime number calculations
module math_bhg =

    open math_common 

    /// Sieve of Atkin
    /// http://blog.hoomla.se/post/The-sieve-of-Atkin-in-F.aspx
    /// Primes up to but not including limit
    let sieveOfAtkin limit =
        // initialize the sieve
        let sieve = Array.create (limit + 1) false
        // put in candidate primes: 
        // integers which have an odd number of
        // representations by certain quadratic forms
        let inline invCand n pred =
            if n < limit && pred then sieve.[n] <- not sieve.[n] 
        let sqrtLimit = int (sqrt (float limit))
        for x = 1 to sqrtLimit do
            for y = 1 to sqrtLimit do
                let xPow2 = x * x
                let yPow2 = y * y
                let n = 4 * xPow2 + yPow2 in invCand n (let m = n % 12 in m = 1 || m = 5)
                let n = 3 * xPow2 + yPow2 in invCand n (n % 12 = 7)
                let n = 3 * xPow2 - yPow2 in invCand n (x > y && n % 12 = 11)
        // eliminate composites by sieving
        let rec eliminate n =
            if n <= sqrtLimit 
            then if sieve.[n]
                 then let nPow2 = n * n
                      for k in nPow2 .. nPow2 .. limit do
                          Array.set sieve k false
                 eliminate (n + 2)
        eliminate 5
        // Generate list from the sieve (backwards)
        let rec generateList acc n =
            if n >= 5 then generateList (if sieve.[n] then n :: acc else acc) (n - 1)
            else acc
        2 :: 3 :: (generateList [] limit)                

    let mutable primeCache : int list option = None
    let primesTill limit =
        match primeCache with
            | Some(l) when l.Length >= limit -> l
            | _ ->
                let l = sieveOfAtkin limit
                primeCache <- Some(l)
                l
    /// Factorize N! into a list of primes ^ n
    let nFactPrimes N =
        if N = 1 then [ (1,1) ] else
            let pList = primesTill (N+1)

            /// Calculate exponent for a single prime number
            let procPrime p = 
                let rec total pj res = if pj > N then res else total (pj * p) (res + N/pj)
                total p 0

            /// Calculate exponents for all primes
            let rec procPrimes primes res =
                match primes with
                    | [] -> List.rev res
                    | hd::_ when hd > N -> res
                    | hd::tl ->
                        let t = procPrime hd
                        procPrimes tl (if t<> 0 then ((hd,t)::res) else res)
            procPrimes pList []
    
    let mergePrimesBase (pList1:(int*int) list) (pList2:(int*int) list) op = 
        let rec merge l1 l2 res =
            match l1,l2 with
                | [] , _ -> (List.rev res)@(l2 |> List.map (fun (p,c) -> p,op 0 c) )
                | _, []  -> (List.rev res)@l1 // (l1 |> List.map (fun (p,c) -> p,op c 0) )
                | (p1,c1)::tl1,(p2,c2)::tl2 when p1 = p2 -> merge tl1 tl2 ((p1,op c1 c2)::res)
                | (p1,c1)::tl1,(p2,_)::_ when p1 < p2 -> merge tl1 l2 ((p1,c1)::res)
                | _,(p2,c2)::tl2  -> merge l1 tl2 ((p2,op 0 c2)::res)
        merge pList1 pList2 []

    let mergePrimes p1 p2 = mergePrimesBase p1 p2 (+)
    let dividePrimes p1 p2 = mergePrimesBase p1 p2 (-)
    let expandPrimesDumb p = p |> List.map (fun (prime,power) -> double(prime)**double(power)) |> List.reduce (*)
    let expandPrimes p = 
        let top,bot = p |> List.partition (fun (_,exp) -> exp >= 0)
        let expand (series:(int*int)list) =
            seq { for p,exp in series do
                    match exp with
                        | 0 -> ()
                        | x when x > 0 ->
                            for _ in {0..exp-1} do
                                yield p
                        | _ (* when y < 0 *) ->
                            for _ in {0..(-exp)-1} do
                                yield p
                } |> List.ofSeq
        let topList = expand top
        let botList = expand bot
        let rec combine a b (t:float) =
            match a,b with
                | [],[] -> t
                | [],(hd::tl) -> combine [] tl (t/float(hd))
                |(hd::tl),[] -> combine tl [] (t*float(hd))
                | _,(h2::t2) when t > 0.0 -> combine a t2 (t/float(h2))
                | (h1::t1),_ when t <= 0.0 -> combine t1 b (t*float(h1))
                | _ -> failwithf "ERROR: unexpected outcome in prime list condensation\n"

        combine topList botList 1.0



    let nCkDumb n k =
        let a = seq {(n-k+1) .. n} |> Seq.map(float) |> Seq.fold (*) (double(1.0))
        let b = seq {1..k} |> Seq.map(float) |> Seq.fold (*) (double(1.0))
        a/b
    let nCkPrimes n k = 
            if k > n then [ (0,1) ] // Can't choose more than we have
            else
                let a = nFactPrimes n
                let b = nFactPrimes k
                let c = nFactPrimes (n-k)
                let d = mergePrimes b c
                dividePrimes a d

    /// hyperGeo N n k m  = prob of drawing k white balls from n draws of an urn with N balls, m of which are white
    let hyperGeo N n k m = 
        if k > n then 0.0 // Can't see more hits than picks
        elif k > m then 0.0 // Can't see more hits than white balls
        else
            let a = nCkPrimes m k
            let b = nCkPrimes (N-m) (n-k)
            let c = nCkPrimes N n
            dividePrimes (mergePrimes a b) c |> expandPrimes
        //dividePrimes (mergePrimes (nCkPrimes m k)  (nCkPrimes (N-m) (n-k)) ) (nCkPrimes N n) |> expandPrimes

    /// hyperGeo N n k m  = prob of drawing k white balls from n draws of an urn with N balls, m of which are white
    let hyperGeoDumb N n k m = ((nCkDumb m k) * (nCkDumb (N-m) (n-k)) / (nCkDumb N n)) 

    /// N,n,k,m : Choosing n balls from an N ball urn, m of which are white, what is the probability that k or less are white.
    let hyperGeoCum N n k m =
        seq { for i in {0..k} do
                yield hyperGeo N n i m
                } |> Seq.sum


    let binomial n k p =  (nCkDumb n k) * (p**float(k)) * ((1.0-p)**float(n-k))

    let binomialCum n k p =
        seq { for i in {0..k} do
                yield binomial n i p
              } |> Seq.sum

    (*
    let z = bignum.FromInt
    let nCkB n k =
        let a = seq {(n-k+1) .. n} |> Seq.map(bignum.FromInt) |> Seq.fold (*) (bignum.FromInt 1)
        let b = seq {1..k} |> Seq.map(bignum.FromInt) |> Seq.fold (*) (bignum.FromInt 1)
        a/b          
    
    
    /// hyperGeo N n k m  = prob of drawing k white balls from n draws of an urn with N balls, m of which are white
    let hyperGeo N n k m = ((nCkB m k) * (nCkB (N-m) (n-k)) / (nCkB N n)) |> bignum.ToDouble


    let hyperGeoCum N n k m =
        seq { for i in {0..k} do
                    yield hyperGeo N n i m
            } |> Seq.sum

    *)