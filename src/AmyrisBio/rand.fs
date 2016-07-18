#light
namespace Amyris
module rand =
    open math
    open Microsoft.FSharp.Math

    let random = new System.Random()
    let sqRatio = ref 0.0
    let q = ref 0.0
    let phase = ref false

    //  below are native implementations of stat distro's using names familiar to any good useR
    // refer back to math.fs for some of the underlying functions

    /// Box-Muller transform for generating random variables from normal (Gaussian) distribution
    let rnorm mean sd =
        if !phase then
            phase:= false
            (!sqRatio * !q * sd)+mean
        else
            let rec tillGoodEnough () =
                let lp = random.NextDouble()*2.0-1.0
                let lq = random.NextDouble()*2.0-1.0
                let v = lp*lp + lq*lq
                if v > 1.0 || v < 0.25 then
                    tillGoodEnough()
                else
                    q := lq // store for next gaussian
                    v,lp
            let v,p = tillGoodEnough()
            let sqratio = sqrt(-2.0*log(random.NextDouble()) / v)
            (sd * sqratio * p)+mean


    //--
    /// test Normal RNG; change it a little for the other distros
    let test() =
        let samples = Array.init 1000 (fun _ -> rnorm 10.0 2.0) 
        let m = mean samples
        let s = sd samples

        printfn "mean=%f sd = %f" m s
        assert(abs(m - 10.0) < 0.1)
        assert(abs(s - 2.0) < 0.1)
        ()
    //--


    /// inversion method for generating random variables from (negative) exponential distribution
    let rexp mean = -mean*log(random.NextDouble());
    

    /// equation for generating random veriables from poisson distribution
    let rpois (mean:float) (k:int) = (mean**float(k)) * exp(-mean) / float(fact k)     


    /// generating random variables from gamma distribution
    /// possible sources for this algorithm(?): R.C.H. Cheng, "The generation of Gamma variables with
    /// non-integral shape parameters", Applied Statistics, (1977), 26, No. 1, p71-74; ALGORITHM GS of 
    /// Statistical Computing - Kennedy & Gentle
    let rgamma (k:float) (sigma:float) =

        (*
                    1. Let m be 1.
                    2. Generate V3m − 2, V3m − 1 and V3m as independent uniformly distributed on (0,1] variables.
                    3. If V_{3m - 2} \le v_0, where v_0 = \frac e {e + \delta}, then go to step 4, else go to step 5.
                    4. Let \xi_m = V_{3m - 1}^{1 / \delta}, \ \eta_m = V_{3m} \xi _m^ {\delta - 1}. Go to step 6.
                    5. Let \xi_m = 1 - \ln {V_{3m - 1}}, \ \eta_m = V_{3m} e^{-\xi_m}.
                    6. If \eta_m > \xi_m^{\delta - 1} e^{-\xi_m}, then increment m and go to step 2.
                    7. Assume ξ = ξm to be the realization of Γ(δ,1)
                    *)

        let rec uniNonZero() = match random.NextDouble() with | x when x>0.0 -> x | _ -> uniNonZero()

        let kInt = floor k |> int // integer part of k
        let delta = k - (floor k) // fractional part of k

        let v0 = E / (E + delta)

        let rec step()=
            let v3m_2 = uniNonZero()
            let v3m_1 = uniNonZero()
            let v3m = uniNonZero()
            let em,nm =
                if v3m_2 <= v0 then
                    let em = v3m_1 ** (1.0/delta)
                    let nm = v3m*(em**(delta-1.0))
                    em,nm
                else
                    let em = 1.0-log (v3m_1)
                    let nm = v3m*(E ** -em)
                    em,nm
            if nm > (em**(delta-1.0)) then step() else em

        let em = step()
        sigma*(em - (seq { for i in {1..kInt} do
                            yield log (uniNonZero()) } |> Seq.sum )
                )


    /// generates random variables from beta distribution; based on algorithm from Knuth Vol. 2 p. 134
    let rbeta alpha beta = (rgamma alpha 1.0) |> ( fun y -> if y=0.0 then 0.0 else y/(y + (rgamma beta 1.0)) )


    /// weibull random variable generator from jain, p. 499
    let rweib alpha beta = alpha*((-log(1.0-random.NextDouble()))**(1.0/beta))


    /// so that Darren's odler code does not break
    let getGamma = rgamma
    let getNorm = rnorm
    let getNegExp = rexp
    let getPois = rpois


    /// Beta Distribution probability gteq x aka cdf(x) := regularized incomplete beta function see ibeta a b x
    /// c++ boost doc http://www.boost.org/doc/libs/1_39_0/libs/math/doc/sf_and_dist/html/math_toolkit/dist/dist_ref/dists/beta_dist.html
    /// general overview of beta distribution identities http://en.wikipedia.org/wiki/Beta_distribution
    let pbeta (x:float) (a:float) (b:float) : float = ibeta a b x


    /// Student's t Distribution probability gteq x aka cdf(x); loosely based on beta distribution
    /// c++ boost doc http://www.boost.org/doc/libs/1_39_0/libs/math/doc/sf_and_dist/html/math_toolkit/dist/dist_ref/dists/students_t_dist.html
    /// general overview of student's t identities at http://en.wikipedia.org/wiki/Student%27s_t-distribution
    let pt (t:float) (dof:float) : float = 
        printf "t=%f dof=%f" t dof
        if (t-0.0) < Epsilon then 0.0
        else if dof < 2.0*(t**2.0) then ibeta (dof/2.0) 0.5 (dof/(dof+(t**2.0))) |> (*) 0.5
        else ibeta 0.5 (dof/2.0) ((t**2.0)/(dof+(t**2.0))) |> (-) 1.0 |> (*) 0.5  // complementary function to keep p<0.9



    /// returns lower-tail probability, used to assess how well one sample fits to a second sample or a pop. mean
    /// takes at least one list of values (sample1) and either a point value or second list of values
    /// optionally specify if this is a one_sample ind't t-test, testing how data conforms to theoretical pop. mean
    /// must specify if we can assume equal variance between the two lists -> if not, becomes Welch's t-test
    /// optionally specify if we want a lower-tail probability calculation; assumes two-tail as default
    /// thanks to (B.E. Cooper, 1968) for the original FORTRAN AS3 Algorithm, "The Integral of Student's t-Distribution"
    /// thanks to (Ajay Shah, 1991) for publishing an f2c translation of the AS3 Algorithm on netlib.org for testing
    let tTest (sample1:float[]) (sample2:float[]) (one_sample:bool) (equal_var:bool) (lower_tail:bool) : float = 
        let count1 = sample1.Length
        let count2 = sample2.Length
        let count1f = float count1
        let count2f = float count2
        let mean1 = mean sample1 // this is simple, general knowledge behind the t-statistic calculation
        let mean2 = mean sample2
        let delta_mean = mean1-mean2 // difference in mean values (numerator of t-statistic)
        let var1 = if equal_var=true && count1=1 && count2<>1 then var sample2 else var sample1
        let var2 = if equal_var=true && count2=1 && count1<>1 then var sample1 else var sample2
        let std_err = // standard error (denominator of t-statistic)
            match one_sample,equal_var,count1>1,count2>1,count1=count2 with   
                                | true,true,_,_,_ -> failwith "only specify equal_var if not one_sample"
                                | true,_,true,true,_ -> failwith "one_sample must specify one true mean"
                                | _,_,false,false,true -> failwith "use a paired difference test to compare two means"
                                | false,false,false,_,_ -> failwith "assume equal variance if one sample is n=1"
                                | false,false,_,false,_ -> failwith "assume equal variance if one sample is n=1"
                                | true,_,true,false,false -> sqrt(var1/count1f) // independent one-sample t-test
                                | true,_,false,true,false -> sqrt(var2/count2f) // same, just if user inputs reverse for some reason
                                | _,true,_,_,true -> sqrt( (var1+var2) / (count1f+count2f) ) // independent two-sample t-test, equal variance, equal sample sizes
                                | _,true,_,_,false -> sqrt((((count1f-1.0)*var1)+((count2f-1.0)*var2))/(count1f+count2f-2.0))*sqrt((1.0/count1f)+(1.0/count2f)) // independent two-sample t-test, equal variance, unequal sample sizes
                                | _,false,_,_,_ -> sqrt((var1/count1f)+(var2/count2f)) // Welch's t-test: actually not a pooled variance, just heuristic when variances are unqual
        let t = delta_mean / std_err // final t-statistic calculation
        let dof = // not an int in the case of W-S equation, and unequal variance is a common real-world situation
            if equal_var=true // if equal variance use typical df definition, else use Welch–Satterthwaite equation
            then count1f+count2f-2.0
            else (((var1/count1f)+(var2/count2f))**2.0)/((((var1/count1f)**2.0)/(count1f-1.0))+(((var2/count2f)**2.0)/(count2f-1.0)))
        if System.Double.IsNaN(t) then failwith "NaN generated; remember test with unequal variance requires sample with n>1" 
        let pi_inv:double = 0.3183098861837906715377675 // one over pi; used in cdf integration heuristic
        if dof < 1.0 then failwith "dof < 1: (%d)\n" dof
        pt t dof


    /// returns list of q-values, min false discovery rates for a list of p-values termed significant under multiple-testing
    /// takes a list of significant p-values and lambda, threshold over which true null hypotheses are uniformly distributed 
    /// uses the lambda estimate and list of p-values to estimate pi0 (proportion of true null hypotheses); 0.5 is a good first try
    /// still need to introduce "smoothing" and "bootstrap" pi0 estimates from (Storey and Tibshirani, 2003)
    /// tested with John Storey's qvalue() algorithm (open source under LGPL) in cran.r-project.org
    /// intentionally retains order of the p-values for efficient array processing
    let qValues (pvals : float[]) (lambda : float) = 
        printf "%A" pvals
        if pvals.Length=0 then [|0.0|] else
        let p_max = Array.max pvals
        let p_min = Array.min pvals
        if p_min < 0.0 || p_max > 1.0 then failwith "p-values are not in correct form"
        if lambda<0.0 || lambda>1.0 then failwith "lambda needs to be in range [0.0,1.0)"
        let pi0 = Array.choose (fun p -> if p>=lambda then Some(p) else None) pvals |> mean |> (fun m -> m/(1.0-lambda)) |> min 1.0
        if pi0<0.0 then failwith "pi0 estimated less than 0.0; please check that your p-values and lambda are correct"
        let m = float pvals.Length
        let q = Array.map (fun x -> x,(Array.filter (fun y -> y<x) pvals).Length) pvals // p comma count p's greater than p
                    |> Array.map (fun (x,c) -> (pi0*m*x/(float c)),c) // map what SHOULD be the minimum FDR (handle exceptions below >>)
        let qPlus = Array.map ( fun (x,c) -> if c=(int m) || c=0 then x,c else (Array.find (fun y -> (snd y)=(c+1)) q) ) q // NOT index plus one; element with exactly count plus of more significant p-values replaces current
        let qMin = Array.map2 (fun x y -> min (fst x) (fst y)) q qPlus
        if (Array.max qMin)>1.0 || (Array.min qMin)<0.0 
        then failwith "error in the calculation: [need to implement printing q-vals here]" // [|for i in q do yield "%f" (fst i)|] |> Seq.ofArray |> String.concat ", " |>
        qMin


    /// Holmes multiple-testing correction - based on Holmes (1979)
    /// order your k p-values, consider smallest p-value a true positive if l.t. fpr/k; now remove from list and repeat with k-1; 
    /// once one false positive is discovered, simply reject all p-values g.t.e. to this value
    let holmesCorrect (pvals_in : float[]) (fpr_in : float) = 
        let acc_in = Array.empty
        let rec hC (pvals : float[]) (fpr : float) (acc : float[]) = 
            let low_now = Array.min pvals
            if low_now < (fpr/(float pvals.Length)) 
            then hC pvals fpr (Array.append acc [|low_now|])
            else acc
        hC pvals_in fpr_in acc_in