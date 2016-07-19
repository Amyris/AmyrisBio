namespace Amyris.Bio
/// beta , gamma distribution implementations
module math_spf = 
    
    open math_common
    open math_bhg

    /// Lanczos approximation deconvolved from the power terms; this allows the approx to 
    /// be ported into gamma, beta, and other gamma-like functions without compounding the 
    /// inaccuracy by, e.g., calling beta = gam(a)gam(b)/gam(a+b); coefficients computed by 
    /// Paul Godfrey according to method detailed below (he has published this all over the internet)
    /// Great overview discussion for c++ boost docs: http://www.boost.org/doc/libs/1_39_0/libs/math/doc/sf_and_dist/html/math_toolkit/backgrounders/lanczos.html
    let lanczos1 z _ (*g*) = 
        // best results when 4<=g<=5; note above: g=607.0/128.0 optimized for complex gamma 
        let c0 = 0.99999999999999709182 // this is the base Lanczos coefficient
        let c = [| 57.156235665862923517; // just the next 14 terms, so total N=15
                -59.597960355475491248;
                14.136097974741747174;
                -0.49191381609762019978;
                0.33994649984811888699e-4;
                0.46523628927048575665e-4;
                -0.98374475304879564677e-4;
                0.15808870322491248884e-3;
                -0.21026444172410488319e-3;
                0.21743961811521264320e-3;
                -0.16431810653676389022e-3;
                0.84418223983852743293e-4;
                -0.26190838401581408670e-4;
                0.36899182659531622704e-5 |]
        let c_index = Array.init c.Length (fun i -> float (i+1)) // needed to foldBack on Lanczos coefficients 1->14
        //let dropone pp k agg = agg+(pp/(z+k-1.0))
        Array.foldBack2 (fun pp k agg -> agg+(pp/(z+k-1.0))) c c_index 0.0 |> (+) c0 // final sum with base coefficient

        // "Lanczos approximation for the complex plane calculated using vpa [var prec] digits(256) 
        // the best set of coeffs was selected from a solution space of g=0 to 32 with 1 to 32 terms
        // these coeffs really give superb performance of 15 sig. digits for 0<=real(z)<=171
        // coeffs should sum to about g*g/2+23/24" (Paul Godfrey)
        //
        // c is calculated from c=D*B*C*F
        // where D is [1  -Sloane's sequence A002457],
        // fact(2n+2)/(2fact(n)fact(n+1)), diagonal
        // and B is rows from the odd  cols of A & Stegun table 24.1 (Binomial),
        // unit Upper triangular
        // and C is cols from the even cols of A & Stegun table 22.3 (Cheby),
        // C(1)=1/2, Lower triangular
        // and F is sqrt(2)/pi*gamma(z+0.5).*exp(z+g+0.5).*(z+g+0.5).^-(z+0.5)
        // gamma(z+0.5) can be computed using factorials,2^z, and sqrt(pi).
        // A & Stegun formulas 6.1.8 and 6.1.12)
        // 
        // for z=0:(g+1) where g=9 for this example. Num. Recipes used g=5
        // gamma functions were made for g=5 to g=13.
        // g=9 gave the best error performance
        // for n=1:171. This accuracy is comparable to Matlab's
        // real only Gamma function.
        //
        // this will generate the Num. Recp. values
        // (if you used vpa for calculating F)
        // notice that (D*B*C) always contains only integers
        // (except for the 1/2 value)
        //
        // Note: 24*sum(c) ~= Integer = 12*g*g+23 if you calculated c correctly
        // for this example g=5 so 24*sum(c) = 322.99... ~= 12*5*5+23 = 323
        //
        // Spouge's approximate formula for the c's can be calculated from applying
        // L'Hopitals rule to the partial fraction expansion of the sum portion of
        // Lanczos's formula. These values are easyer to calculate than D*B*C*F
        // but are not as accurate. (An approx. of an approx.)


    /// Gamma function relies on Lanczos coefficients optimized by Paul Godfrey for Octave (g=607.0/128.0, n=15)
    /// great overview of strategy at http://www.boost.org/doc/libs/1_39_0/libs/math/doc/sf_and_dist/html/math_toolkit/backgrounders/lanczos.html, also http://functions.wolfram.com/GammaBetaErf/Gamma2/10/
    /// this is a reimplementation of Godfrey's gamma function Lanczos approx for complex nums
    /// permission for use in octave detailed here: http://octave.1599824.n4.nabble.com/Better-Gamma-implementation-td1675513.html
    let gamma x = /// this is built for complex numbers too, but cannot figure out how to declare them in F#
        if abs((round x)-x)<Epsilon && x>0.0 then fact (x-1.0)
        else 
            let g=607.0/128.0 // best results when 4<=g<=5 for complex gamma; optimized by Paul Godfrey for this function
            let RT2PI=sqrt(PI+PI) // square root of two pi -- accurate enough?  
            let z = abs(x) - 1.0 // really, we want the abs of the realPart of any comlex numbers, just cannot find F# complex number type
            let zh = z+0.5 // z-0.5
            let zgh = zh+g // z+g=0.5, part of numerator of power term for Lanczos
            let zp = zgh**(zh*0.5) // avoids overflow error for z>141
            // Lanczos calculated to n=15, should be sufficient for our purposes; see lanczos z for details
            let ss = lanczos1 (z+1.0) g // calls lanczos subroutine, need to add 1.0 b/c lanczos already subtracts 1.0 from z
            if abs((round z)-z)<Epsilon && z<=0.0 (* && not imaginary *) then failwithf "gamma(%f) evaluates to infinity!" z 
            else match z=0.0,z=1.0,z>0.0 with   | true,false,false -> 1.0
                                                | false,true,true -> 1.0
                                                | false,false,true -> (RT2PI*ss)*((zp**2.0)*(exp (-zgh))) 
                                                | false,false,false -> (RT2PI*ss)*((zp**2.0)*(exp (-zgh))) |> (fun y -> (-PI)/(z*y*sin(PI*z)))
                                                | _,_,_ -> failwithf "z=%A unexpected" z


    /// Upper incomplete gamma function evaluates integral x->Inf; Lower evaluates 0->x
    /// note: these values are already normalized; to denormalize, just transform with gamma s
    /// note2: these calculations may be suboptimal for integers; better approx for natural s,x is binomial
    /// note3: these calculations may break for large or disparate values s,x; implement boost special cases if necessary
    let rec igamma_up s x = 
        if s<=0.0 then failwith "incomplete gamma function cannot handle s<=0.0 yet"
        else if x<0.0 then failwith "integral approximation unstable for x < zero"
        else if x=0.0 then gamma s // integral 0->Inf is just the regular gamma!!
        else if x>1.1 then // Legendre's continued fraction representation of gamma
            let tiny = 10e-30 // tiny shift for fnot, Cnot, Dnot to prevent instability around zero
            let anot = 1.0 // first term a is one
            let bnot = 0.0 // not leading coefficient; we will add this ourselves after multiplying by power terms
            let fnot = if bnot=0.0 then tiny else bnot // use tiny shift to prevent C and D from hitting zero
            let Cnot = fnot 
            let Dnot = 0.0 
            /// Winitzki simplification by first trasnforming to exp(x)*x**(1-a)*gamma = a1/(b1+(a2/(b2+(a3/(b3+ ...
            /// Lentz approx of continued fraction based on Numerical Recipes ed.2, Ch. 5.2
            let rec gFrac f _(*a why not used?*) _ (*b*) C D m = 
                let a' = if m%2.0=1.0 then (m-s)/x else m/x
                let b' = 1.0
                let Dn = b' + a'*D
                let D' = if Dn=0.0 then tiny else Dn
                let Cn = b' + a'*C
                let C' = if Cn=0.0 then tiny else Cn
                let delta = C'/D'
                let f' = f*delta
                if abs (delta-1.0) < Epsilon then f'
                else gFrac f' a' b' C' D' (m+1.0)
            (x**(s-1.0))*(exp (-x))/(gFrac fnot anot bnot Cnot Dnot 1.0)
        else if x<s then // this is better computed by lower incomplete gamma
                (gamma s) - (igamma_lo s x) 
        else // this uses a weird summation leveraged from boost documentation
                let n = 100 // iterations of the summation
                Seq.init n (fun i -> float i)   |> Seq.fold (fun sum k -> sum + ((-1.0)**k)*(x**k)/(s+k)/(fact k)) 0.0
                                                |> (*) (x**s)
                                                |> (+) (((gamma (s+1.0))-(x**s)-2.0)/s) 
    and igamma_lo s x = 
        if s<=0.0 then failwith "incomplete gamma function cannot handle s<=0.0 yet"
        else if x<0.0 then failwith "lower integration limit is zero; x must be larger than 0.0"
        else if x=0.0 then 0.0 
        else if x>s && x>1.1 then // this is better computed by upper incomplete gamma
            (gamma s) - (igamma_up s x)
        else // this performs best with the series representation of the lower gamma
                let n=100 // how far to carry out the summation? ... only n=15 for regular gamma ... 
                Seq.init n (fun i -> float i)   |> Seq.fold (fun sum k -> sum + ((x**(float n))/(s**(k+1.0)))) 0.0
                                                |> (*) ((x**s)*(exp (-x)))


    /// Lanczos approximation developed by Pugh in http://bh0.physics.ubc.ca/People/matt/Doc/ThesesOthers/Phd/pugh.pdf
    /// for discussion of what makes this method unique from Godfrey's, check out the c++ boost docs: 
    /// http://www.boost.org/doc/libs/1_37_0/libs/math/doc/sf_and_dist/html/math_toolkit/backgrounders/lanczos.html
    /// boost says that Pugh's version is more accurate; we use it here for beta functions only
    let lanczos2 x = 
        let n0 = 56906521.91347156388090791033559122686859
        let num = [| 103794043.1163445451906271053616070238554 ; //11 terms 
                86363131.28813859145546927288977868422342 ;
                43338889.32467613834773723740590533316085 ;
                14605578.08768506808414169982791359218571 ;
                3481712.15498064590882071018964774556468 ;
                601859.6171681098786670226533699352302507 ;
                75999.29304014542649875303443598909137092 ;
                6955.999602515376140356310115515198987526 ;
                449.9445569063168119446858607650988409623 ;
                19.51992788247617482847860966235652136208 ;
                0.5098416655656676188125178644804694509993 |]
        let n12 = 0.006061842346248906525783753964555936883222
        let d0 = 0.0
        let den = [|39916800.0; 
                120543840.0; 
                150917976.0; 
                105258076.0; 
                45995730.0; 
                13339535.0; 
                2637558.0; 
                357423.0; 
                32670.0; 
                1925.0; 
                66.0|] //11 terms
        let d12 = 1.0
        if x<1.0 then // start at end, multiply by z, foldBack
                let ns = Array.foldBack (fun e agg -> (agg*x)+e) num n12
                let ds = Array.foldBack (fun e agg -> (agg*x)+e) den d12
                ((ns*x)+n0)/((ds*x)+d0)
        else // start at front, divide by z, fold
                let ns = Array.fold (fun agg e -> (agg/x)+e) n0 num
                let ds = Array.fold (fun agg e -> (agg/x)+e) d0 den
                ((ns/x)+n12)/((ds/x)+d12)

    
    /// Beta function implemented that uses *some* of boost's optimizations for c++
    /// only recursive because for low alpha,beta <1 bootstrap the alpha,beta up above 1
    /// http://www.boost.org/doc/libs/1_39_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_beta/beta_function.html
    let rec beta a b = 
        let x = if a<b then b else a // we want a>b
        let y = if a<b then a else b
        let z = x+y
        if x<0.0 || y<0.0 then failwith "a and b must be larger than zero" // some special cases
        else if x=1.0 then 1.0/y
        else if x<1.0 then (z*(z+1.0)/(x*y))*(beta (x+1.0) (y+1.0)) // using identity to move to more stable location
        else if y<1.0 then ((x+y)/y)*(beta x (y+1.0))
        else if z-x<Epsilon then gamma y
        else if z-y<Epsilon then gamma x
        else // here is the meat of the function
            let g=6.024680040776729583740234375 // separately optimized by c++ boost guys
            let xg = x+g-0.5
            let yg = y+g-0.5
            let xyg = x+y+g-0.5
            let xmy = x-y-0.5
            let l1 = (lanczos2 x) // boost's implementation of lanczos "lanczos2"
            let l2 = (lanczos2 y)
            let l3 = (lanczos2 z)
            let lterm = l1*l2/l3
            let pterm = if (abs (b*xmy)) < (xyg*100.0) && x > 100.0 then exp (xmy*(log (-b/xyg))) else (xg/xyg)**xmy
            let pterm2 = if xyg > 1e10 then ((xg/xyg)*(yg/xyg))**y else (xg*yg/xyg/xyg)**y // just handling possible overflow case
            sqrt(E/yg)*pterm*pterm2*lterm


    /// Lower incomplete normalized beta function; this is a desperate man's implementation of the algorithms found in 
    /// http://www.boost.org/doc/libs/1_39_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_beta/ibeta_function.html :)
    /// note: to denormalize, all you do is transform with regular beta a b; higher incomplete beta is not very useful ... 
    let rec ibeta_lo a b x = 
        // main rules are that a and b must be larger than zero, x must be in range [0,1] (this is normalized function!)
        if a<0.0 || b<0.0 || x<0.0 || x>1.0 then failwith "remember: ibeta_lo (alpha) (beta) (x); a and b must be larger than zero, x must be in range [0,1]"
        // if we have natural numbers<40, faster and more accurate to use the binomial approximation
        else if a<40.0 && b<40.0 && abs(a-round(a))<Epsilon && abs(b-round(b))<Epsilon then binomial (int (a+b-1.0)) (int (a-1.0)) x
        // this is for min a b >1 in c++ boost docs http://www.boost.org/doc/libs/1_39_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_beta/ibeta_function.html
        // it is an extension of TOMS 708 aka DiDonato and Morris, "Algorithm 708 Significant digit computation of the 
        // incomplete beta function ratios" in ACM 1992 see http://portal.acm.org/citation.cfm?doid=131766.131776
        else 
            let tiny = 10e-30 // tiny shift for fnot, Cnot, Dnot to prevent instability around zero
            let anot = 0.0 // alpha reduces to zero for m=0.0
            let bnot = ((a)*(a-((a+b)*x)+1.0)/(a+1.0)) // beta reduces to this for m=0.0
            let fnot = if bnot<tiny then tiny else bnot // use tiny shift to prevent C and D from hitting zero
            let Cnot = fnot 
            let Dnot = 0.0 
            let rec gFrac f _ (*alpha*) _ (*beta*) C D i = 
                let m = float i
                let alpha' = (a+m-1.0)*(a+b+m-1.0)*m*(b-m)*(x**2.0)/((a+(2.0*m)-1.0)**2.0)
                let beta' = m + (m*(b-m)*x/(a+(2.0*m)-1.0)) + ((a+m)*(a-((a+b)*x)+1.0+(m*(2.0-x)))/(a+(2.0*m)+1.0))
                let D' = alpha'*D |> (+) beta' |> (fun y -> if y<tiny then 1.0/tiny else 1.0/y)
                let C' = alpha'/C |> (+) beta' |> (fun y -> if y<tiny then tiny else y)
                let delta = C'*D'
                let fn = f*delta
                if (abs (delta-1.0)) < Epsilon then fn
                else if i>500 then failwithf "continued fraction has not converged after %i interations" i
                else gFrac fn alpha' beta' C' D' (i+1)
            let pterm = (x**a)*((1.0-x)**b)/(beta a b)
            gFrac fnot anot bnot Cnot Dnot 1 |> (/) pterm // division accounts for real alpha0=1.0; took this trick from c++ boost docs
