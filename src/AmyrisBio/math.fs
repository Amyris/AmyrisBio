namespace Amyris.Bio

/// common high-level math module; this includes all supporting modules
module math = 

    // this is not particularly useful without further debugging; the reassignments 
    // for some reason assert a particular type to each function, based on how they are used within this module :(  EB 4/12/2011
    
    open math_common // common math functions
    open math_spf // special functions (gamma, beta, etc.)
    open math_stat // all other statistics functions
    open math_bhg // bionmial, hypergeometric code
    
    let mean = math_common.mean
    let ssq = math_common.ssq
    let var = math_common.var
    let sd = math_common.sd
    let median = math_common.median
    let fact = math_common.fact
    let hcf = math_common.hcf

    let hyperGeo = math_bhg.hyperGeo
    let hyperGeoCum = math_bhg.hyperGeoCum
    let binomial = math_bhg.binomial
    let binomialCum = math_bhg.binomialCum
    let primesTill = math_bhg.primesTill

    let gamma = math_spf.gamma
    let igamma_up = math_spf.igamma_up
    let igamma_lo = math_spf.igamma_lo
    let beta = math_spf.beta
    let ibeta_lo = math_spf.ibeta_lo

    let rnorm = math_stat.rnorm
    let rpois = math_stat.rpois
    let rexp = math_stat.rexp
    let rgamma = math_stat.rgamma
    let rbeta = math_stat.rbeta
    let rweib = math_stat.rweib

    let pbeta = math_stat.pbeta
    let pt = math_stat.pt
    let ttest = math_stat.tTest
    let qvalues = math_stat.qValues
    let holmes = math_stat.holmesCorrect