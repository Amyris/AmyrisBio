(*** hide ***)
// This block of code is omitted in the generated HTML documentation. Use 
// it to define helpers that you do not want to show in the documentation.
#I "../../bin/AmyrisBio"

(**
Amyris.Bio is a library of useful functions for computational biology,
=====================================================================

All modules are available in the Amyris. namespace.  See individual modules for specific functionality.

*)
#r "AmyrisBio.dll"
open Amyris.Bio.biolib

// generate reverse complement of a DNA sequence array
let reverseComp = revComp ("GATTACA".ToCharArray())

let a = "GATTACA".ToCharArray()
let b = "GATACA".ToCharArray()

// Align two sequences using a seed mer size of 10
let ok,alignA,alignB = Amyris.Bio.smithWaterman.smithWaterman 2 a b

printfn "%s\n%s" alignA alignB
