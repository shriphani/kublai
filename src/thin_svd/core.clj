(ns thin-svd.core
  (:refer-clojure :exclude [* - + == /])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all])
  (:import [com.github.fommil.netlib ARPACK]
           [org.netlib.util intW doubleW]))

(defn eigen-decomposition
  "Args:
    mul: multiply routine
    dimension: dimension of the matrix
    k: number of eigenvalues
    tol: tolerance
    max-iterations: number of iterations to stop after"
  [mul dimension num-evs tol max-iterations]

  (let [arpack (ARPACK/getInstance)

        tolW   (doubleW. tol)
        nev    (intW. num-evs)

        ncv    (min (* 2 num-evs)
                    dimension)

        bmat   "I"
        which  "LM"

        iparam (int-array 11)

        _ (aset-int iparam 0 1)
        _ (aset-int iparam 2 max-iterations)
        _ (aset-int iparam 6 1)

        ido (intW. 0)
        info (intW. 0)
        resid (double-array dimension)
        v     (double-array (* dimension
                               ncv))
        workd (double-array (* ncv
                               3))
        workl (double-array (* ncv
                               (+ 8 ncv)))
        ipntr (int-array 11)]
    
    (.dsaupd arpack
             ido
             bmat
             dimension
             which
             (.val nev)
             tolW
             resid
             ncv
             v
             dimension
             iparam
             ipntr
             workd
             workl
             (count workl)
             info)
    (let [w (vector workd)]
      )))
