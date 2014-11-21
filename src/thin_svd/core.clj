(ns thin-svd.core
  (:refer-clojure :exclude [* - + == /])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all])
  (:import [com.github.fommil.netlib ARPACK]
           [java.util Arrays]
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
    (let [w (array
             (vec workd))]

      (while (not= 99 (.val ido))      
        (if-not (and (not= 1 (.val ido))
                     (not= -1 (.val ido)))
          (let [input-offset (- (aget ipntr 0) 1)
                output-offset (- (aget ipntr 1) 1)

                x (subvector w
                             input-offset
                             (min dimension
                                  (- (-> w shape first)
                                     input-offset)))
                y (subvector w
                             output-offset
                             (min dimension
                                  (- (-> w shape first)
                                     output-offset)))
                
                new-y (mul x)]
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
                     info))))

      (let [d (double-array (.val nev))
            select (boolean-array ncv)
            z (Arrays/copyOfRange v 0 (* (.val nev) dimension))]
        (.dseupd arpack
                 true
                 "A"
                 select
                 d
                 z
                 dimension
                 0.0
                 bmat
                 dimension
                 which
                 nev
                 tol
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

        (let [computed (aget iparam 4)
              eigenpairs (sort-by
                          first
                          (map-indexed
                           (fn [i x]
                             [x
                              (vec
                               (Arrays/copyOfRange z (* i dimension) (+ dimension (* i dimension))))])
                           (Arrays/copyOfRange d 0 computed)))
              eigenvector-matrix (matrix (map second eigenpairs))
              eigenvalue-matrix  (diagonal-matrix (map first eigenpairs))]
          {:Q eigenvector-matrix
           :A eigenvalue-matrix})))))
