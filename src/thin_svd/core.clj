(ns thin-svd.core
  (:refer-clojure :exclude [* - + == /])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]))

(defn eigen-decomposition
  "Args:
    mul: multiply routine
    dimension: dimension of the matrix
    k: number of eigenvalues
    tol: tolerance
    max-iterations: number of iterations to stop after"
  [mul dimension num-evs tol max-iterations]
  '*)
