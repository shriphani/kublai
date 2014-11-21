(ns thin-svd.test-setup
  (:refer-clojure :exclude [* - + == /])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [thin-svd.core :as core]))

(set-current-implementation :vectorz)

(defn test-svd
  "Build a matrix and svd it"
  []
  (let [M (matrix [[1 2 3 4]
                   [5 6 7 8]
                   [9 10 11 12]
                   [13 14 15 16]])
        mul #(* M %)
        dimension 4
        num-evs 2
        tol 1e-6
        max-iterations 100]
    (core/eigen-decomposition mul dimension num-evs tol max-iterations)))

