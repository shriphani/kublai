(ns thin-svd.test-setup
  (:refer-clojure :exclude [* - + == /])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.linear :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [thin-svd.core :as core]))

(set-current-implementation :vectorz)

(defn test-svd
  "Build a matrix and svd it"
  []
  (let [M (matrix [[1 2 3 4]
                   [2 5 6 7]
                   [3 6 8 9]
                   [4 7 9 10]])
        mul #(mmul M %)
        dimension 4
        num-evs 2
        tol 0.0
        max-iterations 300]
    (core/eigen-decomposition mul dimension num-evs tol max-iterations)))

