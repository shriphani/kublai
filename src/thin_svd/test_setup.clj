(ns thin-svd.test-setup
  (:refer-clojure :exclude [* - + == /])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]))

(set-current-implementation :vectorz)

(defn test-svd
  "Build a matrix and svd it"
  []
  (matrix [[1 2 3 4]
           [5 6 7 8]
           [9 10 11 12]
           [13 14 15 16]]))

