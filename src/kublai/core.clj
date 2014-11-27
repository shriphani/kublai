(ns kublai.core
  (:refer-clojure :exclude [* - + == /])
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.operators :refer :all]
            [kublai.arpack :as arpack]))

(defn symmetric-eigen-decomposition
  [M n]
  (arpack/arpack-symmetric-eigen-decomposition #(mmul M %)
                                               (-> M shape first)
                                               n))

(defn non-symmetric-eigen-decomposition
  [M n]
  (arpack/arpack-non-symmetric-eigen-decomposition #(mmul M %)
                                                   (-> M shape first)
                                                   n))

(defn eigs
  ([M n]
     (if (symmetric? M)
       (eigs M n :symmetric)
       (eigs M n :not-symmetric)))
  ([M n option]
     (cond (= option :symmetric)
           (symmetric-eigen-decomposition M n)

           (= option :non-symmetric)
           (non-symmetric-eigen-decomposition M n)

           :else
           (throw
            (Exception. "Supported options are :symmetric and :non-symmetric")))))

(defn svd
  "Args: we expect a matrix and the number
   of singular values you want to retrieve"
  [M n]
  (let [gram-matrix (* (transpose M)
                       M)

        cov-matrix  (* M (transpose M))

        gram-dimension (-> gram-matrix shape first)
        cov-dimension  (-> cov-matrix shape first)
        
        {U :Q lambda :A} (symmetric-eigen-decomposition cov-matrix n)
        {V :Q lambda2 :A} (symmetric-eigen-decomposition gram-matrix n)]
    {:U U
     :V* (transpose V)
     :S (sqrt lambda)}))
