# kublai

A Clojure library that implements truncated matrix decompositions for
<code>[core.matrix](https://github.com/mikera/core.matrix)</code> matrices.

Kublai uses the [ARPACK library](http://www.caam.rice.edu/software/ARPACK/) for implementing these decompositions.

The library is independent of the <code>core.matrix</code> implementation used.

To include it in your projects:

[![Clojars Project](http://clojars.org/kublai/latest-version.svg)](http://clojars.org/kublai)

## Motivation

See [my blog post on the topic]()

## Usage

To run a truncated symmetric-eigendecomposition:

```clojure
user> (def M (matrix [[1 2 3 4]
                   [2 5 6 7]
                   [3 6 8 9]
                   [4 7 9 10]]))
#'user/M
user> (use 'kublai.core :reload-all)
nil
user> (eigs M 2 :symmetric) ;; compute 2 eigenvectors for this matrix
{:Q [[-0.22593827269074584 -0.4432218615090191 -0.5727878807113498 -0.6514754961809404] [0.7253136654558885 0.3184697313242928 0.1424607347013554 -0.5934661371986827]], :A [[24.06253512439672 0.0] [0.0 -0.8054849155764637]]}
user> 
```

To run a truncated svd:

```clojure
user> (def M (matrix [[1 2 3 4]
                   [5 6 7 8]
                   [9 10 11 12]
                   [13 14 15 16]]))
#'user/M
user> (svd M 2)
{:U [[-0.13472212372225584 0.8257420598345273] [-0.3407576960799602 0.4288172018031381] [-0.5467932684376645 0.03189234377176592] [-0.7528288407953688 -0.365032514259624]], :V* [[0.4284123959267892 0.4743725155726848 0.5203326352185806 0.5662927548644766] [0.7186534763126667 0.27380780936493887 -0.17103785758268963 -0.6158835245304229]], :S [[38.62265683187287 0.0] [0.0 2.0713230668787377]]}
user> (clojure.pprint/pprint (svd M 2))
{:U
 [[-0.13472212372225592 -0.825742059834525]
  [-0.34075769607996026 -0.42881720180314464]
  [-0.5467932684376648 -0.03189234377175876]
  [-0.7528288407953688 0.3650325142596215]],
 :V*
 [[-0.4284123959267895
   -0.4743725155726852
   -0.5203326352185804
   -0.5662927548644764]
  [-0.7186534763126535
   -0.27380780936497917
   0.1710378575827312
   0.615883524530409]],
 :S [[38.62265683187287 0.0] [0.0 2.0713230668787403]]}
nil
```

The truncated SVDs perform very well on extremely large matrices. For
example, here's a timing test:

```clojure
(let [M1 (reshape (matrix (range 500000)) [10000 50])
        M2 (reshape (matrix (range 5000000)) [10000 500])
        M3 (reshape (matrix (range 5000000)) [1000 5000])]
    (time (kublai/svd M1 10))
    (time (kublai/svd M2 10))
    (time (kublai/svd M3 10)))
```

So for each of the matrices, we retrieve 10 singular vectors:

```
"Elapsed time: 17372.943 msecs"
"Elapsed time: 78085.404 msecs"
"Elapsed time: 41511.266 msecs"
```

## License

Copyright © 2014 Shriphani Palakodety

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
