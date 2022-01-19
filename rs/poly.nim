{.push raises: [Defect].}

{.deadCodeElim: on.}

import std/sequtils

import ./gf
export gf

template addSubs(p, q: openArray[GFSymbol]): seq[GFSymbol] =
  ## Perform addition/substraction which in GF(2^p) are the same
  ##
  var
    r = newSeq[GFSymbol](max(p.len, q.len))

  for i in 0..<p.len:
    r[i + r.len - p.len] = p[i]

  for i in 0..<q.len:
    r[i + (r.len - q.len)] = r[i + (r.len - q.len)] + q[i]

  r

template `+`*(p, q: sink openArray[GFSymbol]): seq[GFSymbol] =
  move addSubs(p, q)

template `-`*(p, q: openArray[GFSymbol]): seq[GFSymbol] =
  move addSubs(p, q)

proc `*`*(p, q: sink openArray[GFSymbol]): seq[GFSymbol] =
  ## Multiply two polynomials, inside Galois Field
  ##
  ## Optimized function by precomputation of log.
  ##

  var
    # Pre-allocate the result array
    r = newSeq[GFSymbol](p.len + q.len - 1)

    # Precompute the logarithm of p
    lp = p.mapIt(GFLog[it.uint])

  # Compute the polynomial multiplication
  # Just like the outer product of two vectors,
  # we multiply each coefficients of p with all
  # coefficients of q
  for j in 0..<q.len:
    if q[j] != 0: # log(0) is undefined, we need to check that
      let lq = GFLog[q[j].uint] # Optimization: precache the logarithm of the current coefficient of q
      for i in 0..<p.len:
        if p[i] != 0.GFSymbol: # log(0) is undefined, need to check that...
          r[i + j] = (r[i + j] xor GFExp[(lp[i] + lq) mod Degree].GFSymbol)

  return move r

proc mulSimple*(p, q: sink openArray[GFSymbol]): seq[GFSymbol] =
  ## Multiply two polynomials in a Galois Field
  ##
  ## simple equivalent way of multiplying two polynomials
  ## without precomputation, but slower
  ##

  # Pre-allocate the result array
  var r = newSeq[GFSymbol](max(p.len, q.len))
  # Compute the polynomial multiplication
  for j in 0..<q.len:
    for i in 0..<p.len:
      r[i + j] = r[i + j] xor p[i] * q[j]

  return move r

template scale*(p: sink seq[GFSymbol], x: int | GFSymbol): seq[GFSymbol] =
  p.mapIt(it * x)

template neg*[T](poly: openArray[GFSymbol]): openArray[GFSymbol] =
  ## Returns the polynomial with all coefficients negated.
  ## In GF(2^p), negation does not change the coefficient,
  ## so we return the polynomial as-is.
  ##

  poly

proc `div`*(
  dividend,
  divisor: sink openArray[GFSymbol]): (seq[GFSymbol], seq[GFSymbol]) =
  ## Fast polynomial division by using Extended Synthetic Division and
  ## optimized for GF(2^p) computations - doesn't work with standard
  ## polynomials outside of this galois field, see the Wikipedia article
  ## for generic algorithm.
  ##

  # NOTE: this function expects polynomials to follow the opposite
  # convention at decoding: the terms must go from the biggest to lowest
  # degree - while most other functions here expect a list from lowest
  # to biggest degree
  #
  # eg: 1 + 2x + 5x^2 = [5, 2, 1], NOT [1, 2, 5]

  var
    res = @dividend

  for i in 0..res.len - divisor.len:
    if res[i] == 0.GFSymbol:
      continue

    # skip the first coefficient of the divisor -
    # only used to normalize the dividend coefficient
    for j in 1..<divisor.len:
      if divisor[j] == 0.GFSymbol: # log(0) is undefined
        continue

      res[i + j] = res[i + j] xor (divisor[j] * res[i])

  return (
    res[0..res.len - divisor.len],
    res[res.len - divisor.high..res.high]) # return quotient, remainder.

template `/`*(
  dividend,
  divisor: sink openArray[GFSymbol]): (seq[GFSymbol], seq[GFSymbol]) =
  dividend div divisor

proc eval*(poly: sink openArray[GFSymbol], x: GFSymbol | int): GFSymbol =
  ## Evaluates a polynomial in GF(n^p) given the value for x.
  ## This is based on Horner's scheme for maximum efficiency.
  ##
  var y = poly[0]
  for i in 1..<poly.len:
    y = ((y * x.GFSymbol) xor poly[i])

  return move y
