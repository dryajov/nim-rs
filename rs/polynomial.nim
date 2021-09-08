import std/sequtils

import ./gf
export gf

proc `+`*(p, q: seq[GFUint]): seq[GFUint] =
  var
    r = newSeq[GFUint](max(p.len, q.len))

  for i in 0..<p.len:
    r[i + r.len - p.len] = p[i]

  for i in 0..<q.len:
    r[i+r.len - q.len] = r[i+r.len - q.len] + q[i]

  return r

proc `*`*(p, q: seq[GFUint]): seq[GFUint] =
  ## Multiply two polynomials, inside Galois Field
  ## (but the procedure is generic).
  ##
  ## Optimized function by precomputation of log.
  ##

  var
    # Pre-allocate the result array
    r = newSeq[GFUint](p.len + q.len - 1)

    # Precompute the logarithm of p
    lp = p.mapIt(GFLog[it])

  # Compute the polynomial multiplication
  # Just like the outer product of two vectors,
  # we multiply each coefficients of p with all
  # coefficients of q
  for j in 0..<q.len:
    if q[j].uint != 0: # log(0) is undefined, we need to check that
      let lq = GFLog[q[j]] # Optimization: precache the logarithm of the current coefficient of q
      for i in 0..<p.len:
        if p[i] != 0.GFUint: # log(0) is undefined, need to check that...
          r[i + j] = (r[i + j] xor GFExp[lp[i] + lq].GFUint)

  return r

proc mulSimple*(p, q: seq[GFUint]): seq[GFUint] =
  ## Multiply two polynomials in a Galois Field
  ##
  ## simple equivalent way of multiplying two polynomials
  ## without precomputation, but slower
  ##

  # Pre-allocate the result array
  var r = newSeq[GFUint](max(p.len, q.len))
  # Compute the polynomial multiplication
  for j in 0..<q.len:
    for i in 0..<p.len:
      r[i + j] = r[i + j] xor p[i] * q[j]

  return r

proc scale*(p: seq[GFUint], x: int | GFUint): seq[GFUint] =
  p.mapIt(it * x)

proc neg*(poly: seq[GFUint]): seq[GFUint] =
  ## Returns the polynomial with all coefficients negated.
  ## In GF(2^p), negation does not change the coefficient,
  ## so we return the polynomial as-is.
  ##

  # TODO: we support arbitrary GF(x^p) fields
  poly

proc `div`*(
  dividend, divisor: seq[GFUint]): (seq[GFUint], seq[GFUint]) =
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
    res = dividend

  # resize to hold both the result and remainder
  res.setLen(dividend.len + divisor.high)
  for i in 0..<dividend.len:
    let coef = res[i]
    if coef == 0.GFUint:
      continue

    # skip the first coeficient of the divisor -
    # only used to normalize the divident coeficient
    for j in 1..<divisor.len:
      if divisor[j] == 0.GFUint: # log(0) is undefined
        continue

      res[i + j] = res[i + j] xor (divisor[j] * coef)

  return (res[0..dividend.high], res[dividend.len..res.high]) # return quotient, remainder.

proc `/`*(dividend, divisor: seq[GFUint]): (seq[GFUint], seq[GFUint]) =
  dividend div divisor

proc eval*(poly: seq[GFUint], x: GFUint | int): GFUint =
  ## Evaluates a polynomial in GF(2^p) given the value for x.
  ## This is based on Horner's scheme for maximum efficiency.
  ##
  var y = poly[0]
  for i in 1..poly.high:
    y = ((y * x.GFUint) xor poly[i])

  return y
