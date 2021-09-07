import std/sequtils

import ./gf

export gf

type
  GFPolynomial* = distinct seq[GFUint]

proc scale*(p: GFPolynomial, x: int): GFPolynomial =
  p.mapIt(it * x)

proc `+`*(p, q: GFPolynomial): GFPolynomial =
  var r = newSeq[SomeGFuint](max(len(p),len(q)))
  for i in 0..<len(p):
    r[i+len(r)-len(p)] = p[i]

  for i in 0..<len(q):
    r[i+len(r)-len(q)] = r[i+len(r)-len(q)] xor q[i]

  return r

proc `*`*(p, q: GFPolynomial): GFPolynomial =
  ## Multiply two polynomials, inside Galois Field
  ## (but the procedure is generic).
  ##
  ## Optimized function by precomputation of log.
  ##

  var
    # Pre-allocate the result array
    r = newSeq[GFUint](((len(p) + len(q) - 1)))

    # Precompute the logarithm of p
    lp = p.mapIt(GFLog[it])

  # Compute the polynomial multiplication (just like the outer product of two vectors, we multiply each coefficients of p with all coefficients of q)
  for j in 0..<len(q):
    var qj = q[j] # optimization: load the coefficient once
    if qj != 0: # log(0) is undefined, we need to check that
      var lq = GFLog[qj] # Optimization: precache the logarithm of the current coefficient of q
      for i in 0..<len(p):
        if p[i] != 0: # log(0) is undefined, need to check that...
          r[i + j] = (r[i + j] xor GFExp[lp[i] + lq]) # equivalent to: r[i + j] = gf_add(r[i+j], gf_mul(p[i], q[j]))

  return r

proc mulSimple*(p, q: seq[int]): seq[int] =
  ## Multiply two polynomials, inside Galois Field
  ##
  ## simple equivalent way of multiplying two polynomials
  ## without precomputation, but thus it's slower
  ##

  # Pre-allocate the result array
  var r = repeat(0, max(len(p),len(q)))
  # Compute the polynomial multiplication (just like the outer product of two vectors, we multiply each coefficients of p with all coefficients of q)
  for j in 0..<len(q):
    for i in 0..<len(p):
      r[i + j] = r[i + j] xor p[i] * q[j] # equivalent to: r[i + j] = gf_add(r[i+j], gf_mul(p[i], q[j])) -- you can see it's your usual polynomial multiplication

  return r

proc neg*(poly: GFPolynomial): int =
  ## Returns the polynomial with all coefficients negated.
  ## In GF(2^p), negation does not change the coefficient,
  ## so we return the polynomial as-is.
  ##

  # TODO: we support arbitrary GF(x^p) fields
  poly


proc `div`(dividend, divisor: GFPolynomial): (GFPolynomial, GFPolynomial) =
  ## Fast polynomial division by using Extended Synthetic Division and
  ## optimized for GF(2^p) computations (doesn't work with standard
  ## polynomials outside of this galois field, see the Wikipedia article
  ## for generic algorithm).
  ##

  # CAUTION: this function expects polynomials to follow the opposite convention at decoding:
  # the terms must go from the biggest to lowest degree (while most other functions here expect
  # a list from lowest to biggest degree). eg: 1 + 2x + 5x^2 = [5, 2, 1], NOT [1, 2, 5]

  var msg_out = dividend # Copy the dividend list and pad with 0 where the ecc bytes will be computed
  #normalizer = divisor[0] # precomputing for performance
  for i in 0..<(len(dividend) - (len(divisor)-1)):
    #msg_out[i] /= normalizer # for general polynomial division (when polynomials are non-monic), the usual way of using
                              # synthetic division is to divide the divisor g(x) with its leading coefficient, but not needed here.
    var coef = msg_out[i] # precaching
    if coef != 0: # log(0) is undefined, so we need to avoid that case explicitly (and it's also a good optimization).
      for j in 1..<divisor.len: # in synthetic division, we always skip the first coefficient of the divisior,
                                        # because it's only used to normalize the dividend coefficient
        if divisor[j] != 0: # log(0) is undefined
            msg_out[i + j] = msg_out[i + j] xor divisor[j] * coef # equivalent to the more mathematically correct
                                                          # (but xoring directly is faster): msg_out[i + j] += -divisor[j] * coef

  # The resulting msg_out contains both the quotient and the remainder, the remainder being the size of the divisor
  # (the remainder has necessarily the same degree as the divisor -- not length but degree == length-1 -- since it's
  # what we couldn't divide from the dividend), so we compute the index where this separation is, and return the quotient and remainder.
  let separator = (len(divisor) - 1)
  return (msg_out[0..separator], msg_out[separator..msg_out.high]) # return quotient, remainder.

proc Eval*(poly: GFPolynomial, x: int): GFPolynomial =
  ## Evaluates a polynomial in GF(2^p) given the value for x.
  ## This is based on Horner's scheme for maximum efficiency.
  ##

  # writeStackTrace()
  var y = poly[0]
  for i in 1..<len(poly):
    y = ((y * x) xor poly[i])

  return y

