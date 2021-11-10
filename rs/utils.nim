{.push raises: [Defect].}

{.deadCodeElim: on.}

import std/[math, sequtils, tables]
import ./gftype

proc mulNoLUT*(
  x, y, prim: SomeUnsignedInt,
  order = Order,
  carryless = true): SomeUnsignedInt =
  ## Galois Field integer multiplication using Russian
  ## Peasant Multiplication algorithm - faster than the
  ## standard multiplication + modular reduction.
  ##
  ## If `prim` is 0 and `carryless` = false, then the function
  ## produces the result for a standard integers multiplication -
  ## no carry-less arithmetics nor modular reduction.
  ##

  var
    r = 0'u
    # use the smallest term as the leading since
    # it reduces the number of steps
    (y, x) = if y > x: (x, y) else: (y, x)

  while y > 0: # while y is above 0
    if bool(y and 1):
      # y is odd, then add the corresponding x to r
      #
      # Note: that since we're in GF(2), the addition
      # is in fact an XOR - very important because in GF(2)
      # the multiplication and additions are carry-less,
      # thus it changes the result!
      r = if carryless: r xor x else: r + x

    y = y shr 1
    x = x shl 1
    if prim > 0'u and bool(x and order):
      # if x >= field size then apply modular
      # reduction using the primitive polynomial
      x = (x xor prim)

  r

proc rwhPrimes1(n: int): seq[uint] =
  ## Returns  a list of primes < n
  ##

  # TODO: rewrite this!
  let
    size = ceil(n / 2).int
    primeRange = floor((float(n).pow 0.5) + 1).int

  var
    sieve = initTable[int, bool]()
  for i in countup(3, primeRange - 1, 2):
    if sieve.getOrDefault(floor(i/2).int, true):
      for j in countup(floor((i*i)/2).int, (size - 1), i):
        sieve[j] = false

  result = @[2'u]
  for i in 1..<(size - 1):
     if sieve.getOrDefault(floor(i/2).int, true):
       result.add((2 * i + 1).uint)

proc primePolys*(
  generator = Char,
  cExp = Exp,
  fastPrimes = true,
  single = true): seq[uint] =
  ## Prepare the finite field characteristic (2^p - 1),
  ## this also represent the maximum possible value
  ## in this field
  ##

  let
    DegreeNext = ((generator ^ (cExp + 1'u)) - 1'u)
    primCandidates = if fastPrimes:
      # generate maybe prime polynomials and check
      # later if they really are irreducible
      rwhPrimes1(DegreeNext.int).filterIt(
        it > Degree # filter out too small primes
      )
    else:
      # try each possible prime polynomial, but skip even numbers
      # (because divisible by 2 so necessarily not irreducible)
      toSeq(
        countup(Degree + 2, DegreeNext - 1, generator))
        .mapIt(it.uint)

  # Start of the main loop
  var correctPrimes: seq[uint]
  for prim in primCandidates: # try potential candidates primitive irreducible polys
    var
      # memory variable to indicate if a value was
      # already generated in the field (value at
      # index x is set to 1) or not (set to 0 by default)
      seen = newSeq[bool](Degree + 1)
      conflict = false # flag to know if there was at least one conflict

    # Second loop, build the whole Galois Field
    var x = 1'u
    for i in 0..<Degree:
      # Compute the next value in the field
      # (ie, the next power of alpha/generator)
      x = mulNoLUT(x, generator, prim, (Degree + 1'u))

      # Rejection criterion: if the value overflowed
      # (above Degree) or is a duplicate of a
      # previously generated power of alpha, then we
      # reject this polynomial (not prime)
      if x > Degree or seen[x.int]:
        conflict = true
        break
      else:
        # Else we flag this value as seen
        # (to maybe detect future duplicates),
        # and we continue onto the next power
        # of alpha
        seen[x] = true

    # End of the second loop: if there's no conflict
    # (no overflow nor duplicated value), this is a
    # prime polynomial!
    if not conflict:
      correctPrimes.add(prim)
      if single:
        return @[prim]

  # Return the list of all prime polynomials
  return correctPrimes
