{.push raises: [Defect].}

# {.experimental: "codeReordering".}

import std/[sequtils, math, typetraits]

const
  Char* {.intdefine.} = 2'u
  Exp* {.intdefine.} = 8'u
  Order* = (Char ^ Exp)
  Degree* = Order - 1'u

proc multNoLUT*(
  x, y, prim: SomeUnsignedInt,
  fieldSize = Order,
  carryless = true): SomeUnsignedInt =
  ## Galois Field integer multiplication using Russian
  ## Peasant Multiplication algorithm (faster than the
  ## standard multiplication + modular reduction). If
  ## prim is 0 and carryless=False, then the function
  ## produces the result for a standard integers multiplication
  ## (no carry-less arithmetics nor modular reduction).
  ##

  var
    r = 0'u
    y = y
    x = x

  while y > 0: # while y is above 0
    if bool(y and 1):
      # y is odd, then add the corresponding x to r
      #
      # Note: that since we're in GF(2), the addition
      # is in fact an XOR (very important because in GF(2)
      # the multiplication and additions are carry-less,
      # thus it changes the result!).
      r = if carryless: r xor x else: r + x

    y = y shr 1
    x = x shl 1
    if prim > 0'u and bool(x and fieldSize):
      # GF modulo: if x >= field size then apply modular
      # reduction using the primitive polynomial
      x = (x xor prim)

  r

proc rwhPrimes1(n: int): seq[uint] =
  # http://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
  ## Returns  a list of primes < n
  ##

  # TODO: rewrite this!
  let
    size = ceil(n / 2).int
    primeRange = floor((float(n).pow 0.5) + 1).int

  var sieve = repeat(true, size - 1)
  for i in countup(3, primeRange - 1, 2):
    if sieve[floor(i/2).int]:
      for j in countup(floor((i*i)/2).int, sieve.high, i):
        sieve[j] = false

  result = @[2'u]
  for i in 1..<(size - 1):
     if sieve[i]:
       result.add((2 * i + 1).uint)

proc findPrimePolys(
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
      x = multNoLUT(x, generator, prim, (Degree + 1'u))

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

template GFUintOp*(typ: type) {.dirty.} =
  proc `+`*(x: typ, y: uint): typ {.borrow, noSideEffect.}
  proc `+`*(x: uint, y: typ): typ {.borrow, noSideEffect.}

  proc `-`*(x: typ, y: uint): typ {.borrow, noSideEffect.}
  proc `-`*(x: uint, y: typ): typ {.borrow, noSideEffect.}

  # Not closed over type in question (Slot or Epoch)
  proc `mod`*(x: typ, y: typ): typ {.borrow, noSideEffect.}
  proc `mod`*(x: typ, y: uint): uint {.borrow, noSideEffect.}

  proc `xor`*(x: typ, y: typ): typ {.borrow, noSideEffect.}
  proc `xor`*(x: typ, y: uint): uint {.borrow, noSideEffect.}

  proc `div`*(x: typ, y: uint): uint {.borrow, noSideEffect.}
  proc `div`*(x: uint, y: typ): uint {.borrow, noSideEffect.}

  proc `*`*(x: typ, y: uint): uint {.borrow, noSideEffect.}

  proc `+=`*(x: var typ, y: typ) {.borrow, noSideEffect.}
  proc `+=`*(x: var typ, y: uint) {.borrow, noSideEffect.}
  proc `-=`*(x: var typ, y: typ) {.borrow, noSideEffect.}
  proc `-=`*(x: var typ, y: uint) {.borrow, noSideEffect.}

  # Comparison operators
  proc `<`*(x: typ, y: typ): bool {.borrow, noSideEffect.}
  proc `<`*(x: typ, y: uint): bool {.borrow, noSideEffect.}
  proc `<`*(x: uint, y: typ): bool {.borrow, noSideEffect.}

  proc `<=`*(x: typ, y: typ): bool {.borrow, noSideEffect.}
  proc `<=`*(x: typ, y: uint): bool {.borrow, noSideEffect.}
  proc `<=`*(x: uint, y: typ): bool {.borrow, noSideEffect.}

  proc `==`*(x: typ, y: typ): bool {.borrow, noSideEffect.}
  proc `==`*(x: typ, y: uint): bool {.borrow, noSideEffect.}
  proc `==`*(x: uint, y: typ): bool {.borrow, noSideEffect.}

  # Nim integration
  proc `$`*(x: typ): string {.borrow, noSideEffect.}

type
  GFUint* = distinct uint

const
  PrimePoly* {.intdefine.} = findPrimePolys(Char, Exp, true)[0]
  (GFExp*, GFLog*) = block:
    var
      gfExp = newSeq[GFUint](Degree * 2) # Anti-log (exponential) table.
      gfLog = newSeq[GFUint](Order)      # Log table, log[0] is impossible and thus unused

    # For each possible value in the galois field 2^8, we will pre-compute
    # the logarithm and anti-logarithm (exponential) of this value
    #
    # To do that, we generate the Galois Field F(2^p) by building a list
    # starting with the element 0 followed by the (p-1) successive powers of
    # the generator a : 1, a, a^1, a^2, ..., a^(p-1).
    var x = 1'u
    for i in 0..<Degree:
      gfExp[i] = x.GFUint # compute anti-log for this value and store it in a table
      gfLog[x] = i.GFUint # compute log at the same time
      x = multNoLUT(x, Char, PrimePoly, Order)

    # Optimization: double the size of the anti-log table so
    # that we don't need to mod 255 to stay inside the bounds
    # (because we will mainly use this table for the multiplication
    # of two GF numbers, no more).
    gfExp[gfExp.len/2..gfExp.high] = gfExp[0..gfExp.high/2]

    (gfExp, gfLog)

GFUintOp GFUint

proc `[]`*[T](a: seq[T], i: GFUint): T = a[i]

proc `+`*(x, y: GFUint): GFUint =
  x xor y

proc `-`*(x, y: GFUint): GFUint =
  # in binary galois field, substraction
  # is just the same as addition (since we mod 2)
  x xor y

proc `*`*(x, y: GFUint): GFUint =
  if x == 0 or y == 0:
    return 0.GFUint

  GFExp[((GFLog[x] + GFLog[y]) mod Degree)]

proc `div`*(x, y: GFUint): GFUint =
  if y == 0:
    # TODO: use DivByZeroDefect once we drop 1.2.6
    raise newException(DivByZeroError, "Can't divide by 0!")

  if x == 0:
    return 0.GFUint

  GFExp[((GFLog[x] + Degree) - GFLog[y]) mod Degree]

proc `/`*(x, y: GFUint): GFUint =
  x div y

proc `^`*(x: GFUint, power: GFUint | int): GFUint =
  GFExp[(GFLog[x] * power.GFUint) mod (Order - 1)]

func inverse*(x: GFUint): GFUint =
  GFExp[(Degree - GFLog[x])]
