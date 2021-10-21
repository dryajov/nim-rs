{.push raises: [Defect].}

import std/[sequtils, math, typetraits]
import ./utils
import ./gftype

export gftype

proc initTables*(
  primePoly: SomeUnsignedInt,
  chr = Char,
  order = Order,
  degree = Degree): (seq[uint], seq[uint]) =
  var
    gfExp = newSeq[uint](degree * 2) # Anti-log (exponentiation) table.
    gfLog = newSeq[uint](order)      # Log table, log[0] is impossible and thus unused

  # For each possible value in the galois field 2^8, we will pre-compute
  # the logarithm and anti-logarithm (exponentiation) of this value
  #
  # To do that, we generate the Galois Field F(2^p) by building a list
  # starting with the element 0 followed by the (p-1) successive powers of
  # the generator a : 1, a, a^1, a^2, ..., a^(p-1).
  var x = 1'u
  for i in 0..<degree:
    gfExp[i] = x.uint # compute anti-log for this value and store it in a table
    gfLog[x] = i.uint # compute log at the same time
    x = mulNoLUT(x, chr, primePoly, order)

  # double the size of the anti-log table so that we don't
  # need to mod 255 to stay inside the bounds
  gfExp[gfExp.len / 2..gfExp.high] = gfExp[0..gfExp.high / 2]

  (gfExp, gfLog)

const
  PrimePoly* {.intdefine.} = primePolys(Char, Exp, true)[0]
  (GFExp*, GFLog*) = initTables(PrimePoly)

GFUintOp GFUint, bitsToUint(Exp)

proc `[]`*[T](a: seq[T], i: GFUint): T =
  a[i.int]

proc `+`*(x, y: GFUint): GFUint =
  x xor y

proc `-`*(x, y: GFUint): GFUint =
  # in binary galois field, substraction
  # is just the same as addition (since we mod 2)
  x xor y

proc `*`*(x, y: GFUint): GFUint =
  if x == 0 or y == 0:
    return 0.GFUint

  GFExp[((GFLog[x] + GFLog[y]) mod Degree)].GFUint

proc `div`*(x, y: GFUint): GFUint =
  if y == 0:
    # TODO: use DivByZeroDefect once we drop 1.2.6
    raise newException(DivByZeroError, "Can't divide by 0!")

  if x == 0:
    return 0.GFUint

  GFExp[((GFLog[x] + Degree) - GFLog[y]) mod Degree].GFUint

proc `/`*(x, y: GFUint): GFUint =
  x div y

proc `^`*(x: GFUint, power: int): GFUint =
  GFExp[(GFLog[x] * power.uint) mod Degree].GFUint

proc inverse*(x: GFUint): GFUint =
  GFExp[(Degree - GFLog[x])].GFUint
