{.push raises: [Defect].}

{.deadCodeElim: on.}

import std/typetraits
import ./utils
import ./gftype

export gftype

const
  PrimePoly* {.intdefine.} = 0 # this will be initialized either as a
                               # compile time define or by finding a
                               # suitable prime in the block bellow
  (GFExp*, GFLog*) = block:
    var
      gfExp {.noInit.}: array[Degree, uint] # Anti-log (exponentiation) table.
      gfLog {.noInit.}: array[Order, uint]  # Log table, log[0] is impossible and thus unused

    when PrimePoly > 0 and PrimePoly <= Degree:
      {.fatal: "-d:PrimePoly has to be larger or equal to the Degree of the field, set to `0` to source one automatically".}

    const
      primePoly = if PrimePoly == 0:
        primePolys(Char, Exp, true)[0]
      else:
        PrimePoly

    # For each possible value in the galois field 2^p, we will pre-compute
    # the logarithm and anti-logarithm (exponentiation) of this value
    #
    # To do that, we generate the Galois Field F(2^p) by building a list
    # starting with the element 0 followed by the (p-1) successive powers of
    # the generator a: 1, a, a^1, a^2, ..., a^(p-1).
    var x = 1'u
    for i in 0..<Degree:
      gfExp[i] = x.uint # compute anti-log for this value and store it in a table
      gfLog[x] = i.uint # compute log at the same time
      x = mulNoLUT(x = x, y = Char, prim = primePoly, order = Order)

    # double the size of the anti-log table so that we don't
    # need to mod 255 to stay inside the bounds

    (@gfExp & @gfExp, gfLog)

# TODO: This is only practical for GF(2^8)
  # (MulTable*, DivTable*) = block:
  #   var
  #     mulTable {.noInit.}: array[Order, array[Order, uint]] # multiplication table
  #     divTable {.noInit.}: array[Order, array[Order, uint]] # divission table

  #   for x in 0..<Order:
  #     for y in 0..<Order:
  #       mulTable[x][y] = GFExp[((GFLog[x] + GFLog[y]))]
  #       divTable[x][y] = GFExp[((GFLog[x] + Degree) - GFLog[y])]

  #   (mulTable, divTable)

GFUintOp GFUint, bitsToUint(Exp)

proc `[]`*[T](a: seq[T], i: GFSymbol): T {.inline.} =
  a[i.int]

proc `+`*(x, y: GFSymbol): GFSymbol {.inline.} =
  x xor y

proc `-`*(x, y: GFSymbol): GFSymbol {.inline.} =
  # in binary galois field, substraction
  # is just the same as addition (since we mod 2)
  x + y

proc `*`*(x, y: GFSymbol): GFSymbol {.inline.} =
  if x == 0 or y == 0:
    return 0.GFSymbol

  GFExp[((GFLog[x.uint] + GFLog[y.uint]))].GFSymbol
  # MulTable[x.uint][y.uint].GFSymbol

proc `div`*(x, y: GFSymbol): GFSymbol {.raises: [DivByZeroError], inline.} =
  if y == 0:
    # TODO: use DivByZeroDefect once we drop 1.2.6
    raise newException(DivByZeroError, "Can't divide by 0!")

  if x == 0:
    return 0.GFSymbol

  GFExp[((GFLog[x.int] + Degree) - GFLog[y.int])].GFSymbol
  # DivTable[x.uint][y.uint].GFSymbol

proc `/`*(x, y: GFSymbol): GFSymbol {.inline.} =
  x div y

proc `^`*(x: GFSymbol, power: int): GFSymbol {.inline.} =
  GFExp[(GFLog[x.int] * power.uint) mod Degree].GFSymbol

proc inverse*(x: GFSymbol): GFSymbol {.inline.} =
  GFExp[(Degree - GFLog[x.int]) mod Degree].GFSymbol
