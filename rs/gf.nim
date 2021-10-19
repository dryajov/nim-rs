{.push raises: [Defect].}

import std/[sequtils, math, typetraits]
import ./utils
import ./gftype

export gftype

const
  PrimePoly* {.intdefine.} = primePolys(Char, Exp, true)[0]
  (GFExp*, GFLog*) = initTables(PrimePoly)

GFUintOp GFUint

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
