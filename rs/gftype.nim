# import nimprof

{.push raises: [Defect].}

{.deadCodeElim: on.}

import std/math

template bitsToUint*(bits: untyped): untyped =
  when bits >= 2 and bits <= 8: uint8
  elif bits > 8 and bits <= 16: uint16
  elif bits > 16 and bits <= 32: uint32
  elif bits > 32 and bits <= 64: uint64

template GFUintOp*(typ, borowing: type) {.dirty.} =
  proc `+`*(x: typ, y: borowing): typ {.borrow, noSideEffect.}
  proc `+`*(x: borowing, y: typ): typ {.borrow, noSideEffect.}

  proc `-`*(x: typ, y: borowing): typ {.borrow, noSideEffect.}
  proc `-`*(x: borowing, y: typ): typ {.borrow, noSideEffect.}

  # Not closed over type in question (Slot or Epoch)
  proc `mod`*(x: typ, y: typ): typ {.borrow, noSideEffect.}
  proc `mod`*(x: typ, y: borowing): borowing {.borrow, noSideEffect.}

  proc `xor`*(x: typ, y: typ): typ {.borrow, noSideEffect.}
  proc `xor`*(x: typ, y: borowing): borowing {.borrow, noSideEffect.}

  proc `div`*(x: typ, y: borowing): borowing {.borrow, noSideEffect.}
  proc `div`*(x: borowing, y: typ): borowing {.borrow, noSideEffect.}

  proc `*`*(x: typ, y: borowing): borowing {.borrow, noSideEffect.}

  proc `+=`*(x: var typ, y: typ) {.borrow, noSideEffect.}
  proc `+=`*(x: var typ, y: borowing) {.borrow, noSideEffect.}
  proc `-=`*(x: var typ, y: typ) {.borrow, noSideEffect.}
  proc `-=`*(x: var typ, y: borowing) {.borrow, noSideEffect.}

  # Comparison operators
  proc `<`*(x: typ, y: typ): bool {.borrow, noSideEffect.}
  proc `<`*(x: typ, y: borowing): bool {.borrow, noSideEffect.}
  proc `<`*(x: borowing, y: typ): bool {.borrow, noSideEffect.}

  proc `<=`*(x: typ, y: typ): bool {.borrow, noSideEffect.}
  proc `<=`*(x: typ, y: borowing): bool {.borrow, noSideEffect.}
  proc `<=`*(x: borowing, y: typ): bool {.borrow, noSideEffect.}

  proc `==`*(x: typ, y: typ): bool {.borrow, noSideEffect.}
  proc `==`*(x: typ, y: borowing): bool {.borrow, noSideEffect.}
  proc `==`*(x: borowing, y: typ): bool {.borrow, noSideEffect.}

  # Nim integration
  proc `$`*(x: typ): string {.borrow, noSideEffect.}

const
  Exp* {.intdefine.} = 8'u  # can be redefined in powers of two - 8, ... 16...
                            # keep in mind that going above 32 is not very practical
                            # using LUT tables
  Char* = 2'u               # assume GF(2)
  Order* = (Char ^ Exp)
  Degree* = Order - 1'u

type
  GFUint* = distinct bitsToUint(Exp)  # used for GF arithmetic
  GFSymbol* = range[0.GFUint..Degree.GFUint]

  RSError* = object of CatchableError # Base error type

# TODO: maybe there is a way?
# proc `=copy`(dest: var GFSymbol; source: GFSymbol) {.error.} # disable copying
