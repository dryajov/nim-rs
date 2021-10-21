{.push raises: [Defect].}

{.deadCodeElim: on.}

import std/math

const
  Exp* {.intdefine.} = 8'u  # can be redifined in powers of two -2,  4, 8, 16, 32
                            # keep in mind that going above 32 is not very practical
  Char* = 2'u               # assume GF(2)
  Order* = (Char ^ Exp)
  Degree* = Order - 1'u

template bitsToUint*(bits: untyped): untyped =
  when bits == 8: uint8
  elif bits == 16: uint16
  elif bits == 32: uint32
  else: {.fatal: "bits have to be 8, 16 or 32!".}

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

type
  GFUint* = distinct bitsToUint(Exp)
