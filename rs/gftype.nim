{.push raises: [Defect].}

import std/math

const
  Char* = 2'u               # assume GF(2)
  Exp* {.intdefine.} = 8'u  # can be redifined in powers of two -2,  4, 8, 16, 32
                            # keep in mind that going above 32 is not very practical
  Order* = (Char ^ Exp)
  Degree* = Order - 1'u

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
