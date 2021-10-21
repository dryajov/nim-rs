import sequtils
import pkg/unittest2
import pkg/stew/byteutils

import rs

suite "Test Reed-Solomon Coding":
  test "Generator polynomials":
    check:
      generator(10) == @[
        1, 216, 194, 159, 111,
        199, 94, 95, 113, 157,
        193].mapIt( it.GFUint )

      generator(8) == @[
        1, 255, 11, 81, 54, 239,
        173, 200, 24].mapIt( it.GFUint )

  test "Encode":
    let msg = @[
      72, 101, 108, 108, 111, 32, 82, 101,
      101, 100, 45, 83, 111, 108, 111, 109,
      111, 110, 33].mapIt( it.GFUint )

    check encode(msg, 10) == @[
      72, 101, 108, 108, 111, 32, 82, 101,
      101, 100, 45, 83, 111, 108, 111, 109,
      111, 110, 33, 10, 54, 200, 1, 174,
      73, 223, 252, 169, 147].mapIt( it.GFUint )

  test "Decode":
    let msg = @[
      72, 101, 108, 108, 111, 32, 82, 101,
      101, 100, 45, 83, 111, 108, 111, 109,
      111, 110, 33, 10, 54, 200, 1, 174,
      73, 223, 252, 169, 147].mapIt( it.GFUint )

    let (correct, code) = correctMsg(msg, 10)
    check string.fromBytes(correct.mapIt( it.byte )) == "Hello Reed-Solomon!"
