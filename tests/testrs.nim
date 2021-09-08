import sequtils
import pkg/unittest2

import ../rs/rs

suite "Test Reed-Solomon Coding":
  test "Generator polynomials":
    check:
      generator(10) == @[
        1.GFUint, 216.GFUint, 194.GFUint,
        159.GFUint, 111.GFUint, 199.GFUint,
        94.GFUint, 95.GFUint, 113.GFUint,
        157.GFUint, 193.GFUint]

      generator(8) == @[
        1.GFUint, 255.GFUint, 11.GFUint,
        81.GFUint, 54.GFUint, 239.GFUint,
        173.GFUint, 200.GFUint, 24.GFUint]

  test "Encode":
    let msg = @[
      72, 101, 108, 108, 111, 32, 82, 101,
      101, 100, 45, 83, 111, 108, 111, 109,
      111, 110, 33].mapIt(it.GFUint)

    check encode(msg, 10) == @[
      72, 101, 108, 108, 111, 32, 82, 101,
      101, 100, 45, 83, 111, 108, 111, 109,
      111, 110, 33, 10, 54, 200, 1, 174,
      73, 223, 252, 169, 147].mapIt(it.GFUint)
