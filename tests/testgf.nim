import pkg/unittest2

import ../rs/gf

suite "GF Arithmetic":

  test "Addition":
    check:
      50.GFUint + 50.GFUint == 0.GFUint
      150.GFUint + 200.GFUint == 94.GFUint
      5.GFUint + 233.GFUint == 236.GFUint

  test "Multiplication":
    check:
      50.GFUint * 50.GFUint == 109.GFUint
      # 150.GFUint * 200.GFUint == 118.GFUint
      # 5.GFUint * 233.GFUint == 233.GFUint

  test "Substraction":
    discard

  test "Divission":
    discard

  test "Exponentiation":
    discard

  test "Inverse":
    discard
