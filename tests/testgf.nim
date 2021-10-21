import pkg/unittest2

import gf

# TODO: Add tests for fields other than GF(2^8) and generator 0x11D

suite "GF Arithmetic":

  test "Addition":
    check:
      50.GFUint + 50.GFUint == 0.GFUint
      150.GFUint + 200.GFUint == 94.GFUint
      5.GFUint + 233.GFUint == 236.GFUint

  test "Multiplication":
    check:
      50.GFUint * 50.GFUint == 109.GFUint
      150.GFUint * 200.GFUint == 118.GFUint
      5.GFUint * 233.GFUint == 106.GFUint

  test "Substraction":
    check:
      50.GFUint - 50.GFUint == 0.GFUint
      150.GFUint - 200.GFUint == 94.GFUint
      5.GFUint - 233.GFUint == 236.GFUint

  test "Divission":
    check:
      50.GFUint / 50.GFUint == 1.GFUint
      150.GFUint / 200.GFUint == 22.GFUint
      5.GFUint / 233.GFUint == 185.GFUint

  test "Exponentiation":
    check:
      50.GFUint ^ 50 == 116.GFUint
      150.GFUint ^ 200 == 193.GFUint
      5.GFUint ^ 233 == 255.GFUint

  test "Inverse":
    check:
      GFUint(50).inverse == 111.GFUint
      GFUint(150).inverse == 15.GFUint
      GFUint(5).inverse == 167.GFUint
