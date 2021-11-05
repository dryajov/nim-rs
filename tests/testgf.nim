import pkg/unittest2

import gf

# TODO: Add tests for fields other than GF(2^8) and generator 0x11D

suite "GF Arithmetic":

  test "Addition":
    check:
      50.GFSymbol + 50.GFSymbol == 0.GFSymbol
      150.GFSymbol + 200.GFSymbol == 94.GFSymbol
      5.GFSymbol + 233.GFSymbol == 236.GFSymbol

  test "Multiplication":
    check:
      50.GFSymbol * 50.GFSymbol == 109.GFSymbol
      150.GFSymbol * 200.GFSymbol == 118.GFSymbol
      5.GFSymbol * 233.GFSymbol == 106.GFSymbol

  test "Substraction":
    check:
      50.GFSymbol - 50.GFSymbol == 0.GFSymbol
      150.GFSymbol - 200.GFSymbol == 94.GFSymbol
      5.GFSymbol - 233.GFSymbol == 236.GFSymbol

  test "Divission":
    check:
      50.GFSymbol / 50.GFSymbol == 1.GFSymbol
      150.GFSymbol / 200.GFSymbol == 22.GFSymbol
      5.GFSymbol / 233.GFSymbol == 185.GFSymbol

  test "Exponentiation":
    check:
      50.GFSymbol ^ 50 == 116.GFSymbol
      150.GFSymbol ^ 200 == 193.GFSymbol
      5.GFSymbol ^ 233 == 255.GFSymbol

  test "Inverse":
    check:
      GFSymbol(50).inverse == 111.GFSymbol
      GFSymbol(150).inverse == 15.GFSymbol
      GFSymbol(5).inverse == 167.GFSymbol
