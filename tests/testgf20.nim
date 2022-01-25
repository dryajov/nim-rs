import std/unittest

import gf

suite "GF Arithmetic in GF " & $Exp:

  test "Addition":
    check:
      500000.GFSymbol + 500000.GFSymbol == 0.GFSymbol
      150000.GFSymbol + 200000.GFSymbol == 83120.GFSymbol
      5000.GFSymbol + 233000.GFSymbol == 236960.GFSymbol

  test "Multiplication":
    check:
      500000.GFSymbol * 500000.GFSymbol == 355167.GFSymbol
      150000.GFSymbol * 200000.GFSymbol == 411715.GFSymbol
      5000.GFSymbol * 233000.GFSymbol == 724362.GFSymbol

  test "Substraction":
    check:
      500000.GFSymbol - 500000.GFSymbol == 0.GFSymbol
      150000.GFSymbol - 200000.GFSymbol == 83120.GFSymbol
      5000.GFSymbol - 233000.GFSymbol == 236960.GFSymbol

  test "Divission":
    check:
      500000.GFSymbol / 500000.GFSymbol == 1.GFSymbol
      150000.GFSymbol / 200000.GFSymbol == 288686.GFSymbol
      5000.GFSymbol / 233000.GFSymbol == 921553.GFSymbol

  test "Exponentiation":
    check:
      500000.GFSymbol ^ 500000 == 458293.GFSymbol
      150000.GFSymbol ^ 200000 == 52046.GFSymbol
      5000.GFSymbol ^ 233000 == 535174.GFSymbol

  test "Inverse":
    check:
      GFSymbol(500000).inverse == 123976.GFSymbol
      GFSymbol(150000).inverse == 589900.GFSymbol
      GFSymbol(5000).inverse == 317.GFSymbol
