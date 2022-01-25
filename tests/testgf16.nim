import std/unittest

import gf

suite "GF Arithmetic in GF " & $Exp:

  test "Addition":
    check:
      65000.GFSymbol + 65000.GFSymbol == 0.GFSymbol
      596.GFSymbol + 512.GFSymbol == 84.GFSymbol
      5.GFSymbol + 1024.GFSymbol == 1029.GFSymbol

  test "Multiplication":
    check:
      596.GFSymbol * 1000.GFSymbol == 7055.GFSymbol
      1500.GFSymbol * 2000.GFSymbol == 53994.GFSymbol
      50000.GFSymbol * 233.GFSymbol == 49399.GFSymbol

  test "Substraction":
    check:
      50000.GFSymbol - 50000.GFSymbol == 0.GFSymbol
      15000.GFSymbol - 20000.GFSymbol == 29880.GFSymbol
      5.GFSymbol - 23300.GFSymbol == 23297.GFSymbol

  test "Divission":
    check:
      5000.GFSymbol / 50000.GFSymbol == 63916.GFSymbol
      1500.GFSymbol / 2000.GFSymbol == 4192.GFSymbol
      5.GFSymbol / 23300.GFSymbol == 6690.GFSymbol

  test "Exponentiation":
    check:
      50000.GFSymbol ^ 50 == 21829.GFSymbol
      1500.GFSymbol ^ 200 == 11473.GFSymbol
      50.GFSymbol ^ 233 == 5600.GFSymbol

  test "Inverse":
    check:
      GFSymbol(5000).inverse == 21863.GFSymbol
      GFSymbol(15000).inverse == 43376.GFSymbol
      GFSymbol(50).inverse == 18682.GFSymbol
