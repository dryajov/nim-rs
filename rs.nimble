# Package

version       = "0.1.0"
author        = "Dmitriy Ryajov"
description   = "Nim Reed-Solomon erasure and error correcting coder"
license       = "MIT"
srcDir        = "src"

task testgf8, "Test in GF(2^8)":
  exec "nim c -r -d:release --hints:off -d:Exp=8 -d:PrimePoly=285 tests/testall.nim"
  rmFile "tests/testall".toExe

task testgf16, "Test in GF(2^16)":
  exec "nim c -r -d:release --hints:off -d:Exp=16 -d:PrimePoly=65593 tests/testall.nim"
  rmFile "tests/testall".toExe

task testgf18, "Test in GF(2^18)":
  exec "nim c -r -d:release --hints:off -d:Exp=18 -d:PrimePoly=262221 tests/testall.nim"
  rmFile "tests/testall".toExe

task testgf20, "Test in GF(2^20)":
  exec "nim c -r -d:release --hints:off -d:Exp=20 -d:PrimePoly=1048681 tests/testall.nim"
  rmFile "tests/testall".toExe

task test, "Tes Reed-Solomon":
  exec "nimble testgf8"
  exec "nimble testgf16"
  exec "nimble testgf18"
  exec "nimble testgf20"
