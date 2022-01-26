import std/os
import gftype

when Exp == 8:
  import ./testgf8
elif Exp == 16:
  import ./testgf16
elif Exp == 18 and not defined(i386):
  import ./testgf18
elif Exp == 20 and not defined(i386):
  import ./testgf20
