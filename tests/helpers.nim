import std/random
import std/sequtils

# Taken from stew/byteutils
func toBytes*(s: string): seq[byte] =
  ## Convert a string to the corresponding byte sequence - since strings in
  ## nim essentially are byte sequences without any particular encoding, this
  ## simply copies the bytes without a null terminator
  when nimvm:
    var r = newSeq[byte](s.len)
    for i, c in s:
      r[i] = cast[byte](c)
    r
  else:
    @(s.toOpenArrayByte(0, s.high))

func fromBytes*(T: type string, v: openArray[byte]): string =
  if v.len > 0:
    result = newString(v.len)
    when nimvm:
      for i, c in v:
        result[i] = cast[char](c)
    else:
      copyMem(addr result[0], unsafeAddr v[0], v.len)

proc randomAlpha*(len: int): seq[byte] =
  var
    read = 0
    res {.noinit.} = newSeq[byte](len)
    alpha {.noinit.} = toSeq(byte('A')..byte('z'))

  while read < len:
    shuffle(alpha)
    for a in alpha:
      if read >= len:
        break

      res[read] = a
      read.inc

  return res

template benchmark*(benchmarkName: string, code: untyped) =
  block:
    var countBytes {.inject.} = 0.uint64
    let t0 = epochTime()
    code
    let elapsed = epochTime() - t0
    let elapsedStr = elapsed.formatFloat(format = ffDecimal, precision = 3)
    let elapsedMillis = (elapsed * 60 * 1000)
    let mibs = if countBytes > 0 and elapsedMillis > 0: " MiB/s " & $(countBytes div elapsedMillis.uint64) else: ""
    echo "CPU Time [", benchmarkName, "] ", elapsedStr, "s", mibs

template encodedBytes*(): untyped =
  when Exp == 8:
    @[10, 54, 200, 1, 174, 73, 223, 252,
    169, 147].mapIt( it.GFSymbol )
  elif Exp == 16:
    @[36016, 51702, 9536, 44839, 35358, 8384,
    60905, 26494, 59225, 2106].mapIt( it.GFSymbol )
  elif Exp == 18:
    @[168077, 13380, 24621, 141427, 54004,
      88735, 191583, 165529, 230647, 166598]
      .mapIt( it.GFSymbol )
  elif Exp == 20:
    @[608501, 613582, 894931, 464752, 699334,
    57844, 1014134, 84555, 292666, 699046]
    .mapIt( it.GFSymbol )

template genPoly*(): untyped =
  when Exp == 8:
    @[1, 152, 185, 240, 5, 111, 99, 6, 220,
    112, 150, 69, 36, 187, 22, 228, 198, 121,
    121, 165, 174]
    .mapIt( it.GFSymbol )
  elif Exp == 16:
    @[1, 65176, 27223, 6912, 26751, 43073, 33579,
    13334, 45534, 31631, 20966, 31595, 36271, 57445,
    7168, 19847, 45485, 28236, 30650, 56422, 20577]
    .mapIt( it.GFSymbol )
  elif Exp == 18:
    @[1, 261928, 65365, 63074, 93193, 126482, 250408,
    176177, 5059, 249940, 230626, 5039, 261785, 124137,
    218065, 74688, 2109, 96242, 130031, 123299, 216356]
    .mapIt( it.GFSymbol )
  elif Exp == 20:
    @[1, 1048575, 419230, 982646, 251755, 856531, 780318,
    946087, 208868, 852853, 271828, 709278, 157523, 437933,
    393829, 510568, 578081, 185173, 1046911, 334336, 597300]
    .mapIt( it.GFSymbol )
