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
