import std/sequtils
import std/random
import std/unittest

import ./helpers
import rs

template encodedBytes*(): untyped =
  when Exp == 8:
    (@[72, 101, 108, 108, 111, 32, 82, 101,
      101, 100, 45, 83, 111, 108, 111, 109,
      111, 110, 33].mapIt( it.GFSymbol ),
      @[10, 54, 200, 1, 174, 73, 223, 252,
      169, 147].mapIt( it.GFSymbol ))
  elif Exp == 16:
    (@[72, 101, 108, 108, 111, 32, 82, 101,
      101, 100, 45, 83, 111, 108, 111, 109, 111,
      110, 33]
      .mapIt( it.GFSymbol ),
      @[36016, 51702, 9536, 44839, 35358, 8384,
      60905, 26494, 59225, 2106]
      .mapIt( it.GFSymbol ))
  elif Exp == 18:
    (@[72, 101, 108, 108, 111, 32, 82, 101,
      101, 100, 45, 83, 111, 108, 111, 109, 111,
      110, 33]
      .mapIt( it.GFSymbol ),
      @[168077, 13380, 24621, 141427, 54004,
        88735, 191583, 165529, 230647, 166598]
        .mapIt( it.GFSymbol ))
  elif Exp == 20:
    (@[72, 101, 108, 108, 111, 32, 82, 101,
      101, 100, 45, 83, 111, 108, 111, 109, 111,
      110, 33]
      .mapIt( it.GFSymbol ),
      @[608501, 613582, 894931, 464752, 699334,
      57844, 1014134, 84555, 292666, 699046]
      .mapIt( it.GFSymbol ))

suite "Test Reed-Solomon Coding in GF " & $Exp:
  setup:
    randomize()

  test "Encode":
    let msg = "Hello Reed-Solomon!"
      .toBytes()
      .mapIt( it.GFSymbol )

    check encode(msg, 10) == encodedBytes()

  test "Decode with erasures":
    var
      (msg, par) = encodedBytes()
      message = msg & par
      count = 0
      erased: seq[int]

    # Add random `t` erasures
    while count < 10:
      while true:
        let pos = rand(0..<message.len)
        if erased.find(pos) > -1:
          continue

        erased.add(pos)
        message[pos] = 0.GFSymbol
        break

      count += 1

    let
      (corrected, _) = decode(message, 10, erasePos = toSeq(erased), erasures = true)
    check string.fromBytes(corrected.mapIt( it.byte )) == "Hello Reed-Solomon!"

  test "Decode with errors":
    var
      (msg, par) = encodedBytes()
      message = msg & par
      count = 0

    # Add random floor(t/2) errors
    while count < 5:
      let pos = rand(0..<message.len)
      message[pos] = rand(0..Order.int).GFSymbol
      count += 1

    let
      (corrected, _) = decode(message, 10)

    check string.fromBytes(corrected.mapIt( it.byte )) == "Hello Reed-Solomon!"

  test "Decode with erasures and errors":
    var
      (msg, par) = encodedBytes()
      message = msg & par
      count = 0
      erased: seq[int]

    # Add both, random floor(t/2) erasures and errors
    while count < 5:
      if count mod 2 == 0:
        let pos = rand(0..<message.len)
        message[pos] = rand(0..Order.int).GFSymbol
      else:
        while true:
          let pos = rand(0..<message.len)
          if erased.find(pos) > -1:
            continue

          erased.add(pos)
          message[pos] = 0.GFSymbol
          break

      count += 1

    let
      (corrected, _) = decode(message, 10, erasePos = toSeq(erased))
    check string.fromBytes(corrected.mapIt( it.byte )) == "Hello Reed-Solomon!"
