import std/sequtils
import std/random
import std/unittest
import std/times
import std/os
import std/strutils

import ./helpers
import rs

let msg = "Hello Reed-Solomon!"
  .toBytes()
  .mapIt( it.GFSymbol )

suite "Test Reed-Solomon Coding in GF " & $Exp:
  setup:
    randomize()

  test "Encode":
    check msg.encode(10) == encodedBytes()

  test "Decode with erasures":
    var
      par = encodedBytes()
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
      par = encodedBytes()
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
      par = encodedBytes()
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

suite "Test shards":
  setup:
    randomize()

  test "Test shard coding":
    const
      k = 20
      n = 40
      dataSize = 256 * 1024

    var
      data: seq[seq[GFSymbol]]

    benchmark("Filling data"):
      for i in 0..<k:
        data.add(randomAlpha(dataSize).mapIt( it.GFSymbol ))

    var shards = data
    benchmark("Encode Loop"):
      for i in 0..<dataSize:
        var msg {.noinit.} = newSeq[GFSymbol](k)
        for s in 0..<data.len:
          msg[s] = data[s][i]
        countBytes += msg.len.uint64

        let parity = msg.encode(n - k, gen = genPoly())
        check parity.len == n - k

        for m in 0..<parity.len:
          if (k + m) > shards.high:
            shards.add(newSeq[GFSymbol](dataSize))

          shards[k + m][i] = parity[m]

    var
      erased: seq[int]
      count = 0

    while count < k:
      while true:
        let row = rand(0..<k)
        if erased.find(row) > -1:
          continue

        erased.add(row)
        shards[row] = @[]
        break

      count += 1

    var
      rebuilt: seq[seq[GFSymbol]]

    benchmark("Decode Loop"):
      for i in 0..<dataSize:
        var msg {.noinit.} = newSeq[GFSymbol](n)
        for s in 0..<shards.len:
          if shards[s].len > 0:
            msg[s] = shards[s][i]
          else:
            msg[s] = 0.GFSymbol

        countBytes += msg.len.uint64
        let (orig, _) = msg.decode(k, erasePos = erased, erasures = true)

        for o in 0..<orig.len:
          if o > rebuilt.high:
            rebuilt.add(newSeq[GFSymbol](dataSize))

          rebuilt[o][i] = orig[o]

      check data == rebuilt
