import std/sequtils
import std/random
import std/unittest
import std/times
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
    var
      parity = newSeq[GFSymbol](10)

    msg.encode(parity = parity, nsym = 10)
    check parity == encodedBytes()

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

    check:
      message.decode(10, erasePos = toSeq(erased), erasures = true) == erased.len
      message == (msg & par)
      string.fromBytes(message[0..<msg.len].mapIt( it.byte )) == "Hello Reed-Solomon!"

  test "Decode with errors":
    var
      par = encodedBytes()
      message = msg & par
      count = 0

    # Add random floor(t/2) errors
    var errors: seq[int]
    while count < 5:
      let pos = rand(0..<message.len)
      if errors.find(pos) > -1:
        continue

      errors.add(pos)
      message[pos] = rand(0..Order.int).GFSymbol
      count += 1

    check:
      message.decode(10) == errors.len
      message == (msg & par)
      string.fromBytes(message[0..<msg.len].mapIt( it.byte )) == "Hello Reed-Solomon!"

  test "Decode with erasures and errors":
    var
      par = encodedBytes()
      message = msg & par
      count = 0
      erased: seq[int]
      errors: seq[int]

    # Add both, random floor(t/2) erasures and errors
    while count < 6:
      while true:
        let pos = rand(0..<message.len)
        if count mod 2 == 0:
          if errors.find(pos) > -1:
            continue

          errors.add(pos)
          message[pos] = rand(0..Order.int).GFSymbol
        else:
          if erased.find(pos) > -1:
            continue

          erased.add(pos)
          message[pos] = 0.GFSymbol

        break

      count += 1

    check:
      message.decode(10, erasePos = toSeq(erased)) == (erased.len + errors.len)
      message == (msg & par)
      string.fromBytes(message[0..<msg.len].mapIt( it.byte )) == "Hello Reed-Solomon!"

suite "Test shards":
  setup:
    randomize()

  test "Test shard coding":
    const
      k = 20
      n = 40
      dataSize = 256
      genPoly = generator(n - k)

    var
      data: seq[seq[GFSymbol]]

    benchmark("Filling data"):
      for i in 0..<k:
        var d = randomAlpha(dataSize).mapIt( it.GFSymbol )
        countBytes += d.len.uint64
        data.add(d)

    var shards = data
    benchmark("Encode Loop"):
      for i in 0..<dataSize:
        var msg {.noinit.} = newSeq[GFSymbol](k)
        for s in 0..<data.len:
          msg[s] = data[s][i]
        countBytes += msg.len.uint64

        var parity {.noInit.} = newSeq[GFSymbol](n-k)
        msg.encode(nsym = n - k, parity = parity, gen = genPoly)
        # msg.encode(nsym = n - k, parity = parity)
        check parity.len == n - k

        for m in 0..<parity.len:
          if (k + m) > shards.high:
            var d {.noInit.} = newSeq[GFSymbol](dataSize)
            shards.add(d)

          shards[k + m][i] = parity[m]

    var
      erased: seq[int]
      count = 0

    benchmark("Erasures Loop"):
      while count < k:
        while true:
          let row = rand(0..<n)
          if erased.find(row) > -1:
            continue

          erased.add(row)
          shards[row] = @[]
          break

        count += 1

    benchmark("Decode Loop"):
      var msg {.noinit.} = newSeq[GFSymbol](n)
      for i in 0..<dataSize:
        for s in 0..<shards.len:
          if shards[s].len > 0 and not(shards[s].len < dataSize):
            msg[s] = shards[s][i]
          else:
            msg[s] = 0.GFSymbol

        discard msg.decode(k, erasePos = erased, erasures = true)
        countBytes += msg.len.uint64

        for o in erased:
          shards[o].add(msg[o])

    check data == shards[0..<k]
