{.push raises: [Defect].}

{.deadCodeElim: on.}

import std/[algorithm, sequtils]
import ./poly

export poly

func generator*(
  nsym: int,
  alpha = 2.GFSymbol): seq[GFSymbol] =
  ## Generate an irreducible generator polynomial
  ## using consecutive alphas (a_n) that are roots of the
  ## polynomial and consequently its factors (x-a)
  ##

  var g = @[1.GFSymbol]
  for i in 0..<nsym:
    g = g * @[1.GFSymbol, alpha ^ i]

  return g

func encode*(
  msg: openArray[GFSymbol],
  parity: var seq[GFSymbol],
  nsym: int,
  alpha = 2.GFSymbol,
  gen: seq[GFSymbol] = @[]) {.raises: [Defect, RSError].} =
  ## Reed-Solomon encoding
  ##

  if (msg.len + nsym).uint > Order:
    raise newException(
      RSError,
      "Message is too long ( " & $(msg.len + nsym) & " when max is " & $Order & " )")

  if parity.len < nsym:
    raise newException(RSError, "Not enough parity bytes!")

  let
    gen = if gen.len <= 0:
      generator(nsym, alpha)
    else:
      gen

  var res {.noInit.} = newSeq[GFSymbol](msg.len)
  divide(msg, gen, res, parity)

func syndromes*(
  msg: openArray[GFSymbol],
  nsym: int,
  alpha = 2.GFSymbol): seq[GFSymbol] {.inline.} =
  ## Given the received codeword msg and the number of
  ## error correcting symbols (nsym), computes the syndromes
  ## polynomial.
  ##
  ## Mathematically, it's essentially equivalent to a Fourrier
  ## Transform (Chien search being the inverse).
  ##

  var synd {.noinit.} = newSeq[GFSymbol](nsym + 1)
  for i in 0..<nsym:
    synd[i + 1] = msg.eval(alpha ^ i)
    # echo synd

  return synd

func errorLocator*(
  synd: openArray[GFSymbol],
  nsym: int,
  eraseLoc: openArray[GFSymbol] = @[],
  eraseCount = 0): seq[GFSymbol] {.raises: [Defect, RSError].} =
  ## Find error/errata locator and evaluator polynomials
  ## with Berlekamp-Massey algorithm
  ##

  # Init the polynomials

  # if the erasure locator polynomial is supplied,
  # we init with its value, so that we include erasures
  # in the final locator polynomial
  var (errLoc, oldLoc) = if eraseLoc.len > 0:
    (@eraseLoc, @eraseLoc)
  else:
    # This is the main variable we want to fill,
    # also called Sigma in other notations or more
    # formally the errors/errata locator polynomial.
    (@[1.GFSymbol], @[1.GFSymbol])

  let
    syndShift = if synd.len > nsym: synd.len - nsym else: 0

  for i in 0..<(nsym - eraseCount):
    # if an erasures locator polynomial was provided to init the
    # errors locator polynomial, then we must skip the FIRST eraseCount
    # iterations (not the last iterations, this is very important!)
    let K = if eraseLoc.len > 0:
      (eraseCount + i + syndShift)
    else:
      # if erasures locator is not provided, then either there's no
      # erasures to account or we use the Forney syndromes, so we don't
      # need to use eraseCount nor eraseLoc (the erasures have been
      # trimmed out of the Forney syndromes).
      i + syndShift

    # Compute the discrepancy Delta
    var
      delta = synd[K]

    for j in 1..<errLoc.len:
      # delta is also called discrepancy. Here we do a partial polynomial
      # multiplication - ie, we compute the polynomial multiplication only
      # for the term of degree K.

      delta = delta + (errLoc[^(j+1)] * synd[K - j])
    # echo " K ", K, " delta ", delta, " err loc ", errLoc.reversed * synd # debugline

    # Shift polynomials to compute the next degree
    oldLoc = oldLoc & @[0.GFSymbol]

    # Iteratively estimate the errata locator and evaluator polynomials
    if delta != 0: # Update only if there's a discrepancy
      if oldLoc.len > errLoc.len:
        # Computing errata locator polynomial Sigma
        let newLoc = oldLoc.scale(delta)
        oldLoc = errLoc.scale(delta.inverse) # effectively we are doing errLoc * 1/delta = errLoc // delta
        errLoc = newLoc

      # Update with the discrepancy
      errLoc = errLoc + oldLoc.scale(delta)

  # Check if the result is correct, that there's not too many errors to correct
  var drop = 0
  while errLoc.len > 0 and errLoc[0] == 0: drop.inc()
  errLoc = if drop > 0: errLoc[0..<drop] else: errLoc

  let errs = errLoc.len - 1
  if ((errs - eraseCount) * 2) + eraseCount > nsym:
    raise newException(RSError, "Too many errors to correct") # too many errors to correct

  return errLoc

func errataLocator*(
  pos: openArray[int],
  alpha = 2.GFSymbol): seq[GFSymbol] =
  ## Compute the erasures/errors/errata locator polynomial from the
  ## erasures/errors/errata positions - the positions must be relative
  ## to the x coefficient, eg: "hello worldxxxxxxxxx" is tampered to
  ## "h_ll_ worldxxxxxxxxx" with xxxxxxxxx being the ecc of length
  ## n-k=9, here the string positions are [1, 4], but the coefficients
  ## are reversed since the ecc characters are placed as the first
  ## coefficients of the polynomial, thus the coefficients of the erased
  ## characters are n-1 - [1, 4] = [18, 15] = erasures_loc to be specified
  ## as an argument.
  ##

  # just to init because we will multiply, so it must be 1 so that
  # the multiplication starts correctly without nulling any term
  var loc = @[1.GFSymbol]
  # erasures_loc = product(1 - x*alpha**i) for i in erasures_pos and where
  # alpha is the alpha chosen to evaluate polynomials.
  for i in pos:
    loc = loc * (@[1.GFSymbol] + @[alpha ^ i, 0.GFSymbol]) # TODO: not sure why we need to add the 0
                                                           # coefficient and it works without it

  return loc

func errorEvaluator*(
  synd: openArray[GFSymbol],
  errLoc: openArray[GFSymbol],
  nsym: int): seq[GFSymbol] =
  ## Compute the error or erasures if you supply
  ## sigma = erasures locator polynomial, or errata
  ## evaluator polynomial Omega from the syndrome
  ## and the error/erasures/errata locator Sigma.
  ##

  # Omega(x) = [ Synd(x) * Error_loc(x) ] mod x^(n-k+1)

  # first multiply syndromes * errata_locator, then do a
  # polynomial division to truncate the polynomial to the
  # required length
  let remainder = @synd * @errLoc
  return remainder[(remainder.len - (nsym + 1))..remainder.high]

func correctErrata*(
  msg, synd: var openArray[GFSymbol],
  errPos: openArray[int], # errPos is a list of the positions of the errors/erasures/errata
  alpha = 2.GFSymbol): seq[GFSymbol] {.raises: [Defect, RSError].} =
  ## Forney algorithm, computes the values (error magnitude)
  ## to correct the input message.
  ##

  # Calculate errata locator polynomial to correct both errors and erasures.
  # By combining the errors positions given by the error locator polynomial
  # found by BM with the erasures positions given by caller.
  let
    # Need to convert the positions to coefficients degrees for the errata
    # locator algo to work
    #
    # eg: instead of [0, 1, 2] it will become [len(msg)-1, len(msg)-2, len(msg)-3]
    coefPos = errPos.mapIt((msg.len - 1) - it)
    errLoc = errataLocator(coefPos, alpha)
    # calculate errata evaluator polynomial often called Omega or Gamma in academic papers

  var
    evaluator = errorEvaluator(synd.reversed(), errLoc, (errLoc.len - 1))
    errEval = evaluator.reversed()

  # Second part of Chien search to get the error location polynomial X from
  # the error positions in errPos - the roots of the error locator polynomial,
  # ie, where it evaluates to 0
  var X = newSeqOfCap[GFSymbol](coefPos.len) # will store the position of the errors
  for i in 0..<coefPos.len:
    let
      e = Degree - coefPos[i].uint
      p = alpha ^ (GFLog.len - e - 1).int

    X.add(p)

  # Forney algorithm: compute the magnitudes
  var
    # will store the values that need to be corrected
    # (substracted) to the message containing errors.
    # This is sometimes called the error magnitude
    # polynomial.
    E = newSeq[GFSymbol](msg.len)
    Xlength = X.len

  for i in 0..<Xlength:
    let
      Xi_inv = X[i].inverse

    # Compute the formal derivative of the error locator polynomial
    # (see Blahut, Algebraic codes for data transmission, pp 196-197).
    # the formal derivative of the errata locator is used as the
    # denominator of the Forney Algorithm, which simply says that the
    # ith error value is given by
    # error_evaluator(gf_inverse(X[i])) / error_locator_derivative(gf_inverse(X[i])).
    # See Blahut, Algebraic codes for data transmission, pp 196-197.
    var errLocPrimeTmp = newSeqOfCap[GFSymbol](Xlength)
    for j in 0..<Xlength:
      if j != i:
          errLocPrimeTmp.add( 1.GFSymbol - (Xi_inv * X[j]) )

    # compute the product, which is the denominator of the Forney algorithm (errata locator derivative)
    var errLocPrime = 1.GFSymbol
    for coef in errLocPrimeTmp:
      errLocPrime = errLocPrime * coef

    # Compute y (evaluation of the errata evaluator polynomial)
    # This is a more faithful translation of the theoretical
    # equation contrary to the old forney method. Here it is an
    # exact reproduction:
    # Yl = omega(Xl.inverse()) / prod(1 - Xj*Xl.inverse()) for j in len(X)
    var
      y = (X[i] ^ 1) * errEval.reversed().eval(Xi_inv) # adjust to fcr parameter

    # Check: errLocPrime (the divisor) should not be zero.
    if errLocPrime == 0:
      raise newException(RSError, "Could not find error magnitude")

    # Compute the magnitude
    # let
      # magnitude value of the error, calculated by the Forney algorithm
      # (an equation in fact): dividing the errata evaluator with the errata
      # locator derivative gives us the errata magnitude (ie, value to repair)
      # the ith symbol
      # magnitude = y / errLocPrime

    E[errPos[i]] = y / errLocPrime # store the magnitude for this error into the magnitude polynomial

  # Apply the correction of values to get our message corrected!
  # NOTE: that the ecc bytes also gets corrected!
  # this isn't the Forney algorithm, we just apply the result of decoding here

  # equivalent to Ci = Ri - Ei where Ci is the correct message,
  # Ri the received (senseword) message, and Ei the errata magnitudes
  # (minus is replaced by XOR since it's equivalent in GF(2^p)).
  # So in fact here we substract from the received message the errors
  # magnitude, which logically corrects the value to what it should be.
  return (msg - E)

func findErrors(
  errLoc: seq[GFSymbol],
  msgLen: int,
  alpha = 2.GFSymbol): seq[int] {.raises: [Defect, RSError].} =
  ## Find the roots (ie, where evaluation = zero) of error polynomial
  ## by brute-force trial, this is a sort of Chien's search
  ## (but less efficient, Chien's search is a way to evaluate the
  ## polynomial such that each evaluation only takes constant time).
  ##

  var
    errs = errLoc.len - 1
    errPos: seq[int]

  # normally we should try all 2^8 possible values, but here we
  # optimize to just check the interesting symbols
  for i in 0..<msgLen:
    if errLoc.eval(alpha ^ i) == 0: # It's a 0? Bingo, it's a root of the error locator polynomial,
                                    # in other words this is the location of an error

      errPos.add(msgLen - 1 - i)    # since the `errLoc` is the inverted coeficients of the original `errLoc`
                                    # polynomial, the inverse index is the actual position of the error

  # Sanity check: the number of errors/errata positions found should be
  # exactly the same as the length of the errata locator polynomial
  if errPos.len != errs:
    # couldn't find error locations
    raise newException(
      RSError,
      "Too many (or few) errors found by Chien Search for the errata locator polynomial!")

  return errPos

func forneySyndromes*(
  synd: openArray[GFSymbol],
  pos: openArray[int],
  nmess: int,
  alpha = 2.GFSymbol): seq[GFSymbol] =
  ## Compute Forney syndromes, which computes a modified
  ## syndromes to compute only errors (erasures are trimmed out).
  ## Do not confuse this with Forney algorithm, which allows
  ## to correct the message based on the location of errors.
  ##

  let
    # prepare the coefficient degree positions
    # (instead of the erasures positions)
    reversedPos = pos.mapIt(nmess - 1 - it)

  var
    # make a copy and trim the first coefficient which
    # is always 0 by definition
    fsynd = synd[1..synd.high]

  for i in 0..<reversedPos.len:
    let x = alpha ^ reversedPos[i]
    for j in 0..<fsynd.len - 1:
      fsynd[j] = (fsynd[j] * x) + fsynd[j + 1]

  # Theoretical way of computing the modified Forney syndromes: fsynd = (eraseLoc * synd) % x^(n-k)
  # See Shao, H. M., Truong, T. K., Deutsch, L. J., & Reed, I. S. (1986, April). A single chip VLSI Reed-Solomon decoder.
  # In Acoustics, Speech, and Signal Processing, IEEE International Conference on ICASSP'86. (Vol. 11, pp. 2151-2154). IEEE.ISO 690

  # eraseLoc = errataLocator(reversedPos, generator=generator) # computing the erasures locator polynomial
  # fsynd = gf_poly_mul(eraseLoc[::-1], synd[1:]) # then multiply with the syndrome to get the untrimmed forney syndrome
  # fsynd = fsynd[len(pos):] # then trim the first erasePos coefficients which are useless. Seems to be not necessary,
  # but this reduces the computation time later in BM (thus it's an optimization).

  return fsynd

func decode*(
  msgOut: var seq[GFSymbol],
  nsym: int,
  alpha = 2.GFSymbol,
  erasePos: openArray[int] = [],
  erasures = false): int
  {.raises: [Defect, RSError].} =
  ## Reed-Solomon main decoding function
  ##

  if msgOut.len > Order.int:
    raise newException(
      RSError,
      "Message is too long (" & $msgOut.len & " when max is " & $Order & ")")

  # check if there are too many erasures to correct (beyond the Singleton bound)
  if erasePos.len > nsym:
    raise newException(RSError, "Too many erasures to correct")

  # prepare the syndromes using only errors -
  # ie: errors = characters that were either replaced by null byte
  # or changed to another character, but we don't know their positions
  var synd = msgOut.syndromes(nsym, alpha)

  # check if there's any error/erasure in the input codeword.
  # If not (all syndromes coefficients are 0), then just return
  # the message as-is.
  if max(synd) == 0:
    return 0

  # Find errors locations
  var errPos: seq[int]
  if not erasures:
    let
      # compute the Forney syndromes, which hide the erasures from the original
      # syndrome (so that BM will just have to deal with errors, not erasures)
      fsynd = forneySyndromes(synd, erasePos, msgOut.len, alpha)
      # compute the error locator polynomial using Berlekamp-Massey
      errLoc = errorLocator(fsynd, nsym, eraseCount = erasePos.len)
      # locate the message errors using Chien search (or bruteforce search)

    errPos = findErrors(errLoc.reversed(), msgOut.len, alpha)
    if errPos.len <= 0:
      raise newException(RSError, "Could not locate error")

  # Find errors values and apply them to correct the message
  # compute errata evaluator and errata magnitude polynomials,
  # then correct errors and erasures

  # NOTE: we here use the original syndrome, not the forney syndrome
  # (because we will correct both errors and erasures, so we need the full syndrome)
  msgOut = msgOut.correctErrata(
    synd,
    @erasePos & errPos,
    alpha)

  # check if the final message is fully repaired
  if max(msgOut.syndromes(nsym, alpha)) > 0:
    raise newException(RSError, "Could not correct message")

  # return the successfully decoded message
  return (@erasePos & errPos).len # also return the corrected ecc block so that the user can check()

func check(
  msg: openArray[GFSymbol],
  nsym: int,
  fcr = 0,
  alpha = 2.GFSymbol): bool =
  ## Returns true if the message + ecc has no error of false
  ## otherwise (may not always catch a wrong decoding or a wrong
  ## message, particularly if there are too many errors
  ## -- above the Singleton bound --, but it usually does)
  return ( max(syndromes(msg, nsym, alpha)) == 0 )
