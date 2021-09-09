{.push raises: [Defect].}

{.experimental: "codeReordering".}

import std/[strutils, algorithm, sequtils]
import ./polynomial

export polynomial

proc generator*(
  nsym: int,
  fcr = 0,
  generator = Char.GFUint): seq[GFUint] =
  ## Generate an irreducible generator polynomial
  ## (necessary to encode a message into Reed-Solomon)
  ##

  var g = @[1.GFUint]
  for i in 0..<nsym:
    g = g * @[1.GFUint, generator ^ (i + fcr)]

  return g

proc generatorAll*(
  maxSym: int,
  fcr = 0,
  generator = Char.GFUint): seq[seq[GFUint]] =
  ## Generate all irreducible generator polynomials up to maxSym
  ## (usually you can use n, the length of the message+ecc). Very
  ## useful to reduce processing time if you want to encode using
  ## variable schemes and nsym rates.
  ##

  var gen = newSeq[seq[GFUint]](maxSym)
  # TODO: first two are always 1 - maybe in GF(2^8)?
  gen[0] = @[1.GFUint]
  gen[1] = @[1.GFUint]
  for nsym in 0..<maxSym:
    gen[nsym] = generator(nsym, fcr, generator)

  return gen

proc encode*(
  msg: seq[GFUint],
  nsym: int,
  fcr = 0,
  generator = Char.GFUint,
  gen: seq[GFUint] = @[]): seq[GFUint] {.raises: [ValueError].} =
  ## Reed-Solomon encoding
  ##

  if (msg.len + nsym).uint > Order:
    raise newException(
      ValueError,
      "Message is too long ($1 when max is $2)" % [$(msg.len + nsym), $Order])

  var gen = if gen.len <= 0:
      generator(nsym, fcr, generator)
    else:
      gen

  var (res, remainder) = msg / gen

  # TODO: Should not append to the message, return tuple
  # to allow implementations choose how codes are stored
  return msg & remainder # remainder is the RS code

################### REED-SOLOMON DECODING ###################

proc clacSyndromes(
  msg: seq[GFUint],
  nsym: int,
  fcr = 0,
  generator = Char.GFUint): seq[GFUint] =
  ## Given the received codeword msg and the number of
  ## error correcting symbols (nsym), computes the syndromes
  ## polynomial.
  ##
  ## Mathematically, it's essentially equivalent to a Fourrier Transform
  ## (Chien search being the inverse).
  ##

  # Note the "[0] +" : we add a 0 coefficient for the lowest degree (the constant).
  # This effectively shifts the syndrome, and will shift every computations depending
  # on the syndromes - such as the errors locator polynomial, errors evaluator polynomial,
  # etc. but not the errors positions.
  #
  # This is not necessary, you can adapt subsequent computations to start from 0 instead
  # of skipping the first iteration (ie, the often seen range(1, n-k+1)),
  var synd = newSeq[GFUint](nsym)
  for i in 0..<nsym:
    synd[i] = msg.eval(generator ^ (i + fcr))

  return @[0.GFUint] & synd # pad with one 0 for mathematical precision (else we can end up with weird calculations sometimes)

proc findErrorLocator(
  synd: seq[GFUint],
  nsym: int,
  eraseLoc: seq[GFUint] = @[],
  eraseCount = 0): seq[GFUint] {.raises: [CatchableError].} =
  ## Find error/errata locator and evaluator polynomials
  ## with Berlekamp-Massey algorithm
  ##

  # The idea is that BM will iteratively estimate the error locator polynomial.
  # To do this, it will compute a Discrepancy term called Delta, which will tell
  # us if the error locator polynomial needs an update or not - hence why it's called
  # discrepancy: it tells us when we are getting off board from the correct value

  # Init the polynomials

  # if the erasure locator polynomial is supplied,
  # we init with its value, so that we include erasures
  # in the final locator polynomial
  var (errLoc, oldLoc) = if eraseLoc.len > 0:
    (eraseLoc, eraseLoc)
  else:
    # This is the main variable we want to fill,
    # also called Sigma in other notations or more
    # formally the errors/errata locator polynomial.
    (@[1.GFUint], @[1.GFUint])

  # BM is an iterative algorithm, and we need the errata
  # locator polynomial of the previous iteration in order
  # to update other necessary variables.

  # L = 0 # update flag variable, not needed here because we use
  # an alternative equivalent way of checking if update is needed -
  # but using the flag could potentially be faster depending on
  # if using length(list) is taking linear time in your language,
  # here in Python it's constant so it's as fast.

  # Fix the syndrome shifting: when computing the syndrome, some
  # implementations may prepend a 0 coefficient for the lowest
  # degree term (the constant). This is a case of syndrome shifting,
  # thus the syndrome will be bigger than the number of ecc symbols
  # (I don't know what purpose serves this shifting). If that's the
  # case, then we need to account for the syndrome shifting when we
  # use the syndrome such as inside BM, by skipping those prepended
  # coefficients.
  #
  # Another way to detect the shifting is to detect the 0 coefficients:
  # by definition, a syndrome does not contain any 0 coefficient (except
  # if there are no errors/erasures, in this case they are all 0). This
  # however doesn't work with the modified Forney syndrome, which set to
  # 0 the coefficients corresponding to erasures, leaving only the coefficients
  # corresponding to errors.
  var
    syndShift = if len(synd) > nsym: synd.len - nsym else: 0

  # generally: nsym-eraseCount == len(synd), except when you input
  # a partial eraseLoc and using the full syndrome instead of the
  # Forney syndrome, in which case nsym-eraseCount is more correct
  # (len(synd) will fail badly with IndexError).
  for i in 0..<(nsym - eraseCount):
    # if an erasures locator polynomial was provided to init the
    # errors locator polynomial, then we must skip the FIRST eraseCount
    # iterations (not the last iterations, this is very important!)
    var K = if eraseLoc.len > 0:
      (eraseCount + i + syndShift)
    else:
      # if erasures locator is not provided, then either there's no
      # erasures to account or we use the Forney syndromes, so we don't
      # need to use eraseCount nor eraseLoc (the erasures have been
      # trimmed out of the Forney syndromes).
      i + syndShift

    # Compute the discrepancy Delta

    # Here is the close-to-the-books operation to compute the
    # discrepancy Delta: it's a simple polynomial multiplication
    # of error locator with the syndromes, and then we get the Kth
    # element.
    #
    # delta = gf_poly_mul(errLoc[::-1], synd)[K] # theoretically it
    # should be gf_poly_add(synd[::-1], [1])[::-1] instead of just
    # synd, but it seems it's not absolutely necessary to correctly decode.
    #
    # But this can be optimized: since we only need the Kth element,
    # we don't need to compute the polynomial multiplication for any
    # other element but the Kth. Thus to optimize, we compute the polymul
    # only at the item we need, skipping the rest - avoiding a nested loop,
    # thus we are linear time instead of quadratic.
    #
    # This optimization is actually described in several figures of the book
    # "Algebraic codes for data transmission", Blahut, Richard E., 2003,
    # Cambridge university press.

    var
      delta = synd[K]

    for j in 1..<errLoc.len:
      # delta is also called discrepancy. Here we do a partial polynomial
      # multiplication (ie, we compute the polynomial multiplication only
      # for the term of degree K). Should be equivalent to brownanrs.polynomial.mul_at().
      delta = (delta xor errLoc[^(j+1)] * synd[K - j])
    #print "delta", K, delta, list(gf_poly_mul(errLoc[::-1], synd)) # debugline

    # Shift polynomials to compute the next degree
    oldLoc = oldLoc & @[0.GFUint]

    # Iteratively estimate the errata locator and evaluator polynomials
    if delta != 0: # Update only if there's a discrepancy
        if oldLoc.len > errLoc.len: # Rule B (rule A is implicitly defined because rule A just says that we skip any modification for this iteration)
        #if 2*L <= K+eraseCount: # equivalent to len(oldLoc) > len(errLoc), as long as L is correctly computed
            # Computing errata locator polynomial Sigma
            let newLoc = oldLoc.scale(delta)
            oldLoc = errLoc.scale(delta.inverse) # effectively we are doing errLoc * 1/delta = errLoc // delta
            errLoc = newLoc
            # Update the update flag
            #L = K - L # the update flag L is tricky: in Blahut's schema, it's mandatory to use `L = K - L - eraseCount` (and indeed in a previous draft of this function, if you forgot to do `- eraseCount` it would lead to correcting only 2*(errors+erasures) <= (n-k) instead of 2*errors+erasures <= (n-k)), but in this latest draft, this will lead to a wrong decoding in some cases where it should correctly decode! Thus you should try with and without `- eraseCount` to update L on your own implementation and see which one works OK without producing wrong decoding failures.

        # Update with the discrepancy
        errLoc = errLoc + oldLoc.scale(delta)

  # Check if the result is correct, that there's not too many errors to correct
  var drop = 0
  while len(errLoc) > 0 and errLoc[0] == 0: drop.inc()
  errLoc = if drop > 0: errLoc[0..<drop] else: errLoc

  var errs = len(errLoc) - 1
  if ((errs - eraseCount) * 2) + eraseCount > nsym:
    raise newException(CatchableError, "Too many errors to correct")    # too many errors to correct

  return errLoc

proc findErrataLocator*(
  ePos: seq[int],
  generator = Char.GFUint): seq[GFUint] =
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
  var eLoc = @[1.GFUint]
  # erasures_loc = product(1 - x*alpha**i) for i in erasures_pos and where
  # alpha is the alpha chosen to evaluate polynomials.
  for i in ePos:
    eLoc = eLoc * @[1.GFUint] + @[generator ^ i, 0.GFUint]

  return eLoc

proc findErrorEvaluator*(
  synd: seq[GFUint],
  errLoc: seq[GFUint],
  nsym: int): seq[GFUint] =
  ## Compute the error or erasures if you supply
  ## sigma=erasures locator polynomial, or errata
  ## evaluator polynomial Omega from the syndrome
  ## and the error/erasures/errata locator Sigma.
  ##

  # Omega(x) = [ Synd(x) * Error_loc(x) ] mod x^(n-k+1)

  # first multiply syndromes * errata_locator, then do a
  # polynomial division to truncate the polynomial to the
  # required length
  let remainder = synd * errLoc
  return remainder[(remainder.len - (nsym + 1))..remainder.high]

proc correctErrata*(
  msg: seq[GFUint],
  synd: seq[GFUint],
  errPos: seq[int], # errPos is a list of the positions of the errors/erasures/errata
  fcr = 0,
  generator = Char.GFUint): seq[GFUint] {.raises: [CatchableError].} =
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
    coefPos = errPos.mapIt(msg.len - 1 - it)
    errLoc = findErrataLocator(coefPos, generator)
    # calculate errata evaluator polynomial
    # often called Omega or Gamma in academic papers
    errEval = findErrorEvaluator(
      synd.reversed(),
      errLoc,
      (errLoc.len - 1)).reversed()

  # Second part of Chien search to get the error location polynomial X from
  # the error positions in errPos - the roots of the error locator polynomial,
  # ie, where it evaluates to 0
  var X: seq[GFUint] # will store the position of the errors
  for i in 0..<coefPos.len:
    let
      e = Order - coefPos[i].uint
      p = generator ^ (GFLog.len - e - 1).int

    X.add(p)

  # Forney algorithm: compute the magnitudes
  var
    # will store the values that need to be corrected
    # (substracted) to the message containing errors.
    # This is sometimes called the error magnitude
    # polynomial.
    E = newSeq[GFUint](msg.len)
    Xlength = X.len

  for i, Xi in X:
    let
      Xi_inv = Xi.inverse

    # Compute the formal derivative of the error locator polynomial
    # (see Blahut, Algebraic codes for data transmission, pp 196-197).
    # the formal derivative of the errata locator is used as the
    # denominator of the Forney Algorithm, which simply says that the
    # ith error value is given by
    # error_evaluator(gf_inverse(Xi)) / error_locator_derivative(gf_inverse(Xi)).
    # See Blahut, Algebraic codes for data transmission, pp 196-197.
    var errLocPrimeTmp: seq[GFUint]
    for j in 0..<Xlength:
      if j != i:
          errLocPrimeTmp.add( 1.GFUint - (Xi_inv * X[j]) )

    # compute the product, which is the denominator of the Forney algorithm (errata locator derivative)
    var errLocPrime = 1.GFUint
    for coef in errLocPrimeTmp:
      errLocPrime = errLocPrime * coef
    # equivalent to: errLocPrime = functools.reduce(gf_mul, errLocPrimeTmp, 1)

    # Compute y (evaluation of the errata evaluator polynomial)
    # This is a more faithful translation of the theoretical
    # equation contrary to the old forney method. Here it is an
    # exact reproduction:
    # Yl = omega(Xl.inverse()) / prod(1 - Xj*Xl.inverse()) for j in len(X)
    var
      y = errEval.reversed().eval(Xi_inv) # numerator of the Forney algorithm (errata evaluator evaluated)
    y = Xi ^ (1 - fcr) * y # adjust to fcr parameter

    # Check: errLocPrime (the divisor) should not be zero.
    if errLocPrime == 0:
      raise newException(CatchableError, "Could not find error magnitude")

    # Compute the magnitude
    var
      # magnitude value of the error, calculated by the Forney algorithm
      # (an equation in fact): dividing the errata evaluator with the errata
      # locator derivative gives us the errata magnitude (ie, value to repair)
      # the ith symbol
      magnitude = y / errLocPrime

    E[errPos[i]] = magnitude # store the magnitude for this error into the magnitude polynomial

  # Apply the correction of values to get our message corrected!
  # NOTE: that the ecc bytes also gets corrected!
  # this isn't the Forney algorithm, we just apply the result of decoding here

  # equivalent to Ci = Ri - Ei where Ci is the correct message,
  # Ri the received (senseword) message, and Ei the errata magnitudes
  # (minus is replaced by XOR since it's equivalent in GF(2^p)).
  # So in fact here we substract from the received message the errors
  # magnitude, which logically corrects the value to what it should be.
  return (msg + E)

proc findErrors(
  errLoc: seq[GFUint],
  nmess: int,
  generator = 2.GFUint): seq[int] {.raises: [CatchableError].} = # nmess is len(msg)
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
  for i in 0..<nmess:
    if errLoc.eval(generator ^ i) == 0: # It's a 0? Bingo, it's a root of the error locator polynomial,
                                                         # in other terms this is the location of an error
      errPos.add(nmess - 1 - i)

  # Sanity check: the number of errors/errata positions found should be
  # exactly the same as the length of the errata locator polynomial
  if errPos.len != errs:
    # couldn't find error locations
    raise newException(
      CatchableError,
      "Too many (or few) errors found by Chien Search for the errata locator polynomial!")

  return errPos

proc forneySyndromes(
  synd: seq[GFUint],
  pos: seq[int],
  nmess: int,
  generator = 2.GFUint): seq[GFUint] =
  ## Compute Forney syndromes, which computes a modified
  ## syndromes to compute only errors (erasures are trimmed out).
  ## Do not confuse this with Forney algorithm, which allows
  ## to correct the message based on the location of errors.
  ##

  var
    erasePosReversed = pos.mapIt(nmess - 1 - it) # prepare the coefficient degree positions (instead of the erasures positions)

  # Optimized method, all operations are inlined
  var
    fsynd = synd[1..synd.high]      # make a copy and trim the first coefficient which is always 0 by definition
  for i in 0..<pos.len:
    let
      x = generator ^ erasePosReversed[i]
    for j in 0..<(fsynd.len - 1):
      fsynd[j] = (fsynd[j] * x) xor fsynd[j + 1]
    #fsynd.pop() # useless? it doesn't change the results of computations to leave it there

  # Theoretical way of computing the modified Forney syndromes: fsynd = (eraseLoc * synd) % x^(n-k)
  # See Shao, H. M., Truong, T. K., Deutsch, L. J., & Reed, I. S. (1986, April). A single chip VLSI Reed-Solomon decoder. In Acoustics, Speech, and Signal Processing, IEEE International Conference on ICASSP'86. (Vol. 11, pp. 2151-2154). IEEE.ISO 690
  #eraseLoc = findErrataLocator(erasePosReversed, generator=generator) # computing the erasures locator polynomial
  #fsynd = gf_poly_mul(eraseLoc[::-1], synd[1:]) # then multiply with the syndrome to get the untrimmed forney syndrome
  #fsynd = fsynd[len(pos):] # then trim the first erasePos coefficients which are useless. Seems to be not necessary, but this reduces the computation time later in BM (thus it's an optimization).

  return fsynd

proc correctMsg*(
  msg: seq[GFUint],
  nsym: int,
  fcr = 0,
  generator = Char.GFUint,
  erasePos: seq[int] = @[],
  erasures = false):
  (seq[GFUint], seq[GFUint]) {.raises: [ValueError, CatchableError].} =
  ## Reed-Solomon main decoding function
  ##

  if msg.len > Order.int:
    # Note that it is in fact possible to encode/decode messages
    # that are longer than Order, but because this will be above
    # the field, this will generate more error positions during
    # Chien Search than it should, because this will generate
    # duplicate values, which should normally be prevented thank's
    # to the prime polynomial reduction - eg, because it can't
    # discriminate between error at position 1 or 256, both being
    # exactly equal under galois field 2^8. So it's really not
    # advised to do it, but it's possible - but then you're not
    # guaranted to be able to correct any error/erasure on symbols
    # with a position above the length of Order -- if you really
    # need a bigger message without chunking, then you should
    # better enlarge c_exp so that you get a bigger field.
    raise newException(
      ValueError,
      "Message is too long (%1 when max is %2)" % [$len(msg), $Order])

  var
    msgOut = msg  # copy of message

  # check if there are too many erasures to correct (beyond the Singleton bound)
  if erasePos.len > nsym:
    raise newException(CatchableError, "Too many erasures to correct")

  # prepare the syndrome polynomial using only errors
  # (ie: errors = characters that were either replaced by null byte
  # or changed to another character, but we don't know their positions)
  var synd = clacSyndromes(msgOut, nsym, fcr, generator)
  # check if there's any error/erasure in the input codeword.
  # If not (all syndromes coefficients are 0), then just return the message as-is.

  if max(synd) == 0:
    return (
      msgOut[0..<(msgOut.len - nsym)],
      msgOut[(msgOut.len - nsym)..<msgOut.len])  # no errors

  # Find errors locations
  var errPos: seq[int]
  if not erasures:
    let
      # compute the Forney syndromes, which hide the erasures from the original
      # syndrome (so that BM will just have to deal with errors, not erasures)
      fsynd = forneySyndromes(synd, erasePos, len(msgOut), generator)
      # compute the error locator polynomial using Berlekamp-Massey
      errLoc = findErrorLocator(fsynd, nsym, eraseCount=len(erasePos))
      # locate the message errors using Chien search (or bruteforce search)

    errPos = findErrors(errLoc.reversed(), len(msgOut), generator)
    if errPos.len <= 0:
      raise newException(CatchableError, "Could not locate error")

  # Find errors values and apply them to correct the message
  # compute errata evaluator and errata magnitude polynomials,
  # then correct errors and erasures

  # NOTE: we here use the original syndrome, not the forney syndrome
  # (because we will correct both errors and erasures, so we need the full syndrome)
  msgOut = correctErrata(
    msgOut,
    synd,
    erasePos & errPos,
    fcr,
    generator)

  # check if the final message is fully repaired
  synd = clacSyndromes(msgOut, nsym, fcr, generator)
  if max(synd) > 0:
    raise newException(CatchableError, "Could not correct message")

  # return the successfully decoded message
  return (
    msgOut[0..<(msgOut.len-nsym)],
    msgOut[(msgOut.len-nsym)..<msgOut.len]) # also return the corrected ecc block so that the user can check()

proc correctMsgNoFsynd(
  msg: seq[GFUint],
  nsym: int,
  fcr = 0,
  generator = Char.GFUint,
  erasePos: seq[int] = @[],
  erasures = false):
  (seq[GFUint], seq[GFUint]) {.raises: [ValueError, CatchableError].} =
  ## Reed-Solomon main decoding function, without using the modified Forney syndromes
  ## This demonstrates how the decoding process is done without using the Forney syndromes
  ## (this is the most common way nowadays, avoiding Forney syndromes require to use a modified
  ## Berlekamp-Massey that will take care of the erasures by itself, it's a simple matter of
  ## modifying some initialization variables and the loop ranges)
  ##

  if msg.len > Order.int:
    raise newException(
      ValueError,
      "Message is too long ($1 when max is $2)" % [$len(msg), $Order])

  var
    msgOut = msg  # copy of message

  # erasures: set them to null bytes for easier decoding
  # (but this is not necessary, they will be corrected anyway,
  # but debugging will be easier with null bytes because
  # the error locator polynomial values will only depend on
  # the errors locations, not their values)
  if erasePos.len > 0:
    for ePos in erasePos:
      msgOut[ePos] = 0.GFUint

  # check if there are too many erasures
  if len(erasePos) > nsym:
    raise newException(CatchableError, "Too many erasures to correct")

  # prepare the syndrome polynomial using only errors
  # (ie: errors = characters that were either replaced
  # by null byte or changed to another character, but
  # we don't know their positions)
  var
    synd = clacSyndromes(msgOut, nsym, fcr, generator)

  # check if there's any error/erasure in the input codeword.
  # If not (all syndromes coefficients are 0), then just return
  # the codeword as-is.
  if max(synd) == 0:
    return (msgOut[0..<nsym], msgOut[nsym..<msgOut.len])  # no errors

  # prepare erasures locator and evaluator polynomials
  var
    eraseLoc: seq[GFUint]
    #erase_eval = None
    eraseCount = 0

  if erasePos.len > 0:
    eraseCount = erasePos.len
    eraseLoc = findErrataLocator(
      erasePos.mapIt((msgOut.len - 1 - it)),
      generator = generator)
    #erase_eval = findErrorEvaluator(synd[::-1], eraseLoc, len(eraseLoc)-1)

  # prepare errors/errata locator polynomial
  var errLoc: seq[GFUint]
  if erasures:
    errLoc = eraseLoc.reversed()
    #errEval = erase_eval[::-1]
  else:
    errLoc = findErrorLocator(synd, nsym, eraseLoc=eraseLoc, eraseCount=eraseCount)
    errLoc = errLoc.reversed()
    #errEval = findErrorEvaluator(synd[::-1], errLoc[::-1], len(errLoc)-1)[::-1] # find error/errata evaluator polynomial (not really necessary since we already compute it at the same time as the error locator poly in BM)

  # locate the message errors
  var errPos = findErrors(errLoc, len(msgOut), generator) # find the roots of the errata locator polynomial (ie: the positions of the errors/errata)
  if errPos.len <= 0:
    raise newException(CatchableError, "Could not locate error")

  # compute errata evaluator and errata magnitude polynomials, then correct errors and erasures
  msgOut = correctErrata(msgOut, synd, errPos, fcr = fcr, generator = generator)
  # check if the final message is fully repaired
  synd = clacSyndromes(msgOut, nsym, fcr, generator)

  if max(synd) > 0:
    raise newException(CatchableError, "Could not correct message")

  # return the successfully decoded message
  return (msgOut[0..nsym], msgOut[nsym..<msgOut.len]) # also return the corrected ecc block so that the user can check()

proc rs_check(msg: seq[GFUint], nsym: int, fcr = 0, generator = Char.GFUint): bool =
  ## Returns true if the message + ecc has no error of false
  ## otherwise (may not always catch a wrong decoding or a wrong
  ## message, particularly if there are too many errors
  ## -- above the Singleton bound --, but it usually does)
  return ( max(clacSyndromes(msg, nsym, fcr, generator)) == 0 )

when isMainModule:
  import stew/byteutils

  let msg = @[
    0x48, 0x65, 0x6C, 0x6C, 0x6F, 0x20, 0x52,
    0x65, 0x65, 0x64, 0x2D, 0x53, 0x6F, 0x6C,
    0x6F, 0x6D, 0x6F, 0x6E, 0x21 ]
    .mapIt(it.GFUint)

  var msgOut = encode(msg, 10)

  # var strMsg: string
  # for i in 0..<len(msg):
  #   strMsg &= "0x" & uint8(msg[i]).toHex
  #   strMsg &= " "

  echo msgOut

  msgOut[0] = 0.GFUint
  # msgOut[10] = 0.GFUint
  # msgOut[16] = 0.GFUint
  # msgOut[14] = 0.GFUint
  # msgOut[26] = 0.GFUint

  let (correct, code) = correctMsg(msgOut, 10)

  echo correct
  # strMsg = ""
  # for i in 0..<len(correct):
  #   strMsg &= "0x" & uint8(correct[i]).toHex
  #   strMsg &= " "

  # echo "MESSAGE ", strMsg
