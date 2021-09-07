{.push raises: [Defect].}

{.experimental: "codeReordering".}

import ./polynomial

# ################### REED-SOLOMON ENCODING ###################

proc rs_generator_poly(nsym: int, fcr: GFUint = 0, generator: GFUint = 2): GFPolynomial =
  ## Generate an irreducible generator polynomial
  ## (necessary to encode a message into Reed-Solomon)
  ##

  var g = @[1.GFUint]
  for i in 0..<nsym:
    g = g * @[1.GFUint, generator ^ (i + fcr)]

  return g

proc rs_generator_poly_all(max_nsym: int, fcr: GFUint = 0, generator: GFUint = 2): seq[seq[int]] =
  ## Generate all irreducible generator polynomials up to max_nsym
  ## (usually you can use n, the length of the message+ecc). Very
  ## useful to reduce processing time if you want to encode using
  ## variable schemes and nsym rates.
  ##

  var g_all = newSeq[seq[GFUint]](max_nsym)
  g_all[0] = @[1.GFUint]
  g_all[1] = @[1.GFUint]
  for nsym in 0..<max_nsym:
    g_all[nsym] = rs_generator_poly(nsym, fcr, generator)

  return g_all

# proc rs_simple_encode_msg(
#   msg_in: seq[int],
#   nsym: int,
#   fcr = 0,
#   generator = 2,
#   gen: seq[int] = @[]): seq[int] =
#   ## Simple Reed-Solomon encoding (mainly an example for you to understand
#   ## how it works, because it's slower than the inlined function below)
#   ##

#   if (len(msg_in) + nsym) > field_charac:
#     raise newException(
#       ValueError, "Message is too long ($1 when max is $2)" % [$(len(msg_in) + nsym), $field_charac])

#   var gen = if gen.len <= 0: rs_generator_poly(nsym, fcr, generator) else: gen

#   # Pad the message, then divide it by the irreducible generator polynomial
#   var (_, remainder) = gf_poly_div(msg_in[0..<len(gen)-1], gen)

#   # The remainder is our RS code! Just append it to our original
#   # message to get our full codeword (this represents a polynomial
#   # of max 256 terms)
#   return msg_in & remainder

# proc rs_encode_msg(
#   msg_in: seq[int],
#   nsym: int,
#   fcr = 0,
#   generator = 2,
#   gen: seq[int] = @[]): seq[int] =
#   ## Reed-Solomon main encoding function, using polynomial division
#   ## (Extended Synthetic Division, the fastest algorithm available to my knowledge),
#   ## better explained at http://research.swtch.com/field
#   ##

#   if (len(msg_in) + nsym) > field_charac:
#     raise newException(
#       ValueError,
#       "Message is too long ($1 when max is $2)" % [$(len(msg_in) + nsym), $field_charac])

#   var gen = if gen.len <= 0:
#       rs_generator_poly(nsym, fcr, generator)
#     else:
#       gen

#   # Init msg_out with the values inside msg_in and pad with len(gen)-1 bytes (which is the number of ecc symbols).
#   # Initializing the Synthetic Division with the dividend (= input message polynomial)
#   var msg_out = msg_in
#   msg_out.setLen((len(msg_in) + (len(gen)-1)))

#   # Synthetic division main loop
#   for i in 0..<len(msg_in):
#     # Note that it's msg_out here, not msg_in. Thus, we reuse the updated value at each iteration
#     # (this is how Synthetic Division works: instead of storing in a temporary register the intermediate values,
#     # we directly commit them to the output).

#     let coef = msg_out[i]
#     # log(0) is undefined, so we need to manually check for this case.
#     if coef != 0:
#       # in synthetic division, we always skip the first coefficient of the divisior, because it's only used to normalize the dividend coefficient (which is here useless since the divisor, the generator polynomial, is always monic)
#       for j in 1..<len(gen):
#         #if gen[j] != 0: # log(0) is undefined so we need to check that, but it slow things down in fact and it's useless in our case (reed-solomon encoding) since we know that all coefficients in the generator are not 0
#         msg_out[i+j] = (msg_out[i+j] xor gf_mul(gen[j], coef)) # equivalent to msg_out[i+j] += gf_mul(gen[j], coef)

#   # At this point, the Extended Synthetic Divison is done, msg_out contains the quotient in msg_out[:len(msg_in)]
#   # and the remainder in msg_out[len(msg_in):]. Here for RS encoding, we don't need the quotient but only the remainder
#   # (which represents the RS code), so we can just overwrite the quotient with the input message, so that we get
#   # our complete codeword composed of the message + code.
#   return msg_in &  msg_out[msg_in.len..msg_out.high]


# ################### REED-SOLOMON DECODING ###################

# proc rs_calc_syndromes(
#   msg: seq[int],
#   nsym: int,
#   fcr = 0,
#   generator = 2): seq[int] =
#   ## Given the received codeword msg and the number of
#   ## error correcting symbols (nsym), computes the syndromes
#   ## polynomial.
#   ##
#   ## Mathematically, it's essentially equivalent to a Fourrier Transform
#   ## (Chien search being the inverse).
#   ##

#   # Note the "[0] +" : we add a 0 coefficient for the lowest degree (the constant).
#   # This effectively shifts the syndrome, and will shift every computations depending
#   # on the syndromes (such as the errors locator polynomial, errors evaluator polynomial,
#   # etc. but not the errors positions).
#   #
#   # This is not necessary, you can adapt subsequent computations to start from 0 instead
#   # of skipping the first iteration (ie, the often seen range(1, n-k+1)),
#   var synd = repeat(0, nsym)
#   for i in 0..<nsym:
#     synd[i] = gf_poly_eval(msg, gf_pow(generator, (i + fcr)))

#   return @[0] & synd # pad with one 0 for mathematical precision (else we can end up with weird calculations sometimes)

# proc rs_find_error_locator(
#   synd: GFPolynomial,
#   nsym: int,
#   erase_loc: GFPolynomial = @[],
#   erase_count = 0): GFPolynomial =
#   ## Find error/errata locator and evaluator polynomials with Berlekamp-Massey algorithm
#   ##

#   # The idea is that BM will iteratively estimate the error locator polynomial.
#   # To do this, it will compute a Discrepancy term called Delta, which will tell
#   # us if the error locator polynomial needs an update or not (hence why it's called
#   # discrepancy: it tells us when we are getting off board from the correct value).

#   # Init the polynomials
#   var (err_loc, old_loc) = if erase_loc.len > 0: # if the erasure locator polynomial is supplied, we init with its value, so that we include erasures in the final locator polynomial
#     (erase_loc, erase_loc)
#   else:
#     (@[1], @[1]) # This is the main variable we want to fill, also called Sigma in other notations or more formally the errors/errata locator polynomial.
#                  # BM is an iterative algorithm, and we need the errata locator polynomial of the previous iteration in order to update other necessary variables.
#   #L = 0 # update flag variable, not needed here because we use an alternative equivalent way of checking if update is needed (but using the flag could potentially be faster depending on if using length(list) is taking linear time in your language, here in Python it's constant so it's as fast.

#   # Fix the syndrome shifting: when computing the syndrome, some implementations may prepend a 0 coefficient for the lowest degree term (the constant). This is a case of syndrome shifting, thus the syndrome will be bigger than the number of ecc symbols (I don't know what purpose serves this shifting). If that's the case, then we need to account for the syndrome shifting when we use the syndrome such as inside BM, by skipping those prepended coefficients.
#   # Another way to detect the shifting is to detect the 0 coefficients: by definition, a syndrome does not contain any 0 coefficient (except if there are no errors/erasures, in this case they are all 0). This however doesn't work with the modified Forney syndrome, which set to 0 the coefficients corresponding to erasures, leaving only the coefficients corresponding to errors.
#   var synd_shift = if len(synd) > nsym: len(synd) - nsym else: 0

#   for i in 0..<(nsym - erase_count): # generally: nsym-erase_count == len(synd), except when you input a partial erase_loc and using the full syndrome instead of the Forney syndrome, in which case nsym-erase_count is more correct (len(synd) will fail badly with IndexError).
#     var K = if erase_loc.len > 0: # if an erasures locator polynomial was provided to init the errors locator polynomial, then we must skip the FIRST erase_count iterations (not the last iterations, this is very important!)
#       (erase_count + i + synd_shift)
#     else: # if erasures locator is not provided, then either there's no erasures to account or we use the Forney syndromes, so we don't need to use erase_count nor erase_loc (the erasures have been trimmed out of the Forney syndromes).
#       i + synd_shift

#     # Compute the discrepancy Delta
#     # Here is the close-to-the-books operation to compute the discrepancy Delta: it's a simple polynomial multiplication of error locator with the syndromes, and then we get the Kth element.
#     #delta = gf_poly_mul(err_loc[::-1], synd)[K] # theoretically it should be gf_poly_add(synd[::-1], [1])[::-1] instead of just synd, but it seems it's not absolutely necessary to correctly decode.
#     # But this can be optimized: since we only need the Kth element, we don't need to compute the polynomial multiplication for any other element but the Kth. Thus to optimize, we compute the polymul only at the item we need, skipping the rest (avoiding a nested loop, thus we are linear time instead of quadratic).
#     # This optimization is actually described in several figures of the book "Algebraic codes for data transmission", Blahut, Richard E., 2003, Cambridge university press.
#     var delta = synd[K]
#     for j in 1..<len(err_loc):
#       delta = (delta xor err_loc[^(j+1)] * synd[K - j]) # delta is also called discrepancy. Here we do a partial polynomial multiplication (ie, we compute the polynomial multiplication only for the term of degree K). Should be equivalent to brownanrs.polynomial.mul_at().
#     #print "delta", K, delta, list(gf_poly_mul(err_loc[::-1], synd)) # debugline

#     # Shift polynomials to compute the next degree
#     old_loc = old_loc & @[0]

#     # Iteratively estimate the errata locator and evaluator polynomials
#     if delta != 0: # Update only if there's a discrepancy
#         if len(old_loc) > len(err_loc): # Rule B (rule A is implicitly defined because rule A just says that we skip any modification for this iteration)
#         #if 2*L <= K+erase_count: # equivalent to len(old_loc) > len(err_loc), as long as L is correctly computed
#             # Computing errata locator polynomial Sigma
#             var new_loc = old_loc.scale(delta)
#             old_loc = err_loc.scale(delta.inverse) # effectively we are doing err_loc * 1/delta = err_loc // delta
#             err_loc = new_loc
#             # Update the update flag
#             #L = K - L # the update flag L is tricky: in Blahut's schema, it's mandatory to use `L = K - L - erase_count` (and indeed in a previous draft of this function, if you forgot to do `- erase_count` it would lead to correcting only 2*(errors+erasures) <= (n-k) instead of 2*errors+erasures <= (n-k)), but in this latest draft, this will lead to a wrong decoding in some cases where it should correctly decode! Thus you should try with and without `- erase_count` to update L on your own implementation and see which one works OK without producing wrong decoding failures.

#         # Update with the discrepancy
#         err_loc = err_loc + old_loc.scale(delta)

#   # Check if the result is correct, that there's not too many errors to correct
#   var drop = 0
#   while len(err_loc) > 0 and err_loc[0] == 0: drop.inc()
#   err_loc = if drop > 0: err_loc[0..<drop] else: err_loc

#   var errs = len(err_loc) - 1
#   if ((errs - erase_count) * 2) + erase_count > nsym:
#     raise newException(CatchableError, "Too many errors to correct")    # too many errors to correct

#   return err_loc

# proc rs_find_errata_locator(e_pos: seq[int], generator: GFUint = 2): GFPolynomial =
#   ## Compute the erasures/errors/errata locator polynomial from the
#   ## erasures/errors/errata positions (the positions must be relative
#   ## to the x coefficient, eg: "hello worldxxxxxxxxx" is tampered to
#   ## "h_ll_ worldxxxxxxxxx" with xxxxxxxxx being the ecc of length
#   ## n-k=9, here the string positions are [1, 4], but the coefficients
#   ## are reversed since the ecc characters are placed as the first
#   ## coefficients of the polynomial, thus the coefficients of the erased
#   ## characters are n-1 - [1, 4] = [18, 15] = erasures_loc to be specified
#   ## as an argument.
#   ##

#   var e_loc = @[1] # just to init because we will multiply, so it must be 1 so that the multiplication starts correctly without nulling any term
#   # erasures_loc = product(1 - x*alpha**i) for i in erasures_pos and where alpha is the alpha chosen to evaluate polynomials.
#   for i in e_pos:
#     e_loc = e_loc * @[1.GFUint] + @[generator ^ i, 0]

#   return e_loc

# proc rs_find_error_evaluator(synd: GFPolynomial, err_loc: GFPolynomial, nsym: int): GFPolynomial =
#   ## Compute the error (or erasures if you supply
#   ## sigma=erasures locator polynomial, or errata)
#   ## evaluator polynomial Omega from the syndrome
#   ## and the error/erasures/errata locator Sigma.
#   ##

#   # Omega(x) = [ Synd(x) * Error_loc(x) ] mod x^(n-k+1)

#   # first multiply syndromes * errata_locator, then do a
#   # polynomial division to truncate the polynomial to the
#   # required length
#   var (_, remainder) = synd div err_loc * @[1] & newSeq[GFUint](nsym + 1)

#   # Faster way that is equivalent
#   #remainder = gf_poly_mul(synd, err_loc) # first multiply the syndromes with the errata locator polynomial
#   #remainder = remainder[len(remainder)-(nsym+1):] # then slice the list to truncate it (which represents the polynomial), which
#                                                     # is equivalent to dividing by a polynomial of the length we want

#   return remainder

# proc rs_correct_errata(
#   msg_in: seq[int],
#   synd: seq[int],
#   err_pos: seq[int],
#   fcr = 0,
#   generator = 2): seq[int] = # err_pos is a list of the positions of the errors/erasures/errata
#   ## Forney algorithm, computes the values (error magnitude)
#   ## to correct the input message.
#   ##

#   # calculate errata locator polynomial to correct both errors and erasures (by combining the errors positions given by the error locator polynomial found by BM with the erasures positions given by caller)
#   var
#     coef_pos = err_pos.mapIt(len(msg_in) - 1 - it) # need to convert the positions to coefficients degrees for the errata locator algo to work (eg: instead of [0, 1, 2] it will become [len(msg)-1, len(msg)-2, len(msg) -3])
#     err_loc = rs_find_errata_locator(coef_pos, generator)
#     # calculate errata evaluator polynomial (often called Omega or Gamma in academic papers)
#     err_eval = rs_find_error_evaluator(synd.reversed(), err_loc, (len(err_loc)-1)).reversed()

#   # Second part of Chien search to get the error location polynomial X from
#   # the error positions in err_pos (the roots of the error locator polynomial,
#   # ie, where it evaluates to 0)
#   var X: seq[int] # will store the position of the errors
#   for i in 0..<len(coef_pos):
#     X.add(
#       gf_pow(
#         generator,
#         (gf_log.len - (field_charac - coef_pos[i])-1)) )

#   # Forney algorithm: compute the magnitudes
#   var E = repeat(0, len(msg_in)) # will store the values that need to be corrected (substracted) to the message containing errors. This is sometimes called the error magnitude polynomial.
#   var Xlength = len(X)
#   for i, Xi in X:
#     var Xi_inv = gf_inverse(Xi)

#     # Compute the formal derivative of the error locator polynomial (see Blahut, Algebraic codes for data transmission, pp 196-197).
#     # the formal derivative of the errata locator is used as the denominator of the Forney Algorithm, which simply says that the ith error value is given by error_evaluator(gf_inverse(Xi)) / error_locator_derivative(gf_inverse(Xi)). See Blahut, Algebraic codes for data transmission, pp 196-197.
#     var err_loc_prime_tmp: seq[int]
#     for j in 0..<Xlength:
#       if j != i:
#           err_loc_prime_tmp.add( gf_sub(1, gf_mul(Xi_inv, X[j])) )

#     # compute the product, which is the denominator of the Forney algorithm (errata locator derivative)
#     var err_loc_prime = 1
#     for coef in err_loc_prime_tmp:
#       err_loc_prime = gf_mul(err_loc_prime, coef)
#     # equivalent to: err_loc_prime = functools.reduce(gf_mul, err_loc_prime_tmp, 1)

#     # Compute y (evaluation of the errata evaluator polynomial)
#     # This is a more faithful translation of the theoretical equation contrary to the old forney method. Here it is an exact reproduction:
#     # Yl = omega(Xl.inverse()) / prod(1 - Xj*Xl.inverse()) for j in len(X)
#     var y = gf_poly_eval(err_eval.reversed(), Xi_inv) # numerator of the Forney algorithm (errata evaluator evaluated)
#     y = gf_mul(gf_pow(Xi, 1 - fcr), y) # adjust to fcr parameter

#     # Check: err_loc_prime (the divisor) should not be zero.
#     if err_loc_prime == 0:
#       raise newException(CatchableError, "Could not find error magnitude")    # Could not find error magnitude

#     # Compute the magnitude
#     var magnitude = gf_div(y, err_loc_prime) # magnitude value of the error, calculated by the Forney algorithm (an equation in fact): dividing the errata evaluator with the errata locator derivative gives us the errata magnitude (ie, value to repair) the ith symbol
#     E[err_pos[i]] = magnitude # store the magnitude for this error into the magnitude polynomial

#   # Apply the correction of values to get our message corrected! (note that the ecc bytes also gets corrected!)
#   # (this isn't the Forney algorithm, we just apply the result of decoding here)
#   return gf_poly_add(msg_in, E) # equivalent to Ci = Ri - Ei where Ci is the correct message, Ri the received (senseword) message, and Ei the errata magnitudes (minus is replaced by XOR since it's equivalent in GF(2^p)). So in fact here we substract from the received message the errors magnitude, which logically corrects the value to what it should be.

# proc rs_find_errors(err_loc: seq[int], nmess: int, generator = 2): seq[int] = # nmess is len(msg_in)
#   ## Find the roots (ie, where evaluation = zero) of error polynomial
#   ## by brute-force trial, this is a sort of Chien's search
#   ## (but less efficient, Chien's search is a way to evaluate the
#   ## polynomial such that each evaluation only takes constant time).
#   ##

#   var errs = len(err_loc) - 1
#   var err_pos: seq[int]
#   for i in 0..<nmess: # normally we should try all 2^8 possible values, but here we optimize to just check the interesting symbols
#     if gf_poly_eval(err_loc, gf_pow(generator, i)) == 0: # It's a 0? Bingo, it's a root of the error locator polynomial,
#                                                          # in other terms this is the location of an error
#       err_pos.add(nmess - 1 - i)

#   # Sanity check: the number of errors/errata positions found should be exactly the same as the length of the errata locator polynomial
#   if len(err_pos) != errs:
#     # couldn't find error locations
#     raise newException(
#       CatchableError,
#       "Too many (or few) errors found by Chien Search for the errata locator polynomial!")

#   return err_pos

# proc rs_forney_syndromes(synd: seq[int], pos: seq[int], nmess: int, generator = 2): seq[int] =
#   # Compute Forney syndromes, which computes a modified syndromes to compute only errors (erasures are trimmed out). Do not confuse this with Forney algorithm, which allows to correct the message based on the location of errors.
#   var erase_pos_reversed = pos.mapIt(nmess - 1 - it) # prepare the coefficient degree positions (instead of the erasures positions)

#   # Optimized method, all operations are inlined
#   var fsynd = synd[1..synd.high]      # make a copy and trim the first coefficient which is always 0 by definition
#   for i in 0..<len(pos):
#     let x = gf_pow(generator, erase_pos_reversed[i])
#     for j in 0..<(len(fsynd) - 1):
#       fsynd[j] = gf_mul(fsynd[j], x) xor fsynd[j + 1]
#     #fsynd.pop() # useless? it doesn't change the results of computations to leave it there

#   # Theoretical way of computing the modified Forney syndromes: fsynd = (erase_loc * synd) % x^(n-k)
#   # See Shao, H. M., Truong, T. K., Deutsch, L. J., & Reed, I. S. (1986, April). A single chip VLSI Reed-Solomon decoder. In Acoustics, Speech, and Signal Processing, IEEE International Conference on ICASSP'86. (Vol. 11, pp. 2151-2154). IEEE.ISO 690
#   #erase_loc = rs_find_errata_locator(erase_pos_reversed, generator=generator) # computing the erasures locator polynomial
#   #fsynd = gf_poly_mul(erase_loc[::-1], synd[1:]) # then multiply with the syndrome to get the untrimmed forney syndrome
#   #fsynd = fsynd[len(pos):] # then trim the first erase_pos coefficients which are useless. Seems to be not necessary, but this reduces the computation time later in BM (thus it's an optimization).

#   return fsynd

# proc rs_correct_msg(
#   msg_in: seq[int],
#   nsym: int,
#   fcr = 0,
#   generator = 2,
#   erase_pos: seq[int] = @[],
#   only_erasures = false): (seq[int], seq[int]) =
#   ## Reed-Solomon main decoding function
#   ##

#   if len(msg_in) > field_charac:
#     # Note that it is in fact possible to encode/decode messages that are longer than field_charac, but because this will be above the field, this will generate more error positions during Chien Search than it should, because this will generate duplicate values, which should normally be prevented thank's to the prime polynomial reduction (eg, because it can't discriminate between error at position 1 or 256, both being exactly equal under galois field 2^8). So it's really not advised to do it, but it's possible (but then you're not guaranted to be able to correct any error/erasure on symbols with a position above the length of field_charac -- if you really need a bigger message without chunking, then you should better enlarge c_exp so that you get a bigger field).
#     raise newException(ValueError, "Message is too long (%1 when max is %2)" % [$len(msg_in), $field_charac])

#   var msg_out = msg_in     # copy of message
#   # erasures: set them to null bytes for easier decoding (but this is not necessary, they will be corrected anyway, but debugging will be easier with null bytes because the error locator polynomial values will only depend on the errors locations, not their values)
#   if erase_pos.len > 0:
#     for e_pos in erase_pos:
#       msg_out[e_pos] = 0

#   # check if there are too many erasures to correct (beyond the Singleton bound)
#   if len(erase_pos) > nsym:
#     raise newException(CatchableError, "Too many erasures to correct")

#   # prepare the syndrome polynomial using only errors (ie: errors = characters that were either replaced by null byte
#   # or changed to another character, but we don't know their positions)
#   var synd = rs_calc_syndromes(msg_out, nsym, fcr, generator)
#   # check if there's any error/erasure in the input codeword.
#   #mIf not (all syndromes coefficients are 0), then just return the message as-is.

#   if max(synd) == 0:
#     return (
#       msg_out[0..<(msg_out.len-nsym)],
#       msg_out[(msg_out.len-nsym)..<msg_out.len])  # no errors

#   # Find errors locations
#   var err_pos: seq[int]
#   if not only_erasures:
#     let
#       # compute the Forney syndromes, which hide the erasures from the original
#       # syndrome (so that BM will just have to deal with errors, not erasures)
#       fsynd = rs_forney_syndromes(synd, erase_pos, len(msg_out), generator)
#       # compute the error locator polynomial using Berlekamp-Massey
#       err_loc = rs_find_error_locator(fsynd, nsym, erase_count=len(erase_pos))
#       # locate the message errors using Chien search (or bruteforce search)

#     err_pos = rs_find_errors(err_loc.reversed(), len(msg_out), generator)
#     if err_pos.len <= 0:
#       raise newException(CatchableError, "Could not locate error")

#   # Find errors values and apply them to correct the message
#   # compute errata evaluator and errata magnitude polynomials, then correct errors and erasures
#   msg_out = rs_correct_errata(
#     msg_out,
#     synd,
#     erase_pos & err_pos,
#     fcr,
#     generator) # note that we here use the original syndrome, not the forney syndrome
#                                                         # (because we will correct both errors and erasures, so we need the full syndrome)
#   # check if the final message is fully repaired
#   synd = rs_calc_syndromes(msg_out, nsym, fcr, generator)
#   if max(synd) > 0:
#     raise newException(CatchableError, "Could not correct message")

#   # return the successfully decoded message
#   return (
#     msg_out[0..<(msg_out.len-nsym)],
#     msg_out[(msg_out.len-nsym)..<msg_out.len]) # also return the corrected ecc block so that the user can check()

# proc rs_correct_msg_nofsynd(
#   msg_in: seq[int],
#   nsym: int,
#   fcr = 0,
#   generator = 2,
#   erase_pos: seq[int] = @[],
#   only_erasures = false): (seq[int], seq[int]) =
#   ## Reed-Solomon main decoding function, without using the modified Forney syndromes
#   ## This demonstrates how the decoding process is done without using the Forney syndromes
#   ## (this is the most common way nowadays, avoiding Forney syndromes require to use a modified
#   ## Berlekamp-Massey that will take care of the erasures by itself, it's a simple matter of
#   ## modifying some initialization variables and the loop ranges)
#   ##

#   if len(msg_in) > field_charac:
#     raise newException(ValueError, "Message is too long ($1 when max is $2)" % [$len(msg_in), $field_charac])

#   var msg_out = msg_in     # copy of message
#   # erasures: set them to null bytes for easier decoding (but this is not necessary, they will be corrected anyway, but debugging will be easier with null bytes because the error locator polynomial values will only depend on the errors locations, not their values)
#   if erase_pos.len > 0:
#     for e_pos in erase_pos:
#       msg_out[e_pos] = 0

#   # check if there are too many erasures
#   if len(erase_pos) > nsym:
#     raise newException(CatchableError, "Too many erasures to correct")

#   # prepare the syndrome polynomial using only errors (ie: errors = characters that were either replaced by null byte or changed to another character, but we don't know their positions)
#   var synd = rs_calc_syndromes(msg_out, nsym, fcr, generator)
#   # check if there's any error/erasure in the input codeword. If not (all syndromes coefficients are 0), then just return the codeword as-is.
#   if max(synd) == 0:
#     return (msg_out[0..<nsym], msg_out[nsym..<msg_out.len])  # no errors

#   # prepare erasures locator and evaluator polynomials
#   var
#     erase_loc: seq[int]
#     #erase_eval = None
#     erase_count = 0

#   if erase_pos.len > 0:
#     erase_count = len(erase_pos)
#     erase_loc = rs_find_errata_locator(
#       erase_pos.mapIt(len(msg_out)-1-it), generator=generator)
#     #erase_eval = rs_find_error_evaluator(synd[::-1], erase_loc, len(erase_loc)-1)

#   # prepare errors/errata locator polynomial
#   var err_loc: seq[int]
#   if only_erasures:
#     err_loc = erase_loc.reversed()
#     #err_eval = erase_eval[::-1]
#   else:
#     err_loc = rs_find_error_locator(synd, nsym, erase_loc=erase_loc, erase_count=erase_count)
#     err_loc = err_loc.reversed()
#     #err_eval = rs_find_error_evaluator(synd[::-1], err_loc[::-1], len(err_loc)-1)[::-1] # find error/errata evaluator polynomial (not really necessary since we already compute it at the same time as the error locator poly in BM)

#   # locate the message errors
#   var err_pos = rs_find_errors(err_loc, len(msg_out), generator) # find the roots of the errata locator polynomial (ie: the positions of the errors/errata)
#   if err_pos.len <= 0:
#     raise newException(CatchableError, "Could not locate error")

#   # compute errata evaluator and errata magnitude polynomials, then correct errors and erasures
#   msg_out = rs_correct_errata(msg_out, synd, err_pos, fcr = fcr, generator = generator)
#   # check if the final message is fully repaired
#   synd = rs_calc_syndromes(msg_out, nsym, fcr, generator)

#   if max(synd) > 0:
#     raise newException(CatchableError, "Could not correct message")

#   # return the successfully decoded message
#   return (msg_out[0..nsym], msg_out[nsym..<msg_out.len]) # also return the corrected ecc block so that the user can check()

# proc rs_check(msg: seq[int], nsym: int, fcr = 0, generator = 2): bool =
#   ## Returns true if the message + ecc has no error of false
#   ## otherwise (may not always catch a wrong decoding or a wrong
#   ## message, particularly if there are too many errors
#   ## -- above the Singleton bound --, but it usually does)
#   return ( max(rs_calc_syndromes(msg, nsym, fcr, generator)) == 0 )

# when isMainModule:
#   import stew/byteutils

#   initTables()

#   # echo "GF_LOG ", gf_log
#   # echo "GF_EXP ", gf_exp
#   # echo "FIELD ", field_charac

#   let msg_in = @[
#     0x48, 0x65, 0x6C, 0x6C, 0x6F, 0x20, 0x52,
#     0x65, 0x65, 0x64, 0x2D, 0x53, 0x6F, 0x6C,
#     0x6F, 0x6D, 0x6F, 0x6E, 0x21 ]

#   var msg = rs_encode_msg(msg_in, 10)

#   var strMsg: string
#   for i in 0..<len(msg):
#     strMsg &= "0x" & uint8(msg[i]).toHex
#     strMsg &= " "

#   echo strMsg

#   msg[0] = 0
#   msg[10] = 0
#   msg[16] = 0
#   msg[14] = 0
#   msg[26] = 0

#   let (correct, code) = rs_correct_msg(msg, 10)

#   strMsg = ""
#   for i in 0..<len(correct):
#     strMsg &= "0x" & uint8(correct[i]).toHex
#     strMsg &= " "

#   echo "MESSAGE ", strMsg

