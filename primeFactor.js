import BigInteger from './BigInteger.js';
import nthRoot from './nthRoot.js';

function min(a, b) {
  return BigInteger.lessThan(a, b) ? a : b;
}

function modPow(base, exponent, modulus, accumulator) {
  return !BigInteger.lessThan(BigInteger.BigInt(0), exponent) ? accumulator : modPow(BigInteger.remainder(BigInteger.multiply(base, base), modulus), BigInteger.divide(exponent, BigInteger.BigInt(2)), modulus, !BigInteger.lessThan(BigInteger.remainder(exponent, BigInteger.BigInt(2)), BigInteger.BigInt(1)) ? BigInteger.remainder(BigInteger.multiply(accumulator, base), modulus) : accumulator);
}

function bitLength(n) {
  var x = BigInteger.toNumber(n);
  return x < 1 / 0 ? Math.floor(Math.log(x) / Math.log(2)) : 1024 + bitLength(BigInteger.divide(n, BigInteger.exponentiate(2, 1024)));
}

function log(n) {
  //n = f * 2**e
  //Math.log(n) = Math.log(f) + Math.log(2) * e;
  var x = BigInteger.toNumber(n);
  return x < 1 / 0 ? Math.log(x) : Math.log(BigInteger.toNumber(BigInteger.divide(n, BigInteger.exponentiate(2, bitLength(n) - 53))) + Math.log(2) * (bitLength(n) - 53));
}

// isPrime implementation is stolen from:
// https://github.com/peterolson/BigInteger.js
// https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Deterministic_variants
function isPrime(n) {
  if (BigInteger.lessThan(n, BigInteger.BigInt(2))) {
    throw new RangeError();
  }
  if (!BigInteger.lessThan(BigInteger.BigInt(2), n)) {
    return true;
  }
  if (BigInteger.lessThan(BigInteger.remainder(n, BigInteger.BigInt(2)), BigInteger.BigInt(1))) {
    return false;
  }
  var s = 0;
  var d = BigInteger.subtract(n, BigInteger.BigInt(1));
  while (BigInteger.lessThan(BigInteger.remainder(d, BigInteger.BigInt(2)), BigInteger.BigInt(1))) {
    d = BigInteger.divide(d, BigInteger.BigInt(2));
    s += 1;
  }
  for (var a = min(BigInteger.subtract(n, BigInteger.BigInt(1)), BigInteger.BigInt(Math.floor(2 * Math.pow(log(n), 2)))); !BigInteger.lessThan(a, BigInteger.BigInt(2)); a = BigInteger.subtract(a, BigInteger.BigInt(1))) {
    var adn = modPow(a, d, n, BigInteger.BigInt(1));
    if (BigInteger.lessThan(adn, BigInteger.BigInt(1)) || BigInteger.lessThan(BigInteger.BigInt(1), adn)) {
      for (var r = 0, x = adn; BigInteger.lessThan(x, BigInteger.subtract(n, BigInteger.BigInt(1))); r += 1, x = BigInteger.remainder(BigInteger.multiply(x, x), n)) {
        if (r === s - 1) {
          return false;
        }
      }
    }
  }
  return true;
}

function abs(a) {
  return BigInteger.lessThan(a, BigInteger.BigInt(0)) ? BigInteger.unaryMinus(a) : a;
}

function gcd(a, b) {
  while (BigInteger.lessThan(BigInteger.BigInt(0), b)) {
    var t = BigInteger.remainder(a, b);
    a = b;
    b = t;
  }
  return a;
}

// Pollard's rho implementation is stolen from:
// https://github.com/jiggzson/nerdamer/blob/master/nerdamer.core.js
// https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#C_code_sample
function primeFactorByPollardRho(n) {
  if (BigInteger.lessThan(n, BigInteger.BigInt(2))) {
    throw new RangeError();
  }
  if (isPrime(n)) {
    return n;
  }
  var factor = n;
  for (var x0 = BigInteger.BigInt(2); !BigInteger.lessThan(factor, n); x0 = BigInteger.add(x0, BigInteger.BigInt(1))) {
    if (BigInteger.lessThan(BigInteger.remainder(n, x0), BigInteger.BigInt(1))) {
      //?
      return x0;
    }
    var xFixed = x0;
    var cycleSize = 2;
    var x = x0;
    factor = BigInteger.BigInt(1);
    while (BigInteger.lessThan(factor, BigInteger.BigInt(2))) {
      var test = BigInteger.BigInt(1);
      var testStart = x;
      var found = false;
      for (var count = 1; count <= cycleSize && BigInteger.lessThan(factor, BigInteger.BigInt(2)); count += 1) {
        x = BigInteger.remainder(BigInteger.add(BigInteger.multiply(x, x), BigInteger.BigInt(1)), n);
        //factor = gcd(abs(x - xFixed), n);
        test = BigInteger.remainder(BigInteger.multiply(test, abs(BigInteger.subtract(x, xFixed))), n);
        if (found || count === cycleSize || count % 16 === 0) {
          // https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Variants
          factor = gcd(test, n);
          if (!found && !BigInteger.lessThan(factor, BigInteger.BigInt(2))) {
            cycleSize *= 2;
            factor = BigInteger.BigInt(1);
            x = testStart;
            found = true;
          }
          test = BigInteger.BigInt(1);
          testStart = x;
        }
      }
      cycleSize *= 2;
      xFixed = x;
    }
  }
  var a = primeFactor(factor);
  var b = primeFactor(BigInteger.divide(n, factor));
  return BigInteger.lessThan(a, b) ? a : b;
}

var WHEEL3 = [1, 2, 2, 4, 2, 4, 2, 4, 6, 2, 6];

function primeFactorUsingWheel(n) {
  var steps = WHEEL3;
  var cycle = 3;
  var sn = Math.floor(Math.sqrt(n + 0.5));
  var i = 2;
  var s = 0;
  while (i <= sn) {
    if (n % i === 0) {
      return i;
    }
    i += steps[s];
    s += 1;
    if (s === steps.length) {
      s = cycle;
    }
  }
  return n;
}

function primeFactor(n) {
  var x = BigInteger.toNumber(n);
  if (x < 1) {
    throw new TypeError("primeFactor of a negative integer");
  }

  //! optimize n = f**2
  var f = nthRoot(n, 2);
  if (BigInteger.equal(BigInteger.exponentiate(f, BigInteger.BigInt(2)), n)) {
    return primeFactor(f);
  }
  //! optimize n = f**3
  var f = nthRoot(n, 3);
  if (BigInteger.equal(BigInteger.exponentiate(f, BigInteger.BigInt(3)), n)) {
    return primeFactor(f);
  }

  if (x <= 9007199254740991) {
    return BigInteger.BigInt(primeFactorUsingWheel(x));
  }
  var s = BigInteger.toNumber(gcd(n, BigInteger.BigInt(304250263527210))); // a primorial
  if (s !== 1) {
    //TODO: use-cases - ?
    return BigInteger.BigInt(primeFactorUsingWheel(s));
  }
  return primeFactorByPollardRho(n);
}

export default primeFactor;
