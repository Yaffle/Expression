import BigInteger from './BigInteger.js';
import nthRoot from './nthRoot.compiled.js';

function min(a, b) {
  return BigInteger.lessThan(a, b) ? a : b;
}

function modPow(base, exponent, modulus) {
  // exponent can be huge, use non-recursive variant
  var accumulator = BigInteger.BigInt(1);

  while (BigInteger.notEqual(exponent, BigInteger.BigInt(0))) {
    if (BigInteger.equal(BigInteger.remainder(exponent, BigInteger.BigInt(2)), BigInteger.BigInt(0))) {
      exponent = BigInteger.divide(exponent, BigInteger.BigInt(2));
      base = BigInteger.remainder(BigInteger.multiply(base, base), modulus);
    } else {
      exponent = BigInteger.subtract(exponent, BigInteger.BigInt(1));
      accumulator = BigInteger.remainder(BigInteger.multiply(accumulator, base), modulus);
    }
  }

  return accumulator;
}

function naturalLogarithm(n) {
  var number = BigInteger.toNumber(n);

  if (number < 1 / 0) {
    return Math.log(number);
  } // https://github.com/tc39/proposal-bigint/issues/205


  var s = n.toString(16);
  var p = Math.floor(Math.log((Number.MAX_SAFE_INTEGER + 1) / 32 + 0.5) / Math.log(2));
  var l = Math.floor(p / 4);
  return Math.log(Number('0x' + s.slice(0, l)) / Math.pow(2, 4 * l)) + 4 * Math.log(2) * s.length;
} // isPrime implementation is stolen from:
// https://github.com/peterolson/BigInteger.js
// https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Deterministic_variants


function isPrime(n) {
  if (BigInteger.lessThan(n, BigInteger.BigInt(2))) {
    throw new RangeError();
  }

  if (BigInteger.equal(n, BigInteger.BigInt(2))) {
    return true;
  }

  if (BigInteger.equal(BigInteger.remainder(n, BigInteger.BigInt(2)), BigInteger.BigInt(0))) {
    return false;
  }

  var r = 0;
  var d = BigInteger.subtract(n, BigInteger.BigInt(1));

  while (BigInteger.equal(BigInteger.remainder(d, BigInteger.BigInt(2)), BigInteger.BigInt(0))) {
    d = BigInteger.divide(d, BigInteger.BigInt(2));
    r += 1;
  }

  for (var a = BigInteger.BigInt(2), to = min(BigInteger.subtract(n, BigInteger.BigInt(2)), BigInteger.BigInt(Math.floor(2 * Math.pow(naturalLogarithm(n), 2)))); BigInteger.lessThanOrEqual(a, to); a = BigInteger.add(a, BigInteger.BigInt(1))) {
    var adn = modPow(a, d, n);

    if (BigInteger.notEqual(adn, BigInteger.BigInt(1))) {
      for (var i = 0, x = adn; BigInteger.notEqual(x, BigInteger.subtract(n, BigInteger.BigInt(1))); i += 1, x = BigInteger.remainder(BigInteger.multiply(x, x), n)) {
        if (i === r - 1) {
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
  while (BigInteger.notEqual(b, BigInteger.BigInt(0))) {
    let r1 = BigInteger.remainder(a, b);
    let r2 = BigInteger.subtract(b, r1);
    let r = BigInteger.lessThan(r1, r2) ? r1 : r2;
    a = b;
    b = r;
  }

  return a;
}

function f(x, c, mod) {
  //return ((x * x) % mod + c) % mod;
  return BigInteger.remainder(BigInteger.add(BigInteger.multiply(x, x), c), mod);
} // https://cp-algorithms.com/algebra/factorization.html


function brent(n, x0 = BigInteger.BigInt(2), c = BigInteger.BigInt(1)) {
  var x = x0;
  var g = 1;
  var q = BigInteger.BigInt(1);
  var xs, y;
  var m = 128;
  var l = 1;

  while (g == 1) {
    y = x;

    for (var i = 1; i < l; _x = i, i = BigInteger.add(i, BigInteger.BigInt(1)), _x) {
      var _x;

      x = f(x, c, n);
    }

    var k = 0;

    while (k < l && g == 1) {
      xs = x;

      for (var i = 0; i < m && i < l - k; _x2 = i, i = BigInteger.add(i, BigInteger.BigInt(1)), _x2) {
        var _x2;

        x = f(x, c, n);
        q = BigInteger.remainder(BigInteger.multiply(q, abs(BigInteger.subtract(y, x))), n);
      }

      g = gcd(q, n);
      k += m;
    }

    l *= 2;
  }

  if (g == n) {
    do {
      xs = f(xs, c, n);
      g = gcd(abs(BigInteger.subtract(xs, y)), n);
    } while (g == 1);
  }

  return g;
} // Pollard's rho implementation is stolen from:
// https://github.com/jiggzson/nerdamer/blob/master/nerdamer.core.js
// https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#C_code_sample

/*
function factorByPollardRho(n, x0 = 2n, c = 1n) {
    var factor = n;
    if (n % x0 === 0n) {//?
      return x0;
    }
    var xFixed = x0;
    var cycleSize = 2;
    var x = x0;
    factor = 1n;
    while (factor === 1n) {
      var test = 1n;
      var testStart = x;
      var found = false;
      for (var count = 1; count <= cycleSize && factor === 1n; count += 1) {
        x = (x * x + c) % n;
        //factor = gcd(abs(x - xFixed), n);
        test = (test * abs(x - xFixed)) % n;
        if (found || count === cycleSize || count % 16 === 0) {
          // https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Variants
          factor = gcd(test, n);
          if (!found && factor !== 1n) {
            cycleSize *= 2;
            factor = 1n;
            x = testStart;
            found = true;
          }
          test = 1n;
          testStart = x;
        }
      }
      cycleSize *= 2;
      xFixed = x;
    }
    return factor;
}
*/


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

function primeFactorUsingWheelBig(n, max) {
  var steps = WHEEL3;
  var cycle = 3;
  var i = BigInteger.BigInt(2);
  var s = 0;

  while (BigInteger.lessThanOrEqual(i, max)) {
    if (BigInteger.equal(BigInteger.remainder(n, i), BigInteger.BigInt(0))) {
      return i;
    }

    i = BigInteger.add(i, BigInteger.BigInt(steps[s]));
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
  } //! optimize n = f**2


  var squareRoot = nthRoot(n, 2);

  if (BigInteger.equal(BigInteger.exponentiate(squareRoot, BigInteger.BigInt(2)), n)) {
    return primeFactor(squareRoot);
  } //! optimize n = f**3


  var cubicRoot = nthRoot(n, 3);

  if (BigInteger.equal(BigInteger.exponentiate(cubicRoot, BigInteger.BigInt(3)), n)) {
    return primeFactor(cubicRoot);
  }

  var s = BigInteger.toNumber(gcd(n, BigInteger.BigInt(304250263527210))); // a primorial - https://en.wikipedia.org/wiki/Primorial

  if (s !== 1) {
    //TODO: use-cases - ?
    return BigInteger.BigInt(primeFactorUsingWheel(s));
  }

  if (x <= 9007199254740991) {
    return BigInteger.BigInt(primeFactorUsingWheel(x));
  }

  if (isPrime(n)) {
    return n;
  }

  var x0 = 2 - 1;
  var g = n;

  do {
    x0 += 1;
    g = brent(n, BigInteger.BigInt(x0));
  } while (g == n);

  var factor = g;

  if (BigInteger.lessThan(BigInteger.divide(n, factor), factor)) {
    factor = BigInteger.divide(n, factor);
  }

  var a = primeFactor(factor);
  var b = BigInteger.greaterThan(a, nthRoot(nthRoot(BigInteger.divide(n, factor), 2), 2)) ? primeFactor(BigInteger.divide(n, factor)) : primeFactorUsingWheelBig(BigInteger.divide(n, factor), a);
  return min(a, b);
}

primeFactor._isPrime = isPrime;
export default primeFactor;

