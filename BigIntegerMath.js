/*global BigInteger*/

(function () {
"use strict";

// floor(log(n) / log(b)), n >= 1, base >= 2, https://programmingpraxis.com/contents/standard-prelude/
function ilogb(n, base) {
  var i = 0;
  while (BigInteger.compareTo(n, BigInteger.pow(base, BigInteger.pow(2, i))) >= 0) {
    i = BigInteger.add(i, 1);
  }
  var e = 0;
  var t = 1;
  while (BigInteger.compareTo(i, 0) >= 0) {
    var b = BigInteger.pow(2, i);
    if (BigInteger.compareTo(n, BigInteger.multiply(t, BigInteger.pow(base, b))) >= 0) {
      t = BigInteger.multiply(t, BigInteger.pow(base, b));
      e = BigInteger.add(e, b);
    }
    i = BigInteger.subtract(i, 1);
  }
  return e;
}

function ilog2(n) {
  return ilogb(n, 2);
}

// https://stackoverflow.com/a/15979957/839199y
function nthRoot(S, n) {
  if (n === 2) {
    var t = BigInteger.toNumber(S);
    if (t < 4503599627370496) { // 2**52
      return BigInteger.fromNumber(Math.floor(Math.sqrt(t + 0.5)));
    }
    if (t < 4503599627370496 * 4503599627370496) { // 2**104
      var y = BigInteger.fromNumber(Math.floor(Math.sqrt(t) + 0.5));
      if (BigInteger.compareTo(BigInteger.multiply(y, y), S) > 0) {
        y = BigInteger.subtract(y, 1);
      }
      return y;
    }
  }
  var e = ilog2(S);
  if (e < n) {
    return 1;
  }
  var f = BigInteger.divide(BigInteger.add(e, n), BigInteger.multiply(2, n));
  var x = BigInteger.multiply(BigInteger.add(nthRoot(BigInteger.divide(S, BigInteger.pow(2, f * n)), n), 1), BigInteger.pow(2, f));
  var xprev = BigInteger.add(x, 1);
  while (BigInteger.compareTo(xprev, x) > 0) {
    xprev = x;
    x = BigInteger.divide(BigInteger.add(BigInteger.multiply(x, BigInteger.subtract(n, 1)), BigInteger.divide(S, BigInteger.pow(x, n - 1))), n);
  }
  return xprev;
}

function powBig(x, count, accumulator) {
  accumulator = accumulator == undefined ? 1 : accumulator;
  if (count < 0) {
    throw new RangeError();
  }
  if (count > 9007199254740991) {
    throw new RangeError();
  }
  return (count < 1 ? accumulator : (2 * Math.floor(count / 2) !== count ? powBig(x, count - 1, BigInteger.multiply(accumulator, x)) : powBig(BigInteger.multiply(x, x), Math.floor(count / 2), accumulator)));
}

function pow(x, count) {
  if (typeof x === "number" && count >= 0 && count < 53) {
    var value = 0 + Math.pow(x, count);
    if (value >= -9007199254740991 && value <= 9007199254740991) {
      return value;
    }
  }
  return powBig(x, count, 1);
}

var WHEEL2 = [1, 2, 2, 4];
function primeFactorUsingWheel2(n) {
  var i = 2;
  var s = 0;
  var r = Math.floor(Math.sqrt(n + 0.5));
  while (i <= r) {
    if (n % i === 0) {
      return i;
    }
    i += WHEEL2[s];
    s += 1;
    if (s === WHEEL2.length) {
      s = 2;
    }
  }
  return n;
}

function gcd(a, b) {
  while (BigInteger.compareTo(b, 0) !== 0) {
    var t = BigInteger.remainder(a, b);
    a = b;
    b = t;
  }
  return a;
}

function modPow(a, n, m, accumulator) {
  accumulator = accumulator == undefined ? 1 : accumulator;
  return BigInteger.compareTo(n, 0) === 0 ? accumulator : modPow(BigInteger.remainder(BigInteger.multiply(a, a), m), BigInteger.divide(n, 2), m, BigInteger.compareTo(BigInteger.remainder(n, 2), 1) === 0 ? BigInteger.remainder(BigInteger.multiply(accumulator, a), m) : accumulator);
}

function min(a, b) {
  return BigInteger.compareTo(a, b) < 0 ? a : b;
}

function isPrime(n) {
  // isPrime implementation is from https://github.com/peterolson/BigInteger.js
  // https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Deterministic_variants
  if (BigInteger.compareTo(n, 2) < 0) {
    throw new RangeError();
  }
  if (BigInteger.compareTo(n, 2) === 0) {
    return true;
  }
  if (BigInteger.compareTo(BigInteger.remainder(n, 2), 0) === 0) {
    return false;
  }
  var s = 0;
  var d = BigInteger.subtract(n, 1);
  while (BigInteger.compareTo(BigInteger.remainder(d, 2), 0) === 0) {
    d = BigInteger.divide(d, 2);
    s += 1;
  }
  function test(x0, n, s) {
    for (var r = 0, x = x0; r <= s - 1; r += 1, x = BigInteger.remainder(BigInteger.multiply(x, x), n)) {
      if (BigInteger.compareTo(x, BigInteger.subtract(n, 1)) === 0) {
        return false;
      }
    }
    return true;
  }
  for (var a = min(BigInteger.subtract(n, 1), BigInteger.pow(BigInteger.add(ilog2(n), 1), 2)); BigInteger.compareTo(a, 2) >= 0; a = BigInteger.subtract(a, 1)) {
    var x = modPow(a, d, n);
    if (BigInteger.compareTo(x, 1) !== 0 && test(x, n, s)) {
      return false;
    }
  }
  return true;
}

function abs(a) {
  return BigInteger.compareTo(a, 0) < 0 ? BigInteger.negate(a) : a;
}

// Pollard's rho stolen from https://github.com/jiggzson/nerdamer/blob/master/nerdamer.core.js
function ifactor(n) {
  if (BigInteger.compareTo(n, 2) < 0) {
    throw new RangeError();
  }
  if (isPrime(n)) {
    return n;
  }
  // https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#C++_code_sample
  var factor = n;
  for (var x0 = 2; BigInteger.compareTo(factor, n) === 0; x0 = BigInteger.add(x0, 1)) {
    if (BigInteger.compareTo(BigInteger.remainder(n, x0), 0) === 0) {
      //?
      return x0;
    }
    var xFixed = x0;
    var cycleSize = 2;
    var x = x0;
    factor = 1;
    while (BigInteger.compareTo(factor, 1) === 0) {
      var test = 1;
      var testStart = x;
      var found = false;
      for (var count = 1; count <= cycleSize && BigInteger.compareTo(factor, 1) === 0; count += 1) {
        x = BigInteger.remainder(BigInteger.add(BigInteger.multiply(x, x), 1), n);
        test = BigInteger.remainder(BigInteger.multiply(test, abs(BigInteger.subtract(x, xFixed))), n);
        if (found || count === cycleSize || count % 16 === 0) {
          // https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Variants
          factor = gcd(test, n);
          if (!found && BigInteger.compareTo(factor, 1) !== 0) {
            cycleSize *= 2;
            factor = 1;
            x = testStart;
            found = true;
          }
          test = 1;
          testStart = x;
        }
      }
      cycleSize *= 2;
      xFixed = x;
    }
  }
  var a = ifactor(factor);
  var b = ifactor(BigInteger.divide(n, factor));
  return min(a, b);
}

function primeFactor(n) {
  if (BigInteger.compareTo(n, BigInteger.fromNumber(9007199254740991)) <= 0) {
    return BigInteger.fromNumber(primeFactorUsingWheel2(n));
  }
  return ifactor(n);
}

BigInteger.nthRoot = nthRoot;
BigInteger.pow = pow;
BigInteger.gcd = gcd;
BigInteger.primeFactor = primeFactor;

}());
