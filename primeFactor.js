import nthRoot from './nthRoot.js';

function min(a, b) {
  return a < b ? a : b;
}

function modPow(base, exponent, modulus) {
  // exponent can be huge, use non-recursive variant
  var accumulator = 1n;
  while (exponent !== 0n) {
    if (exponent % 2n === 0n) {
      exponent /= 2n;
      base = (base * base) % modulus;
    } else {
      exponent -= 1n;
      accumulator = (accumulator * base) % modulus;
    }
  }
  return accumulator;
}

function naturalLogarithm(n) {
  var number = Number(n);
  if (number < 1 / 0) {
    return Math.log(number);
  }
  // https://github.com/tc39/proposal-bigint/issues/205
  var s = n.toString(16);
  var p = Math.floor(Math.log((Number.MAX_SAFE_INTEGER + 1) / 32 + 0.5) / Math.log(2));
  var l = Math.floor(p / 4);
  return Math.log(Number('0x' + s.slice(0, l)) / Math.pow(2, 4 * l)) + 4 * Math.log(2) * s.length;
}

// isPrime implementation is stolen from:
// https://github.com/peterolson/BigInteger.js
// https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Deterministic_variants
function isPrime(n) {
  if (n < 2n) {
    throw new RangeError();
  }
  if (n === 2n) {
    return true;
  }
  if (n % 2n === 0n) {
    return false;
  }
  var r = 0;
  var d = n - 1n;
  while (d % 2n === 0n) {
    d /= 2n;
    r += 1;
  }
  for (var a = 2n, to = min(n - 2n, BigInt(Math.floor(2 * Math.pow(naturalLogarithm(n), 2)))); a <= to; a += 1n) {
    var adn = modPow(a, d, n);
    if (adn !== 1n) {
      for (var i = 0, x = adn; x !== n - 1n; i += 1, x = (x * x) % n) {
        if (i === r - 1) {
          return false;
        }
      }
    }
  }
  return true;
}

function abs(a) {
  return a < 0n ? -a : a;
}

function gcd(a, b) {
  while (b !== 0n) {
    let r1 = a % b;
    let r2 = b - r1;
    let r = r1 < r2 ? r1 : r2;
    a = b;
    b = r;
  }
  return a;
}

function f(x, c, mod) {
  //return ((x * x) % mod + c) % mod;
  return (x * x + c) % mod;
}

// https://cp-algorithms.com/algebra/factorization.html
function brent(n, x0 = 2n, c = 1n) {
    var x = x0;
    var g = 1;
    var q = 1n;
    var xs, y;

    var m = 128;
    var l = 1;
    while (g == 1) {
        y = x;
        for (var i = 1; i < l; i++)
            x = f(x, c, n);
        var k = 0;
        while (k < l && g == 1) {
            xs = x;
            for (var i = 0; i < m && i < l - k; i++) {
                x = f(x, c, n);
                q = (q * abs(y - x)) % n;
            }
            g = gcd(q, n);
            k += m;
        }
        l *= 2;
    }
    if (g == n) {
        do {
            xs = f(xs, c, n);
            g = gcd(abs(xs - y), n);
        } while (g == 1);
    }
    return g;
}

// Pollard's rho implementation is stolen from:
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
  var i = 2n;
  var s = 0;
  while (i <= max) {
    if (n % i === 0n) {
      return i;
    }
    i += BigInt(steps[s]);
    s += 1;
    if (s === steps.length) {
      s = cycle;
    }
  }
  return n;
}

function primeFactor(n) {
  var x = Number(n);
  if (x < 1) {
    throw new TypeError("primeFactor of a negative integer");
  }

  //! optimize n = f**2
  var squareRoot = nthRoot(n, 2);
  if (squareRoot**2n === n) {
    return primeFactor(squareRoot);
  }
  //! optimize n = f**3
  var cubicRoot = nthRoot(n, 3);
  if (cubicRoot**3n === n) {
    return primeFactor(cubicRoot);
  }

  var s = Number(gcd(n, BigInt(304250263527210))); // a primorial - https://en.wikipedia.org/wiki/Primorial
  if (s !== 1) {
    //TODO: use-cases - ?
    return BigInt(primeFactorUsingWheel(s));
  }
  if (x <= 9007199254740991) {
    return BigInt(primeFactorUsingWheel(x));
  }
  if (isPrime(n)) {
    return n;
  }
  var x0 = 2 - 1;
  var g = n;
  do {
    x0 += 1;
    g = brent(n, BigInt(x0));
  } while (g == n);
  var factor = g;
  if (n / factor < factor) {
    factor = n / factor;
  }
  var a = primeFactor(factor);
  var b = a > nthRoot(nthRoot(n / factor, 2), 2) ? primeFactor(n / factor) : primeFactorUsingWheelBig(n / factor, a);
  return min(a, b);
}

primeFactor._isPrime = isPrime;

export default primeFactor;
