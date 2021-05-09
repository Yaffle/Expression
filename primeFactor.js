import nthRoot from './nthRoot.js';
//import bitLength from './bitLength.js';


// https://github.com/tc39/proposal-bigint/issues/205
// https://github.com/tc39/ecma262/issues/1729
// bitLength(a) = floor(log2(a)) + 1 if a > 0
function bitLength(a) {
  const s = a.toString(16);
  const c = s.charCodeAt(0) - '0'.charCodeAt(0);
  if (c <= 0) {
    throw new RangeError();
  }
  return (s.length - 1) * 4 + (32 - Math.clz32(Math.min(c, 8)));
}

//export default bitLength;

//function min(a, b) {
//  return a < b ? a : b;
//}

function modPow(base, exponent, modulus) {
  // exponent can be huge, use non-recursive variant
  let accumulator = 1n;
  while (exponent !== 0n) {
    let q = exponent >> 1n;
    if (exponent !== q + q) {
      accumulator = (accumulator * base) % modulus;
    }
    exponent = q;
    base = (base * base) % modulus;
  }
  return accumulator;
}

function range(start, end) {
  let i = start - 1;
  const tmp = {
    next: function () {
      i += 1;
      return {value: i > end ? undefined : i, done: i > end};
    }
  };
  tmp[globalThis.Symbol.iterator] = function () {
    return this;
  };
  return tmp;
}

// isPrime implementation is stolen from:
// https://github.com/peterolson/BigInteger.js
// https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Deterministic_variants
function isPrime(n) {
  const number = Number(n);
  if (number <= Number.MAX_SAFE_INTEGER) {
    return primeFactorUsingWheel(number) === number;
  }
  if (n < 2n) {
    throw new RangeError();
  }
  if (n === 2n) {
    return true;
  }
  if (n % 2n === 0n) {
    return false;
  }
  let r = 0;
  let d = n - 1n;
  while (d % 2n === 0n) {
    d /= 2n;
    r += 1;
  }
  // https://en.wikipedia.org/wiki/Miller–Rabin_primality_test#Testing_against_small_sets_of_bases
  const getRange = function (n) {
    const lnN = bitLength(n) * Math.log(2);
    const max = Math.min(Number(n - 2n), Math.floor(2 * lnN * Math.log(lnN)));
    return range(2, max);
  };
  let bases = null;
  if (n < 3317044064679887385961981n) {
    bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41].filter(x => x < n);
  } else {
    bases = getRange(n);
  }
  for (const a of bases) {
    const adn = modPow(BigInt(a), d, n);
    if (adn !== 1n) {
      for (let i = 0, x = adn; x !== n - 1n; i += 1, x = (x * x) % n) {
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
    var r1 = a % b;
    var r2 = b - r1;
    var r = r1 < r2 ? r1 : r2;
    a = b;
    b = r;
  }
  return a;
}

function f(x, c, mod) {
  return ((x * x) % mod + c) % mod;
}

// https://cp-algorithms.com/algebra/factorization.html#toc-tgt-9
function brent(n, x0 = 2n, c = 1n) {
    var x = x0;
    var g = 1n;
    var q = 1n;
    var xs, y;

    var m = 128;
    var l = 1;
    while (g === 1n) {
        y = x;
        for (var i = 1; i < l; i++)
            x = f(x, c, n);
        var k = 0;
        while (k < l && g === 1n) {
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
    if (g === n) {
        do {
            xs = f(xs, c, n);
            g = gcd(abs(xs - y), n);
        } while (g === 1n);
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

function remainder(n, i) {
  // `n % i` is slower in Chrome and Firefox ...
  return n - Math.floor(n / i) * i;
}

function primeFactorUsingWheel(n) {
  var steps = WHEEL3;
  var cycle = 3;
  var sn = Math.floor(Math.sqrt(n + 0.5));
  var i = 2;
  var s = 0;
  while (i <= sn) {
    if (remainder(n, i) === 0) {
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

function someFactor(n) {
  var x = Number(n);
  if (x < 1) {
    throw new TypeError("primeFactor cannot be called for numbers less than 2");
  }
  if (x === 1) {
    return n;
  }
  var s = gcd(BigInt(n), BigInt(304250263527210)); // a primorial - https://en.wikipedia.org/wiki/Primorial
  if (s > 1n) {
    //TODO: use-cases - ?
    return s === n ? BigInt(primeFactorUsingWheel(Number(s))) : s;
  }

  //! optimize n = f**2
  var squareRoot = nthRoot(n, 2);
  if (squareRoot**2n === n) {
    return squareRoot;
  }
  // https://en.wikipedia.org/wiki/Fermat%27s_factorization_method
  var a = squareRoot + 1n;
  var b2 = a*a - BigInt(n);
  var b = nthRoot(b2, 2);
  if (b*b === b2) {
    //console.debug("Fermat's method", n, a - b);
    return a - b;
  }
  //! optimize n = f**3
  var cubicRoot = nthRoot(n, 3);
  if (cubicRoot**3n === n) {
    return cubicRoot;
  }

  if (x <= Number.MAX_SAFE_INTEGER) {
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
  } while (g === n);
  return g;
}

  // https://en.wikipedia.org/wiki/Find_first_set#CTZ
  function countTrailingZeros(x, base) {
    //console.log(x, base);
    x = BigInt(x);
    base = BigInt(base);
    if (x < 0n || base < 0n) {
      throw new RangeError();
    }
    if (x === base) {
      return 1;
    }
    //if (base == 2 && typeof x === "bigint") { return bitLength(x & -x) - 1; } //TODO: ?
    if (x === 0n) {
      throw new TypeError();
    }
    var k = 1;
    while (x % base**BigInt(k) === 0n) {
      k *= 2;
    }
    var n = 0;
    for (var i = k / 2; i >= 1; i /= 2) {
      var v = base**BigInt(i);
      if (x % v === 0n) {
        n += i;
        x = x / v;
      }
    }
    return n;
  }

    primeFactor._rationalNumberToDecimalString = function (n, d, rounding) {
      // 1 / denominator
      const getPeriodOfRepeatingDecimalSegment = function (denominator, limit) {
        // https://softwareengineering.stackexchange.com/a/192081
        // https://en.wikipedia.org/wiki/Repeating_decimal#Other_properties_of_repetend_lengths
        // "If k = 2**a*5**b*n where n > 1 and n is not divisible by 2 or 5, then the length of the transient of 1/k is max(a, b), and the period equals r, where r is the smallest integer such that 10r ≡ 1 (mod n)."
        denominator = BigInt(denominator);
        if (denominator % 2n === 0n || denominator % 5n === 0n) {
          throw new RangeError("should be called with denominator not divible by 2 or 5");
        }
        const n = denominator;
        let period = 0;
        if (n !== 1n) {
          let z = 1n;
          do {
            period += 1;
            z = (10n * z) % n;
          } while (period <= limit && z !== 1n);
        }
        return period;
      };
      const bigDecimalToPlainString = function (significand, exponent, minFraction, minSignificant) {
        let e = exponent + significand.length - 1;
        significand = significand.replace(/0+$/g, '');
        const zeros = Math.max(0, Math.max(e + 1, minSignificant) - significand.length);
        if (e <= -1) {
          significand = "0".repeat(0 - e) + significand;
          e = 0;
        }
        significand += "0".repeat(zeros);
        significand += "0".repeat(Math.max(minFraction - (significand.length - (e + 1)), 0));
        return significand.slice(0, e + 1) + (significand.length > e + 1 ? "." + significand.slice(e + 1) : "");
      };
      // Something like Number#toPrecision: when value is between 10**-6 and 10**p? - to fixed, otherwise - to exponential:
      const toPrecision = function (significand, exponent, minSignificant) {
        const e = exponent + significand.length - 1;
        if (e < -6 || e >= minSignificant) {
          return bigDecimalToPlainString(significand, -(significand.length - 1), 0, minSignificant) + 'e' + (e < 0 ? '-' : '') + Math.abs(e).toString();
        }
        return bigDecimalToPlainString(significand, exponent, 0, minSignificant);
      };
      const toFixed = function (significand, exponent, minFraction) {
        return bigDecimalToPlainString(significand, exponent, minFraction, 0);
      };
      const digitsToDecimalNumber = function (significand, exponent, rounding) {
        // significand * 10**exponent
        if (rounding.significantDigits != undefined) {
          return toPrecision(significand, exponent, rounding.significantDigits);
        }
        return toFixed(significand, exponent, rounding.fractionDigits);
      };
      n = BigInt(n);
      d = BigInt(d);
      let sign = +1;
      if (d < 0n) {
        d = -d;
        sign = -sign;
      }
      if (n < 0n) {
        n = -n;
        sign = -sign;
      }
      const floorOfLog10 = function (n, d) {
        //TODO: optimize - ?
        let guess = Math.floor((bitLength(n) - 1 - bitLength(d)) / Math.log2(10));
        while ((guess < 0 ? 10n**BigInt(-guess) * n : n) >= (guess > 0 ? 10n**BigInt(guess) * d : d)) {
          guess += 1;
        }
        return guess - 1;
      };
      const a = primeFactor._countTrailingZeros(d, 2);
      const b = primeFactor._countTrailingZeros(d, 5);
      if (a > 0 && n % 2n === 0n) {
        throw new RangeError("not implemented");
      }
      if (b > 0 && n % 5n === 0n) {
        throw new RangeError("not implemented");
      }
      const d1 = d / (2n**BigInt(a) * 5n**BigInt(b));
      if ((rounding.fractionDigits != undefined && Math.max(a, b) <= rounding.fractionDigits || rounding.significantDigits != undefined && Math.max(a, b) + (floorOfLog10(n, d) + 1) - primeFactor._countTrailingZeros(n, 10) <= rounding.significantDigits) && n % d1 === 0n) { // exact result
        const scaling = Math.max(a, b);
        const result = ((10n**BigInt(scaling) * n) / d).toString(); //TODO: optimize - ?
        const minRounding = rounding.fractionDigits != undefined ? {fractionDigits: Math.max(a, b)} : {significantDigits: Math.max(a, b) + (floorOfLog10(n, d) + 1) <= rounding.significantDigits ? Math.max(a, b) + (floorOfLog10(n, d) + 1) : Math.max(a, b) + (floorOfLog10(n, d) + 1) - primeFactor._countTrailingZeros(n, 10)};
        const f = (sign < 0 ? '-' : '') + digitsToDecimalNumber(result, -scaling, minRounding);
        return f;
      } else {
        const scaling = rounding.fractionDigits != undefined ? rounding.fractionDigits : rounding.significantDigits - (floorOfLog10(n, d) + 1);
        const sn = scaling > 0 ? 10n**BigInt(scaling) * n : n;
        const sd = scaling < 0 ? 10n**BigInt(-scaling) * d : d;
        const result = ((sn + sd / 2n) / sd).toString();
        let f = (sign < 0 ? '-' : '') + digitsToDecimalNumber(result, -scaling, rounding);
        const lengthOfTransient = Math.max(a, b);
        const d1 = d / (2n**BigInt(a) * 5n**BigInt(b));
        const period = getPeriodOfRepeatingDecimalSegment(d1, f.length);
        ///^0\.(\d+?)(\d*?)(?:\1\2)*\1$/.exec('0.123123')
        if (period !== 0 && period <= f.length) { // a repeating decimal, the result is not exact
          const j = f.indexOf('.'); //?
          let offset = j + 1 + lengthOfTransient - (f.indexOf('e') !== -1 ? -Number(f.slice(f.indexOf('e') + 1)) - 0 : 0) - primeFactor._countTrailingZeros(n, 10);
          if (offset < j + 1) {
            //TODO: fix
            offset = j + 1;//!?
          }
          const lastFractionDigit = f.indexOf('e') !== -1 ? f.indexOf('e') : f.length;
          if (j !== -1 && (offset + period < lastFractionDigit || offset + period === lastFractionDigit && f.charCodeAt(offset) < '5'.charCodeAt(0))) {
            f = f.slice(0, offset) + '(' + f.slice(offset, offset + period) + ')' + f.slice(offset + period);
          }
        }
        if (!/[^0\.]/.test(f) && sign >= 0 && n !== 0n) {
          f = '+' + f;
        }
        return f;
      }
    };

function primeFactor(n) {
  n = BigInt(n);
  var factor = someFactor(n);
  if (factor === n) {
    return factor;
  }
  //var otherFactor = n / factor**BigInt(countTrailingZeros(n, factor));
  //if (otherFactor === 1n) {
  //  return primeFactor(factor);
  //}
  //TODO: divide by gcd (?)
  //if (otherFactor < factor) {
  //  var tmp = factor;
  //  factor = otherFactor;
  //  otherFactor = tmp;
  //}
  var a = primeFactor(factor);
  //var b = a > nthRoot(nthRoot(otherFactor, 2), 2) || a > 1e9 ? primeFactor(otherFactor) : primeFactorUsingWheelBig(otherFactor, a);
  //return min(a, b);
  return a;
}

// from https://github.com/juanelas/bigint-mod-arith/blob/master/lib/index.browser.mod.js :
/**
 * @typedef {Object} egcdReturn A triple (g, x, y), such that ax + by = g = gcd(a, b).
 * @property {bigint} g
 * @property {bigint} x
 * @property {bigint} y
 */
/**
 * An iterative implementation of the extended euclidean algorithm or extended greatest common divisor algorithm.
 * Take positive integers a, b as input, and return a triple (g, x, y), such that ax + by = g = gcd(a, b).
 *
 * @param {number|bigint} a
 * @param {number|bigint} b
 *
 * @returns {egcdReturn} A triple (g, x, y), such that ax + by = g = gcd(a, b).
 */
function eGcd (a, b) {
  if (Number(a) <= Number.MAX_SAFE_INTEGER && Number(b) <= Number.MAX_SAFE_INTEGER) {
    return eGcdSmall(a, b);
  }
  a = BigInt(a)
  b = BigInt(b)
  if (a <= 0n || b <= 0n) throw new RangeError('a and b MUST be > 0') // a and b MUST be positive

  let x = 0n
  let y = 1n
  let u = 1n
  let v = 0n

  while (a !== 0n) {
    const q = b / a
    const r = b % a
    const m = x - (u * q)
    const n = y - (v * q)
    b = a
    a = r
    x = u
    y = v
    u = m
    v = n
  }
  return {
    g: b,
    x: x,
    y: y
  }
}

function eGcdSmall (a, b) {
  a = Number(a)
  b = Number(b)
  if (a <= 0 || b <= 0) throw new RangeError('a and b MUST be > 0') // a and b MUST be positive

  let x = 0
  let y = 1
  let u = 1
  let v = 0

  while (a !== 0) {
    const q = Math.floor(b / a)
    const r = b % a
    const m = x - (u * q)
    const n = y - (v * q)
    b = a
    a = r
    x = u
    y = v
    u = m
    v = n
  }
  return {
    g: b,
    x: x,
    y: y
  }
}

function modInverse(a, m) {
  a = BigInt(a);
  m = BigInt(m);
  console.assert(a >= 0n);
  console.assert(m > 0n);
  if (a > m) {
    a = a % m;
  }
  let inv = BigInt(eGcd(a, m).x);
  inv = inv < 0n ? BigInt(inv) + m : inv;
  console.assert(inv >= 0n && inv < m);
  return inv;
}

function nextPrime(n) {
  if (n < Number.MAX_SAFE_INTEGER) {
    n = Number(n);
    n += 1;
    while (primeFactorUsingWheel(n) !== n) {
      n += 1;
    }
    return n;
  }
  n = BigInt(n);
  if (n < 2n) {
    return 2n;
  }
  n += 1n;
  if (n % 2n === 0n) {
    n += 1n;
  }
  while (!isPrime(n)) {
    n += 2n;
  }
  return n;
}

primeFactor._bitLength = bitLength;
primeFactor._isPrime = isPrime;
primeFactor._countTrailingZeros = countTrailingZeros;
primeFactor._someFactor = someFactor;
primeFactor._modInverse = modInverse;
primeFactor._nextPrime = nextPrime;
export default primeFactor;
