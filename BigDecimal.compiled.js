import BigInteger from './BigInteger.js';

/*jslint bigint: true, vars: true, indent: 2*/
// Usage:
// BigDecimal.BigDecimal(bigint)
// BigDecimal.BigDecimal(string)
// BigDecimal.BigDecimal(number) (only integers)
// BigDecimal.toBigInt(a) (not in the spec)
// BigDecimal.toNumber(a) (not in the spec, only integers)
// BigDecimal.unaryMinus(a)
// BigDecimal.add(a, b[, rounding])
// BigDecimal.subtract(a, b[, rounding])
// BigDecimal.multiply(a, b[, rounding])
// BigDecimal.divide(a, b, rounding)
// BigDecimal.lessThan(a, b)
// BigDecimal.greaterThan(a, b)
// BigDecimal.equal(a, b)
// BigDecimal.round(a, rounding)
// a.toString()
// Math: (not in the spec)
// BigDecimal.log(a, rounding)
// BigDecimal.exp(a, rounding)
// BigDecimal.sin(a, rounding)
// BigDecimal.cos(a, rounding)
// BigDecimal.atan(a, rounding)
function BigDecimal(significand, exponent) {
  this.significand = significand;
  this.exponent = exponent;
}

BigDecimal.BigDecimal = function (value) {
  if (value instanceof BigDecimal) {
    return value;
  }

  if (typeof value === "string") {
    //throw new TypeError();
    var match = /^\s*([+\-])?(\d+)?\.?(\d+)?(?:e([+\-]?\d+))?\s*$/.exec(value);

    if (match == null) {
      throw new RangeError();
    }

    return create(BigInteger.BigInt((match[1] || "") + (match[2] || "") + (match[3] || "")), BigInteger.subtract(BigInteger.BigInt(match[4] || "0"), BigInteger.BigInt((match[3] || "").length)));
  }

  var a = create(BigInteger.BigInt(value), BigInteger.BigInt(0));
  var b = normalize(a);
  while (a !== b) {
    a = b;
    b = normalize(a);
  }
  return a;
};

BigDecimal.toNumber = function (a) {
  return BigInteger.toNumber(BigDecimal.toBigInt(a));
};

BigDecimal.toBigInt = function (a) {
  if (BigInteger.lessThan(a.exponent, BigInteger.BigInt(0)) && BigInteger.notEqual(BigInteger.remainder(a.significand, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.unaryMinus(a.exponent))), BigInteger.BigInt(0))) {
    throw new RangeError("The bigdecimal " + a.toString() + " cannot be converted to a BigInt because it is not an integer");
  }

  return BigInteger.lessThan(a.exponent, BigInteger.BigInt(0)) ? BigInteger.divide(a.significand, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.unaryMinus(a.exponent))) : BigInteger.multiply(a.significand, BigInteger.exponentiate(BigInteger.BigInt(10), a.exponent));
};

function create(significand, exponent) {
  return new BigDecimal(significand, exponent);
}

function bigIntMax(a, b) {
  return BigInteger.lessThan(a, b) ? b : a;
}

function bigIntMin(a, b) {
  return BigInteger.lessThan(a, b) ? a : b;
}

function bigIntSign(a) {
  return BigInteger.lessThan(a, BigInteger.BigInt(0)) ? BigInteger.unaryMinus(BigInteger.BigInt(1)) : BigInteger.greaterThan(a, BigInteger.BigInt(0)) ? BigInteger.BigInt(1) : BigInteger.BigInt(0);
}

function bigIntAbs(a) {
  return BigInteger.lessThan(a, BigInteger.BigInt(0)) ? BigInteger.unaryMinus(a) : a;
}

function digits(a) {
  // floor(log10(abs(a.significand))) + 1
  function naturalLogarithm(n) {
    var number = BigInteger.toNumber(n);
    if (number < 1 / 0) {
      return Math.log(number);
    }
    // https://github.com/tc39/proposal-bigint/issues/205
    var s = n.toString(16);
    var p = Math.floor(Math.log((Number.MAX_SAFE_INTEGER + 1) / 32 + 0.5) / Math.log(2));
    var l = Math.floor(p / 4);
    return Math.log(Number('0x' + s.slice(0, l)) / Math.pow(2, 4 * l)) + 4 * Math.log(2) * s.length;
  }

  var n = bigIntMax(bigIntAbs(a.significand), BigInteger.BigInt(1));
  var number = BigInteger.toNumber(n);

  if (number < (9007199254740991 + 1) / 16) {
    return Math.floor(Math.log(number + 0.5) / Math.log(10)) + 1;
  }

  if (number < 1 / 0) {
    var e = Math.log(number) / Math.log(10);

    if (Math.floor(e * (1 + 2 / (9007199254740991 + 1))) < e) {
      return Math.floor(e) + 1;
    }

    var i = Math.floor(e + 0.5);
    return BigInteger.greaterThanOrEqual(n, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(i))) ? i + 1 : i;
  }

  var e = naturalLogarithm(n) / Math.log(10);

  if (Math.floor(e * (1 - 32 / (9007199254740991 + 1))) === Math.floor(e) && Math.floor(e * (1 + 32 / (9007199254740991 + 1))) === Math.floor(e)) {
    return Math.floor(e) + 1;
  }

  var i = Math.floor(e + 0.5);
  return BigInteger.greaterThanOrEqual(n, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(i))) ? i + 1 : i;
}

function normalize(a) {
  if (BigInteger.greaterThanOrEqual(bigIntAbs(a.significand), BigInteger.BigInt(67108864))) {
    // Math.sqrt((Number.MAX_SAFE_INTEGER + 1) / 2)
    var dividend = a.significand;
    var divisor = BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(15));
    var e = 15;

    if (BigInteger.equal(BigInteger.remainder(dividend, divisor), BigInteger.BigInt(0))) {
      while (BigInteger.equal(BigInteger.remainder(dividend, BigInteger.multiply(divisor, divisor)), BigInteger.BigInt(0))) {
        divisor = BigInteger.multiply(divisor, divisor);
        e *= 2;
      }

      var quotient = BigInteger.divide(dividend, divisor);
      return create(quotient, BigInteger.add(a.exponent, BigInteger.BigInt(e)));
    }
  }

  if (BigInteger.equal(a.significand, BigInteger.BigInt(0)) && BigInteger.notEqual(a.exponent, BigInteger.BigInt(0))) {
    return create(BigInteger.BigInt(0), BigInteger.BigInt(0));
  }

  return a;
}

function round(a, d, r, rounding) {
  if (rounding != null) {
    var k = BigInteger.BigInt(0);

    if (rounding.maximumSignificantDigits != null) {
      console.assert(rounding.maximumSignificantDigits > 0);
      k = bigIntMax(BigInteger.BigInt(digits(a) - rounding.maximumSignificantDigits), BigInteger.BigInt(0));
    }

    if (rounding.maximumFractionDigits != null) {
      console.assert(rounding.maximumFractionDigits >= 0);
      k = bigIntMax(BigInteger.subtract(BigInteger.BigInt(0), BigInteger.add(a.exponent, BigInteger.BigInt(rounding.maximumFractionDigits))), BigInteger.BigInt(0));
    }

    if (BigInteger.greaterThan(k, BigInteger.BigInt(0)) || BigInteger.notEqual(r, BigInteger.BigInt(0))) {
      var divisor = BigInteger.exponentiate(BigInteger.BigInt(10), k);
      var dividend = a.significand;
      var quotient = BigInteger.divide(dividend, divisor);
      var remainder = BigInteger.subtract(dividend, BigInteger.multiply(quotient, divisor));
      divisor = BigInteger.multiply(divisor, d);
      remainder = BigInteger.add(BigInteger.multiply(remainder, d), r);

      if (BigInteger.notEqual(remainder, BigInteger.BigInt(0))) {
        if (rounding.roundingMode === "floor") {
          if (BigInteger.lessThan(remainder, BigInteger.BigInt(0))) {
            quotient = BigInteger.subtract(quotient, BigInteger.BigInt(1));
          }
        } else if (rounding.roundingMode === "ceil") {
          if (BigInteger.greaterThan(remainder, BigInteger.BigInt(0))) {
            quotient = BigInteger.add(quotient, BigInteger.BigInt(1));
          }
        } else if (rounding.roundingMode === "half-up") {
          if (BigInteger.greaterThanOrEqual(BigInteger.multiply(remainder, BigInteger.BigInt(2)), divisor)) {
            quotient = BigInteger.add(quotient, BigInteger.BigInt(1));
          }

          if (BigInteger.greaterThanOrEqual(BigInteger.multiply(BigInteger.unaryMinus(remainder), BigInteger.BigInt(2)), divisor)) {
            quotient = BigInteger.subtract(quotient, BigInteger.BigInt(1));
          }
        } else if (rounding.roundingMode === "half-down") {
          if (BigInteger.greaterThan(BigInteger.multiply(remainder, BigInteger.BigInt(2)), divisor)) {
            quotient = BigInteger.add(quotient, BigInteger.BigInt(1));
          }

          if (BigInteger.greaterThan(BigInteger.multiply(BigInteger.unaryMinus(remainder), BigInteger.BigInt(2)), divisor)) {
            quotient = BigInteger.subtract(quotient, BigInteger.BigInt(1));
          }
        } else if (rounding.roundingMode === "half-even") {
          if (BigInteger.greaterThan(BigInteger.multiply(remainder, BigInteger.BigInt(2)), divisor) || BigInteger.equal(BigInteger.multiply(remainder, BigInteger.BigInt(2)), divisor) && BigInteger.notEqual(BigInteger.remainder(quotient, BigInteger.BigInt(2)), BigInteger.BigInt(0))) {
            quotient = BigInteger.add(quotient, BigInteger.BigInt(1));
          }

          if (BigInteger.greaterThan(BigInteger.multiply(BigInteger.unaryMinus(remainder), BigInteger.BigInt(2)), divisor) || BigInteger.equal(BigInteger.multiply(BigInteger.unaryMinus(remainder), BigInteger.BigInt(2)), divisor) && BigInteger.notEqual(BigInteger.remainder(quotient, BigInteger.BigInt(2)), BigInteger.BigInt(0))) {
            quotient = BigInteger.subtract(quotient, BigInteger.BigInt(1));
          }
        } else {
          throw new RangeError("supported roundingMode (floor/ceil/half-even/half-up/half-down) is not given");
        }
      }

      return create(quotient, BigInteger.add(a.exponent, k));
    }
  }

  if (BigInteger.notEqual(r, BigInteger.BigInt(0))) {
    throw new RangeError("rounding is not given for inexact operation");
  }

  return a;
}

BigDecimal.unaryMinus = function (a) {
  return create(BigInteger.unaryMinus(a.significand), a.exponent);
};

BigDecimal.add = function (a, b, rounding = null) {
  if (BigInteger.equal(b.significand, BigInteger.BigInt(0))) {
    return round(a, BigInteger.BigInt(1), BigInteger.BigInt(0), rounding);
  }

  if (BigInteger.equal(a.significand, BigInteger.BigInt(0))) {
    return round(b, BigInteger.BigInt(1), BigInteger.BigInt(0), rounding);
  }

  if (rounding != null && rounding.maximumSignificantDigits != null && BigInteger.greaterThan(BigInteger.subtract(a.exponent, b.exponent), BigInteger.BigInt(digits(b) + (rounding.maximumSignificantDigits + 1)))) {
    b = create(bigIntSign(b.significand), BigInteger.subtract(a.exponent, BigInteger.BigInt(rounding.maximumSignificantDigits + 1)));
  }

  if (rounding != null && rounding.maximumSignificantDigits != null && BigInteger.greaterThan(BigInteger.subtract(b.exponent, a.exponent), BigInteger.BigInt(digits(a) + (rounding.maximumSignificantDigits + 1)))) {
    a = create(bigIntSign(a.significand), BigInteger.subtract(b.exponent, BigInteger.BigInt(rounding.maximumSignificantDigits + 1)));
  }

  var exponent = bigIntMax(a.exponent, b.exponent);
  return round(create(BigInteger.add(BigInteger.multiply(a.significand, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.subtract(exponent, b.exponent))), BigInteger.multiply(b.significand, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.subtract(exponent, a.exponent)))), bigIntMin(a.exponent, b.exponent)), BigInteger.BigInt(1), BigInteger.BigInt(0), rounding);
};

BigDecimal.subtract = function (a, b, rounding = null) {
  return BigDecimal.add(a, BigDecimal.unaryMinus(b), rounding);
};

BigDecimal.multiply = function (a, b, rounding = null) {
  return normalize(round(create(BigInteger.multiply(a.significand, b.significand), BigInteger.add(a.exponent, b.exponent)), BigInteger.BigInt(1), BigInteger.BigInt(0), rounding));
};

BigDecimal.divide = function (a, b, rounding) {
  if (BigInteger.equal(a.significand, BigInteger.BigInt(0))) {
    return a;
  }

  var exponent = BigInteger.subtract(a.exponent, b.exponent);
  var scaling = BigInteger.BigInt(0);

  if (rounding != null && rounding.maximumSignificantDigits != null) {
    scaling = BigInteger.BigInt(rounding.maximumSignificantDigits + digits(b) - digits(a));
  } else if (rounding != null && rounding.maximumFractionDigits != null) {
    //scaling = BigInt(rounding.maximumFractionDigits) + bigIntMax(a.exponent, 0n) + bigIntMax(0n - b.exponent, 0n) - bigIntMin(a.exponent - b.exponent + BigInt(digits(a) - digits(b)), 0n);
    scaling = BigInteger.add(BigInteger.BigInt(rounding.maximumFractionDigits), exponent);
  } else {
    // Try to do exact division:
    scaling = BigInteger.BigInt(Math.ceil(digits(b) / (Math.log(2) / Math.log(10)) + 1));
  }

  var dividend = BigInteger.multiply(a.significand, BigInteger.greaterThan(scaling, BigInteger.BigInt(0)) ? BigInteger.exponentiate(BigInteger.BigInt(10), scaling) : BigInteger.BigInt(1));
  var divisor = BigInteger.multiply(b.significand, BigInteger.lessThan(scaling, BigInteger.BigInt(0)) ? BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.unaryMinus(scaling)) : BigInteger.BigInt(1));
  var quotient = BigInteger.divide(dividend, divisor);
  var remainder = BigInteger.subtract(dividend, BigInteger.multiply(quotient, divisor));
  return round(create(quotient, BigInteger.subtract(exponent, scaling)), BigInteger.lessThan(divisor, BigInteger.BigInt(0)) ? BigInteger.unaryMinus(divisor) : divisor, BigInteger.lessThan(divisor, BigInteger.BigInt(0)) ? BigInteger.unaryMinus(remainder) : remainder, rounding);
};

function compare(a, b) {
  if (BigInteger.lessThanOrEqual(a.significand, BigInteger.BigInt(0)) && BigInteger.greaterThanOrEqual(b.significand, BigInteger.BigInt(0))) {
    return !(BigInteger.equal(a.significand, BigInteger.BigInt(0)) && BigInteger.equal(b.significand, BigInteger.BigInt(0))) ? -1 : 0;
  }

  if (BigInteger.greaterThanOrEqual(a.significand, BigInteger.BigInt(0)) && BigInteger.lessThanOrEqual(b.significand, BigInteger.BigInt(0))) {
    return BigInteger.equal(a.significand, BigInteger.BigInt(0)) && BigInteger.equal(b.significand, BigInteger.BigInt(0)) ? 0 : +1;
  }

  var differenceOfLogarithms = BigInteger.add(BigInteger.subtract(a.exponent, b.exponent), BigInteger.BigInt(digits(a) - digits(b)));

  if (BigInteger.notEqual(differenceOfLogarithms, BigInteger.BigInt(0))) {
    return BigInteger.lessThan(a.significand, BigInteger.BigInt(0)) && BigInteger.lessThan(b.significand, BigInteger.BigInt(0)) ? BigInteger.greaterThan(differenceOfLogarithms, BigInteger.BigInt(0)) ? -1 : +1 : BigInteger.lessThan(differenceOfLogarithms, BigInteger.BigInt(0)) ? -1 : +1;
  }

  var exponent = bigIntMax(a.exponent, b.exponent);
  var difference = BigInteger.subtract(BigInteger.multiply(a.significand, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.subtract(exponent, b.exponent))), BigInteger.multiply(b.significand, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.subtract(exponent, a.exponent))));
  return BigInteger.lessThan(difference, BigInteger.BigInt(0)) ? -1 : BigInteger.greaterThan(difference, BigInteger.BigInt(0)) ? +1 : 0;
}

BigDecimal.lessThan = function (a, b) {
  return compare(a, b) < 0;
};

BigDecimal.greaterThan = function (a, b) {
  return compare(a, b) > 0;
};

BigDecimal.equal = function (a, b) {
  return compare(a, b) === 0;
};

BigDecimal.round = function (a, rounding) {
  return round(a, BigInteger.BigInt(1), BigInteger.BigInt(0), rounding);
};

BigDecimal.prototype.toString = function () {
  //! https://tc39.es/ecma262/#sec-number.prototype.tostring
  if (arguments.length !== 0) {
    throw new RangeError("not implemented");
  }

  var x = BigDecimal.BigDecimal(this); //! https://tc39.es/ecma262/#sec-numeric-types-number-tostring

  if (BigDecimal.equal(x, BigDecimal.BigDecimal(0))) {
    return "0";
  }

  var sign = "";

  if (BigDecimal.lessThan(x, BigDecimal.BigDecimal(0))) {
    x = BigDecimal.unaryMinus(x);
    sign = "-";
  }

  var getSignificand = function (a, log10) {
    var s = BigDecimal.divide(a, exponentiate(BigDecimal.BigDecimal(10), log10));

    var m = BigDecimal.BigDecimal(Math.pow(10, 15));
    while (!BigDecimal.equal(BigDecimal.round(BigDecimal.multiply(s, m), {
      maximumFractionDigits: 0,
      roundingMode: "half-even"
    }), BigDecimal.multiply(s, m))) {
      m = BigDecimal.multiply(m, m);
    }

    return BigDecimal.toBigInt(BigDecimal.multiply(s, m)).toString().replace(/0+$/g, "") || "0";
  };

  var e = getCountOfDigits(x);
  var significand = getSignificand(x, e);

  if (!BigDecimal.greaterThan(BigDecimal.divide(BigDecimal.BigDecimal(1), BigDecimal.BigDecimal(BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(6)))), x) && BigDecimal.lessThan(x, BigDecimal.BigDecimal(BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(21))))) {
    e = BigInteger.toNumber(e);
    var zeros = Math.max(0, BigInteger.subtract(e, significand.length));

    if (e <= 0) {
      significand = "0".repeat(1 - e) + significand;
      e = 1;
    }

    significand = significand + "0".repeat(zeros);
    return sign + significand.slice(0, e) + (BigInteger.greaterThan(significand.length, e) ? "." + significand.slice(e) : "");
  }

  return sign + (significand.length === 1 ? significand : significand.slice(0, 1) + "." + significand.slice(1)) + "e" + (BigInteger.greaterThanOrEqual(BigInteger.subtract(e, BigInteger.BigInt(1)), BigInteger.BigInt(0)) ? "+" : "") + BigInteger.subtract(e, BigInteger.BigInt(1)).toString();
};

function exponentiate(a, n) {
  if (BigInteger.lessThan(n, BigInteger.BigInt(0))) {
    return BigDecimal.divide(BigDecimal.BigDecimal(1), exponentiate(a, BigInteger.unaryMinus(n)));
  }

  console.assert(BigInteger.greaterThanOrEqual(n, BigInteger.BigInt(0)));
  var accumulator = BigDecimal.BigDecimal(1);
  var x = a;

  while (BigInteger.greaterThan(n, BigInteger.BigInt(0))) {
    if (BigInteger.notEqual(BigInteger.remainder(n, BigInteger.BigInt(2)), BigInteger.BigInt(0))) {
      accumulator = BigDecimal.multiply(accumulator, x);
      n = BigInteger.subtract(n, BigInteger.BigInt(1));
    } else {
      n = BigInteger.divide(n, BigInteger.BigInt(2));
      x = BigDecimal.multiply(x, x);
    }
  }

  return accumulator;
}

function getCountOfDigits(a) {
  // floor(log10(abs(a))) + 1
  if (BigInteger.equal(a.significand, BigInteger.BigInt(0))) {
    throw new RangeError();
  }

  return BigInteger.add(BigInteger.BigInt(digits(a)), a.exponent);
}

function abs(a) {
  return BigDecimal.lessThan(a, BigDecimal.BigDecimal(0)) ? BigDecimal.unaryMinus(a) : a;
}

function sign(x) {
  return BigDecimal.lessThan(x, BigDecimal.BigDecimal(0)) ? BigDecimal.BigDecimal(-1) : BigDecimal.greaterThan(x, BigDecimal.BigDecimal(0)) ? BigDecimal.BigDecimal(1) : BigDecimal.BigDecimal(0);
}

function significandDigits(a) {
  var maximumSignificantDigits = 1;

  while (!BigDecimal.equal(BigDecimal.round(a, {
    maximumSignificantDigits: maximumSignificantDigits,
    roundingMode: "half-even"
  }), a)) {
    maximumSignificantDigits *= 2;
  }

  return maximumSignificantDigits;
}

function tryToMakeCorrectlyRounded(specialValue, f, name) {
  function getExpectedResultIntegerDigits(x) {
    if (name === "exp") {
      // e**x <= 10**k
      // k >= x / log(10)
      return Math.ceil(BigInteger.toNumber(BigDecimal.toBigInt(BigDecimal.round(x, {
        maximumFractionDigits: 0,
        roundingMode: "half-even"
      }))) / Math.log(10));
    }

    if (name === "log") {
      // log(x) <= 10**k
      // log10(log10(x)*log(10)) <= k
      return Math.ceil(Math.log(Math.ceil(Math.max(BigInteger.toNumber(getCountOfDigits(x)), 1) * Math.log(10))) / Math.log(10));
    }

    return 1;
  } // (?) https://en.wikipedia.org/wiki/Rounding#Table-maker's_dilemma


  return function (x, rounding) {
    if (BigDecimal.equal(x, BigDecimal.BigDecimal(specialValue))) {
      return f(x, rounding);
    }

    var result = BigDecimal.BigDecimal(0);
    var i = 0;
    var error = BigDecimal.BigDecimal(0);

    do {
      i += 1;
      var internalRounding = {
        maximumSignificantDigits: Math.ceil(Math.max(rounding.maximumSignificantDigits || rounding.maximumFractionDigits + 1 + getExpectedResultIntegerDigits(x) - 1, significandDigits(x)) * Math.cbrt(Math.pow(2, i - 1))) + 2,
        roundingMode: "half-even"
      };
      result = f(x, internalRounding); // round(result - error) === round(result + error)

      error = BigDecimal.divide(abs(result), BigDecimal.BigDecimal(BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(internalRounding.maximumSignificantDigits)))); //if (i > 0) {
      //console.log(i, f.name, x + "", result + "", error + "", BigDecimal.round(BigDecimal.subtract(result, error), rounding) + "", BigDecimal.round(BigDecimal.add(result, error), rounding) + "");
      //}

      if (i > 10 && rounding.maximumSignificantDigits != null) {
        throw new Error();
      }
    } while (!BigDecimal.equal(BigDecimal.round(BigDecimal.subtract(result, error), rounding), BigDecimal.round(BigDecimal.add(result, error), rounding)));

    if (i > 1) {
      console.log(i);
    }

    return BigDecimal.round(result, rounding);
  };
}

function sqrt(x, rounding) {
  // from https://en.wikipedia.org/wiki/Square_root#Computation
  var lastResult = x;
  var result = BigDecimal.divide(x, BigDecimal.BigDecimal(2));

  while (BigDecimal.lessThan(result, lastResult)) {
    lastResult = result;
    result = BigDecimal.divide(BigDecimal.add(BigDecimal.divide(BigDecimal.BigDecimal(x), result, rounding), result), BigDecimal.BigDecimal(2), rounding);
  }

  return result;
}

BigDecimal.log = tryToMakeCorrectlyRounded(1, function log(x, rounding) {
  if (!BigDecimal.greaterThan(x, BigDecimal.BigDecimal(0))) {
    throw new RangeError();
  } // https://ru.wikipedia.org/wiki/Логарифм#Разложение_в_ряд_и_вычисление_натурального_логарифма


  var internalRounding = {
    maximumSignificantDigits: rounding.maximumSignificantDigits + Math.ceil(Math.log(rounding.maximumSignificantDigits + 0.5) / Math.log(10)),
    roundingMode: "half-even"
  };

  if (true) {
    //! ln(f * 10**k) = ln(f) + k * ln(10), where 0.1 <= f <= 10
    var k = BigInteger.subtract(getCountOfDigits(x), BigInteger.BigInt(1));
    var f = BigDecimal.divide(x, exponentiate(BigDecimal.BigDecimal(10), k));
    var ff = BigDecimal.round(BigDecimal.multiply(f, f), {
      maximumSignificantDigits: 3,
      roundingMode: "half-even"
    });

    if (BigDecimal.greaterThan(ff, BigDecimal.BigDecimal(10))) {
      k = BigInteger.add(k, BigInteger.BigInt(1));
      f = BigDecimal.divide(f, BigDecimal.BigDecimal(10));
    }

    if (BigDecimal.lessThan(ff, BigDecimal.divide(BigDecimal.BigDecimal(1), BigDecimal.BigDecimal(10)))) {
      k = BigInteger.subtract(k, BigInteger.BigInt(1));
      f = BigDecimal.multiply(f, BigDecimal.BigDecimal(10));
    }

    if (BigInteger.notEqual(k, BigInteger.BigInt(0))) {
      return BigDecimal.add(BigDecimal.log(f, internalRounding), BigDecimal.multiply(BigDecimal.BigDecimal(BigInteger.multiply(BigInteger.BigInt(2), k)), BigDecimal.log(BigDecimal.BigDecimal(sqrt(BigDecimal.BigDecimal(10), internalRounding)), internalRounding)));
    }
  } //! log(x) = log((1 + g) / (1 - g)) = 2*(g + g**3/3 + g**5/5 + ...)


  var g = BigDecimal.divide(BigDecimal.subtract(x, BigDecimal.BigDecimal(1)), BigDecimal.add(x, BigDecimal.BigDecimal(1)), internalRounding);
  var n = 1;
  var term = BigDecimal.BigDecimal(1);
  var sum = term;
  var lastSum = BigDecimal.BigDecimal(0);
  var gg = BigDecimal.multiply(g, g, internalRounding);

  while (!BigDecimal.equal(lastSum, sum)) {
    n += 2;
    term = BigDecimal.multiply(term, BigDecimal.BigDecimal(n - 2));
    term = BigDecimal.multiply(term, gg);
    term = BigDecimal.divide(term, BigDecimal.BigDecimal(n), internalRounding);
    lastSum = sum;
    sum = BigDecimal.add(sum, term, internalRounding);
  }

  return BigDecimal.multiply(BigDecimal.multiply(BigDecimal.BigDecimal(2), g), sum);
}, "log");
BigDecimal.exp = tryToMakeCorrectlyRounded(0, function exp(x, rounding) {
  //! k = round(x / ln(10));
  //! exp(x) = exp(x - k * ln(10) + k * ln(10)) = exp(x - k * ln(10)) * 10**k
  var k = BigDecimal.divide(x, BigDecimal.BigDecimal(BigDecimal.divide(BigDecimal.BigDecimal(2302585092994046), BigDecimal.BigDecimal(1000000000000000))), {
    maximumFractionDigits: 0,
    roundingMode: "half-even"
  });
  var internalRounding = {
    maximumSignificantDigits: rounding.maximumSignificantDigits + Math.ceil(Math.log(rounding.maximumSignificantDigits + 0.5) / Math.log(10)),
    roundingMode: "half-even"
  };

  if (!BigDecimal.equal(k, BigDecimal.BigDecimal(0))) {
    var r = BigDecimal.subtract(x, BigDecimal.multiply(k, BigDecimal.log(BigDecimal.BigDecimal(10), {
      maximumSignificantDigits: internalRounding.maximumSignificantDigits + BigInteger.toNumber(getCountOfDigits(k)),
      roundingMode: "half-even"
    })));
    return BigDecimal.multiply(BigDecimal.exp(r, internalRounding), exponentiate(BigDecimal.BigDecimal(10), BigDecimal.toBigInt(k)));
  } // https://en.wikipedia.org/wiki/Exponential_function#Computation


  var n = 0;
  var term = BigDecimal.BigDecimal(1);
  var sum = term;
  var lastSum = BigDecimal.BigDecimal(0);

  while (!BigDecimal.equal(lastSum, sum)) {
    n += 1;
    term = BigDecimal.multiply(term, x);
    term = BigDecimal.divide(term, BigDecimal.BigDecimal(n), internalRounding);
    lastSum = sum;
    sum = BigDecimal.add(sum, term, internalRounding);
  }

  return sum;
}, "exp");

function divideByHalfOfPI(x, rounding) {
  // x = k*pi/2 + r + 2*pi*n, where |r| < pi/4
  if (BigDecimal.lessThan(x, BigDecimal.BigDecimal(0))) {
    throw new RangeError();
  }

  if (BigDecimal.greaterThan(x, BigDecimal.divide(BigDecimal.BigDecimal(785398163397448), BigDecimal.BigDecimal(1000000000000000)))) {
    var halfOfPi = BigDecimal.multiply(BigDecimal.BigDecimal(2), BigDecimal.atan(BigDecimal.BigDecimal(1), {
      maximumSignificantDigits: rounding.maximumSignificantDigits + BigInteger.toNumber(getCountOfDigits(x)),
      roundingMode: "half-even"
    }));
    var i = BigDecimal.divide(x, halfOfPi, {
      maximumFractionDigits: 0,
      roundingMode: "half-even"
    });
    var remainder = BigDecimal.subtract(x, BigDecimal.multiply(i, halfOfPi));
    return {
      remainder: remainder,
      k: (BigInteger.toNumber(BigInteger.remainder(BigDecimal.toBigInt(i), BigInteger.BigInt(4))) + 4) % 4
    };
  }

  return {
    remainder: x,
    k: 0
  };
}

BigDecimal.sin = tryToMakeCorrectlyRounded(0, function (x, rounding) {
  if (BigDecimal.lessThan(x, BigDecimal.BigDecimal(0))) {
    return BigDecimal.unaryMinus(BigDecimal.sin(BigDecimal.unaryMinus(x), rounding));
  }

  var tmp = divideByHalfOfPI(x, rounding);
  var a = tmp.remainder;
  var k = tmp.k;

  if (k === 1) {
    return BigDecimal.cos(a, rounding);
  }

  if (k === 2) {
    return BigDecimal.unaryMinus(BigDecimal.sin(a, rounding));
  }

  if (k === 3) {
    return BigDecimal.unaryMinus(BigDecimal.cos(a, rounding));
  } // https://en.wikipedia.org/wiki/Lookup_table#Computing_sines


  var internalRounding = {
    maximumSignificantDigits: rounding.maximumSignificantDigits + Math.ceil(Math.log(rounding.maximumSignificantDigits + 0.5) / Math.log(10)),
    roundingMode: "half-even"
  };
  var n = 1;
  var term = BigDecimal.BigDecimal(1);
  var sum = term;
  var lastSum = BigDecimal.BigDecimal(0);

  while (!BigDecimal.equal(lastSum, sum)) {
    n += 2;
    term = BigDecimal.multiply(term, BigDecimal.multiply(a, a));
    term = BigDecimal.divide(term, BigDecimal.BigDecimal(-n * (n - 1)), internalRounding);
    lastSum = sum;
    sum = BigDecimal.add(sum, term, internalRounding);
  }

  return BigDecimal.multiply(a, sum);
});
BigDecimal.cos = tryToMakeCorrectlyRounded(0, function (x, rounding) {
  if (BigDecimal.lessThan(x, BigDecimal.BigDecimal(0))) {
    return BigDecimal.cos(BigDecimal.unaryMinus(x), rounding);
  }

  var tmp = divideByHalfOfPI(x, rounding);
  var a = tmp.remainder;
  var k = tmp.k;

  if (k === 1) {
    return BigDecimal.unaryMinus(BigDecimal.sin(a, rounding));
  }

  if (k === 2) {
    return BigDecimal.unaryMinus(BigDecimal.cos(a, rounding));
  }

  if (k === 3) {
    return BigDecimal.sin(a, rounding);
  } // https://en.wikipedia.org/wiki/Trigonometric_functions#Power_series_expansion


  var internalRounding = {
    maximumSignificantDigits: rounding.maximumSignificantDigits + Math.ceil(Math.log(rounding.maximumSignificantDigits + 0.5) / Math.log(10)),
    roundingMode: "half-even"
  };
  var n = 0;
  var term = BigDecimal.BigDecimal(1);
  var sum = term;
  var lastSum = BigDecimal.BigDecimal(0);

  while (!BigDecimal.equal(lastSum, sum)) {
    n += 2;
    term = BigDecimal.multiply(term, BigDecimal.multiply(a, a));
    term = BigDecimal.divide(term, BigDecimal.BigDecimal(-n * (n - 1)), internalRounding);
    lastSum = sum;
    sum = BigDecimal.add(sum, term, internalRounding);
  }

  return sum;
});
BigDecimal.atan = tryToMakeCorrectlyRounded(0, function (x, rounding) {
  if (BigDecimal.greaterThan(abs(x), BigDecimal.BigDecimal(1))) {
    var halfOfPi = BigDecimal.multiply(BigDecimal.atan(BigDecimal.BigDecimal(1), rounding), BigDecimal.BigDecimal(2));
    return BigDecimal.multiply(sign(x), BigDecimal.subtract(halfOfPi, BigDecimal.atan(BigDecimal.divide(BigDecimal.BigDecimal(1), abs(x), rounding), rounding)));
  } // https://en.wikipedia.org/wiki/Inverse_trigonometric_functions#Infinite_series


  var internalRounding = {
    maximumSignificantDigits: rounding.maximumSignificantDigits + Math.ceil(Math.log(rounding.maximumSignificantDigits + 0.5) / Math.log(10)),
    roundingMode: "half-even"
  };
  var n = 0;
  var term = BigDecimal.divide(BigDecimal.BigDecimal(1), BigDecimal.add(BigDecimal.BigDecimal(1), BigDecimal.multiply(x, x)), internalRounding);
  var sum = term;
  var lastSum = BigDecimal.BigDecimal(0);

  while (!BigDecimal.equal(lastSum, sum)) {
    n += 1;
    term = BigDecimal.multiply(term, BigDecimal.BigDecimal(2 * n * (2 * n)));
    term = BigDecimal.divide(term, BigDecimal.BigDecimal(2 * n * (2 * n + 1)), internalRounding);
    term = BigDecimal.multiply(term, BigDecimal.multiply(x, x));
    term = BigDecimal.divide(term, BigDecimal.add(BigDecimal.BigDecimal(1), BigDecimal.multiply(x, x)), internalRounding);
    lastSum = sum;
    sum = BigDecimal.add(sum, term, internalRounding);
  }

  return BigDecimal.multiply(x, sum);
});
export default BigDecimal;

