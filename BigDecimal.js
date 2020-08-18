
import BigInteger from './BigInteger.js';

// https://en.wikipedia.org/wiki/Fixed-point_arithmetic
// https://github.com/tc39/proposal-decimal
function BigDecimal(significand, exponent) {
  this.significand = significand;
  this.exponent = exponent;
}
BigDecimal.toBigInt = function (a) {
  var scalingCoefficient = BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(0 - a.exponent));
  return BigInteger.divide(a.significand, scalingCoefficient);
};
BigDecimal.BigDecimal = function (n) {
  return new BigDecimal(n, BigInteger.BigInt(0));
};
BigDecimal.lessThan = function (a, b) {
  var exponent = Math.max(a.exponent, b.exponent);
  var x = BigInteger.multiply(a.significand, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(exponent - b.exponent)));
  var y = BigInteger.multiply(b.significand, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(exponent - a.exponent)));
  return BigInteger.lessThan(x, y);
};
BigDecimal.greaterThan = function (x, y) {
  return BigDecimal.lessThan(y, x);
};
BigDecimal.equal = function (x, y) {
  return !BigDecimal.lessThan(x, y) && !BigDecimal.lessThan(y, x);
};
BigDecimal.unaryMinus = function (x) {
  return new BigDecimal(BigInteger.unaryMinus(x.significand), x.exponent);
};
BigDecimal.add = function (a, b, rounding) {
  var exponent = Math.max(a.exponent, b.exponent);
  var x = BigInteger.multiply(a.significand, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(exponent - b.exponent)));
  var y = BigInteger.multiply(b.significand, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(exponent - a.exponent)));
  return new BigDecimal(BigInteger.add(x, y), Math.min(a.exponent, b.exponent));
};
BigDecimal.subtract = function (x, y, rounding) {
  return BigDecimal.add(x, BigDecimal.unaryMinus(y), rounding);
};
var divide = function (x, y, roundingMode) {
  if (BigInteger.lessThan(y, BigInteger.BigInt(0))) {
    x = BigInteger.unaryMinus(x);
    y = BigInteger.unaryMinus(y);
  }
  var q = BigInteger.divide(x, y);
  var r = BigInteger.subtract(x, BigInteger.multiply(q, y));
  if (roundingMode === 'ceil') {
    var ceil = BigInteger.lessThan(BigInteger.BigInt(0), r) ? BigInteger.add(q, BigInteger.BigInt(1)) : q;
    return ceil;
  } else if (roundingMode === 'floor') {
    var floor = BigInteger.lessThan(r, BigInteger.BigInt(0)) ? BigInteger.subtract(q, BigInteger.BigInt(1)) : q;
    return floor;
  } else if (roundingMode === 'half-up') {
    if (!BigInteger.lessThan(BigInteger.add(r, r), y)) {
      q = BigInteger.add(q, BigInteger.BigInt(1));
    }
    if (!BigInteger.greaterThan(BigInteger.add(r, r), BigInteger.unaryMinus(y))) {
      q = BigInteger.subtract(q, BigInteger.BigInt(1));
    }
    return q;
  } else if (roundingMode === 'half-down') {
    if (BigInteger.greaterThan(BigInteger.add(r, r), y)) {
      q = BigInteger.add(q, BigInteger.BigInt(1));
    }
    if (BigInteger.lessThan(BigInteger.add(r, r), BigInteger.unaryMinus(y))) {
      q = BigInteger.subtract(q, BigInteger.BigInt(1));
    }
    return q;
  } else if (roundingMode === 'half-even') {
    if (BigInteger.greaterThan(BigInteger.add(r, r), y) || (BigInteger.equal(BigInteger.add(r, r), y) && !BigInteger.equal(BigInteger.remainder(q, BigInteger.BigInt(2)), BigInteger.BigInt(0)))) {
      q = BigInteger.add(q, BigInteger.BigInt(1));
    }
    if (BigInteger.lessThan(BigInteger.add(r, r), BigInteger.unaryMinus(y)) || (BigInteger.equal(BigInteger.add(r, r), BigInteger.unaryMinus(y)) && !BigInteger.equal(BigInteger.remainder(q, BigInteger.BigInt(2)), BigInteger.BigInt(0)))) {
      q = BigInteger.subtract(q, BigInteger.BigInt(1));
    }
    return q;
  }
  throw new RangeError("unknown rounding mode: " + roundingMode);
};
BigDecimal.round = function (bd, rounding) {
  if (rounding != null && bd.exponent + rounding.maximumFractionDigits < 0) {
    var scalingCoefficient = BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(0 - (bd.exponent + rounding.maximumFractionDigits)));
    return new BigDecimal(divide(bd.significand, scalingCoefficient, rounding.roundingMode), 0 - rounding.maximumFractionDigits);
  }
  return bd;
};
BigDecimal.multiply = function (a, b, rounding) {
  return BigDecimal.round(new BigDecimal(BigInteger.multiply(a.significand, b.significand), a.exponent + b.exponent), rounding);
};
BigDecimal.divide = function (a, b, rounding) {
  console.assert(a.exponent <= 0 && b.exponent <= 0);
  var x = BigInteger.multiply(a.significand, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(0 - b.exponent)));
  var y = BigInteger.multiply(b.significand, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(0 - a.exponent)));
  var scalingCoefficient = BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(rounding.maximumFractionDigits));
  return new BigDecimal(divide(BigInteger.multiply(x, scalingCoefficient), y, rounding.roundingMode), -rounding.maximumFractionDigits);
};

var toFixedPoint = function (a, rounding) {
  var s = BigDecimal.BigDecimal(BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(rounding.maximumSignificantDigits || rounding.maximumFractionDigits)));
  return BigDecimal.toBigInt(BigDecimal.multiply(a, s));
};
var toBigDecimal = function (x, rounding) {
  var s = BigDecimal.BigDecimal(BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(rounding.maximumSignificantDigits || rounding.maximumFractionDigits)));
  return BigDecimal.divide(BigDecimal.BigDecimal(x), s, rounding);
};
BigDecimal.atan = function (x, rounding) {
  if (!BigDecimal.equal(x, BigDecimal.BigDecimal(BigInteger.BigInt(1)))) {
    throw new RangeError("not supported");
  }
  // https://en.wikipedia.org/wiki/Approximations_of_π#Arctangent
  function factorial(n) {
    var result = BigInteger.BigInt(1);
    for (var i = 1; i <= n; i += 1) {
      result = BigInteger.multiply(result, BigInteger.BigInt(i));
    }
    return result;
  }
  function calculatePI(precision) {
    var iterations = Math.floor(10 * (precision + 1) / 3) + 1;
    var scale = factorial(2 * iterations + 1);
    var pi = BigInteger.BigInt(0);
    var n = 0;
    var nf = BigInteger.BigInt(1);
    var _2np1f = BigInteger.BigInt(1);
    while (n <= iterations) {
      var x = BigInteger.multiply(BigInteger.multiply(BigInteger.divide(scale, _2np1f), BigInteger.multiply(nf, nf)), BigInteger.exponentiate(BigInteger.BigInt(2), BigInteger.BigInt(n + 1)));
      pi = BigInteger.add(pi, x);
      n += 1;
      nf = BigInteger.multiply(nf, BigInteger.BigInt(n));
      _2np1f = BigInteger.multiply(_2np1f, BigInteger.BigInt(2 * n));
      _2np1f = BigInteger.multiply(_2np1f, BigInteger.BigInt(2 * n + 1));
    }
    return BigInteger.divide(BigInteger.multiply(pi, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(precision))), scale);
  }
  if (!BigDecimal.equal(x, BigDecimal.BigDecimal(BigInteger.BigInt(1)))) {
    throw new TypeError();
  }
  var a = toBigDecimal(calculatePI(rounding.maximumSignificantDigits || rounding.maximumFractionDigits), rounding);
  return BigDecimal.divide(a, BigDecimal.BigDecimal(BigInteger.BigInt(4)), rounding);
};
BigDecimal.exp = function (x, rounding) {

  function expp(a, b, scaling) {
    // a > 0, b > 0, scaling > 0
    var y = BigInteger.BigInt(1);
    var t = BigInteger.BigInt(1);
    var u = BigInteger.BigInt(1);
    var k = BigInteger.BigInt(0);
    var denominator = BigInteger.BigInt(1);
    while (BigInteger.greaterThan(BigInteger.multiply(BigInteger.multiply(t, scaling), BigInteger.BigInt(2)), denominator) ||
           BigInteger.greaterThan(BigInteger.multiply(a, BigInteger.BigInt(2)), BigInteger.multiply(b, k))) {
      t = BigInteger.multiply(t, a);
      k = BigInteger.add(k, BigInteger.BigInt(1));
      y = BigInteger.add(BigInteger.multiply(y, BigInteger.multiply(k, b)), t);
      denominator = BigInteger.multiply(denominator, BigInteger.multiply(k, b));
    }
    return {
      numerator: BigInteger.multiply(y, scaling),
      denominator: denominator
    };
  }
  function exp(a, b, scaling) {
    if (BigInteger.lessThan(a, BigInteger.BigInt(0))) {
      var tmp = expp(BigInteger.unaryMinus(a), b, scaling);
      return BigInteger.divide(BigInteger.multiply(BigInteger.multiply(scaling, scaling), tmp.denominator), tmp.numerator);
    }
    var tmp = expp(a, b, scaling);
    return BigInteger.divide(tmp.numerator, tmp.denominator);
  }

  var scalingCoefficient = BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(rounding.maximumSignificantDigits || rounding.maximumFractionDigits));
  return toBigDecimal(exp(toFixedPoint(x, rounding), scalingCoefficient, scalingCoefficient), rounding);
};
BigDecimal.log = function (x, rounding) {
  //TODO: verify
  function factorial(n) {
    var result = BigInteger.BigInt(1);
    for (var i = 1; i <= n; i += 1) {
      result = BigInteger.multiply(result, BigInteger.BigInt(i));
    }
    return result;
  }
  var logarithm = function (n, scalingCoefficient, precision) { // ln(z) = ln(n / scalingCoefficient)
    if (!BigInteger.equal(scalingCoefficient, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(precision)))) {
      throw new RangeError();
    }
    if (BigInteger.lessThan(n, BigInteger.BigInt(1))) {
      throw new TypeError("NotSupportedError");
    }
    // https://en.wikipedia.org/wiki/Logarithm#:~:text=More%20efficient%20series
    var s = BigInteger.BigInt(0);
    var ps = BigInteger.BigInt(1);
    //TODO: fix
    var scale = BigInteger.multiply(factorial(precision), scalingCoefficient);
    for (var i = BigInteger.BigInt(1); ps > BigInteger.BigInt(0); i = BigInteger.add(i, BigInteger.BigInt(2))) {
      ps = BigInteger.divide(BigInteger.multiply(BigInteger.divide(scale, i), BigInteger.exponentiate(BigInteger.subtract(n, scalingCoefficient), i)), BigInteger.exponentiate(BigInteger.add(n, scalingCoefficient), i));
      s = BigInteger.add(s, ps);
    }
    return BigInteger.divide(BigInteger.multiply(BigInteger.multiply(BigInteger.BigInt(2), s), scalingCoefficient), scale);
  };
  var scalingCoefficient = BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(rounding.maximumSignificantDigits || rounding.maximumFractionDigits));
  return toBigDecimal(logarithm(toFixedPoint(x, rounding), scalingCoefficient, rounding.maximumSignificantDigits || rounding.maximumFractionDigits), rounding);
};
var trigonometry = function (x, start, rounding) {
  function factorial(n) {
    var result = BigInteger.BigInt(1);
    for (var i = 1; i <= n; i += 1) {
      result = BigInteger.multiply(result, BigInteger.BigInt(i));
    }
    return result;
  }
  //if (precision > 10) throw new Error();
  //throw new Error();
  function trigonometry(xn, xd, start, precision) {
    var iterations = Math.floor(10 * (precision + 1) / 3) + 1;//TODO: ?
    //iterations += bitLength(x);//?
    //TODO: optimize
    var x = BigInteger.add(BigInteger.divide(xn, xd), BigInteger.BigInt(1));
    iterations += Math.floor(Math.log(Number(x) + 0.5) / Math.log(2)) * (Math.floor(Number(x)) + 1);
    if (iterations > 1000) {
      throw new TypeError("NotSupportedError");//!TODO: fix
    }
    var scale = factorial(iterations);
    scale = BigInteger.multiply(scale, BigInteger.exponentiate(xd, BigInteger.BigInt(iterations)));
    var y = BigInteger.BigInt(0);
    var k = start;
    while (k < iterations) {
      var s = BigInteger.multiply(BigInteger.BigInt(Math.pow(-1, (k - start) / 2)), BigInteger.divide(scale, factorial(k)));
      s = BigInteger.divide(s, BigInteger.exponentiate(xd, BigInteger.BigInt(k)));
      s = BigInteger.multiply(s, BigInteger.exponentiate(xn, BigInteger.BigInt(k)));
      y = BigInteger.add(y, s);
      k += 2;
    }
    y = BigInteger.divide(BigInteger.multiply(y, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(precision))), scale);
    return y;
  }
  var scalingCoefficient = BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(rounding.maximumSignificantDigits || rounding.maximumFractionDigits));
  return toBigDecimal(trigonometry(toFixedPoint(x, rounding), scalingCoefficient, start, rounding.maximumSignificantDigits || rounding.maximumFractionDigits), rounding);
};
BigDecimal.sin = function (x, rounding) {
  return trigonometry(x, 1, rounding);
};
BigDecimal.cos = function (x, rounding) {
  return trigonometry(x, 0, rounding);
};

export default BigDecimal;
