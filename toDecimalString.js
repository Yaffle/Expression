import Expression from './Expression.js';
import BigInteger from './BigInteger.js';
import nthRoot from './nthRoot.js';

//TODO: ???

// https://en.wikipedia.org/wiki/Fixed-point_arithmetic
function FixedPointContext(roundingMode, scalingCoefficient) {
  if (roundingMode !== "FLOOR" && roundingMode !== "CEIL") {
    throw new RangeError();
  }
  this.roundingMode = roundingMode;
  this.scalingCoefficient = scalingCoefficient;
}
FixedPointContext.prototype.sign = function (x) { // returns -1 or +1
  return BigInteger.lessThan(x, BigInteger.BigInt(0)) ? -1 : +1;
};
FixedPointContext.prototype.round = function (x) { // returns BigInteger
  // rounding to closest, half - away from zero
  var y = this.scalingCoefficient;
  // division of x by y with rounding
  var q = BigInteger.divide(x, y);
  var r = BigInteger.subtract(x, BigInteger.multiply(q, y));
  var r2 = BigInteger.add(r, r);
  return !BigInteger.lessThan(BigInteger.unaryMinus(r2), y) ? BigInteger.subtract(q, BigInteger.BigInt(1)) : (!BigInteger.lessThan(r2, y) ? BigInteger.add(q, BigInteger.BigInt(1)) : q);
};
FixedPointContext.prototype.fromInteger = function (n) {
  return BigInteger.multiply(n, this.scalingCoefficient);
};
FixedPointContext.prototype.compareTo = function (x, y) {
  return BigInteger.lessThan(x, y) ? -1 : (BigInteger.lessThan(y, x) ? +1 : 0);
};
FixedPointContext.prototype.negate = function (x) {
  return BigInteger.unaryMinus(x);
};
FixedPointContext.prototype.add = function (x, y) {
  return BigInteger.add(x, y);
};
FixedPointContext.prototype.subtract = function (x, y) {
  return BigInteger.subtract(x, y);
};
FixedPointContext.prototype._divide = function (x, y) {
  var q = BigInteger.divide(x, y);
  var r = BigInteger.subtract(x, BigInteger.multiply(q, y));
  if (this.roundingMode === "FLOOR") {
    return BigInteger.lessThan(r, BigInteger.BigInt(0)) ? BigInteger.subtract(q, BigInteger.BigInt(1)) : q;
  }
  if (this.roundingMode === "CEIL") {
    return BigInteger.lessThan(BigInteger.BigInt(0), r) ? BigInteger.add(q, BigInteger.BigInt(1)) : q;
  }
  throw new RangeError();
};
FixedPointContext.prototype.multiply = function (x, y) {
  return this._divide(BigInteger.multiply(x, y), this.scalingCoefficient);
};
FixedPointContext.prototype.divide = function (x, y) {
  return this._divide(BigInteger.multiply(x, this.scalingCoefficient), y);
};
FixedPointContext.prototype.min = function (x, y) {
  return this.compareTo(x, y) > 0 ? y : x;
};
FixedPointContext.prototype.max = function (x, y) {
  return this.compareTo(x, y) < 0 ? y : x;
};

//function DegenerateInterval(a) {
//  this.a = a;
//}

// https://en.wikipedia.org/wiki/Interval_arithmetic
function Interval(a, b) {
  this.a = a;
  this.b = b;
}
Interval.Context = function (precision) {
  this.scalingCoefficient = BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(precision));
  this.precision = precision;
  this.down = new FixedPointContext("FLOOR", this.scalingCoefficient);
  this.up = new FixedPointContext("CEIL", this.scalingCoefficient);
};
Interval.Context.prototype.negate = function (x) {
  return new Interval(this.down.negate(x.b), this.up.negate(x.a));
};
Interval.Context.prototype.add = function (x, y) {
  return new Interval(this.down.add(x.a, y.a), this.up.add(x.b, y.b));
};
Interval.Context.prototype.subtract = function (x, y) {
  return this.add(x, this.negate(y));
};
Interval.Context.prototype.multiply = function (x, y) {
  return new Interval(this.down.min(this.down.min(this.down.multiply(x.a, y.a), this.down.multiply(x.a, y.b)), this.down.min(this.down.multiply(x.b, y.a), this.down.multiply(x.b, y.b))),
                      this.up.max(this.up.max(this.up.multiply(x.a, y.a), this.up.multiply(x.a, y.b)), this.up.max(this.up.multiply(x.b, y.a), this.up.multiply(x.b, y.b))));
};
Interval.Context.prototype.divide = function (x, y) {
  if (this.down.compareTo(y.a, this.up.fromInteger(BigInteger.BigInt(0))) <= 0 && this.up.compareTo(y.b, this.down.fromInteger(BigInteger.BigInt(0))) >= 0) {
    //throw new RangeError();
    return "CANNOT_DIVIDE";//TODO: FIX
  }
  return new Interval(this.down.min(this.down.min(this.down.divide(x.a, y.a), this.down.divide(x.a, y.b)), this.down.min(this.down.divide(x.b, y.a), this.down.divide(x.b, y.b))),
                      this.up.max(this.up.max(this.up.divide(x.a, y.a), this.up.divide(x.a, y.b)), this.up.max(this.up.divide(x.b, y.a), this.up.divide(x.b, y.b))));
};
Interval.Context.prototype.nthRoot = function (A, n) {
  var c = this.scalingCoefficient;
  var sA = BigInteger.multiply(A, BigInteger.exponentiate(c, BigInteger.BigInt(n)));

  var x0 = nthRoot(sA, BigInteger.BigInt(n));
  var x1 = BigInteger.lessThan(BigInteger.exponentiate(x0, BigInteger.BigInt(n)), sA) ? BigInteger.add(x0, BigInteger.BigInt(1)) : x0;

  var a = this.down.divide(this.down.fromInteger(x0), this.down.fromInteger(c));
  var b = this.up.divide(this.up.fromInteger(x1), this.up.fromInteger(c));
  return new Interval(a, b);
};
Interval.Context.prototype.fromInteger = function (a, b) {
  return new Interval(this.down.fromInteger(a), this.up.fromInteger(a));
};
Interval.Context.prototype.fromIntegers = function (a, b) {
  return new Interval(this.down.fromInteger(a), this.up.fromInteger(b));
};
//?
Interval.Context.prototype.toInteger = function (x) {
  var signA = this.down.sign(x.a);
  var signB = this.up.sign(x.b);
  if (signA === signB) {
    var candidateA = this.down.round(x.a);
    var candidateB = this.up.round(x.b);
    if (!BigInteger.lessThan(candidateA, candidateB)) {
      return {sign: signA, integer: candidateA};
    }
  }
  return {sign: undefined, integer: undefined};
};
Interval.prototype.toString = function () {
  return "[" + this.a.toString() + ";" + this.b.toString() + "]";
};

var evaluateExpression = function (e, context) {
  if (e instanceof Expression.Integer) {
    var n = e.value;
    return context.fromInteger(n);
  } else if (e instanceof Expression.NthRoot) {
    var a = e.a;
    var n = e.n;
    if (a instanceof Expression.Integer) {
      return context.nthRoot(a.value, n);
    }
  } else if (e instanceof Expression.BinaryOperation) {
    var a = evaluateExpression(e.a, context);
    var b = evaluateExpression(e.b, context);
    if (a === "CANNOT_DIVIDE") {
      return a;
    }
    if (b === "CANNOT_DIVIDE") {
      return b;
    }
    if (a != undefined && b != undefined) {
      var operator = e.getS();
      if (operator === "+") {
        return context.add(a, b);
      } else if (operator === "-") {
        return context.subtract(a, b);
      } else if (operator === "*") {
        return context.multiply(a, b);
      } else if (operator === "/") {
        return context.divide(a, b);
      } else if (operator === "^") { // Expression.PolynomialRoot^3
        var result = a;
        var n = Number.parseInt(e.b.value.toString(), 10);//TODO: FIX!
        for (var i = 1; i < n; i += 1) {
          result = context.multiply(result, a);
        }
        return result;
      }
    }
  } else if (e instanceof Expression.PolynomialRoot) {
    var i = e.polynomial.getZero(e.interval, context.precision);
    var cd = i.a.getDenominator().lcm(i.b.getDenominator());
    return context.divide(context.fromIntegers(i.a.getNumerator().multiply(cd.divide(i.a.getDenominator())).value,
                                               i.b.getNumerator().multiply(cd.divide(i.b.getDenominator())).value),
                         context.fromInteger(cd.value));
  }
  return undefined;
};

var decimalToString = function (sign, number) {
  return (sign < 0 ? "-" : "") + number;
};

var complexToString = function (real, imaginary, imaginarySign) {
  return real + (imaginarySign >= 0 ? "+" : "") + imaginary + "i";
};

var digitsToDecimalNumber = function (sign, value, fractionDigits, decimalToStringCallback) {
  // TODO: fix
  // new Intl.NumberFormat().format(1.1)
  // "<math decimalpoint=\"" + decimalSeparator + "\"></math>" -?
  var digits = (BigInteger.lessThan(value, BigInteger.BigInt(0)) ? BigInteger.unaryMinus(value) : value).toString();
  var decimalSeparator = ".";
  var zeros = "";
  for (var i = 0; i < fractionDigits; i += 1) {
    zeros += "0";
  }
  var number = (fractionDigits === 0 ? digits : (digits.slice(0, -fractionDigits) || "0") + decimalSeparator + (zeros + digits).slice(-fractionDigits));
  return decimalToStringCallback(sign, number);
};

var toDecimalNumberOld = function (numerator, denominator, fractionDigits, decimalToStringCallback) {
  var sign = (BigInteger.lessThan(numerator, BigInteger.BigInt(0)) ? -1 : +1) * (BigInteger.lessThan(denominator, BigInteger.BigInt(0)) ? -1 : +1);
  var abs = function (x) {
    return BigInteger.lessThan(x, BigInteger.BigInt(0)) ? BigInteger.unaryMinus(x) : x;
  };
  var n = abs(numerator);
  var d = abs(denominator);

  //? ((n * 10**(fractionDigits + 1)) ~/ d + 5) ~/ 10

  var x = BigInteger.multiply(n, BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(fractionDigits)));
  var q = BigInteger.divide(x, d);
  var r = BigInteger.subtract(x, BigInteger.multiply(q, d));
  if (!BigInteger.lessThan(BigInteger.add(r, r), d)) {
    q = BigInteger.add(q, BigInteger.BigInt(1));
  }
  var result = q;
  return digitsToDecimalNumber(sign, result, fractionDigits, decimalToStringCallback);
};

var toDecimalStringInternal = function (expression, fractionDigits, decimalToStringCallback, complexToStringCallback) {
  decimalToStringCallback = decimalToStringCallback || decimalToString;
  complexToStringCallback = complexToStringCallback || complexToString;
  //TODO: remove
  if (expression instanceof Expression.Integer) {
    return toDecimalNumberOld(expression.value, BigInteger.BigInt(1), fractionDigits, decimalToStringCallback);
  }
  if (expression instanceof Expression.Division) {
    var numerator = expression.getNumerator();//.unwrap();
    var denominator = expression.getDenominator();//.unwrap();
    if (numerator instanceof Expression.Integer && denominator instanceof Expression.Integer) {
      return toDecimalNumberOld(numerator.value, denominator.value, fractionDigits, decimalToStringCallback);
    }
    if (numerator instanceof Expression.Complex && denominator instanceof Expression.Integer) {
      var real = toDecimalNumberOld(numerator.real.value, denominator.value, fractionDigits, decimalToStringCallback);
      var imaginary = toDecimalNumberOld(numerator.imaginary.value, denominator.value, fractionDigits, decimalToStringCallback);
      var imaginarySign = numerator.imaginary.compareTo(Expression.ZERO);
      return complexToStringCallback(real, imaginary, imaginarySign);
    }
  }
  if (expression instanceof Expression.NthRoot) {
    var a = expression.a;//.unwrap();
    if (a instanceof Expression.Integer) {
      var A = a.value;
      var n = expression.n;
      var c = BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(fractionDigits));
      var sA = BigInteger.multiply(A, BigInteger.exponentiate(c, BigInteger.BigInt(n)));

      var x0 = nthRoot(sA, BigInteger.BigInt(n));
      var x1 = BigInteger.lessThan(BigInteger.exponentiate(x0, BigInteger.BigInt(n)), sA) ? BigInteger.add(x0, BigInteger.BigInt(1)) : x0;

      // root - x0 < x1 - root
      // 2root < x0 + x1
      // 2**n * A < (x0 + x1)**n
      var nearest = BigInteger.lessThan(BigInteger.multiply(BigInteger.exponentiate(BigInteger.BigInt(2), BigInteger.BigInt(n)), sA), BigInteger.exponentiate(BigInteger.add(x0, x1), BigInteger.BigInt(n))) ? x0 : x1;
      return toDecimalNumberOld(nearest, c, fractionDigits, decimalToStringCallback);
    }
  }
  //---
  if (!Expression.has(expression, Expression.NthRoot) &&
      !Expression.has(expression, Expression.PolynomialRoot)) {
    throw new TypeError("toDecimalString:" + fractionDigits + ":" + expression.toString({}));
  }
  if (fractionDigits < 0 || fractionDigits > 9007199254740991) {
    throw new RangeError();
  }
  var sign = 0;
  var result = undefined;
  var guessedPrecision = 0;
  var scale = BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(fractionDigits));
  while (result == undefined) {
    var context = new Interval.Context(guessedPrecision);
    var x = evaluateExpression(expression, context);
    if (x == undefined) {
      return undefined;
    }
    if (x !== "CANNOT_DIVIDE") { // continue the loop otherwise
      x = context.multiply(context.fromInteger(scale), x);
      var tmp = context.toInteger(x);
      sign = tmp.sign;
      result = tmp.integer;
    }
    guessedPrecision = (guessedPrecision === 0 ? 1 : guessedPrecision * 2);
  }
  return digitsToDecimalNumber(sign, result, fractionDigits, decimalToStringCallback);
};

export default toDecimalStringInternal;
