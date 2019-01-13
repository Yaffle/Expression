/*global Expression, BigInteger, self*/

(function (global) {
"use strict";

BigInteger.ZERO = BigInteger.fromNumber(0);
BigInteger.ONE = BigInteger.fromNumber(1);
BigInteger.TWO = BigInteger.fromNumber(2);
BigInteger.TEN = BigInteger.fromNumber(10);

function nthRoot(A, n) {
  var x = BigInteger.nthRoot(A, n);
  return {x0: x, x1: BigInteger.compareTo(A, BigInteger.pow(x, n)) === 0 ? x : BigInteger.add(x, BigInteger.ONE)};
}

//TODO: ???

function FP() {
}
FP.Context = function (roundingMode, precision) {
  this.roundingMode = roundingMode;
  this.precision = precision;
  this.scalingCoefficient = BigInteger.pow(BigInteger.TEN, this.precision);
};
FP.Context.prototype.sign = function (x) { // returns -1 or +1
  return BigInteger.compareTo(x, BigInteger.ZERO) < 0 ? -1 : +1;
};
FP.Context.prototype.round = function (x) { // returns BigInteger
  // rounding to closest, half - away from zero
  var roundDivision = function (a, b) {
    var q = BigInteger.divide(a, b);
    var r = BigInteger.remainder(a, b);
    var r2 = BigInteger.add(r, r);
    return BigInteger.compareTo(BigInteger.negate(r2), b) >= 0 ? BigInteger.subtract(q, BigInteger.ONE) : (BigInteger.compareTo(r2, b) >= 0 ? BigInteger.add(q, BigInteger.ONE) : q);
  };
  return roundDivision(x, this.scalingCoefficient);
};
FP.Context.prototype.fromInteger = function (n) {
  return BigInteger.multiply(n, this.scalingCoefficient);
};
FP.Context.prototype.compareTo = function (x, y) {
  return BigInteger.compareTo(x, y);
};
FP.Context.prototype.negate = function (x) {
  return BigInteger.negate(x);
};
FP.Context.prototype.add = function (x, y) {
  return BigInteger.add(x, y);
};
FP.Context.prototype.subtract = function (x, y) {
  return BigInteger.subtract(x, y);
};
FP.Context.prototype._divide = function (x, y) {
  var q = BigInteger.divide(x, y);
  var r = BigInteger.remainder(x, y);
  if (this.roundingMode === "FLOOR") {
    return BigInteger.compareTo(r, BigInteger.ZERO) < 0 ? BigInteger.subtract(q, BigInteger.ONE) : q;
  }
  if (this.roundingMode === "CEIL") {
    return BigInteger.compareTo(r, BigInteger.ZERO) > 0 ? BigInteger.add(q, BigInteger.ONE) : q;
  }
  throw new RangeError();
};
FP.Context.prototype.multiply = function (x, y) {
  return this._divide(BigInteger.multiply(x, y), this.scalingCoefficient);
};
FP.Context.prototype.divide = function (x, y) {
  return this._divide(BigInteger.multiply(x, this.scalingCoefficient), y);
};
FP.Context.prototype.min = function (x, y) {
  return this.compareTo(x, y) > 0 ? y : x;
};
FP.Context.prototype.max = function (x, y) {
  return this.compareTo(x, y) < 0 ? y : x;
};

function Interval(a, b) {
  this.a = a;
  this.b = b;
}
Interval.Context = function (precision) {
  this.down = new FP.Context("FLOOR", precision);
  this.up = new FP.Context("CEIL", precision);
  this.precision = precision;
  this.zero = this.down.fromInteger(BigInteger.ZERO);
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
  if (this.down.compareTo(y.a, this.zero) <= 0 && this.down.compareTo(y.b, this.zero) >= 0) {
    //throw new RangeError();
    return "CANNOT_DIVIDE";//TODO: FIX
  }
  return new Interval(this.down.min(this.down.min(this.down.divide(x.a, y.a), this.down.divide(x.a, y.b)), this.down.min(this.down.divide(x.b, y.a), this.down.divide(x.b, y.b))),
                      this.up.max(this.up.max(this.up.divide(x.a, y.a), this.up.divide(x.a, y.b)), this.up.max(this.up.divide(x.b, y.a), this.up.divide(x.b, y.b))));
};
Interval.Context.prototype.nthRoot = function (A, n) {
  var c = BigInteger.pow(BigInteger.TEN, this.precision);
  var sA = BigInteger.multiply(A, BigInteger.pow(c, n));
  var tmp = nthRoot(sA, n);
  var a = this.down.divide(this.down.fromInteger(tmp.x0), this.down.fromInteger(c));
  var b = this.up.divide(this.up.fromInteger(tmp.x1), this.up.fromInteger(c));
  return new Interval(a, b);
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
    if (BigInteger.compareTo(candidateA, candidateB) === 0) {
      return {sign: signA, integer: candidateA};
    }
  }
  return {sign: undefined, integer: undefined};
};
Interval.prototype.toString = function () {
  return "[" + this.a.toString() + ";" + this.b.toString() + "]";
};

Expression.prototype.evaluate = function (context) {
  var e = this;
  if (e instanceof Expression.Integer) {
    var n = e.value;
    return context.fromIntegers(n, n);
  } else if (e instanceof Expression.NthRoot) {
    var a = e.a;
    var n = e.n;
    if (a instanceof Expression.Integer) {
      return context.nthRoot(a.value, n);
    }
  } else if (e instanceof Expression.BinaryOperation) {
    var a = e.a.evaluate(context);
    var b = e.b.evaluate(context);
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
    return context.divide(context.fromIntegers(i.a.getNumerator().multiply(cd.divide(i.a.getDenominator())).value, i.b.getNumerator().multiply(cd.divide(i.b.getDenominator())).value), context.fromIntegers(cd.value, cd.value));
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
  var digits = (BigInteger.compareTo(value, BigInteger.ZERO) < 0 ? BigInteger.negate(value) : value).toString();
  var decimalSeparator = ".";
  var zeros = "";
  for (var i = 0; i < fractionDigits; i += 1) {
    zeros += "0";
  }
  var number = (fractionDigits === 0 ? digits : (digits.slice(0, -fractionDigits) || "0") + decimalSeparator + (zeros + digits).slice(-fractionDigits));
  return decimalToStringCallback(sign, number);
};

var toDecimalNumberOld = function (numerator, denominator, fractionDigits, decimalToStringCallback) {
  var sign = (BigInteger.compareTo(numerator, BigInteger.ZERO) < 0 ? -1 : +1) * (BigInteger.compareTo(denominator, BigInteger.ZERO) < 0 ? -1 : +1);
  var abs = function (x) {
    return BigInteger.compareTo(x, BigInteger.ZERO) < 0 ? BigInteger.negate(x) : x;
  };
  var n = abs(numerator);
  var d = abs(denominator);

  //? ((n * 10**(fractionDigits + 1)) ~/ d + 5) ~/ 10

  var x = BigInteger.multiply(n, BigInteger.pow(BigInteger.TEN, fractionDigits));
  var q = BigInteger.divide(x, d);
  var r = BigInteger.subtract(x, BigInteger.multiply(q, d));
  if (BigInteger.compareTo(BigInteger.add(r, r), d) >= 0) {
    q = BigInteger.add(q, 1);
    r = BigInteger.subtract(r, d);
  }
  return digitsToDecimalNumber(sign, q, fractionDigits, decimalToStringCallback);
};

Expression.toDecimalStringInternal = function (expression, fractionDigits, decimalToStringCallback, complexToStringCallback) {
  decimalToStringCallback = decimalToStringCallback || decimalToString;
  complexToStringCallback = complexToStringCallback || complexToString;
  if (expression instanceof Expression.Integer) {
    return toDecimalNumberOld(expression.value, BigInteger.ONE, fractionDigits, decimalToStringCallback);
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
      var c = BigInteger.pow(BigInteger.TEN, fractionDigits);
      var sA = BigInteger.multiply(A, BigInteger.pow(c, n));
      var tmp = nthRoot(sA, n);
      var x0 = tmp.x0;
      var x1 = tmp.x1;
      // root - x0 < x1 - root
      // 2root < x0 + x1
      // 2**n * A < (x0 + x1)**n
      var nearest = BigInteger.compareTo(BigInteger.multiply(BigInteger.pow(BigInteger.TWO, n), sA), BigInteger.pow(BigInteger.add(x0, x1), n)) < 0 ? x0 : x1;
      return toDecimalNumberOld(nearest, c, fractionDigits, decimalToStringCallback);
    }
  }
  //---
  if (!Expression.has(expression, Expression.NthRoot) &&
      !Expression.has(expression, Expression.PolynomialRoot)) {
    self.setTimeout(function () {
      throw new TypeError("toDecimalString:" + fractionDigits + ":" + expression.toString({}));
    }, 0);
  }
  if (fractionDigits < 0 || fractionDigits > 9007199254740991) {
    throw new RangeError();
  }
  var sign = 0;
  var result = undefined;
  var guessedPrecision = fractionDigits + 1;
  while (result == undefined) {
    var context = new Interval.Context(guessedPrecision);
    var x = new Expression.Multiplication(expression, new Expression.Integer(BigInteger.pow(BigInteger.TEN, fractionDigits))).evaluate(context);
    if (x === "CANNOT_DIVIDE") {
      x = context.fromIntegers(BigInteger.negate(BigInteger.ONE), BigInteger.ONE); // to continue the loop
    }
    if (x == undefined) {
      return undefined;
    }
    var tmp = context.toInteger(x);
    sign = tmp.sign;
    result = tmp.integer;
    guessedPrecision *= 2;
  }
  return digitsToDecimalNumber(sign, result, fractionDigits, decimalToStringCallback);
};

}(this));
