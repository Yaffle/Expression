import Expression from './Expression.js';
import BigInteger from './BigInteger.js';
import nthRoot from './nthRoot.js';

//TODO: ???

// https://en.wikipedia.org/wiki/Fixed-point_arithmetic
function FixedPointContext(scalingCoefficient) {
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
  return Interval.degenerate(BigInteger.multiply(n, this.scalingCoefficient));
};
FixedPointContext.prototype.compareTo = function (x, y) {
  return BigInteger.lessThan(x, y) ? -1 : (BigInteger.lessThan(y, x) ? +1 : 0);
};
FixedPointContext.prototype.negate = function (x) {
  return Interval.degenerate(BigInteger.unaryMinus(x));
};
FixedPointContext.prototype.add = function (x, y) {
  return Interval.degenerate(BigInteger.add(x, y));
};
FixedPointContext.prototype.subtract = function (x, y) {
  return BigInteger.subtract(x, y);
};
FixedPointContext.prototype._divide = function (x, y) {
  var q = BigInteger.divide(x, y);
  var r = BigInteger.subtract(x, BigInteger.multiply(q, y));
  var floor = BigInteger.lessThan(r, BigInteger.BigInt(0)) ? BigInteger.subtract(q, BigInteger.BigInt(1)) : q;
  var ceil = BigInteger.lessThan(BigInteger.BigInt(0), r) ? BigInteger.add(q, BigInteger.BigInt(1)) : q;
  return new Interval(floor, ceil);
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


// https://en.wikipedia.org/wiki/Interval_arithmetic
function Interval(a, b) {
  if (a instanceof Interval || b instanceof Interval) {
    throw new TypeError();// to help with debugging
  }
  this.a = a;
  this.b = b;
}
Interval.degenerate = function (a) {
  return new Interval(a, a);
};

Interval.Context = function (precision) {
  this.scalingCoefficient = BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(precision));
  this.precision = precision;
  this.c = new FixedPointContext(this.scalingCoefficient);
};
Interval.Context.prototype.negate = function (x) {
  if (this.c.compareTo(x.a, x.b) === 0) {
    return this.c.negate(x.a);
  }
  return new Interval(this.c.negate(x.b).a, this.c.negate(x.a).b);
};
Interval.Context.prototype.add = function (x, y) {
  if (this.c.compareTo(x.a, x.b) === 0 && this.c.compareTo(y.a, y.b) === 0) {
    return this.c.add(x.a, y.a);
  }
  return new Interval(this.c.add(x.a, y.a).a, this.c.add(x.b, y.b).b);
};
Interval.Context.prototype.subtract = function (x, y) {
  return this.add(x, this.negate(y));
};
Interval.Context.prototype.multiply = function (x, y) {
  if (this.c.compareTo(x.a, x.b) === 0 && this.c.compareTo(y.a, y.b) === 0) {
    return this.c.multiply(x.a, y.a);
  }
  var a = this.c.multiply(x.a, y.a);
  var b = this.c.multiply(x.a, y.b);
  var c = this.c.multiply(x.b, y.a);
  var d = this.c.multiply(x.b, y.b);
  return new Interval(this.c.min(this.c.min(a.a, b.a),
                                 this.c.min(c.a, d.a)),
                      this.c.max(this.c.max(a.b, b.b),
                                 this.c.max(c.b, d.b)));
};
Interval.Context.prototype.divide = function (x, y) {
  var zero = this.c.fromInteger(BigInteger.BigInt(0));
  if (this.c.compareTo(y.a, zero.b) <= 0 && this.c.compareTo(y.b, zero.a) >= 0) {
    //throw new RangeError();
    return "CANNOT_DIVIDE";//TODO: FIX
  }
  if (this.c.compareTo(x.a, x.b) === 0 && this.c.compareTo(y.a, y.b) === 0) {
    return this.c.divide(x.a, y.a);
  }
  var a = this.c.divide(x.a, y.a);
  var b = this.c.divide(x.a, y.b);
  var c = this.c.divide(x.b, y.a);
  var d = this.c.divide(x.b, y.b);
  return new Interval(this.c.min(this.c.min(a.a, b.a),
                                 this.c.min(c.a, d.a)),
                      this.c.max(this.c.max(a.b, b.b),
                                 this.c.max(c.b, d.b)));
};
Interval.Context.prototype.nthRoot = function (A, n) {
  var c = this.scalingCoefficient;
  var sA = BigInteger.multiply(A, BigInteger.exponentiate(c, BigInteger.BigInt(n)));

  var x0 = nthRoot(sA, n);
  var x1 = BigInteger.lessThan(BigInteger.exponentiate(x0, BigInteger.BigInt(n)), sA) ? BigInteger.add(x0, BigInteger.BigInt(1)) : x0;

  var t = this.c.fromInteger(c);
  var a = this.c.divide(this.c.fromInteger(x0).a, t.b).a;
  var b = this.c.divide(this.c.fromInteger(x1).b, t.a).b;
  return new Interval(a, b);
};
Interval.Context.prototype.exp = function (x) {

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

  var a = exp(x.a, this.scalingCoefficient, this.scalingCoefficient);
  var b = exp(x.b, this.scalingCoefficient, this.scalingCoefficient);
  b = BigInteger.add(b, BigInteger.BigInt(1));
  return new Interval(a, b);
};
Interval.Context.prototype.pi = function () {
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
      var x = BigInteger.multiply(BigInteger.multiply(BigInteger.divide(scale, _2np1f), BigInteger.multiply(nf, nf)), BigInteger.exponentiate(2, n + 1));
      pi = BigInteger.add(pi, x);
      n += 1;
      nf = BigInteger.multiply(nf, n);
      _2np1f = BigInteger.multiply(_2np1f, BigInteger.BigInt(2 * n));
      _2np1f = BigInteger.multiply(_2np1f, BigInteger.BigInt(2 * n + 1));
    }
    return BigInteger.divide(BigInteger.multiply(pi, BigInteger.exponentiate(10, precision)), scale);
  }
  var a = calculatePI(this.precision);
  var b = BigInteger.add(a, BigInteger.BigInt(1));
  return new Interval(a, b);
};
Interval.Context.prototype._trigonometry = function (x, start) {
  function factorial(n) {
    var result = BigInteger.BigInt(1);
    for (var i = 1; i <= n; i += 1) {
      result = BigInteger.multiply(result, BigInteger.BigInt(i));
    }
    return result;
  }
  var precision = this.precision;
  //if (precision > 10) throw new Error();
  //throw new Error();
  function trigonometry(xn, xd) {
    var iterations = Math.floor(10 * (precision + 1) / 3) + 1;//TODO: ?
    //iterations += bitLength(x);//?
    //TODO: optimize
    var x = BigInteger.divide(xn, xd) + 1;
    iterations += Math.floor(Math.log(x) / Math.log(2) + 0.5) * (Math.floor(x) + 1);
    if (iterations > 380) {
      throw new TypeError("NotSupportedError");//!TODO: fix
    }
    var scale = factorial(iterations);
    scale = BigInteger.multiply(scale, BigInteger.exponentiate(xd, iterations));
    var y = 0;
    var k = start;
    while (k < iterations) {
      var s = BigInteger.multiply(BigInteger.BigInt(Math.pow(-1, (k - start) / 2)), BigInteger.divide(scale, factorial(k)));
      s = BigInteger.divide(s, BigInteger.exponentiate(xd, k));
      s = BigInteger.multiply(s, BigInteger.exponentiate(xn, k));
      y = BigInteger.add(y, s);
      k += 2;
    }
    y = BigInteger.divide(BigInteger.multiply(y, BigInteger.exponentiate(10, precision)), scale);
    var a = BigInteger.subtract(y, BigInteger.BigInt(1));
    var b = BigInteger.add(y, BigInteger.BigInt(1));
    return new Interval(a, b);
  }
  if (this.c.compareTo(this.c.subtract(x.b, x.a), this.pi().a) >= 0) {
    return new Interval(this.fromInteger(-1).a, this.fromInteger(1).a);
  }
  var a = trigonometry(x.a, this.c.scalingCoefficient, start).a;
  var b = trigonometry(x.b, this.c.scalingCoefficient, start).b;//TODO: ?
  var middle = this.c.divide(this.c.add(x.a, x.b).a, this.c.fromInteger(2).a).a;
  var c = trigonometry(middle, this.c.scalingCoefficient, start);
  if (this.c.compareTo(a, c.b) >= 0 && this.c.compareTo(c.a, b) >= 0) {
    return new Interval(b, a);
  }
  if (this.c.compareTo(a, c.b) <= 0 && this.c.compareTo(c.a, b) <= 0) {
    return new Interval(a, b);
  }
  return new Interval(this.c.fromInteger(-1).a, this.c.fromInteger(1).b);
};
Interval.Context.prototype.sin = function (x) {
  return this._trigonometry(x, 1);
};
Interval.Context.prototype.cos = function (x) {
  return this._trigonometry(x, 0);
};
Interval.Context.prototype.logarithm = function (x) {
  function factorial(n) {
    var result = BigInt(1);
    for (var i = 1; i <= n; i += 1) {
      result = result * BigInt(i);
    }
    return result;
  }
  var logarithm = function (n, scalingCoefficient, precision) { // ln(z) = ln(n / scalingCoefficient)
    if (scalingCoefficient !== BigInt(BigInteger.exponentiate(BigInt(10), BigInt(precision)))) {
      throw new RangeError();
    }
    if (n < BigInt(1)) {
      throw new TypeError("NotSupportedError");
    }
    // https://en.wikipedia.org/wiki/Logarithm#:~:text=More%20efficient%20series
    var s = BigInt(0);
    var ps = BigInt(1);
    var scale = factorial(precision);
    for (var i = BigInt(1); ps > BigInt(0); i += BigInt(2)) {
      ps = (scale / i) * BigInt(BigInteger.exponentiate(n - scalingCoefficient, i)) / BigInt(BigInteger.exponentiate(n + scalingCoefficient, i));
      s += ps;
    }
    return BigInt(2) * s * scalingCoefficient / scale;
  };
  //TODO: verify
  var a = logarithm(BigInt(x.a), BigInt(this.scalingCoefficient), this.precision);
  var b = logarithm(BigInt(x.b), BigInt(this.scalingCoefficient), this.precision);
  return new Interval(a, b + BigInt(1));
};
Interval.Context.prototype.fromInteger = function (a) {
  return this.c.fromInteger(a);
};
Interval.Context.prototype.fromIntegers = function (a, b) {
  return new Interval(this.c.fromInteger(a).a, this.c.fromInteger(b).b);
};
//?
Interval.Context.prototype.toInteger = function (x) {
  if (this.c.compareTo(x.a, x.b) === 0) {
    return {sign: this.c.sign(x.a), integer: this.c.round(x.a)};
  }
  var signA = this.c.sign(x.a);
  var signB = this.c.sign(x.b);
  if (signA === signB) {
    var candidateA = this.c.round(x.a);
    var candidateB = this.c.round(x.b);
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
    var y = evaluateExpression(e.a, context);
    if (y === "CANNOT_DIVIDE" || y == null) {
      return y;
    }
    //TODO: debug
    var yy = new Interval(context.nthRoot(y.a, n).a, context.nthRoot(y.b, n).b);
    var s = context.nthRoot(context.scalingCoefficient, n);
    yy = context.divide(yy, context.multiply(Interval.degenerate(context.scalingCoefficient), s));
    return yy;
  } else if (e instanceof Expression.BinaryOperation) {
    if (e.a === Expression.E && e.getS() === "^") {
      var b = evaluateExpression(e.b, context);
      if (b === "CANNOT_DIVIDE") {
        return b;
      }
      return context.exp(b);
    }

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
        var n = e.b.toNumber();//TODO: FIX!
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
  } else if (e === Expression.E) {
    return context.exp(context.fromInteger(BigInteger.BigInt(1)));
  } else if (e === Expression.PI) {
    return context.pi();
  } else if (e instanceof Expression.Sin || e instanceof Expression.Cos) {
    //TODO: ?
    if (e.a instanceof Expression.Radians) {
      var x = evaluateExpression(e.a.value, context);
      if (e instanceof Expression.Sin) {
        return context.sin(x);
      }
      if (e instanceof Expression.Cos) {
        return context.cos(x);
      }
    }

  } else if (e instanceof Expression.Logarithm) {
    var x = evaluateExpression(e.a, context);
    return context.logarithm(x);
  }

  return undefined;
};

var decimalToString = function (sign, number) {
  return (sign < 0 ? "-" : "") + number;
};

var complexToString = function (real, imaginary) {
  return real + (imaginary.indexOf('-') === -1 ? '+' : '') + imaginary + 'i';
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

  //? ((n * 10**(fractionDigits + 1)) ~/ d + 5) ~/ 10

var toDecimalStringInternal = function (expression, fractionDigits, decimalToStringCallback, complexToStringCallback) {
  decimalToStringCallback = decimalToStringCallback || decimalToString;
  complexToStringCallback = complexToStringCallback || complexToString;

  if (expression instanceof Expression.Division || expression instanceof Expression.Addition) {
    var numerator = expression.getNumerator();//.unwrap();
    var denominator = expression.getDenominator();//.unwrap();
    if (denominator instanceof Expression.Integer) {
      if (numerator instanceof Expression.Addition || numerator instanceof Expression.Multiplication || numerator instanceof Expression.Complex) {
        var realValue = Expression.ZERO;
        var imaginaryValue = Expression.ZERO;
        var ok = true;
        var e = numerator;
        for (var x of e.summands()) {
          var c = null;
          var r = Expression.ONE;
          for (var y of x.factors()) {
            if (c == null && y instanceof Expression.Complex) {
              c = y;
            } else if (y instanceof Expression.NthRoot || y instanceof Expression.Integer) {//TODO: ?
              r = r.multiply(y);
            } else {
              ok = false;
            }
          }
          realValue = realValue.add(r.multiply(c == null ? Expression.ONE : c.real));
          imaginaryValue = imaginaryValue.add(c != null ? r.multiply(c.imaginary) : Expression.ZERO);
        }
        if (ok && !imaginaryValue.equals(Expression.ZERO)) {
          realValue = realValue.divide(denominator);
          imaginaryValue = imaginaryValue.divide(denominator);
          var real = toDecimalStringInternal(realValue, fractionDigits, decimalToStringCallback, complexToStringCallback);
          var imaginary = toDecimalStringInternal(imaginaryValue, fractionDigits, decimalToStringCallback, complexToStringCallback);
          return complexToStringCallback(realValue.equals(Expression.ZERO) ? '' : real, imaginaryValue.equals(Expression.ONE) ? '' : imaginary);
        }
      }
    }
  }
  //TODO: remove - ?
  if (expression instanceof Expression.NthRoot) {
    var a = expression.a;//.unwrap();
    if (a instanceof Expression.Integer) {
      var A = a.value;
      var n = expression.n;
      var scale = BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(fractionDigits));
      var sA = BigInteger.multiply(A, BigInteger.exponentiate(scale, BigInteger.BigInt(n)));

      var x0 = nthRoot(sA, n);
      var x1 = BigInteger.lessThan(BigInteger.exponentiate(x0, BigInteger.BigInt(n)), sA) ? BigInteger.add(x0, BigInteger.BigInt(1)) : x0;

      // root - x0 < x1 - root
      // 2root < x0 + x1
      // 2**n * A < (x0 + x1)**n
      var nearest = BigInteger.lessThan(BigInteger.multiply(BigInteger.exponentiate(BigInteger.BigInt(2), BigInteger.BigInt(n)), sA), BigInteger.exponentiate(BigInteger.add(x0, x1), BigInteger.BigInt(n))) ? x0 : x1;
      return toDecimalStringInternal(new Expression.Division(new Expression.Integer(nearest), new Expression.Integer(scale)), fractionDigits, decimalToStringCallback, complexToStringCallback);
    }
  }
  //---
  if (!Expression.has(expression, Expression.Symbol) &&
      !Expression.has(expression, Expression.NthRoot) &&
      !Expression.has(expression, Expression.PolynomialRoot) &&
      !(expression instanceof Expression.Integer) &&
      !(expression instanceof Expression.Division) &&
      !Expression.has(expression, Expression.Radians) &&
      !Expression.has(expression, Expression.Logarithm)) {
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
