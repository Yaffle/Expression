import Expression from './Expression.js';
import primeFactor from './primeFactor.js';
import {BigDecimal, BigFloat} from './node_modules/@yaffle/bigdecimal/BigDecimal.js';

const BASE = 2;
//const BASE = 10;

function MakeMath(BigDecimal, BASE) {
  function BigDecimalMath() {
  }
  var BIG_DECIMAL_BASE = BigDecimal.round(BigDecimal.BigDecimal(BASE), {maximumSignificantDigits: 1, roundingMode: 'half-even'});
  BigDecimalMath.nextAfter = function (a, rounding) {
    if (rounding == undefined) {
      throw new RangeError();
    }
    if (rounding.roundingMode !== 'floor' && rounding.roundingMode !== 'ceil') {
      throw new RangeError();
    }
    var t = BigDecimal.round(a, rounding);
    if (!BigDecimal.equal(a, t)) {
      return t;
    }
    var _nextAfter = function (a, k, v, rounding) {
      var small = BigDecimal.multiply(BigDecimal.BigDecimal(rounding.roundingMode === 'floor' ? -1 : 1), exponentiate(BIG_DECIMAL_BASE, -k));
      var aim = BigDecimal.multiply(BigDecimal.abs(v), small);
      return BigDecimal.add(a, aim, rounding);
    };
    if (rounding.maximumFractionDigits != undefined) {
      return _nextAfter(a, rounding.maximumFractionDigits, BigDecimal.BigDecimal(1), rounding);
    }
    if (rounding.maximumSignificantDigits != undefined) {
      // a * (1 + 1 / 2**maximumSignificantDigits)
      return _nextAfter(a, rounding.maximumSignificantDigits, a, rounding);
    }
    throw new RangeError();
  };
  //TODO: remove
  //BigDecimalMath.fma = function (a, b, c, rounding) { // a * b + c
  //  return BigDecimal.round(BigDecimal.add(BigDecimal.multiply(a, b), c), rounding);
  //};
  var exponentiate = function (a, n) {
    //console.assert(a == 2);
    if (n < 0) {
      return BigDecimal.divide(BigDecimal.BigDecimal(1), exponentiate(a, -n), null);
    }
    var y = BigDecimal.BigDecimal(1);
    while (n >= 1) {
      if (0 === n % 2) {
        a = BigDecimal.multiply(a, a);
        n = n / 2;
      } else {
        y = y == undefined ? a : BigDecimal.multiply(a, y);
        n = n - 1;
      }
    }
    return y;
  };
  return BigDecimalMath;
};

const BigDecimalMath = MakeMath(BigDecimal, 10);
const BigFloatMath = MakeMath(BigFloat, 2);

//BigDecimalMath.nthRoot(BigDecimal.BigDecimal(2), 2, 3);




// https://en.wikipedia.org/wiki/Interval_arithmetic
function Interval(a, b) {
  if (a instanceof Interval || b instanceof Interval) {
    throw new TypeError();// to help with debugging
  }
  if (BigFloat.greaterThan(a, b)) {
    throw new TypeError();
  }
  this.a = a;
  this.b = b;
}

Interval.Context = function (precision, flag0) {
  //this.precision = precision;
  //this.anyRounding = {maximumSignificantDigits: precision, roundingMode: 'half-even'};
  var anyRounding = flag0 ? {maximumSignificantDigits: precision, roundingMode: 'half-even'} : {maximumFractionDigits: precision - 1, roundingMode: 'half-even'};
  this.anyRounding = anyRounding;
  this.floorRounding = Object.assign({}, anyRounding, {roundingMode: 'floor'});
  this.ceilRounding = Object.assign({}, anyRounding, {roundingMode: 'ceil'});
};
Interval.Context.prototype.unaryMinus = function (x) {
  if (BigFloat.equal(x.a, x.b)) {
    var t = BigFloat.unaryMinus(x.a);
    return new Interval(t, t);
  }
  return new Interval(BigFloat.unaryMinus(x.b), BigFloat.unaryMinus(x.a));
};
Interval.Context.prototype.add = function (x, y) {
  return new Interval(BigFloat.add(x.a, y.a, this.floorRounding), BigFloat.add(x.b, y.b, this.ceilRounding));
};
Interval.Context.prototype.subtract = function (x, y) {
  return new Interval(BigFloat.add(x.a, BigFloat.unaryMinus(y.b), this.floorRounding), BigFloat.add(x.b, BigFloat.unaryMinus(y.a), this.ceilRounding));
};
Interval.Context.prototype._multiply = function (x1, x2, y1, y2, f) {
  var sign = BigFloat.sign;
  if (sign(x1) >= 0) {
    if (sign(y1) >= 0) {
      return f(x1, y1, x2, y2);
    }
    if (sign(y2) <= 0) {
      return f(x2, y1, x1, y2);
    }
    // y1 < 0 && y2 > 0
    return f(x2, y1, x2, y2);
  }
  if (sign(x2) <= 0) {
    if (sign(y2) <= 0) {
      return f(x2, y2, x1, y1);
    }
    if (sign(y1) >= 0) {
      return f(x1, y2, x2, y1);
    }
    // y1 < 0 && y2 > 0
    return f(x1, y2, x1, y1);
  }
  if (sign(y1) >= 0) {// x1 < 0 && x2 > 0
    return f(x1, y2, x2, y2);
  }
  if (sign(y2) <= 0) {// x1 < 0 && x2 > 0
    return f(x2, y1, x1, y1);
  }
  // x1 < 0 && x2 > 0 && y1 < 0 && y2 > 0
  var interval1 = f(x1, y2, x1, y1);
  var interval2 = f(x2, y1, x2, y2);
  return new Interval(BigFloat.min(interval1.a, interval2.a),
                      BigFloat.max(interval1.b, interval2.b));
};
Interval.Context.prototype.multiply = function (x, y) {
  var floorRounding = this.floorRounding;
  var ceilRounding = this.ceilRounding;
  if (BigFloat.equal(x.a, x.b) && BigFloat.equal(y.a, y.b)) {
    var product = BigFloat.multiply(x.a, y.a);
    return new Interval(BigFloat.round(product, floorRounding), BigFloat.round(product, ceilRounding));
  }
  var f = function (a, b, c, d) {
    return new Interval(BigFloat.multiply(a, b, floorRounding), BigFloat.multiply(c, d, ceilRounding));
  };
  return this._multiply(x.a, x.b, y.a, y.b, f);
};
Interval.Context.prototype.divide = function (x, y) {
  var floorRounding = this.floorRounding;
  var ceilRounding = this.ceilRounding;
  if (BigFloat.sign(y.a) <= 0 && BigFloat.sign(y.b) >= 0) {
    if (BigFloat.equal(y.a, y.b)) {
      throw new RangeError();
    }
    return "CANNOT_DIVIDE";//TODO: FIX
  }
  //if (BigFloat.equal(x.a, x.b) && BigFloat.equal(y.a, y.b)) {
  //  var q = BigFloat.divide(x.a, y.a, this.anyRounding);
  //  var r = BigFloat.subtract(x.a, BigFloat.multiply(y.a, q));
  //  return new Interval(floorDivide(x.a, y.a), ceilDivide(x.a, y.a)); //TODO: single division
  //}
  var f = function (a, d, c, b) {
    //Note: b and d are swapped
    return new Interval(BigFloat.divide(a, b, floorRounding), BigFloat.divide(c, d, ceilRounding));
  };
  return this._multiply(x.a, x.b, y.a, y.b, f);
};
Interval.Context.prototype.sqrt = function (x) {
  if (BigFloat.sign(x.a) < 0 && BigFloat.sign(x.b) >= 0) {
    return "CANNOT_DIVIDE";//TODO: FIX
  }
  if (BigFloat.equal(x.a, x.b)) {
    var ya = BigFloat.sqrt(x.a, this.floorRounding);
    var yb = BigFloat.equal(BigFloat.multiply(ya, ya), x.b) ? ya : BigFloatMath.nextAfter(ya, this.ceilRounding);
    return new Interval(ya, yb);
  }
  return this._map(x, function (x, rounding) {
    return BigFloat.sqrt(x, rounding);
  });
};
Interval.Context.prototype.cbrt = function (x) {
  if (BigFloat.equal(x.a, x.b)) {
    var ya = BigFloat.cbrt(x.a, this.floorRounding);
    var yb = BigFloat.equal(BigFloat.multiply(ya, BigFloat.multiply(ya, ya)), x.b) ? ya : BigFloatMath.nextAfter(ya, this.ceilRounding);
    return new Interval(ya, yb);
  }
  return this._map(x, function (x, rounding) {
    return BigFloat.cbrt(x, rounding);
  });
};
/*Interval.Context.prototype.nthRoot = function (x, n) {
  if (n % 2 === 0 && BigFloat.sign(x.a) < 0 && BigFloat.sign(x.b) >= 0) {
    return "CANNOT_DIVIDE";//TODO: FIX
  }
  if (BigFloat.equal(x.a, x.b)) {
    var exponentiate = function (x, n) {
      return n === 1 ? x : (n % 2 === 0 ? exponentiate(BigFloat.multiply(x, x), n / 2) : BigFloat.multiply(x, exponentiate(x, n - 1)));
    };
    var round = function (x, rounding) {
      if (rounding.maximumFractionDigits != undefined) {
        return BigFloat.round(x, {roundingMode: 'half-even', maximumFractionDigits: Math.ceil(rounding.maximumFractionDigits / n)});
      }
      if (rounding.maximumSignificantDigits != undefined) {
        return BigFloat.round(x, {roundingMode: 'half-even', maximumSignificantDigits: Math.ceil(rounding.maximumSignificantDigits / n)});
      }
      throw new RangeError();
    };
    var ya = BigFloatMath.nthRoot(x.a, n, this.floorRounding);
    var yb = BigFloat.equal(exponentiate(round(ya, this.anyRounding), n), x.b) ? ya : BigFloatMath.nextAfter(ya, this.ceilRounding);
    return new Interval(ya, yb);
  }
  return this._map(x, function (x, rounding) {
    return BigFloatMath.nthRoot(x, n, rounding);
  });
};*/
Interval.Context.prototype.exp = function (x) {
  return this._map(x, BigFloat.exp);
};
Interval.Context.prototype.log = function (x, precision) {
  if (BigFloat.sign(x.a) <= 0 && BigFloat.sign(x.b) > 0) {
    return "CANNOT_DIVIDE";//TODO: FIX
  }
  if (BigFloat.equal(x.a, x.b) && BigFloat.equal(x.b, BigFloat.BigFloat(1))) {
    return new Interval(BigFloat.BigFloat(0), BigFloat.BigFloat(0));
  }
  return this._map(x, BigFloat.log);
};
Interval.Context.prototype.atan = function (x, precision) {
  return this._map(x, BigFloat.atan);
};
Interval.Context.prototype._mapValue = function (value, callback) {
  var c = callback(value, this.floorRounding); // TODO: ?
  //var a = BigFloatMath.nextAfter(c, this.floorRounding);
  var a = c;
  var b = c;
  if (!BigFloat.equal(value, BigFloat.BigFloat(0))) {
    b = BigFloatMath.nextAfter(c, this.ceilRounding);
  }
  return new Interval(a, b);
};
Interval.Context.prototype._map = function (x, callback) {
  if (BigFloat.equal(x.a, x.b)) {
    return this._mapValue(x.a, callback);
  }
  var a = callback(x.a, this.floorRounding);
  var b = callback(x.b, this.ceilRounding);
  return new Interval(a, b);
};
Interval.Context.prototype._trigonometry = function (x, which) {
  if (BigFloat.equal(x.a, x.b)) {
    return this._mapValue(x.a, which === 'sin' ? BigFloat.sin : BigFloat.cos);
  }
  var tau = BigFloat.multiply(BigFloat.BigFloat(8), BigFloat.atan(BigFloat.BigFloat(1), this.anyRounding));
  if (!BigFloat.lessThan(BigFloat.subtract(x.b, x.a), tau)) {
    return new Interval(BigFloat.BigFloat(-1), BigFloat.BigFloat(+1));
  }
  var f = function (x, rounding) {
    return which === 'sin' ? BigFloat.sin(x, rounding) : BigFloat.cos(x, rounding);
  };
  var middle = BigFloat.divide(BigFloat.add(x.a, x.b), BigFloat.BigFloat(2), this.anyRounding); // with rounding it works better in case the interval has huge significant digits difference
  var anyRounding = this.anyRounding;
  var extremumPoint = function (q) {
    var shift = BigFloat.multiply(BigFloat.divide(BigFloat.BigFloat(q), BigFloat.BigFloat(4), anyRounding), tau);
    var k = BigFloat.round(BigFloat.divide(BigFloat.subtract(middle, shift, anyRounding), tau, anyRounding), {
      maximumFractionDigits: 0,
      roundingMode: 'half-even'
    });
    return BigFloat.add(BigFloat.multiply(tau, k), shift);
  };
  var minimumPoint = extremumPoint(which === 'sin' ? 3 : 2);
  var maximumPoint = extremumPoint(which === 'sin' ? 1 : 0);
  var fmin = BigFloat.lessThan(minimumPoint, x.a) ? f(x.a, this.floorRounding) : (BigFloat.greaterThan(minimumPoint, x.b) ? f(x.b, this.floorRounding) : BigFloat.BigFloat(-1));
  var fmax = BigFloat.lessThan(maximumPoint, x.a) ? f(x.a, this.ceilRounding) : (BigFloat.greaterThan(maximumPoint, x.b) ? f(x.b, this.ceilRounding) : BigFloat.BigFloat(+1));
  /**/
  return new Interval(fmin, fmax);
};
Interval.Context.prototype.sin = function (x) {
  return this._trigonometry(x, 'sin');
};
Interval.Context.prototype.cos = function (x) {
  return this._trigonometry(x, 'cos');
};
Interval.Context.prototype.fromInteger = function (a) {
  if (BASE !== 2) {
    var abs = function (a) {
      return a < 0 ? BigInteger.unaryMinus(a) : a;
    };
    var k = a != 0 ? bitLength(abs(a)) - Math.ceil(Math.log2(10) * this.anyRounding.maximumSignificantDigits) : 0;
    if (k > 42) {
      //TODO: move to BigFloat.round - ?
      // for performance
      var p2k = BigInteger.exponentiate(BigInteger.BigInt(2), BigInteger.BigInt(k));
      var q = BigInteger.divide(a, p2k);
      var from = a < 0 ? BigInteger.subtract(q, BigInteger.BigInt(1)) : q;
      var to = a < 0 ? q : BigInteger.add(q, BigInteger.BigInt(1));
      return this.multiply(this.fromIntegers(from, to), this.exponentiate(this.fromInteger(BigInteger.BigInt(2)), k));
    }
  }
  //var x = BigFloat.BigFloat(a);
  //if (BigInteger.equal(BigInteger.BigInt(BigFloat.toBigInt(x)), a)) { // TODO: ?
  //  return new Interval(x, x);
  //}
  a = BigFloat.BigFloat(a);
  return new Interval(BigFloat.round(a, this.floorRounding),
                      BigFloat.round(a, this.ceilRounding));
  //return this.fromIntegers(a, a);
};
Interval.Context.prototype.fromIntegers = function (a, b) {
  var a = BigFloat.BigFloat(a);
  var b = BigFloat.BigFloat(b);
  //TODO: test case (!!!)
  console.assert(!BigFloat.lessThan(b, a));
  return new Interval(BigFloat.round(a, this.floorRounding),
                      BigFloat.round(b, this.ceilRounding));
};
Interval.Context.prototype.abs = function (x) {
  if (BigFloat.lessThan(x.a, BigFloat.BigFloat(0))) {
    if (BigFloat.lessThan(x.b, BigFloat.BigFloat(0))) {
      return new Interval(BigFloat.unaryMinus(x.b), BigFloat.unaryMinus(x.a));
    } else {
      return new Interval(BigFloat.BigFloat(0), BigFloat.max(BigFloat.unaryMinus(x.a), x.b));
    }
  }
  return x;
};
Interval.Context.prototype.exponentiate = function (x, n) {
  var y = undefined;
  while (n >= 1) {
    if (n === 2 * Math.floor(n / 2)) {
      x = this.multiply(x, x);
      n = Math.floor(n / 2);
    } else {
      y = y == undefined ? x : this.multiply(x, y);
      n -= 1;
    }
  }
  return y;
};

//?

      // todo: for exact notation (?):
      // 10000 -> 10**4
      // 15000 -> 15*10**3
      // 0.00015 -> 15*10**(-5)



Interval.Context.prototype.formatToDecimal = function (x, rounding) {
  // assume, that the value is not exact
  var signA = BigFloat.sign(x.a);
  var signB = BigFloat.sign(x.b);
  var sign = (signA || signB) === (signB || signA) ? (signA || signB) : 0;
  x = this.abs(x);
  var stringify = function (a, roundingMode) {
    if (rounding.fractionDigits != undefined) {
      return a.toFixed(rounding.fractionDigits, roundingMode);
    }
    return a.toPrecision(rounding.significantDigits, roundingMode);
  };
  var a = stringify(x.a, "half-up");
  var b = BigFloat.equal(x.a, x.b) ? a : stringify(x.b, "half-down");
  var isZero = function (a) {
    return !/[^0\.]/.test(a);
  };
  if (a === b && (sign !== 0 || isZero(a) && isZero(b))) {
    return (sign < 0 ? '-' : (sign > 0 && isZero(a) && isZero(b) ? '+' : '')) + a;
  }
  return undefined;
};
Interval.prototype.toString = function () {
  return "[" + this.a.toString() + ";" + this.b.toString() + "]";
};

var calcAt = function (polynomial, x, context) {
  var result = evaluateExpression(Expression.ZERO, context);
  for (var i = polynomial.getDegree(); i >= 0; i--) {
    result = context.multiply(result, x);
    var tmp = evaluateExpression(polynomial.getCoefficient(i), context);
    if (tmp === "CANNOT_DIVIDE" || tmp == undefined) {
      return tmp;
    }
    result = context.add(result, tmp);
  }
  return result;
};

var evaluateExpression = function (e, context) {
  if (e instanceof Expression.Integer) {
    var n = e.value;
    return context.fromInteger(n);
  } else if (e instanceof Expression.NthRoot) {
    var a = e.a;
    var n = e.n;
    var y = evaluateExpression(a, context);
    if (y === "CANNOT_DIVIDE" || y == undefined) {
      return y;
    }
    if (n == 2) {
      return context.sqrt(y);
    }
    if (n == 3) {
      return context.cbrt(y);
    }
    return context.exp(context.divide(context.log(y), context.fromInteger(n)));
  } else if (e instanceof Expression.BinaryOperation) {
    // slow for some cases:
    if (e instanceof Expression.Addition && Expression.has(e, Expression.PolynomialRootSymbol)) {
      var root = Expression.getVariable(e);//?
      var p = Polynomial.toPolynomial(e, root);
      if (p.hasIntegerCoefficients()) {// trying to avoid slow cases (?)
        //TODO: https://en.wikipedia.org/wiki/Horner%27s_method - ?
        var zero = evaluateExpression(root, context);
        //return evaluateExpression(p.calcAt(), context);
        return calcAt(p, zero, context);
      }
    }
    if (e.a === Expression.E && e.getS() === "^") {
      var b = evaluateExpression(e.b, context);
      if (b === "CANNOT_DIVIDE" || b == undefined) {
        return b;
      }
      return context.exp(b);
    }

    var a = evaluateExpression(e.a, context);
    if (a === "CANNOT_DIVIDE" || a == undefined) {
      return a;
    }
    var b = evaluateExpression(e.b, context);
    if (b === "CANNOT_DIVIDE" || b == undefined) {
      return b;
    }
    var operator = e.getS();
    if (operator === "+") {
      return context.add(a, b);
    } else if (operator === "-") {
      return context.subtract(a, b);
    } else if (operator === "*") {
      return context.multiply(a, b);
    } else if (operator === "/") {
      return context.divide(a, b);
    } else if (operator === "^") { // Expression.PolynomialRootSymbol^3, pi^2, 2**(sqrt(3)), (log(2))^2
      //TODO: to polynomial
      if (!(e.b instanceof Expression.Integer) || e.b.toNumber() <= 0 || e.b.toNumber() > Number.MAX_SAFE_INTEGER) {
        //throw new TypeError();
        var log = context.log(a);
        if (log === "CANNOT_DIVIDE") {
          return log;
        }
        return context.exp(context.multiply(log, b));
      }
      var n = e.b.toNumber();//TODO: FIX!
      return context.exponentiate(a, n);
    }
  } else if (e instanceof Expression.PolynomialRootSymbol) {
    var i = e.toDecimal(context.anyRounding.maximumSignificantDigits || context.anyRounding.maximumFractionDigits);
    // "lcm" is too slow to compute (?)
    /*if (true) {
      var a = BigFloat.divide(BigFloat.BigFloat(i.a.getNumerator().value), BigFloat.BigFloat(i.a.getDenominator().value), context.floorRounding);
      var b = BigFloat.divide(BigFloat.BigFloat(i.b.getNumerator().value), BigFloat.BigFloat(i.b.getDenominator().value), context.ceilRounding);
      return new Interval(a, b);
    }*/
    return context.divide(context.fromIntegers(i.b.getDenominator().multiply(i.a.getNumerator()).value,
                                               i.a.getDenominator().multiply(i.b.getNumerator()).value),
                          context.fromInteger(i.a.getDenominator().multiply(i.b.getDenominator()).value));
  } else if (e === Expression.E) {
    return context.exp(context.fromInteger(1));
  } else if (e === Expression.PI) {
    return context.multiply(context.fromInteger(4), context.atan(context.fromInteger(1)));
  } else if (e instanceof Expression.Function) {
    var x = evaluateExpression((e instanceof Expression.Sin || e instanceof Expression.Cos) && e.a instanceof Expression.Radians ? e.a.value : e.a, context);
    if (x === "CANNOT_DIVIDE" || x == undefined) {
      return x;
    }
    if (e instanceof Expression.Sin) {
      return context.sin(x);
    }
    if (e instanceof Expression.Cos) {
      return context.cos(x);
    }
    if (e instanceof Expression.Logarithm) {
      return context.log(x);
    }
    if (e instanceof Expression.Arctan) {
      return context.atan(x);
    }
  } else if (e instanceof Expression.ExpressionWithPolynomialRoot) {
    return evaluateExpression(e.e, context);
  } else if (e instanceof Expression.ExpressionPolynomialRoot) {
    var i = e.root.toDecimal(context.anyRounding.maximumSignificantDigits || context.anyRounding.maximumFractionDigits);
    // "lcm" is too slow to compute (?)
    /*if (true) {
      var a = BigFloat.divide(BigFloat.BigFloat(i.a.getNumerator().value), BigFloat.BigFloat(i.a.getDenominator().value), context.floorRounding);
      var b = BigFloat.divide(BigFloat.BigFloat(i.b.getNumerator().value), BigFloat.BigFloat(i.b.getDenominator().value), context.ceilRounding);
      return new Interval(a, b);
    }*/
    var root = context.divide(context.fromIntegers(i.b.getDenominator().multiply(i.a.getNumerator()).value,
                                                   i.a.getDenominator().multiply(i.b.getNumerator()).value),
                              context.fromInteger(i.a.getDenominator().multiply(i.b.getDenominator()).value));
    return root;
  }

  return undefined;
};

var decimalToString = function (decimal) {
  return decimal.replace(/[eE]/g, '*10^');
};

var complexToString = function (real, imaginary) {
  return real + (/^[\-\+]/.test(imaginary) ? imaginary.replace(/^([\-\+])[\s\S]+/g, '$1') : (real !== '' ? '+' : '')) + (imaginary !== '1' && imaginary !== '-1' ? imaginary.replace(/^[\-\+]/g, '') + '*' + 'i' : 'i');
};

  //? ((n * 10**(fractionDigits + 1)) ~/ d + 5) ~/ 10

var toDecimalStringInternal = function (expression, rounding, decimalToStringCallback, complexToStringCallback) {
  decimalToStringCallback = decimalToStringCallback || decimalToString;
  complexToStringCallback = complexToStringCallback || complexToString;
  if (rounding.fractionDigits == undefined && rounding.significantDigits == undefined ||
      rounding.fractionDigits != undefined && rounding.significantDigits != undefined) {//?
    throw new RangeError();
  }
  if (rounding.fractionDigits != undefined && (rounding.fractionDigits < 0 || rounding.fractionDigits > Number.MAX_SAFE_INTEGER) ||
      rounding.significantDigits != undefined && (rounding.significantDigits < 1 || rounding.significantDigits > Number.MAX_SAFE_INTEGER) ||
      rounding.roundingMode != undefined) {
    throw new RangeError();
  }

  if (expression instanceof Expression.Complex || Expression.has(expression, Expression.Complex)) {//?TODO: ?
    var numerator = expression.getNumerator();//.unwrap();
    var denominator = expression.getDenominator();//.unwrap();
    if (denominator instanceof Expression.Integer ||
        Expression.has(denominator, Expression.PolynomialRootSymbol) ||
        Expression.has(denominator, Expression.ExpressionPolynomialRoot) ||
        Expression._isPositive(denominator)) { // e^2
      if (numerator instanceof Expression.Addition || numerator instanceof Expression.Multiplication || numerator instanceof Expression.Complex) {
        const tmp = Expression.getComplexNumberParts(numerator);
        let realValue = tmp.real;
        let imaginaryValue = tmp.imaginary;
        if (!imaginaryValue.equals(Expression.ZERO)) {
          realValue = realValue.divide(denominator);
          imaginaryValue = imaginaryValue.divide(denominator);
          var real = realValue.equals(Expression.ZERO) ? '' : toDecimalStringInternal(realValue, rounding, decimalToStringCallback, complexToStringCallback);
          var imaginary = toDecimalStringInternal(imaginaryValue, rounding, decimalToStringCallback, complexToStringCallback);
          return complexToStringCallback(real, imaginary);
        }
      }
    }
  }
   
  if (expression instanceof Expression.Integer || expression instanceof Expression.Division && expression.a instanceof Expression.Integer && expression.b instanceof Expression.Integer) {
    //TODO: ?
    if (true) {//TODO: ?
      return decimalToStringCallback(primeFactor._rationalNumberToDecimalString(expression.getNumerator().toBigInt(), expression.getDenominator().toBigInt(), rounding));
    }
  }
  //TODO: remove - ?
  /*TODO: enable
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
      //return toDecimalStringInternal(new Expression.Division(new Expression.Integer(nearest), new Expression.Integer(scale)), fractionDigits, decimalToStringCallback, complexToStringCallback);
      var context = new Interval.Context(fractionDigits + 1);
      var a = BigFloat.divide(BigFloat.BigFloat(nearest), BigFloat.BigFloat(scale));
      var result = context.formatToDecimal(new Interval(a, BigFloatMath.nextAfter(a, {maximumFractionDigits: fractionDigits + 1, roundingMode: 'ceil'})), fractionDigits);
      return decimalToStringCallback(result);
    }
  }
  */
  //---
  /*if (!Expression.has(expression, Expression.Symbol) &&
      !Expression.has(expression, Expression.NthRoot) &&
      !Expression.has(expression, Expression.PolynomialRootSymbol) &&
      !Expression.has(expression, Expression.ExpressionPolynomialRoot) &&
      !(expression instanceof Expression.Integer) &&
      !(expression instanceof Expression.Division) &&
      !Expression.has(expression, Expression.Radians) &&
      !Expression.has(expression, Expression.Logarithm) &&
      !Expression.has(expression, Expression.Exponentiation)) {
    throw new TypeError("toDecimalString:" + JSON.stringify(rounding) + ":" + expression.toString({}));
  }*/
  console.assert(BASE % 2 === 0);
  var result = undefined;
  var guessedPrecision = 1; //TODO: ?
  //!new
  //TODO: !!!
  //if (expression instanceof Expression.Multiplication && expression.a instanceof Expression.Integer && rounding.fractionDigits != undefined) {
    //guessedPrecision = 2 * bitLength(expression.a.toBigInt());
    //guessedPrecision  = 128;
  //}
  //!
  var flag0 = Expression.has(expression, Expression.Function) || Expression.has(expression, Expression.Exponentiation);
  while (result == undefined) {
    if (guessedPrecision > 60000 && guessedPrecision > (rounding.fractionDigits || rounding.significantDigits) * 4 * Math.log2(10)) {
      debugger;
      throw new TypeError();
    }
    //if (guessedPrecision > 1024) throw new Error();
    var context = new Interval.Context(guessedPrecision, flag0);
    var x = evaluateExpression(expression, context);
    if (x == undefined) {
      return undefined;
    }
    if (x !== "CANNOT_DIVIDE") { // continue the loop otherwise
      result = context.formatToDecimal(x, rounding);
    }
    if (guessedPrecision > 1  && result == undefined) {
      //console.count('guessedPrecision > 1  && result == undefined');
    }
    if (x !== "CANNOT_DIVIDE" && result == undefined && rounding.fractionDigits != undefined && guessedPrecision === 1) {//TODO: ?
      var log10OfValue = BigFloat.max(BigFloat.abs(x.a), BigFloat.abs(x.b)).toFixed(0).length;
      guessedPrecision = Math.ceil(Math.max(log10OfValue * Math.log2(10), 2) / 2 + Math.log2(10) * rounding.fractionDigits / 2 * 2);
    }
    if (x !== "CANNOT_DIVIDE" && result == undefined && rounding.significantDigits != undefined && guessedPrecision === 1) {//TODO: ?
      if (BigFloat.sign(x.a) === BigFloat.sign(x.b)) { // zero is not part of the interval, so we know the minimal value
        var tmp = BigFloat.log(BigFloat.min(BigFloat.abs(x.a), BigFloat.abs(x.b)), {maximumSignificantDigits: 1, roundingMode: 'half-even'}).toFixed(0).length;
        guessedPrecision = Math.max(Math.ceil(tmp * Math.log2(10) / 2), 1);
      }
    }
    guessedPrecision *= 2;
  }
  if (guessedPrecision !== 256) {
    //console.log(guessedPrecision);
  }
  return decimalToStringCallback(result);
};

export default toDecimalStringInternal;

toDecimalStringInternal.testables = {
  BigDecimalMath: BigDecimalMath,
  BigDecimal: BigDecimal,
  BigFloat: BigFloat,
  Interval: Interval
};
