import Expression from './Expression.js';
import BigInteger from './BigInteger.js';
import nthRoot from './nthRoot.compiled.js';
//import BigDecimal from './BigDecimal.js';

import BigDecimal from './BigDecimal/BigDecimal.compiled.js';

//TODO: ???

//TODO: remove roundingMode - ?

function BigDecimalMath() {
}
BigDecimalMath.sign = function (a) { // returns -1, 0, or +1
  return BigDecimal.lessThan(a, BigDecimal.BigDecimal(0)) ? -1 : (BigDecimal.greaterThan(a, BigDecimal.BigDecimal(0)) ? +1 : 0);
};
BigDecimalMath.abs = function (a) {
  return BigDecimal.lessThan(a, BigDecimal.BigDecimal(0)) ? BigDecimal.unaryMinus(a) : (BigDecimal.greaterThan(a, BigDecimal.BigDecimal(0)) ? a : BigDecimal.BigDecimal(0));
};
BigDecimalMath.nthRoot = function (x, n, rounding) {
  if ((true || n > 100) && rounding.maximumSignificantDigits != undefined && rounding.maximumSignificantDigits < n) {//?
    // RPN('2**(1/10000)').toMathML({fractionDigits: 10});
    var internalRounding = {
      maximumSignificantDigits: rounding.maximumSignificantDigits + 1 + Math.max(Math.ceil(Math.log(BigDecimal.log(BigDecimal.round(x, {maximumSignificantDigits: 1, roundingMode: 'ceil'}), {maximumSignificantDigits: 1, roundingMode: 'half-even'}) / n) / Math.log(10)), 0),
      roundingMode: "half-even"
    };
    return BigDecimal.exp(BigDecimal.divide(BigDecimal.log(x, internalRounding), BigDecimal.BigDecimal(n), internalRounding), internalRounding);
  }
  // (f * 10**e)**(1/n) = (f * 10**(e - round(e / n) * n))**(1/n) * 10**round(e / n)
  var toFixedPoint = function (a, scalingCoefficient) {
    return BigDecimal.toBigInt(BigDecimal.multiply(a, BigDecimal.BigDecimal(scalingCoefficient)));
  };
  var toBigDecimal = function (x, scalingCoefficient, rounding) {
    return BigDecimal.divide(BigDecimal.BigDecimal(x), BigDecimal.BigDecimal(scalingCoefficient), rounding);
  };
  function calculateNthRoot(x, n, scalingCoefficient) {
    var sx = BigInteger.multiply(x, BigInteger.exponentiate(scalingCoefficient, BigInteger.BigInt(n - 1)));
    var x0 = nthRoot(sx, n);
    //var x1 = BigInteger.lessThan(BigInteger.exponentiate(x1, BigInteger.BigInt(n)), sA) ? BigInteger.add(x1, BigInteger.BigInt(1)) : x1;
    return x0;
  }
  var scalingCoefficient = BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(rounding.maximumSignificantDigits || rounding.maximumFractionDigits));
  return toBigDecimal(calculateNthRoot(toFixedPoint(x, scalingCoefficient), n, scalingCoefficient), scalingCoefficient, rounding);
};
BigDecimalMath.pi = function (rounding) {
  return BigDecimal.multiply(BigDecimal.BigDecimal(4), BigDecimal.atan(BigDecimal.BigDecimal(1), rounding));
};
BigDecimalMath._nextToward = function (a, rounding, which) {
  // rounding.roundingMode is ignored
  console.assert(rounding != undefined);
  console.assert(which === 'floor' || which === 'ceil');
  if (BigDecimal.lessThan(a, BigDecimal.BigDecimal(0))) {
    return BigDecimal.unaryMinus(BigDecimalMath._nextToward(BigDecimal.unaryMinus(a), rounding, which === 'floor' ? 'ceil' : 'floor'));
  }
  if (rounding.maximumFractionDigits != undefined) {
    var scalingCoefficient = BigDecimal.BigDecimal(BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(rounding.maximumFractionDigits)));
    var small = BigDecimal.divide(BigDecimal.BigDecimal(which === 'floor' ? -1 : 1), scalingCoefficient);
    return BigDecimal.add(a, small, {maximumFractionDigits: rounding.maximumFractionDigits, roundingMode: which === 'floor' ? 'ceil' : 'floor'});
  }
  if (BigDecimal.equal(a, BigDecimal.BigDecimal(0))) {
    return a;
  }
  var t = BigDecimal.round(a, {maximumSignificantDigits: rounding.maximumSignificantDigits, roundingMode: which === 'floor' ? 'floor' : 'ceil'});
  if (!BigDecimal.equal(t, a)) {
    return t;
  }
  var m = BigDecimal.BigDecimal(BigInteger.exponentiate(BigInteger.BigInt(10), BigInteger.BigInt(rounding.maximumSignificantDigits)));
  var c = BigDecimal.add(BigDecimal.BigDecimal(1), BigDecimal.divide(BigDecimal.BigDecimal(which === 'floor' ? -1 : 1), m));
  return BigDecimal.multiply(a, c, {maximumSignificantDigits: rounding.maximumSignificantDigits, roundingMode: which === 'floor' ? 'floor' : 'ceil'});
};
BigDecimalMath.nextFloor = function (a, rounding) {
  return BigDecimalMath._nextToward(a, rounding, 'floor');
};
BigDecimalMath.nextCeil = function (a, rounding) {
  return BigDecimalMath._nextToward(a, rounding, 'ceil');
};

//BigDecimalMath.nthRoot(BigDecimal.BigDecimal(2), 2, 3);

if (false) {
  console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0.995'), {maximumSignificantDigits: 4}).toString() === '0.9951');
  console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0.995'), {maximumSignificantDigits: 3}).toString() === '0.996');
  console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0.995'), {maximumSignificantDigits: 2}).toString() === '1');
  console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0'), {maximumSignificantDigits: 1}).toString() === '0');
  console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0.995'), {maximumFractionDigits: 4}).toString() === '0.9951');
  console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0.995'), {maximumFractionDigits: 3}).toString() === '0.996');
  console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0.995'), {maximumFractionDigits: 2}).toString() === '1');
  console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0'), {maximumFractionDigits: 1}).toString() === '0.1');
  console.assert(BigDecimalMath.nextFloor(BigDecimal.BigDecimal('0'), {maximumFractionDigits: 1}).toString() === '-0.1');
  console.time();
  console.assert(BigDecimalMath.nextFloor(BigDecimal.BigDecimal('1e10000000'), {maximumSignificantDigits: 100000}).toString() === '?');
  console.timeEnd();
  // default: 2289.925048828125ms
  // default: 1335.22216796875ms
}

// https://en.wikipedia.org/wiki/Interval_arithmetic
function Interval(a, b) {
  if (a instanceof Interval || b instanceof Interval) {
    throw new TypeError();// to help with debugging
  }
  if (BigDecimal.greaterThan(a, b)) {
    throw new TypeError();
  }
  this.a = a;
  this.b = b;
}
Interval.degenerate = function (a) {
  return new Interval(a, a);
};

Interval.Context = function (precision) {
  this.precision = precision;
  this.floorRounding = {maximumSignificantDigits: precision, roundingMode: 'floor'};
  this.ceilRounding = {maximumSignificantDigits: precision, roundingMode: 'ceil'};
  this.anyRounding = {maximumSignificantDigits: precision, roundingMode: 'half-even'};
};
Interval.Context.prototype.unaryMinus = function (x) {
  if (BigDecimal.equal(x.a, x.b)) {
    return Interval.degenerate(BigDecimal.unaryMinus(x.a));
  }
  return new Interval(BigDecimal.unaryMinus(x.b), BigDecimal.unaryMinus(x.a));
};
Interval.Context.prototype.add = function (x, y) {
  return new Interval(BigDecimal.add(x.a, y.a, this.floorRounding), BigDecimal.add(x.b, y.b, this.ceilRounding));
};
Interval.Context.prototype.subtract = function (x, y) {
  return new Interval(BigDecimal.subtract(x.a, y.a, this.floorRounding), BigDecimal.subtract(x.b, y.b, this.ceilRounding));
};
Interval.Context.prototype.multiply = function (x, y) {
  var floorRounding = this.floorRounding;
  var ceilRounding = this.ceilRounding;
  var floorMultiply = function (x, y) {
    //TODO: test for wrong roundingMode
    return BigDecimal.multiply(x, y, floorRounding);
  };
  var ceilMultiply = function (x, y) {
    //TODO: test for wrong roundingMode
    return BigDecimal.multiply(x, y, ceilRounding);
  };
  if (BigDecimal.equal(x.a, x.b) && BigDecimal.equal(y.a, y.b)) {
    return new Interval(floorMultiply(x.a, y.a), ceilMultiply(x.a, y.a)); //TODO: single multiplication
  }
  var x1 = x.a;
  var x2 = x.b;
  var y1 = y.a;
  var y2 = y.b;
  if (BigDecimalMath.sign(x1) >= 0) {
    if (BigDecimalMath.sign(y1) >= 0) {
      return new Interval(floorMultiply(x1, y1), ceilMultiply(x2, y2));
    }
    if (BigDecimalMath.sign(y2) <= 0) {
      return new Interval(floorMultiply(x2, y1), ceilMultiply(x1, y2));
    }
    // y1 < 0 && y2 > 0
    return new Interval(floorMultiply(x2, y1), ceilMultiply(x2, y2));
  }
  if (BigDecimalMath.sign(x2) <= 0) {
    if (BigDecimalMath.sign(y2) <= 0) {
      return new Interval(floorMultiply(x2, y2), ceilMultiply(x1, y1));
    }
    if (BigDecimalMath.sign(y1) >= 0) {
      return new Interval(floorMultiply(x1, y2), ceilMultiply(x2, y1));
    }
    // y1 < 0 && y2 > 0
    return new Interval(floorMultiply(x1, y2), ceilMultiply(x1, y1));
  }
  if (BigDecimalMath.sign(y1) >= 0) {// x1 < 0 && x2 > 0
    return new Interval(floorMultiply(x1, y2), ceilMultiply(x2, y2));
  }
  if (BigDecimalMath.sign(y2) <= 0) {// x1 < 0 && x2 > 0
    return new Interval(floorMultiply(x2, y1), ceilMultiply(x1, y1));
  }
  var min = function (a, b) {
    return BigDecimal.greaterThan(a, b) ? b : a;
  };
  var max = function (a, b) {
    return BigDecimal.lessThan(a, b) ? b : a;
  };
  // x1 < 0 && x2 > 0 && y1 < 0 && y2 > 0
  return new Interval(min(floorMultiply(x1, y2), floorMultiply(x2, y1)),
                      max(ceilMultiply(x1, y1), ceilMultiply(x2, y2)));
};
Interval.Context.prototype.divide = function (x, y) {
  var floorRounding = this.floorRounding;
  var ceilRounding = this.ceilRounding;
  var floorDivide = function (x, y) {
    return BigDecimal.divide(x, y, floorRounding);
  };
  var ceilDivide = function (x, y) {
    return BigDecimal.divide(x, y, ceilRounding);
  };
  if (BigDecimalMath.sign(y.a) <= 0 && BigDecimalMath.sign(y.b) >= 0) {
    //throw new RangeError();
    return "CANNOT_DIVIDE";//TODO: FIX
  }
  if (BigDecimal.equal(x.a, x.b) && BigDecimal.equal(y.a, y.b)) {
    return new Interval(floorDivide(x.a, y.a), ceilDivide(x.a, y.a)); //TODO: single division
  }
  var x1 = x.a;
  var x2 = x.b;
  var y1 = y.a;
  var y2 = y.b;
  if (BigDecimalMath.sign(x1) >= 0) {
    if (BigDecimalMath.sign(y1) >= 0) {
      return new Interval(floorDivide(x1, y2), ceilDivide(x2, y1));
    }
    if (BigDecimalMath.sign(y2) <= 0) {
      return new Interval(floorDivide(x2, y2), ceilDivide(x1, y1));
    }
    // y1 < 0 && y2 > 0
    return "CANNOT_DIVIDE";//TODO: FIX
  }
  if (BigDecimalMath.sign(x2) <= 0) {
    if (BigDecimalMath.sign(y2) <= 0) {
      return new Interval(floorDivide(x2, y1), ceilDivide(x1, y2));
    }
    if (BigDecimalMath.sign(y1) >= 0) {
      return new Interval(floorDivide(x1, y1), ceilDivide(x2, y2));
    }
    // y1 < 0 && y2 > 0
    return "CANNOT_DIVIDE";//TODO: FIX
  }
  if (BigDecimalMath.sign(y1) >= 0) {// x1 < 0 && x2 > 0
    return new Interval(floorDivide(x1, y1), ceilDivide(x2, y1));
  }
  if (BigDecimalMath.sign(y2) <= 0) {// x1 < 0 && x2 > 0
    return new Interval(floorDivide(x2, y2), ceilDivide(x1, y2));
  }
  // x1 < 0 && x2 > 0 && y1 < 0 && y2 > 0
  return "CANNOT_DIVIDE";//TODO: FIX
};
Interval.Context.prototype.nthRoot = function (x, n) {
  if (BigDecimalMath.sign(x.a) < 0 && BigDecimalMath.sign(x.b) >= 0) {
    return "CANNOT_DIVIDE";//TODO: FIX
  }
  if (BigDecimal.equal(x.a, x.b)) {
    var c = BigDecimalMath.nthRoot(x.a, n, this.anyRounding);
    var a = BigDecimalMath.nextFloor(c, this.anyRounding);
    var b = BigDecimalMath.nextCeil(c, this.anyRounding);
    return new Interval(a, b);
  }
  var a = BigDecimalMath.nextFloor(BigDecimalMath.nthRoot(x.a, n, this.anyRounding), this.anyRounding);
  var b = BigDecimalMath.nextCeil(BigDecimalMath.nthRoot(x.b, n, this.anyRounding), this.anyRounding);
  return new Interval(a, b);
};
Interval.Context.prototype.exp = function (x) {
  if (BigDecimal.equal(x.a, x.b)) {
    var c = BigDecimal.exp(x.a, this.anyRounding);
    var a = BigDecimalMath.nextFloor(c, this.anyRounding);
    var b = BigDecimalMath.nextCeil(c, this.anyRounding);
    return new Interval(a, b);
  }
  var a = BigDecimalMath.nextFloor(BigDecimal.exp(x.a, this.anyRounding), this.anyRounding);
  var b = BigDecimalMath.nextCeil(BigDecimal.exp(x.b, this.anyRounding), this.anyRounding);
  return new Interval(a, b);
};
Interval.Context.prototype.log = function (x, precision) {
  if (BigDecimalMath.sign(x.a) < 0 && BigDecimalMath.sign(x.b) >= 0) {
    return "CANNOT_DIVIDE";//TODO: FIX
  }
  if (BigDecimal.equal(x.a, x.b)) {
    var c = BigDecimal.log(x.a, this.anyRounding);
    var a = BigDecimalMath.nextFloor(c, this.anyRounding);
    var b = BigDecimalMath.nextCeil(c, this.anyRounding);
    return new Interval(a, b);
  }
  //TODO: test coverage
  var a = BigDecimalMath.nextFloor(BigDecimal.log(x.a, this.anyRounding), this.anyRounding);
  var b = BigDecimalMath.nextCeil(BigDecimal.log(x.b, this.anyRounding), this.anyRounding);
  return new Interval(a, b);
};
Interval.Context.prototype.pi = function () {
  var c = BigDecimalMath.pi(this.anyRounding);
  var a = BigDecimalMath.nextFloor(c, this.anyRounding);
  var b = BigDecimalMath.nextCeil(c, this.anyRounding);
  return new Interval(a, b);
};
Interval.Context.prototype._trigonometry = function (x, which) {
  var min = function (a, b) {
    return BigDecimal.greaterThan(a, b) ? b : a;
  };
  var max = function (a, b) {
    return BigDecimal.lessThan(a, b) ? b : a;
  };
  var clamp = function (interval) {
    var minusOne = BigDecimal.BigDecimal(-1);
    var one = BigDecimal.BigDecimal(+1);
    return new Interval(min(max(interval.a, minusOne), one), min(max(interval.b, minusOne), one));
  };
  var anyRounding = this.anyRounding;
  var f = function (x) {
    return which === 'sin' ? BigDecimal.sin(x, anyRounding) : BigDecimal.cos(x, anyRounding);
  };
  if (BigDecimal.equal(x.a, x.b)) {
    var fc = f(x.a);
    return clamp(new Interval(BigDecimalMath.nextFloor(fc, this.anyRounding), BigDecimalMath.nextCeil(fc, this.anyRounding)));
  }
  var pi = BigDecimalMath.pi(this.anyRounding);
  var tau = BigDecimal.multiply(BigDecimal.BigDecimal(2), pi);
  if (!BigDecimal.lessThan(BigDecimal.subtract(x.b, x.a), tau)) {
    return new Interval(BigDecimal.BigDecimal(-1), BigDecimal.BigDecimal(+1));
  }
  var middle = BigDecimal.divide(BigDecimal.add(x.a, x.b), BigDecimal.BigDecimal(2));
  var extremumPoint = function (q) {
    var shift = BigDecimal.multiply(BigDecimal.divide(BigDecimal.BigDecimal(q), BigDecimal.BigDecimal(2)), pi);
    var k = BigDecimal.round(BigDecimal.divide(BigDecimal.subtract(middle, shift), tau, anyRounding), {
      maximumFractionDigits: 0,
      roundingMode: 'half-even'
    });
    return BigDecimal.add(BigDecimal.multiply(tau, k), shift);
  };
  var minimumPoint = extremumPoint(which === 'sin' ? 3 : 2);
  var maximumPoint = extremumPoint(which === 'sin' ? 1 : 0);
  var fmin = BigDecimal.lessThan(minimumPoint, x.a) ? f(x.a) : (BigDecimal.greaterThan(minimumPoint, x.b) ? f(x.b) : BigDecimal.BigDecimal(-1));
  var fmax = BigDecimal.lessThan(maximumPoint, x.a) ? f(x.a) : (BigDecimal.greaterThan(maximumPoint, x.b) ? f(x.b) : BigDecimal.BigDecimal(+1));
  /**/
  return clamp(new Interval(BigDecimalMath.nextFloor(fmin, this.anyRounding), BigDecimalMath.nextCeil(fmax, this.anyRounding)));
};
Interval.Context.prototype.sin = function (x) {
  return this._trigonometry(x, 'sin');
};
Interval.Context.prototype.cos = function (x) {
  return this._trigonometry(x, 'cos');
};
Interval.Context.prototype.fromInteger = function (a) {
  var x = BigDecimal.BigDecimal(a);
  if (BigInteger.equal(BigDecimal.toBigInt(x), a)) { // TODO: ?
    return Interval.degenerate(x);
  }
  return new Interval(BigDecimal.round(x, this.floorRounding),
                      BigDecimal.round(x, this.ceilRounding));
};
Interval.Context.prototype.fromIntegers = function (a, b) {
  //TODO: test case (!!!)
  console.assert(!BigInteger.greaterThan(a, b));
  var a = BigDecimal.BigDecimal(a);
  var b = BigDecimal.BigDecimal(b);
  return new Interval(BigDecimal.round(a, this.floorRounding),
                      BigDecimal.round(b, this.ceilRounding));
};

//?

  var bigDecimalToPlainString = function (significand, exponent, minFraction, minSignificant) {
    var e = Number(exponent) + 1;
    var zeros = Math.max(0, Math.max(e, minSignificant) - significand.length);
    if (e <= 0) {
      significand = "0".repeat(1 - e) + significand;
      e = 1;
    }
    significand += "0".repeat(zeros);
    significand += "0".repeat(Math.max(minFraction - (significand.length - e), 0));
    return significand.slice(0, e) + (significand.length > e ? "." + significand.slice(e) : "");
  };
  var digitsToDecimalNumber = function (tmp, rounding) {
    var sign = tmp.sign;
    var decimal = tmp.decimal;
    var exact = tmp.exact;
    var match = /^(\d+)?\.?(\d+)?[eE]?([\-\+]?\d+)?$/.exec(decimal);
    // significand * 10**exponent
    var significand = (match[1] || "") + (match[2] || "");
    var z = significand.length;
    significand = significand.replace(/^0+/g, '');
    z -= significand.length;
    var exponent = BigInteger.add(BigInteger.BigInt((match[3] || "")), BigInteger.BigInt((match[1] || "").length - 1 - z));
    significand = significand.replace(/0+$/g, '');

    var plusSign = sign > 0 && !exact && significand === '' ? '+' : '';
    var minusSign = sign < 0 ? '-' : '';
    // Something like Number#toPrecision: when value is between 10**-6 and 10**p? - to fixed, otherwise - to exponential:
    if (rounding.significantDigits != undefined) {
      var e = BigInteger.toNumber(exponent);
      if (e < -6 || e >= rounding.significantDigits) {
        //exponent = BigInteger.subtract(exponent, BigInteger.BigInt(1));
        exponent = e < 0 ? BigInteger.unaryMinus(exponent) : exponent;
        // todo: for exact notation (?):
        // 10000 -> 10**4
        // 15000 -> 15*10**3
        // 0.00015 -> 15*10**(-5)
        var s = bigDecimalToPlainString(significand, 0, 0, exact ? 0 : rounding.significantDigits);
        return {
          plusSign: plusSign,
          minusSign: minusSign,
          number: s,
          exponentMinusSign: e < 0 ? '-' : '',
          exponentInteger: exponent.toString()
        };
      }
      var number = bigDecimalToPlainString(significand, exponent, 0, exact ? 0 : rounding.significantDigits);
      return {plusSign: plusSign, minusSign: minusSign, number: number, exponentMinusSign: '', exponentInteger: ''};
    }
    var number = bigDecimalToPlainString(significand, exponent, exact ? 0 : rounding.fractionDigits, 0);
    //var number = (maximumFractionDigits === 0 ? digits : (digits.slice(0, -maximumFractionDigits) || "0") + "." + ("0".repeat(maximumFractionDigits) + digits).slice(-fractionDigits));
    return {plusSign: plusSign, minusSign: minusSign, number: number, exponentMinusSign: '', exponentInteger: ''};
  };

Interval.Context.prototype.formatToParts = function (x, rounding) {
  if (BigDecimal.equal(x.a, x.b)) {
    //TODO: ?
    var value = BigDecimal.round(x.a, rounding.significantDigits != undefined ? {maximumSignificantDigits: rounding.significantDigits, roundingMode: "half-up"} : {maximumFractionDigits: rounding.fractionDigits, roundingMode: "half-up"});
    return digitsToDecimalNumber({sign: BigDecimalMath.sign(x.a), decimal: BigDecimalMath.abs(value).toString(), exact: BigDecimal.equal(value, x.a)}, rounding); // "half-up" rounding
  }
  var signA = BigDecimalMath.sign(x.a);
  var signB = BigDecimalMath.sign(x.b);
  // assume, that the value is not exact zero
  //if ((signA || signB) === (signB || signA)) {
    var candidateA = BigDecimal.round(x.a, rounding.significantDigits != undefined ? {maximumSignificantDigits: rounding.significantDigits, roundingMode: signA >= 0 ? "half-up" : "half-down"} : {maximumFractionDigits: rounding.fractionDigits, roundingMode: signA >= 0 ? "half-up" : "half-down"}); // "half-ceil" rounding
    var candidateB = BigDecimal.round(x.b, rounding.significantDigits != undefined ? {maximumSignificantDigits: rounding.significantDigits, roundingMode: signB <= 0 ? "half-up" : "half-down"} : {maximumFractionDigits: rounding.fractionDigits, roundingMode: signA >= 0 ? "half-up" : "half-down"}); // "half-floor" rounding
    if (BigDecimal.equal(candidateA, candidateB)) {
      return digitsToDecimalNumber({sign: (signA || signB) === (signB || signA) ? (signA || signB) : 0, decimal: BigDecimalMath.abs(candidateA).toString(), exact: false}, rounding);
    }
  //}
  return undefined;
};
Interval.prototype.toString = function () {
  return "[" + this.a.toString() + ";" + this.b.toString() + "]";
};

if (false) {
  //TODO: how to test it?
  var x = new Interval.Context(2).cos(new Interval(BigDecimal.BigDecimal((-3 / 4 * Math.PI).toString()), BigDecimal.BigDecimal((1 / 4 * Math.PI).toString()))).toString();
  console.assert(x === '[-0.72;-1]');
}

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
    return context.nthRoot(y, n);
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
      } else if (operator === "^") { // Expression.PolynomialRoot^3, pi^2
        if (!(e.b instanceof Expression.Integer) || e.b.toNumber() <= 0 || e.b.toNumber() > 9007199254740991) {
          //throw new TypeError();
          return context.exp(context.multiply(context.log(a), b));
        }
        var n = e.b.toNumber();//TODO: FIX!
        var exponentiate = function (x, n) {
          var y = undefined;
          while (n >= 1) {
            if (n === 2 * Math.floor(n / 2)) {
              x = context.multiply(x, x);
              n = Math.floor(n / 2);
            } else {
              y = y == undefined ? x : context.multiply(x, y);
              n -= 1;
            }
          }
          return y;
        };
        return exponentiate(a, n);
      }
    }
  } else if (e instanceof Expression.PolynomialRoot) {
    var i = e.toDecimal(context.precision);
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
    return context.log(x);
  }

  return undefined;
};

var decimalToString = function (parts) {
  return parts.plusSign + parts.minusSign + parts.number + (parts.exponentInteger !== '' ? '*' + (10) + '^' + parts.exponentMinusSign + parts.exponentInteger : '');
};

var complexToString = function (real, imaginary) {
  return real + (/^[\-\+]/.test(imaginary) ? imaginary.replace(/^([\-\+])[\s\S]+/g, '$1') : (real !== '' ? '+' : '')) + (imaginary !== '1' && imaginary !== '-1' ? imaginary.replace(/^[\-\+]/g, '') + '*' + 'i' : 'i');
};

  //? ((n * 10**(fractionDigits + 1)) ~/ d + 5) ~/ 10

var toDecimalStringInternal = function (expression, rounding, decimalToStringCallback, complexToStringCallback) {
  decimalToStringCallback = decimalToStringCallback || decimalToString;
  complexToStringCallback = complexToStringCallback || complexToString;
  if (rounding.fractionDigits == undefined && rounding.significantDigits == undefined) {
    throw new RangeError();
  }

  if (expression instanceof Expression.Complex || expression instanceof Expression.Division || expression instanceof Expression.Addition) {
    var numerator = expression.getNumerator();//.unwrap();
    var denominator = expression.getDenominator();//.unwrap();
    if (denominator instanceof Expression.Integer || Expression.has(denominator, Expression.PolynomialRoot)) {
      if (numerator instanceof Expression.Addition || numerator instanceof Expression.Multiplication || numerator instanceof Expression.Complex) {
        var realValue = Expression.ZERO;
        var imaginaryValue = Expression.ZERO;
        var ok = true;
        var e = numerator;
        for (var x of e.summands()) {
          var c = undefined;
          var r = Expression.ONE;
          for (var y of x.factors()) {
            if (c == undefined && y instanceof Expression.Complex) {
              c = y;
            } else if (y instanceof Expression.NthRoot || y instanceof Expression.Integer) {//TODO: ?
              r = r.multiply(y);
            } else if ((y instanceof Expression.Logarithm || y instanceof Expression.Sin || y instanceof Expression.Cos) && Expression.isConstant(y.a)) {
              r = r.multiply(y);
            } else if (y instanceof Expression.PolynomialRoot) {//TODO: ?
              r = r.multiply(y);
            } else if (y instanceof Expression.Exponentiation && y.a instanceof Expression.PolynomialRoot && y.b instanceof Expression.Integer) {//TODO: ?
              r = r.multiply(y);
            } else {
              ok = false;
            }
          }
          realValue = realValue.add(r.multiply(c == undefined ? Expression.ONE : c.real));
          imaginaryValue = imaginaryValue.add(c != undefined ? r.multiply(c.imaginary) : Expression.ZERO);
        }
        if (ok && !imaginaryValue.equals(Expression.ZERO)) {
          realValue = realValue.divide(denominator);
          imaginaryValue = imaginaryValue.divide(denominator);
          var real = toDecimalStringInternal(realValue, rounding, decimalToStringCallback, complexToStringCallback);
          var imaginary = toDecimalStringInternal(imaginaryValue, rounding, decimalToStringCallback, complexToStringCallback);
          return complexToStringCallback(realValue.equals(Expression.ZERO) ? '' : real, imaginary);
        }
      }
    }
  }
  if (expression instanceof Expression.Division && expression.a instanceof Expression.Integer && expression.b instanceof Expression.Integer) {
    //TODO: ?
    if (expression.b.compareTo(Expression.ZERO) > 0) {//TODO: ?
      var n = expression.a.toBigInt();
      var d = expression.b.toBigInt();
      //debugger;
      ///^0\.(\d+?)(\d*?)(?:\1\2)*\1$/.exec('0.123123')
      // https://softwareengineering.stackexchange.com/a/192081
      // https://en.wikipedia.org/wiki/Repeating_decimal#Other_properties_of_repetend_lengths
      // "If k = 2**a*5**b*n where n > 1 and n is not divisible by 2 or 5, then the length of the transient of 1/k is max(a, b), and the period equals r, where r is the smallest integer such that 10r ≡ 1 (mod n)."
      var a = nthRoot.countTrailingZeros(d, 2);
      var b = nthRoot.countTrailingZeros(d, 5);
      var tmp = BigInteger.divide(d, BigInteger.multiply(BigInteger.exponentiate(2, a), BigInteger.exponentiate(5, b)));
      if (!BigInteger.equal(tmp, BigInteger.BigInt(1))) { // a repeating decimal, the result is not exact
        var e = (nthRoot.naturalLogarithm(n) - nthRoot.naturalLogarithm(d)) / Math.log(10);
        var maximumFractionDigits = rounding.fractionDigits != undefined ? rounding.fractionDigits : rounding.significantDigits - (e);
        var z = 1;
        var period = 0;
        do {
          period += 1;
          z = BigInteger.remainder(BigInteger.multiply(z, 10), tmp);
        } while (period <= maximumFractionDigits - Math.max(a, b) && z !== 1);
        var result = BigDecimal.divide(BigDecimal.BigDecimal(n), BigDecimal.BigDecimal(d), rounding.significantDigits != undefined ? {maximumSignificantDigits: rounding.significantDigits, roundingMode: 'half-up'} : {maximumFractionDigits: rounding.fractionDigits, roundingMode: 'half-up'});
        var f = digitsToDecimalNumber({sign: BigDecimalMath.sign(result), decimal: BigDecimalMath.abs(result).toString(), exact: false}, rounding);
        var j = f.number.indexOf('.'); //?
        var offset = j + 1 + Math.max(a, b) - (f.exponentInteger !== '' && f.exponentMinusSign === '-' ? Number(f.exponentInteger) - 1 : 0);
        if (j !== -1 && (offset + period < f.number.length || offset + period === f.number.length && f.number.charCodeAt(offset) < '5'.charCodeAt(0))) {
          var s = f.number;
          var number = s.slice(0, offset) + '<span style="text-decoration:overline;">' + s.slice(offset, offset + period) + '</span>' + s.slice(offset + period);
          f = Object.assign({}, f, {number: number});
        }
        return decimalToStringCallback(f);
      }
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
      var a = BigDecimal.divide(BigDecimal.BigDecimal(nearest), BigDecimal.BigDecimal(scale));
      var result = context.formatToParts(new Interval(a, BigDecimalMath.nextCeil(a, {maximumFractionDigits: fractionDigits + 1})), fractionDigits);
      return decimalToStringCallback(result);
    }
  }
  */
  //---
  if (!Expression.has(expression, Expression.Symbol) &&
      !Expression.has(expression, Expression.NthRoot) &&
      !Expression.has(expression, Expression.PolynomialRoot) &&
      !(expression instanceof Expression.Integer) &&
      !(expression instanceof Expression.Division) &&
      !Expression.has(expression, Expression.Radians) &&
      !Expression.has(expression, Expression.Logarithm)) {
    throw new TypeError("toDecimalString:" + JSON.stringify(rounding) + ":" + expression.toString({}));
  }
  if (rounding.fractionDigits != undefined && (rounding.fractionDigits < 0 || rounding.fractionDigits > 9007199254740991) ||
      rounding.significantDigits != undefined && (rounding.significantDigits < 1 || rounding.significantDigits > 9007199254740991) ||
      rounding.roundingMode != undefined) {
    throw new RangeError();
  }
  var result = undefined;
  var guessedPrecision = 1; //TODO: ?
  while (result == undefined) {
    //if (guessedPrecision > 1024) throw new Error();
    var context = new Interval.Context(guessedPrecision);
    var x = evaluateExpression(expression, context);
    if (x == undefined) {
      return undefined;
    }
    if (x !== "CANNOT_DIVIDE") { // continue the loop otherwise
      result = context.formatToParts(x, rounding);
    }
    guessedPrecision *= 2;
  }
  return decimalToStringCallback(result);
};

export default toDecimalStringInternal;
