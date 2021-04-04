import Expression from './Expression.js';
//import BigInteger from './BigInteger.js';
import nthRoot from './nthRoot.js';
//import bitLength from './bitLength.js';
//import BigDecimal from './BigDecimal.js';
import primeFactor from './primeFactor.js';

import {BigFloat} from './BigDecimal/BigDecimal.js';
//import BigDecimal from './BigDecimal/BigDecimalByDecimal.js.js';
//import BigDecimal from './BigDecimal/BigFloat.js';
//BigDecimal.BigDecimal = BigDecimal.BigFloat;

const BigDecimal = BigFloat;
const BASE = 2;

//TODO: ???

//TODO: rename BigDecimal to BigFloat, BigDecimalMath to BigFloatMath

// BigFloat requirements:
// only uses roundingMode 'half-even'
// and rounding to maximumSignificantDigits

function BigDecimalMath() {
}
BigDecimalMath.ZERO = BigDecimal.BigDecimal(0);
BigDecimalMath.BASE = BigDecimal.round(BigDecimal.BigDecimal(BASE), {maximumSignificantDigits: 1, roundingMode: 'half-even'});
BigDecimalMath.sign = function (a) { // returns -1, 0, or +1
  return BigDecimal.lessThan(a, BigDecimalMath.ZERO) ? -1 : (BigDecimal.greaterThan(a, BigDecimalMath.ZERO) ? +1 : 0);
};
BigDecimalMath.abs = function (a) {
  return BigDecimal.lessThan(a, BigDecimalMath.ZERO) ? BigDecimal.unaryMinus(a) : (BigDecimal.greaterThan(a, BigDecimalMath.ZERO) ? a : BigDecimalMath.ZERO);
};
BigDecimalMath.max = function (a, b) {
  return BigDecimal.lessThan(a, b) ? b : a;
};
BigDecimalMath.min = function (a, b) {
  return BigDecimal.greaterThan(a, b) ? b : a;
};
BigDecimalMath.nthRoot = function (x, n, rounding) {
  if ((true || n > 100) && !BigDecimal.lessThan(x, BigDecimal.BigDecimal(2)) && rounding.maximumSignificantDigits != undefined && 2 + rounding.maximumSignificantDigits < n) {//?
    // ExpressionParser.parse('2**(1/10000)').toMathML({fractionDigits: 10});
    var roundingToInteger = {
      maximumFractionDigits: 0,
      roundingMode: 'half-even'
    };
    var internalRounding = {
      maximumSignificantDigits: rounding.maximumSignificantDigits + 1 + Math.max(Math.ceil(Math.log2(Number(BigDecimal.toBigInt(BigDecimal.round(BigDecimal.log(BigDecimal.round(x, {maximumSignificantDigits: 1, roundingMode: 'half-even'}), {maximumSignificantDigits: 1, roundingMode: 'half-even'}), roundingToInteger)).toString()) / n) / Math.log2(BASE)), 0),
      roundingMode: "half-even"
    };
    return BigDecimal.exp(BigDecimal.divide(BigDecimal.log(x, internalRounding), BigDecimal.BigDecimal(n), internalRounding), internalRounding);
  }
  // (f * 10**e)**(1/n) = (f * 10**(e - round(e / n) * n))**(1/n) * 10**round(e / n)
  var fractionDigits = 2;
  while (!BigDecimal.equal(BigDecimal.round(x, {maximumFractionDigits: fractionDigits, roundingMode: 'half-even'}), x)) {
    fractionDigits *= 2;
  }
  var scaling = Math.max(rounding.maximumSignificantDigits || rounding.maximumFractionDigits, fractionDigits / 2);
  var scalingCoefficient = BigDecimalMath.exponentiate(BigDecimalMath.BASE, scaling);
  var sx = BigDecimal.toBigInt(BigDecimal.multiply(x, BigDecimalMath.exponentiate(BigDecimalMath.BASE, scaling * n)));
  var x0 = nthRoot(sx, n);
  //var x1 = BigInteger.lessThan(BigInteger.exponentiate(x1, BigInteger.BigInt(n)), sA) ? BigInteger.add(x1, BigInteger.BigInt(1)) : x1;
  return BigDecimal.divide(BigDecimal.BigDecimal(x0), BigDecimal.BigDecimal(scalingCoefficient));
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
    var scalingCoefficient = BigDecimal.BigDecimal(BigDecimalMath.exponentiate(BigDecimalMath.BASE, rounding.maximumFractionDigits));
    var small = BigDecimal.divide(BigDecimal.BigDecimal(which === 'floor' ? -1 : 1), scalingCoefficient);
    return BigDecimal.add(a, small, {maximumFractionDigits: rounding.maximumFractionDigits, roundingMode: which === 'floor' ? 'ceil' : 'floor'});
  }
  if (BigDecimal.equal(a, BigDecimal.BigDecimal(0))) {
    return a;
  }
  var t = BigDecimal.round(a, {maximumSignificantDigits: rounding.maximumSignificantDigits, roundingMode: 'half-even'});
  if (!BigDecimal.equal(a, t)) {
    if (which === 'floor' && BigDecimal.lessThan(t, a) || which === 'ceil' && BigDecimal.greaterThan(t, a)) {
      return t;
    }
    return BigDecimalMath._nextToward(t, rounding, which);
  }
  var k = 1;
  while (true) {
    var m = BigDecimalMath.exponentiate(BigDecimalMath.BASE, rounding.maximumSignificantDigits);
    var c = BigDecimal.add(BigDecimal.BigDecimal(1), BigDecimal.divide(BigDecimal.BigDecimal(which === 'floor' ? -k : k), m));
    var guess = BigDecimal.multiply(a, c, {maximumSignificantDigits: rounding.maximumSignificantDigits, roundingMode: 'half-even'});
    if (!BigDecimal.equal(a, guess)) {
      return guess;
    }
    ++k;
  }
};
BigDecimalMath.nextFloor = function (a, rounding) {
  return BigDecimalMath._nextToward(a, rounding, 'floor');
};
BigDecimalMath.nextCeil = function (a, rounding) {
  return BigDecimalMath._nextToward(a, rounding, 'ceil');
};
//TODO: remove
BigDecimalMath.fma = function (a, b, c, rounding) { // a * b + c
  return BigDecimal.round(BigDecimal.add(BigDecimal.multiply(a, b), c), rounding);
};
BigDecimalMath.exponentiate = function (a, n) {
  //console.assert(a == 2);
  if (n < 0) {
    return BigDecimal.divide(BigDecimal.BigDecimal(1), BigDecimalMath.exponentiate(a, -n));
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

//BigDecimalMath.nthRoot(BigDecimal.BigDecimal(2), 2, 3);

if (false) {
  console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0.995'), {maximumSignificantDigits: 4}).toString() === '0.9951');
  console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0.995'), {maximumSignificantDigits: 3}).toString() === '0.996');
  console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0.995'), {maximumSignificantDigits: 2}).toString() === '1');
  console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0'), {maximumSignificantDigits: 1}).toString() === '0');
  console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('1'), {maximumSignificantDigits: 1}).toString() === '2');
  //console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0.995'), {maximumFractionDigits: 4}).toString() === '0.9951');
  //console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0.995'), {maximumFractionDigits: 3}).toString() === '0.996');
  //console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0.995'), {maximumFractionDigits: 2}).toString() === '1');
  //console.assert(BigDecimalMath.nextCeil(BigDecimal.BigDecimal('0'), {maximumFractionDigits: 1}).toString() === '0.1');
  //console.assert(BigDecimalMath.nextFloor(BigDecimal.BigDecimal('0'), {maximumFractionDigits: 1}).toString() === '-0.1');
  console.time();
  console.assert(BigDecimalMath.nextFloor(BigDecimal.BigDecimal('1e10000000'), {maximumSignificantDigits: 100000}).toString() === '?');
  console.timeEnd();
  // default: 2289.925048828125ms
  // default: 1335.22216796875ms
}

const USE_DIRECTED_ROUNDING_MODES = true;

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

Interval.Context = function (precision, flag0) {
  //this.precision = precision;
  //this.anyRounding = {maximumSignificantDigits: precision, roundingMode: 'half-even'};
  var anyRounding = flag0 ? {maximumSignificantDigits: precision, roundingMode: 'half-even'} : {maximumFractionDigits: precision - 1, roundingMode: 'half-even'};
  this.anyRounding = anyRounding;
  this.floorRounding = Object.assign({}, anyRounding, {roundingMode: 'floor'});
  this.ceilRounding = Object.assign({}, anyRounding, {roundingMode: 'ceil'});
};
Interval.Context.prototype.unaryMinus = function (x) {
  if (BigDecimal.equal(x.a, x.b)) {
    return Interval.degenerate(BigDecimal.unaryMinus(x.a));
  }
  return new Interval(BigDecimal.unaryMinus(x.b), BigDecimal.unaryMinus(x.a));
};
Interval.Context.prototype._add = function (x, y, which) {
  if (USE_DIRECTED_ROUNDING_MODES) {
    if (which === 'floor') {
      var floorRounding = this.floorRounding;
      return BigDecimal.add(x, y, floorRounding);
    }
    if (which === 'floor') {
      var ceilRounding = this.ceilRounding;
      return BigDecimal.add(x, y, ceilRounding);
    }
  }
  if (BigDecimal.lessThan(BigDecimalMath.abs(x), BigDecimalMath.abs(y))) {
    return this._add(y, x, which);
  }
  // |x| <= |y|
  //TODO: check that round(x) === x and round(y) === y - ?
  var s = BigDecimal.add(x, y, this.anyRounding);
  var e = BigDecimal.subtract(BigDecimal.subtract(s, x), y);
  if (which === 'floor' && BigDecimal.greaterThan(e, BigDecimal.BigDecimal(0))) {
    return BigDecimalMath.nextFloor(s, this.anyRounding);
  }
  if (which === 'ceil' && BigDecimal.lessThan(e, BigDecimal.BigDecimal(0))) {
    return BigDecimalMath.nextCeil(s, this.anyRounding);
  }
  return s;
};
Interval.Context.prototype.add = function (x, y) {
  return new Interval(this._add(x.a, y.a, 'floor'), this._add(x.b, y.b, 'ceil'));
};
Interval.Context.prototype.subtract = function (x, y) {
  return new Interval(this._add(x.a, BigDecimal.unaryMinus(y.b), 'floor'), this._add(x.b, BigDecimal.unaryMinus(y.a), 'ceil'));
};
Interval.Context.prototype.multiply = function (x, y) {
  var anyRounding = this.anyRounding;
  var floorRounding = this.floorRounding;
  var ceilRounding = this.ceilRounding;
  var floorMultiply = function (x, y) {
    if (USE_DIRECTED_ROUNDING_MODES) {
      return BigDecimal.multiply(x, y, floorRounding);
    }
    //TODO: test for wrong roundingMode
    var p = BigDecimal.multiply(x, y, anyRounding);
    var e = BigDecimalMath.fma(x, y, BigDecimal.unaryMinus(p), anyRounding);
    if (BigDecimal.lessThan(e, BigDecimal.BigDecimal(0))) {
      return BigDecimalMath.nextFloor(p, anyRounding);
    }
    return p;
  };
  var ceilMultiply = function (x, y) {
    if (USE_DIRECTED_ROUNDING_MODES) {
      return BigDecimal.multiply(x, y, ceilRounding);
    }
    //TODO: test for wrong roundingMode
    var p = BigDecimal.multiply(x, y, anyRounding);
    var e = BigDecimalMath.fma(x, y, BigDecimal.unaryMinus(p), anyRounding);
    if (BigDecimal.greaterThan(e, BigDecimal.BigDecimal(0))) {
      return BigDecimalMath.nextCeil(p, anyRounding);
    }
    return p;
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
  var anyRounding = this.anyRounding;
  var floorRounding = this.floorRounding;
  var ceilRounding = this.ceilRounding;
  var floorDivide = function (x, y) {
    if (USE_DIRECTED_ROUNDING_MODES) {
      return BigDecimal.divide(x, y, floorRounding);
    }
    var q = BigDecimal.divide(x, y, anyRounding);
    var e = BigDecimalMath.fma(q, y, BigDecimal.unaryMinus(x), anyRounding);
    if (BigDecimal.greaterThan(BigDecimalMath.sign(y) < 0 ? BigDecimal.unaryMinus(e) : e, BigDecimal.BigDecimal(0))) {
      return BigDecimalMath.nextFloor(q, anyRounding);
    }
    return q;
  };
  var ceilDivide = function (x, y) {
    if (USE_DIRECTED_ROUNDING_MODES) {
      return BigDecimal.divide(x, y, ceilRounding);
    }
    var q = BigDecimal.divide(x, y, anyRounding);
    var e = BigDecimalMath.fma(q, y, BigDecimal.unaryMinus(x), anyRounding);
    if (BigDecimal.lessThan(BigDecimalMath.sign(y) < 0 ? BigDecimal.unaryMinus(e) : e, BigDecimal.BigDecimal(0))) {
      return BigDecimalMath.nextCeil(q, anyRounding);
    }
    return q;
  };
  if (BigDecimalMath.sign(y.a) <= 0 && BigDecimalMath.sign(y.b) >= 0) {
    if (BigDecimal.equal(y.a, y.b)) {
      throw new RangeError();
    }
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
  if (BigDecimalMath.sign(x.a) <= 0 && BigDecimalMath.sign(x.b) > 0) {
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
  var middle = BigDecimal.divide(BigDecimal.add(x.a, x.b), BigDecimal.BigDecimal(2), this.anyRounding); // with rounding it works better in case the interval has huge significant digits difference
  var extremumPoint = function (q) {
    var shift = BigDecimal.multiply(BigDecimal.divide(BigDecimal.BigDecimal(q), BigDecimal.BigDecimal(2), anyRounding), pi);
    var k = BigDecimal.round(BigDecimal.divide(BigDecimal.subtract(middle, shift, anyRounding), tau, anyRounding), {
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
  if (BASE !== 2) {
    var abs = function (a) {
      return a < 0 ? BigInteger.unaryMinus(a) : a;
    };
    var k = a != 0 ? bitLength(abs(a)) - Math.ceil(Math.log2(10) * this.anyRounding.maximumSignificantDigits) : 0;
    if (k > 42) {
      //TODO: move to BigDecimal.round - ?
      // for performance
      var p2k = BigInteger.exponentiate(BigInteger.BigInt(2), BigInteger.BigInt(k));
      var q = BigInteger.divide(a, p2k);
      var from = a < 0 ? BigInteger.subtract(q, BigInteger.BigInt(1)) : q;
      var to = a < 0 ? q : BigInteger.add(q, BigInteger.BigInt(1));
      return this.multiply(this.fromIntegers(from, to), this.exponentiate(this.fromInteger(BigInteger.BigInt(2)), k));
    }
  }
  //var x = BigDecimal.BigDecimal(a);
  //if (BigInteger.equal(BigInteger.BigInt(BigDecimal.toBigInt(x)), a)) { // TODO: ?
  //  return Interval.degenerate(x);
  //}
  a = BigDecimal.BigDecimal(a);
  if (USE_DIRECTED_ROUNDING_MODES) {
    var anyRounding = this.anyRounding;
    var floorRounding = this.floorRounding;
    var ceilRounding = this.ceilRounding;
    return new Interval(BigDecimal.round(a, floorRounding), BigDecimal.round(a, ceilRounding));
  }
  var guess = BigDecimal.round(a, this.anyRounding);
  if (BigDecimal.greaterThan(guess, a)) {
    return new Interval(BigDecimalMath.nextFloor(guess, this.anyRounding), guess);
  } else if (BigDecimal.lessThan(guess, a)) {
    return new Interval(guess, BigDecimalMath.nextCeil(guess, this.anyRounding));
  }
  return Interval.degenerate(guess);
  //return this.fromIntegers(a, a);
};
Interval.Context.prototype.fromIntegers = function (a, b) {
  var anyRounding = this.anyRounding;
  var floorRounding = this.floorRounding;
  var ceilRounding = this.ceilRounding;
  var floorRound = function (x) {
    if (USE_DIRECTED_ROUNDING_MODES) {
      return BigDecimal.round(x, floorRounding);
    }
    var guess = BigDecimal.round(x, anyRounding);
    if (BigDecimal.greaterThan(guess, x)) {
      return BigDecimalMath.nextFloor(guess, anyRounding);
    }
    return guess;
  };
  var ceilRound = function (x) {
    if (USE_DIRECTED_ROUNDING_MODES) {
      return BigDecimal.round(x, ceilRounding);
    }
    var guess = BigDecimal.round(x, anyRounding);
    if (BigDecimal.lessThan(guess, x)) {
      return BigDecimalMath.nextCeil(guess, anyRounding);
    }
    return guess;
  };
  var a = BigDecimal.BigDecimal(a);
  var b = BigDecimal.BigDecimal(b);
  //TODO: test case (!!!)
  console.assert(!BigDecimal.greaterThan(a, b));
  return new Interval(floorRound(a),
                      ceilRound(b));
};
Interval.Context.prototype.abs = function (x) {
  if (BigDecimal.lessThan(x.a, BigDecimal.BigDecimal(0))) {
    if (BigDecimal.lessThan(x.b, BigDecimal.BigDecimal(0))) {
      return new Interval(BigDecimal.unaryMinus(x.b), BigDecimal.unaryMinus(x.a));
    } else {
      return new Interval(BigDecimal.BigDecimal(0), BigDecimalMath.max(BigDecimal.unaryMinus(x.a), x.b));
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
  var signA = BigDecimalMath.sign(x.a);
  var signB = BigDecimalMath.sign(x.b);
  var sign = (signA || signB) === (signB || signA) ? (signA || signB) : 0;
  x = this.abs(x);
  var anyRounding = this.anyRounding;
  var stringify = function (a, roundingMode) {
    //TODO: roundingMode - ?
    /*if (roundingMode === "half-down" && BASE === 2) {
      var string = stringify(a, "half-up");
      var match = /^(\d+)\.?(\d*)e?\+?(\-?\d*)$/.exec(string);
      var significand = match[1] + (match[2] || '') + '5';
      var exponent = BigInteger.subtract(BigInteger.BigInt(match[3]), BigInteger.BigInt((match[2] || '').length + 1));
      var a1 = BigDecimal.divide(a, BigDecimalMath.exponentiate(BigDecimal.BigDecimal(2), exponent));
      if (BigDecimal.equal(BigDecimal.round(a1, {maximumFractionDigits: 0, roundingMode: "half-even"}), a1)) {
        var a2 = BigDecimal.divide(a1, BigDecimal.BigDecimal(2));
        if (!BigDecimal.equal(BigDecimal.round(a2, {maximumFractionDigits: 0, roundingMode: "half-even"}), a2)) {
          // a1 should be "small"
          var u = BigInteger.subtract(BigInteger.BigInt(significand), BigInteger.BigInt(10));
          var s = BigInteger.BigInt(BigDecimal.toBigInt(a1));
          if (exponent < 0) {
            s = BigInteger.multiply(s, BigInteger.exponentiate(BigInteger.BigInt(5), BigInteger.unaryMinus(exponent)));
          } else {
            u = BigInteger.multiply(u, BigInteger.exponentiate(BigInteger.BigInt(5), exponent));
          }
          if (BigInteger.equal(s, u)) {
            return stringify(BigDecimalMath.nextFloor(a, {maximumSignificantDigits: Math.max(anyRounding.maximumSignificantDigits, Math.ceil(Math.log2(10) * significand.length)), roundingMode: "half-even"}), "half-up");
          }
        }
      }
      return string;
    }*/
    if (rounding.fractionDigits != undefined) {
      return a.toFixed(rounding.fractionDigits, roundingMode);
    }
    return a.toPrecision(rounding.significantDigits, roundingMode);
  };
  var a = stringify(x.a, "half-up");
  var b = BigDecimal.equal(x.a, x.b) ? a : stringify(x.b, "half-down");
  if (a !== b && (this.anyRounding.maximumSignificantDigits > 60000 || this.anyRounding.maximumFractionDigits > 60000)) {
    debugger;
    throw new Error();
  }
  var isZero = function (a) {
    return !/[^0\.]/.test(a);
  };
  if (a === b && (sign !== 0 || isZero(a) && isZero(b))) {
    return (sign < 0 ? '-' : (sign > 0 && isZero(a) && isZero(b) ? '+' : '')) + a;
  }
  //var a = BigDecimal.toBigInt(BigDecimalMath.round(x.a));
  //var b = BigDecimal.equal(x.a, x.b) ? a : BigDecimal.toBigInt(BigDecimal.unaryMinus(BigDecimalMath.round(BigDecimal.unaryMinus(x.b))));
  //if (sign !== 0 && BigInteger.equal(a, b) || sign === 0 && a == 0 && b == 0) {
    //TODO: ?
    //var value = BigDecimal.round(x.a, rounding.significantDigits != undefined ? {maximumSignificantDigits: rounding.significantDigits, roundingMode: "half-up"} : {maximumFractionDigits: rounding.fractionDigits, roundingMode: "half-up"});
    //var decimal = BigDecimalMath.abs(value).toString();
    //return (sign < 0 ? '-' : (sign > 0 && a == 0 && b == 0 ? '+' : '')) + digitsToDecimalNumber(a.toString(), -fd, rounding); // "half-up" rounding
  //}
    //var candidateA = BigDecimal.round(x.a, rounding.significantDigits != undefined ? {maximumSignificantDigits: rounding.significantDigits, roundingMode: signA >= 0 ? "half-up" : "half-down"} : {maximumFractionDigits: rounding.fractionDigits, roundingMode: signA >= 0 ? "half-up" : "half-down"}); // "half-ceil" rounding
    //var candidateB = BigDecimal.round(x.b, rounding.significantDigits != undefined ? {maximumSignificantDigits: rounding.significantDigits, roundingMode: signB <= 0 ? "half-up" : "half-down"} : {maximumFractionDigits: rounding.fractionDigits, roundingMode: signB >= 0 ? "half-up" : "half-down"}); // "half-floor" rounding
    //if (BigDecimal.equal(candidateA, candidateB)) {
      //return (sign < 0 ? '-' : (sign > 0 && a == 0 && b == 0 ? '+' : '')) + digitsToDecimalNumber({decimal: BigDecimal.toBigInt(BigDecimalMath.abs(BigDecimal.multiply(candidateA, BigDecimal.BigDecimal('1e+' + e)))).toString() + 'e-' + e}, rounding);
    //}
  //return undefined;
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
    // slow for some cases:
    if (e instanceof Expression.Addition && Expression.has(e, Expression.PolynomialRoot)) {
      var root = Expression.getVariable(e);//?
      var p = Polynomial.toPolynomial(e, root);
      if (p.hasIntegerCoefficients()) {// trying to avoid slow cases (?)
        //TODO: https://en.wikipedia.org/wiki/Horner%27s_method - ?
        var zero = evaluateExpression(root, context);
        //return evaluateExpression(p.calcAt(), context);
        var result = context.fromInteger(Expression.ZERO.value);
        for (var i = p.getDegree(); i >= 0; i--) {
          result = context.multiply(result, zero);
          var tmp = evaluateExpression(p.getCoefficient(i), context);
          if (tmp === "CANNOT_DIVIDE" || tmp == undefined) {
            return tmp;
          }
          result = context.add(result, tmp);
        }
        return result;
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
    } else if (operator === "^") { // Expression.PolynomialRoot^3, pi^2
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
  } else if (e instanceof Expression.PolynomialRoot) {
    var i = e.toDecimal(context.anyRounding.maximumSignificantDigits || context.anyRounding.maximumFractionDigits);
    // "lcm" is too slow to compute (?)
    /*if (true) {
      var a = BigDecimal.divide(BigDecimal.BigDecimal(i.a.getNumerator().value), BigDecimal.BigDecimal(i.a.getDenominator().value), context.floorRounding);
      var b = BigDecimal.divide(BigDecimal.BigDecimal(i.b.getNumerator().value), BigDecimal.BigDecimal(i.b.getDenominator().value), context.ceilRounding);
      return new Interval(a, b);
    }*/
    return context.divide(context.fromIntegers(i.b.getDenominator().multiply(i.a.getNumerator()).value,
                                               i.a.getDenominator().multiply(i.b.getNumerator()).value),
                          context.fromInteger(i.a.getDenominator().multiply(i.b.getDenominator()).value));
  } else if (e === Expression.E) {
    return context.exp(context.fromInteger(1));
  } else if (e === Expression.PI) {
    return context.pi();
  } else if (e instanceof Expression.Sin || e instanceof Expression.Cos) {
    var x = evaluateExpression(e.a instanceof Expression.Radians ? e.a.value : e.a, context);
    if (x === "CANNOT_DIVIDE" || x == undefined) {
      return x;
    }
    if (e instanceof Expression.Sin) {
      return context.sin(x);
    }
    if (e instanceof Expression.Cos) {
      return context.cos(x);
    }
  } else if (e instanceof Expression.Logarithm) {
    var x = evaluateExpression(e.a, context);
    if (x === "CANNOT_DIVIDE" || x == undefined) {
      return x;
    }
    return context.log(x);
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
        Expression.has(denominator, Expression.PolynomialRoot) ||
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
      var a = BigDecimal.divide(BigDecimal.BigDecimal(nearest), BigDecimal.BigDecimal(scale));
      var result = context.formatToDecimal(new Interval(a, BigDecimalMath.nextCeil(a, {maximumFractionDigits: fractionDigits + 1})), fractionDigits);
      return decimalToStringCallback(result);
    }
  }
  */
  //---
  /*if (!Expression.has(expression, Expression.Symbol) &&
      !Expression.has(expression, Expression.NthRoot) &&
      !Expression.has(expression, Expression.PolynomialRoot) &&
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
      var tmp = Number(BigDecimalMath.max(BigDecimalMath.abs(x.a), BigDecimalMath.abs(x.b)).toFixed(0));
      guessedPrecision = Math.ceil(Math.log2(Math.min(Math.max(tmp, 2), Number.MAX_VALUE)) / 2 + Math.log2(10) * rounding.fractionDigits / 2 * 2);
    }
    guessedPrecision *= 2;
  }
  if (guessedPrecision !== 256) {
    //console.log(guessedPrecision);
  }
  return decimalToStringCallback(result);
};

export default toDecimalStringInternal;
