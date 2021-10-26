import primeFactor from './primeFactor.js';
import nthRoot from './nthRoot.js';

// API:

// a class for real algebraic numbers
// Operations are implemented as described at https://en.wikipedia.org/wiki/Resultant#Number_theory

// "AbstractPolynomialRoot":
//  #toDecimal(precision)
//  #scale(k), k is an "algebraic expression constant"
//  #translate(k), k is an "algebraic expression constant"
//  #multiply(otherPolynomialRoot)
//  #add(otherPolynomialRoot)
//  #negate()
//  #inverse()
//  #sign()
//  #_pow(n), n is an integer
//  #_nthRoot(n), n is an integer
//  #equals(otherPolynomialRoot)

// PolynomialRoot implements AbstractPolynomialRoot - a basic class to represent real algebraic number exactly
//   .polynomial
//   .interval
// LazyPolynomialRoot implements AbstractPolynomialRoot - a class to represent real algebraic number as a rational expression (see https://en.wikipedia.org/wiki/Algebraic_expression )
//   .e
//   ._root

//TODO: remove references to Expression.ExpressionPolynomialRoot (?)

  var isRational = function (e) {
    return e.getNumerator() instanceof Expression.Integer && e.getDenominator() instanceof Expression.Integer;
  };
  var isPowerOf2 = function (i) {
    return Expression.TWO._pow(primeFactor.bitLength(i.toBigInt()) - 1).equals(i);
  };

function SimpleInterval(a, b) {
  console.assert(isRational(a));
  console.assert(isRational(b));
  console.assert(a.getDenominator().multiply(b.getNumerator()).compareTo(b.getDenominator().multiply(a.getNumerator())) >= 0);
  //TODO: a.getDenominator() is a power of two and b.getDenominator() is a power of two
  //console.assert(isPowerOf2(a.getDenominator()));
  //console.assert(isPowerOf2(b.getDenominator()));
  this.a = a;
  this.b = b;
}
SimpleInterval.prototype.negate = function () {
  return new SimpleInterval(this.b.negate(), this.a.negate());
};
SimpleInterval.prototype.add = function (other) {
  return new SimpleInterval(this.a.add(other.a), this.b.add(other.b));
};
SimpleInterval.prototype.multiply = function (other) {
  var sign = function (e) {
    return e.getNumerator().compareTo(Expression.ZERO);
  };
  if (sign(this.a) === sign(this.b) && sign(other.a) === sign(other.b)) {
    if (sign(this.a) < 0) {
      return this.negate().multiply(other).negate();
    }
    if (sign(other.a) < 0) {
      return this.multiply(other.negate()).negate();
    }
    return new SimpleInterval(this.a.multiply(other.a), this.b.multiply(other.b));
  }
  //TODO: add a test
  var a = this.a.multiply(other.a);
  var b = this.b.multiply(other.a);
  var c = this.a.multiply(other.b);
  var d = this.b.multiply(other.b);
  var tmp = [a, b, c, d];
  tmp.sort((x1, x2) => x1.subtract(x2).getNumerator().compareTo(Expression.ZERO));
  return new SimpleInterval(tmp[0], tmp[3]);
};
SimpleInterval.prototype.scale = function (s) {
  return this.multiply(new SimpleInterval(s, s));//TODO: ?
};
SimpleInterval.prototype.inverse = function () {
  var sign = function (e) {
    return e.getNumerator().compareTo(Expression.ZERO);
  };
  if (sign(this.a) < 0 && sign(this.b) > 0) {
    throw new TypeError();
  }
  return new SimpleInterval(this.b.inverse(), this.a.inverse());
};
SimpleInterval.prototype._pow = function (n) {
  if (n % 2 === 0) {
    if (n === 0) {
      return new SimpleInterval(Expression.ONE, Expression.ONE);
    }
    return this.multiply(this)._pow(n / 2);
  }
  return new SimpleInterval(this.a._pow(n), this.b._pow(n));
};
SimpleInterval.prototype.toString = function () {
  return '[' + this.a + ';' + this.b + ']';
};
SimpleInterval.intersection = function (a, b) {
  var cmp = function (x1, x2) {
    return x1.subtract(x2).getNumerator().compareTo(Expression.ZERO);
  };
  var max = function (x1, x2) {
    return cmp(x1, x2) < 0 ? x2 : x1;
  };
  var min = function (x1, x2) {
    return cmp(x1, x2) < 0 ? x1 : x2;
  };
  // https://scicomp.stackexchange.com/a/26260
  if (cmp(b.a, a.b) > 0 || cmp(a.a, b.b) > 0) {
    return null;
  }
  return {
    a: max(a.a, b.a),
    b: min(a.b, b.b)
  };
};

var toSimpleIntervalOld = function (e, precision) {
  if (isRational(e)) {//TODO: REMOVE
    return new SimpleInterval(e.getNumerator(), e.getNumerator()).scale(e.getDenominator().inverse());
  }
  var tmp = ExpressionParser.parse(toDecimalStringInternal(e, {significantDigits: precision}));
  var epsilonInterval = new SimpleInterval(ExpressionParser.parse('1-5*10**-' + precision), ExpressionParser.parse('1+5*10**-' + precision));
  var result = epsilonInterval.scale(tmp);
  return result;
};

// TODO:
var toSimpleIntervalNew = function (e, precision) { // precision - precision of the computation (?)
  if (e instanceof Expression.Integer) {
    return new SimpleInterval(e, e);
  } else if (e instanceof Expression.BinaryOperation) {
    var a = toSimpleIntervalNew(e.a, precision);
    var b = toSimpleIntervalNew(e.b, precision);
    var s = e.getS();
    if (s === "+") {
      return a.add(b);
    } else if (s === "-") {
      return a.add(b.negate());
    } else if (s === "*") {
      return a.multiply(b);
    } else if (s === "/") {//- why was it commented out - ?
      return a.multiply(b.inverse());
    } else if (s === "^") {
      if (e.b instanceof Expression.Integer) {
        var n = e.b.toBigInt();
        return a._pow(n);
      }
    } else {
      debugger;
    }
  } else if (e instanceof Expression.NthRoot) {
    const n = e.n;
    console.assert(n >= 2 && n % 1 === 0);
    if (e.a instanceof Expression.Integer && e.a.sign() > 0) {
      var a = e.a;
      var scale = Expression.TWO._pow(precision);
      var numerator = Expression.Integer.fromBigInt(nthRoot(a.multiply(scale._pow(n)).toBigInt(), n));
      return new SimpleInterval(numerator.divide(scale), numerator.add(Expression.ONE).divide(scale));
    }
    var a = toSimpleIntervalNew(e.a, precision);
    if (n % 2 === 0) {
      //TODO: !?!?!
      var i = 2;
      while (a.a.getNumerator().sign() < 0 && a.b.getNumerator().sign() > 0) {
        a = toSimpleIntervalNew(e.a, precision * i);
        i *= 2;
      }
    }
    var f = function (x, n, rounding) {
      var scale = Expression.TWO._pow(precision);
      var e = rounding === 'floor' ? (x.getNumerator().sign() >= 0 ? Expression.ZERO : Expression.ONE.negate()) : (x.getNumerator().sign() < 0 ? Expression.ZERO : Expression.ONE);
      return Expression.Integer.fromBigInt(nthRoot(x.getNumerator().multiply(x.getDenominator()._pow(n - 1)).multiply(scale._pow(n)).toBigInt(), n)).add(e).divide(scale.multiply(x.getDenominator()));
    };
    return new SimpleInterval(f(a.a, n, 'floor'), f(a.b, n, 'ceil'));
  } else {
  }
  //TODO: REMOVE(?)
  if (e instanceof Expression.PolynomialRootSymbol) {
    return e.toDecimal(precision);
  }
  if (e instanceof Expression.ExpressionPolynomialRoot) {
    return e.root.toDecimal(precision);
  }
  debugger;
  throw new TypeError("?");
};

const toSimpleInterval = function (e, precision) {
  //return toSimpleIntervalOld(e, precision);
  return toSimpleIntervalNew(e, precision);
  //var result = toSimpleIntervalOld(e, precision);
  //var result2 = toSimpleIntervalNew(e, precision);
  //if (SimpleInterval.intersection(result, result2) == null) {
  //  debugger;
  //}
  //return result2;
};

function Helper(polynomial) {
  this.squareFreeFactors = [];
  var tmp = null;
  const ONE = Polynomial.of(Expression.ONE);
  do {
    tmp = polynomial.squareFreeFactors();
    if (!tmp.a1.equals(ONE)) {//TODO: ?
      this.squareFreeFactors.push(tmp.a1);
    }
    if (tmp.a0.getDegree() !== 0) {
      polynomial = tmp.a0;
    } else {
      polynomial = null;
    }
  } while (polynomial != null);
}
Helper.prototype.calcAt = function (point) {
  var result = Expression.ONE;
  for (var factor of this.squareFreeFactors) {
    result = result.multiply(factor.calcAt(point));
  }
  return result;
};
Helper.prototype.numberOfRoots = function (interval) {
  var result = 0;
  var newFactors = [];
  for (var factor of this.squareFreeFactors) {
    var n = factor.numberOfRoots(interval);
    result += n;
    if (n > 0) {
      newFactors.push(factor);
    }
  }
  this.squareFreeFactors = newFactors;
  return result;
};
Helper.get = function (that, interval) {
  //TODO: do not call twice
  //for (var factor of that.squareFreeFactors) {
  //  if (factor.numberOfRoots(interval) > 0) {
  //    return factor;
  //  }
  //}
  //return null;
  return that.squareFreeFactors.length === 1 ? that.squareFreeFactors[0] : null;
};

  var calculateNewInterval = function (newPolynomial, zeroFunction) {
    if (!newPolynomial.hasIntegerCoefficients()) {
      throw new RangeError("just a check");
    }
    newPolynomial = new Helper(newPolynomial);//!?!?TODO: REMOVE
    var precision = 1;
    var guess = zeroFunction(precision);
    while ((guess.a.getNumerator().sign() !== guess.b.getNumerator().sign() && !newPolynomial.calcAt(Expression.ZERO).equals(Expression.ZERO)) || newPolynomial.numberOfRoots(guess) > 1) {
      precision *= 2;
      guess = zeroFunction(precision);
      if (precision > 1024) throw new Error();//TODO: ?
    }
    var newInterval = guess;
    newPolynomial = Helper.get(newPolynomial, guess);
    return new PolynomialRoot(newPolynomial, newInterval);
  };


function PolynomialRoot(polynomial, interval, options = {}) {
  if (!polynomial.hasIntegerCoefficients()) {
    throw new TypeError();
  }
  if (polynomial.getLeadingCoefficient().compareTo(Expression.ZERO) < 0) {
    return new PolynomialRoot(polynomial.negate(), interval);
  }
  var content = polynomial.getContent();
  if (!content.equals(Expression.ONE)) {
    return new PolynomialRoot(polynomial.scale(content.inverse()), interval);
  }
if (!options.skipFactorization) {//!
  var factor = polynomial.factorize();
  //TODO: pass the zero to help the factorization to return the correct factor:
  // var factor = polynomial.factorize({zero: new PolynomialRoot(polynomial, interval)});
  if (factor != null && !factor.equals(polynomial)) {
    if (factor.numberOfRoots(interval) !== 0) {
      return new PolynomialRoot(factor, interval);
    } else {
      var otherFactor = polynomial.divideAndRemainder(factor, "throw").quotient;
      return new PolynomialRoot(otherFactor, interval);
    }
  }
}
  if (!(interval instanceof SimpleInterval)) {
    throw new TypeError();
  }
  if (interval.a.subtract(interval.b).getNumerator().compareTo(Expression.ZERO) > 0) {
    throw new TypeError();
  }
  // how to represent zero - ?
  //if ((interval.a.getNumerator().sign() || interval.b.getNumerator().sign()) !== (interval.b.getNumerator().sign() || interval.a.getNumerator().sign())) {
  //  throw new TypeError();
  //}
  if (polynomial.numberOfRoots(interval) !== 1) {
    throw new TypeError();
  }
  if (!polynomial.getContent().equals(Expression.ONE)) {
    throw new TypeError();
  }
  //TODO: factorization
  this.polynomial = polynomial;
  //TODO: https://www.wolframalpha.com/input/?i=x**5%2B7x**3%2Bx**2%2Bx%2B1%3D0
  this.interval = interval;
}

PolynomialRoot.prototype.toDecimal = function (precision) {
  var tmp = this.polynomial.getZero(this.interval, precision);
  return new SimpleInterval(tmp.a, tmp.b);
};

PolynomialRoot.prototype.toString = function () {
  // for debugging (?)
  return "[root of " + this.polynomial + " near " + this.interval.a.add(this.interval.b).divide(Expression.TWO).toString() + "]";
};

//TODO: remove (?)
PolynomialRoot.prototype.scale = function (k) {
  //console.assert(k instanceof Expression.Integer || isRational(k));
  // z = k * x, x = z / k
  var newPolynomial = this.polynomial._scaleRoots(k).primitivePart();
  if (!(isRational(k))) { // TODO: remove
    var root = this;
    newPolynomial = toPolynomialWithIntegerCoefficients(newPolynomial);
    return calculateNewInterval(newPolynomial, function (precision) {
      return root.toDecimal(precision).multiply(toSimpleInterval(k, precision));
    });
  }
  var newInterval = this.interval.scale(k);
  return new PolynomialRoot(newPolynomial, newInterval);
};

//TODO: remove (?)
PolynomialRoot.prototype.translate = function (k) {
  //console.assert(k instanceof Expression.Integer || isRational(k));//TODO: ???
  // z = x + k, x = z - k
  var newPolynomial = this.polynomial._translateRoots(k).primitivePart();
  // to avoid intervals, which include zero
  var root = this;
  var newInterval = null;
  if (!(isRational(k))) { // TODO: remove
    newPolynomial = toPolynomialWithIntegerCoefficients(newPolynomial);
    return calculateNewInterval(newPolynomial, function (precision) {
      return root.toDecimal(precision).add(toSimpleInterval(k, precision));
    });
  }
  // to avoid intervals, which include zero
  return calculateNewInterval(newPolynomial, function (precision) {
    return root.toDecimal(precision).add(toSimpleInterval(k, precision));
  });
};

PolynomialRoot.prototype.multiply = function (other) {
  //TODO: remove
if (true) {
  var that = this;
  let g = Math.gcd(that.polynomial.getGCDOfTermDegrees(), other.polynomial.getGCDOfTermDegrees());
  if (g > 1) {
    var tmp = that._pow(g).multiply(other._pow(g))._nthRoot(g);
    //TODO: TEST!!!
    if (g % 2 === 0 && that.sign() * other.sign() < 0) {
      tmp = tmp.negate();
    }
    return tmp;
  }
}
  // z = x * y, y = z / x
  //TODO: variable names
  const $z = new Expression.Polynomial(Polynomial.of(Expression.ONE).shift(1));
  const toPInZ = c => new Expression.Polynomial(Polynomial.of(c));
  var second = other.polynomial._exponentiateRoots(-1).map(toPInZ)._scaleRoots($z);
  var newPolynomial = Polynomial.resultant(that.polynomial.map(toPInZ), second).polynomial.primitivePart();
  return calculateNewInterval(newPolynomial, function (precision) {
    return that.toDecimal(precision).multiply(other.toDecimal(precision));
  });
};

PolynomialRoot.prototype.add = function (other) {
  var that = this;
  if (that.polynomial.isEven() && that.polynomial.equals(other.polynomial) && that.interval.a.equals(other.interval.b.negate()) && that.interval.b.equals(other.interval.a.negate())) {
    return new PolynomialRoot(Polynomial.of(Expression.ONE).shift(1), new SimpleInterval(Expression.ONE.negate(), Expression.ONE));
  }
  // z = x + y, y = z - x
  //TODO: variable names
  const $z = new Expression.Polynomial(Polynomial.of(Expression.ONE).shift(1));
  const toPInZ = c => new Expression.Polynomial(Polynomial.of(c));
  var second = other.polynomial._scaleRoots(Expression.ONE.negate()).map(toPInZ)._translateRoots($z);
  var newPolynomial = Polynomial.resultant(that.polynomial.map(toPInZ), second).polynomial.primitivePart();
  return calculateNewInterval(newPolynomial, function (precision) {
    return that.toDecimal(precision).add(other.toDecimal(precision));
  });
};

//TODO: remove (?)
PolynomialRoot.prototype.negate = function () {
  return new PolynomialRoot(this.polynomial._scaleRoots(Expression.ONE.negate()), this.interval.negate());
};

//TODO: remove (?)
PolynomialRoot.prototype.inverse = function () {
  // z = 1/y, y = 1/z
  var newPolynomial = this.polynomial._exponentiateRoots(-1);
  console.assert(this.interval.a.getNumerator().sign() === this.interval.b.getNumerator().sign());
  var newInterval = new SimpleInterval(this.interval.b.inverse(), this.interval.a.inverse());
  return new PolynomialRoot(newPolynomial, newInterval);
};

PolynomialRoot.prototype.sign = function () {
  if (this.polynomial.getCoefficient(0).equals(Expression.ZERO)) {
    if (this.interval.a.getNumerator().sign() <= 0 && this.interval.b.getNumerator().sign() >= 0) {
      return 0;
    }
  }
  if (this.interval.a.getNumerator().sign() >= 0) {
    return +1;
  }
  if (this.interval.b.getNumerator().sign() <= 0) {
    return -1;
  }
  throw new TypeError("should not happen");
};

PolynomialRoot.prototype._pow = function (n) {
  var pow = function (x, count, accumulator) {
    if (!(count >= 0)) {
      throw new RangeError();
    }
    if (count > Number.MAX_SAFE_INTEGER) {
      throw new RangeError("NotSupportedError");
    }
    return (count < 1 ? accumulator : (2 * Math.floor(count / 2) !== count ? pow(x, count - 1, accumulator.multiply(x)) : pow(x._pow(2), Math.floor(count / 2), accumulator)));
  };
  if (n === 0) {
    return new PolynomialRoot(Polynomial.of(Expression.ONE.negate(), Expression.ONE), new SimpleInterval(Expression.ZERO, Expression.TWO)); // x-1=0
  }
  var g = Math.gcd(n, this.polynomial.getGCDOfTermDegrees());
  if (g === 1) {
    //return Expression.prototype._pow.call(this, n);//TODO: ?
    return pow(this, n - 1, this);
  }
  if (g < n) {
    return this._pow(g)._pow(n / g);
  }
  //TODO: faster method
  var newInterval = undefined;
  if (n % 2 === 0 && this.interval.b.getNumerator().compareTo(Expression.ZERO) <= 0) {
    newInterval = {a: this.interval.b._pow(n), b: this.interval.a._pow(n)};
  } else {
    newInterval = {a: this.interval.a._pow(n), b: this.interval.b._pow(n)}
  }
  return new PolynomialRoot(this.polynomial._exponentiateRoots(n), new SimpleInterval(newInterval.a, newInterval.b));
};

const $α = function () {
  return new Expression.Symbol('α');
  //return new Expression.Polynomial(Polynomial.of(Expression.ONE).shift(1));
};

PolynomialRoot.prototype._nthRoot = function (n) {
  var newPolynomial = this.polynomial._exponentiateRoots(1 / n);
  var root = this;
  return calculateNewInterval(newPolynomial, function (precision) {
    //TODO: 
    //return root.toDecimal(precision).nthRoot(n);
    return toSimpleInterval(Expression.NthRoot.makeRoot(new Expression.ExpressionPolynomialRoot(new LazyPolynomialRoot($α(), root)), n), precision);
  });
};

PolynomialRoot.prototype.equals = function (other) {
  if (this === other) {
    return true;
  }
  if (this.polynomial.getDegree() !== other.polynomial.getDegree()) {
    return false;
  }
  if (this.polynomial.equals(other.polynomial) && this.interval.a.equals(other.interval.a) && this.interval.b.equals(other.interval.b)) {
    return true;
  }
  if (SimpleInterval.intersection(this.interval, other.interval) == null) {
    return false;
  }
  //TODO: ?
  //return this.polynomial.equals(other.polynomial) && intersection(this.interval, other.interval) != null && this.add(other.negate()).equals(Expression.ZERO);
  var interval = this.add(other.negate()).interval;
  return interval.a.getNumerator().compareTo(Expression.ZERO) <= 0 && interval.b.getNumerator().compareTo(Expression.ZERO) >= 0;
};

PolynomialRoot._calculateNewInterval = calculateNewInterval;//TODO: remove
LazyPolynomialRoot._calculateNewInterval = calculateNewInterval;//TODO: remove

PolynomialRoot._toSimpleInterval = toSimpleInterval;//TODO: remove
LazyPolynomialRoot._toSimpleInterval = toSimpleInterval;//TODO: remove

function LazyPolynomialRoot(e, _root) {
  console.assert(e instanceof Expression);//TODO: Expression.Polynomial (?)
  console.assert(_root instanceof PolynomialRoot);
  //TODO:
  //console.assert(Expression.isRealAlgebraicNumber(e));
  this.e = e; // internal symbolic expression with a "root" as a symbol
  this._root = _root;
}

LazyPolynomialRoot.prototype.toDecimal = function (precision) {
  var calcAt = function (polynomial, x, precision) {
    var result = toSimpleInterval(Expression.ZERO, precision);
    for (var i = polynomial.getDegree(); i >= 0; i--) {
      result = result.multiply(x);
      var tmp = toSimpleInterval(polynomial.getCoefficient(i), Math.max(1, precision));
      //TODO: ?
      //if (tmp === "CANNOT_DIVIDE" || tmp == undefined) {
      //  return tmp;
      //}
      result = result.add(tmp);
    }
    return result;
  };
  var alphaValue = this._root.toDecimal(precision);
  var alphaExpression = this.e;
  var p1 = Polynomial.toPolynomial(alphaExpression.getNumerator(), new Expression.Symbol('α'));
  var p2 = Polynomial.toPolynomial(alphaExpression.getDenominator(), new Expression.Symbol('α'));
  var a = calcAt(p1, alphaValue, precision);
  //if (a === "CANNOT_DIVIDE" || a == undefined) {
  //  return a;
  //}
  var b = calcAt(p2, alphaValue, precision);
  //if (b === "CANNOT_DIVIDE" || b == undefined) {
  //  return b;
  //}
  var result = a.multiply(b.inverse());
  //TODO: precision !!!
  return result;
};

LazyPolynomialRoot.prototype.toString = function () {
  return "[" + this.e + ", where " + this._root + "]"; // for debugging
};

function makeExpressionWithPolynomialRoot(e, root, variable) {
  //TODO: use cases - ?
  if (e instanceof Expression.Integer) {
    return e;
  }
  if (e instanceof Expression.Division && e.a instanceof Expression.Integer && e.b instanceof Expression.Integer) {
    return e;
  }
  //!

  if (e.getNumerator().equals(Expression.ZERO)) {
    return Expression.ZERO;
  }
  var v = root;

if (true) {
  if (!(e.getDenominator() instanceof Expression.Integer)) {
    var p = Polynomial.toPolynomial(e.getDenominator(), variable);
    if (p.getDegree() > 0) {
      return makeExpressionWithPolynomialRoot(e.getNumerator().multiply(p.modularInverse(root.polynomial).calcAt(variable)), root, variable);
    }
  }
}

  var c = function (x) {
    //!new 2020-08-27
    //TODO: remove
    //TODO: optimize
    if (true && !(x instanceof Expression.Integer) &&
        !(x instanceof Expression.Multiplication && x.a === Expression.I && x.b === v) &&
        !(x instanceof Expression.Exponentiation)) {
      var p1 = Polynomial.toPolynomial(x.subtract(new Expression.Symbol('$n')), variable);
      //var test = v.polynomial.divideAndRemainder(p1).remainder;
      if (p1.getDegree() >= v.polynomial.getDegree()) {
        p1 = Polynomial.pseudoRemainder(p1, v.polynomial);
      }
      var test = v.polynomial.getDegree() >= p1.getDegree() ? Polynomial.pseudoRemainder(v.polynomial, p1) : v.polynomial;
      if (test.getDegree() === 0) {
        //(x**2-2)(x**2+x-1) = 0
        var pn0 = Polynomial.toPolynomial(test.calcAt(Expression.ZERO).getNumerator(), new Expression.Symbol('$n'));
        var pn = Polynomial.toPolynomial(Expression.getConjugateExpression(test.calcAt(Expression.ZERO).getNumerator()), new Expression.Symbol('$n'));
        //pn = pn.scale(pn.getLeadingCoefficient().inverse());
        pn = pn.primitivePart();
        var tmp = pn.squareFreeFactors();
        var f = tmp.a0;
        if (tmp.a0.getDegree() === 0) {
          f = tmp.a1;
        }
        if (f.getDegree() <= 2) {//TODO: ?
          const roots = f.getroots();
          let c = [];
          let fractionDigits = 3;
          do {
            c.splice(0, c.length);
            for (const root of roots) {
              if (Expression.has(root, Expression.Complex) === Expression.has(x, Expression.Complex)) {
                if (root.getNumerator().toMathML({rounding: {fractionDigits: fractionDigits}}) === new Expression.Multiplication(x, root.getDenominator()).toMathML({rounding: {fractionDigits: fractionDigits}})) {
                  c.push(root);
                }
              }
            }
            fractionDigits *= 2;
          } while (c.length > 2);
          if (c.length === 1) {
            return c[0];
          }
        }
        var lc = pn.getLeadingCoefficient();
        if (lc instanceof Expression.Integer && Expression.isConstant(x)) {//TODO: ? - this is a filter to avoid infinite computation
          var s = toDecimalStringInternal(x.multiply(lc), {fractionDigits: 0});
          var n = Number(s);//TODO: complex - ?
          if (!Number.isNaN(n)) {
            var q = Expression.Integer.fromString(s).divide(lc);
            if (pn.calcAt(q).equals(Expression.ZERO) && pn0.calcAt(q).equals(Expression.ZERO)) {
              const tmp = Expression.ONE.divide(Expression.TWO).divide(lc);
              const interval = {a: q.subtract(tmp), b: q.add(tmp)};
              if (pn.numberOfRoots(interval) === 1) {
                console.debug('q.toString():', q.toString());
                return q;
              }
              //TODO: ?
              //if (pn0.numberOfRoots(interval) === 1) {
              //  debugger;
              //  console.debug('q.toString():', q.toString());
              //  return q;
              //}
            }
          }
        }
      }
    }
    var px = Polynomial.toPolynomial(x, variable);
    if (false) {
      if (v.polynomial.getDegree() === 6 && v.polynomial.getCoefficient(1).equals(Expression.ZERO) && v.polynomial.getCoefficient(2).equals(Expression.ZERO) && v.polynomial.getCoefficient(4).equals(Expression.ZERO) && v.polynomial.getCoefficient(5).equals(Expression.ZERO)) {
        if (px.getDegree() >= 3) {
          var alpha = c(v._pow(3));
          var dv = Polynomial.of(alpha.negate(), Expression.ZERO, Expression.ZERO, Expression.ONE);
          return px.divideAndRemainder(dv).remainder.calcAt(variable);
        }
      }
    }
    //!
    return px.divideAndRemainder(v.polynomial).remainder.calcAt(variable);
  };
  var oldE = e;
  e = c(e.getNumerator()).divide(c(e.getDenominator()));
  if (!oldE.equals(e)) {
    var tmp = (!(oldE.getDenominator() instanceof Expression.Integer) || !(e.getDenominator() instanceof Expression.Integer));
    if (tmp) {
    //console.log('tmp', tmp);
    var oldE1 = e;
    e = c(e.getNumerator()).divide(c(e.getDenominator())); // something may change after the previous step
    if (!oldE1.equals(e)) {
      //debugger;
    }
    }
  }


  //TODO: use polynomial from the start - ?
  if (Polynomial.toPolynomial(e.getNumerator(), variable).hasRoot(v)) {//Note: slow
    return Expression.ZERO;
  }

  return e;
}

PolynomialRoot._makeExpressionWithPolynomialRoot = makeExpressionWithPolynomialRoot;
LazyPolynomialRoot._makeExpressionWithPolynomialRoot = makeExpressionWithPolynomialRoot;

function simplifyExpressionWithPolynomialRoot(e, root) {
  return new LazyPolynomialRoot(makeExpressionWithPolynomialRoot(e, root, new Expression.Symbol('α')), root);
}
PolynomialRoot.SimpleInterval = SimpleInterval;
LazyPolynomialRoot.SimpleInterval = SimpleInterval;
PolynomialRoot.create = function (polynomial, interval, options) {
  return new PolynomialRoot(polynomial, interval, options);
};
LazyPolynomialRoot.create = function (polynomial, interval, options) {
  return new LazyPolynomialRoot(new Expression.Symbol('α'), new PolynomialRoot(polynomial, interval, options));
};
function fromRoot(root) {
  return new LazyPolynomialRoot(new Expression.Symbol('α'), root);
}
LazyPolynomialRoot.prototype.scale = function (k) {
  console.assert(Expression.isRealAlgebraicNumber(k));
  return simplifyExpressionWithPolynomialRoot(this.e.multiply(k), this._root);
};
LazyPolynomialRoot.prototype.translate = function (k) {
  console.assert(Expression.isRealAlgebraicNumber(k));
  return simplifyExpressionWithPolynomialRoot(this.e.add(k), this._root);
};

var toPolynomialWithIntegerCoefficients = function (polynomial) {
  if (!polynomial.hasIntegerCoefficients()) {
    var variable = new Expression.Symbol('$$');
    var e = polynomial.calcAt(variable);
    var c = Expression.getConjugateExpression(e);
    if (c != null && !c.equals(e)) {
      //TODO: what if multiple (?) - ?
      return Polynomial.toPolynomial(c, variable);
    }
  }
  return polynomial;
};

function upgrade(lazyRoot) {
  var root = lazyRoot._root;
  var e = lazyRoot.e;
  var variable = new Expression.Symbol('α');
  if (e.equals(variable)) { // short path
    return root;
  }
  if (e.getNumerator().equals(variable)) { // short path 2
    return root.scale(e.getDenominator().inverse());
  }
  //TODO: modInverse or upgrade numerator and denominator separately - ?
  var p = Polynomial.toPolynomial(e.getNumerator(), variable);
  var p2 = Polynomial.toPolynomial(e.getDenominator(), variable);
  var scale = Expression.ONE;
  if (p2.getDegree() === 0 && p2.hasIntegerCoefficients()) {
    scale = p2.getLeadingCoefficient();
    p2 = Polynomial.of(Expression.ONE);
  }
  var polynomial = p.subtract(Polynomial.of(new Expression.Symbol('β')).multiply(p2));
  polynomial = toPolynomialWithIntegerCoefficients(polynomial);//TODO: ???
  const toPInBeta = c => new Expression.Polynomial(Polynomial.of(c));
  polynomial = polynomial.map(c => new Expression.Polynomial(Polynomial.toPolynomial(c, new Expression.Symbol('β'))));//TODO: ?
  var newPolynomial = Polynomial.resultant(polynomial, root.polynomial.map(toPInBeta)).polynomial.primitivePart();
  if (scale !== Expression.ONE) {
    // "unscale"
    newPolynomial = newPolynomial._scaleRoots(scale.inverse()).primitivePart();
  }
  return PolynomialRoot._calculateNewInterval(newPolynomial, function (precision) {
    return toSimpleInterval(new Expression.ExpressionPolynomialRoot(lazyRoot), precision);
  });
}
LazyPolynomialRoot.prototype.multiply = function (other) {
  if (this._root.equals(other._root)) {
    return simplifyExpressionWithPolynomialRoot(this.e.multiply(other.e), this._root);
  }
  var root = upgrade(this).multiply(upgrade(other));
  return fromRoot(root);
};
LazyPolynomialRoot.prototype.add = function (other) {
  if (this._root.equals(other._root)) {
    return simplifyExpressionWithPolynomialRoot(this.e.add(other.e), this._root);
  }
  var root = upgrade(this).add(upgrade(other));
  return fromRoot(root);
};
LazyPolynomialRoot.prototype.negate = function () {
  return new LazyPolynomialRoot(this.e.negate(), this._root);
};
LazyPolynomialRoot.prototype.inverse = function () {
  return simplifyExpressionWithPolynomialRoot(this.e.inverse(), this._root);
};
LazyPolynomialRoot.prototype.sign = function () {
  if (this.e.equals(Expression.ZERO)) {
    return 0;
  }
  //return this.e;
  //?
  //TODO: ???
  //var s = toDecimalStringInternal(new Expression.ExpressionPolynomialRoot(this), {significantDigits: 1});
  //return s.startsWith('-') ? -1 : +1;
  let precision = 1;
  while (true) {
    var interval = this.toDecimal(precision);
    if (interval.a.getNumerator().sign() >= 0) {
      return +1;
    }
    if (interval.b.getNumerator().sign() <= 0) {
      return -1;
    }
    precision *= 2;
    console.info('hm...');
  }
};
LazyPolynomialRoot.prototype._pow = function (n) {
  //TODO: modular exponentiation (?)
  return simplifyExpressionWithPolynomialRoot(this.e.getNumerator()._pow(n), this._root).multiply(simplifyExpressionWithPolynomialRoot(this.e.getDenominator()._pow(n), this._root).inverse());
};
LazyPolynomialRoot.prototype._nthRoot = function (n) {//?
  if (this.e.equals(new Expression.Symbol('α'))) {//TODO: ?
    if (this._root.interval.a.getNumerator().compareTo(Expression.ZERO) < 0) {
      throw new RangeError();
    }
    var newRoot = this._root._nthRoot(n);
    return new LazyPolynomialRoot(new Expression.Symbol('α'), newRoot);
  }
  if (!(this.e instanceof Expression.Exponentiation)) {
    if (true && n === 2) {
      return this.upgrade()._nthRoot(n);
    }
  }
  return simplifyExpressionWithPolynomialRoot(this.e.a.pow(this.e.b.divide(Expression.Integer.fromNumber(n))), this._root);
};
LazyPolynomialRoot.prototype.equals = function (other) {
  if (this._root.equals(other._root)) {
    //TODO:? ?
    return this.e.equals(other.e) || simplifyExpressionWithPolynomialRoot(this.e.subtract(other.e), this._root).e.equals(Expression.ZERO);
  }
  //!TODO: remove (hack to avoid error)
  if (toDecimalStringInternal(new Expression.ExpressionPolynomialRoot(this), {significantDigits: 3}) !== toDecimalStringInternal(new Expression.ExpressionPolynomialRoot(other), {significantDigits: 3})) {
    return false;
  }
  //!
  var result = upgrade(this).equals(upgrade(other));
  return result;
};

PolynomialRoot.prototype.toPolynomialRoot = function () {//TODO: rename - ?
  return this;
};
LazyPolynomialRoot.prototype.toPolynomialRoot = function () {//TODO: rename - ?
  return upgrade(this);
};
PolynomialRoot.prototype.upgrade = function () {
  return this;
};
LazyPolynomialRoot.prototype.upgrade = function () {
  return fromRoot(upgrade(this));
};

//LazyPolynomialRoot.PolynomialRoot = PolynomialRoot;//TODO: REMOVE!!!

globalThis.testables = globalThis.testables || {};
globalThis.testables.LazyPolynomialRoot = LazyPolynomialRoot;
globalThis.testables.PolynomialRoot = PolynomialRoot;
globalThis.testables.toSimpleIntervalOld = toSimpleIntervalOld;
globalThis.testables.toSimpleIntervalNew = toSimpleIntervalNew;

if (false) {//TODO: move to tests
  console.assert(Object.keys(PolynomialRoot).join(' ') === Object.keys(LazyPolynomialRoot).join(' '));
  console.assert(Object.keys(PolynomialRoot.prototype).join(' ') === Object.keys(LazyPolynomialRoot.prototype).join(' '));
  console.assert(PolynomialRoot.prototype.__proto__ === LazyPolynomialRoot.prototype.__proto__);
}

export default LazyPolynomialRoot;
