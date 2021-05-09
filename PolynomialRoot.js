import primeFactor from './primeFactor.js';

// a class for real algebraic numbers
// Operations are implemented as described at https://en.wikipedia.org/wiki/Resultant#Number_theory

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
SimpleInterval.prototype.toString = function () {
  return '[' + this.a + ';' + this.b + ']';
};

var toSimpleInterval = function (e, precision) {
  if (e instanceof Expression.Integer) {
    return new SimpleInterval(e, e);
  } else if (e instanceof Expression.BinaryOperation) {
    var a = toSimpleInterval(e.a, precision);
    var b = toSimpleInterval(e.b, precision);
    var s = e.getS();
    if (s === "+") {
      return a.add(b);
    } else if (s === "-") {
      return a.add(b.negate());
    } else if (s === "*") {
      return a.multiply(b);
    //} else if (s === "/") {
      //return a.multiply(b.inverse());
    } else {
      debugger;
    }
  } else if (e instanceof Expression.NthRoot) {
    if (e.a instanceof Expression.Integer) {
      var a = e.a;
      var n = e.n;
      var scale = Expression.TWO._pow(precision);
      return new SimpleInterval(nthRoot(a.multiply(scale._pow(n)).toBigInt(), n), scale);
    }
    //var a = toSimpleInterval(e.a, precision);
    //return a._nthRoot(e.n);
    debugger;
  } else {
    debugger;
  }
  throw new TypeError("?");
};

  var calculateNewInterval = function (newPolynomial, zeroFunction) {
    if (!newPolynomial.hasIntegerCoefficients()) {
      throw new RangeError("just a check");
    }
    var precision = 1;
    var guess = zeroFunction(precision);
    var sturmSequence = new SturmSequence(newPolynomial);
    while ((guess.a.getNumerator().sign() !== guess.b.getNumerator().sign() && !newPolynomial.calcAt(Expression.ZERO).equals(Expression.ZERO)) || sturmSequence.numberOfRoots(guess) > 1) {
      precision *= 2;
      guess = zeroFunction(precision);
      if (precision > 1024) throw new Error();//TODO: ?
    }
    return guess;
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
  if (interval.a.subtract(interval.b).getNumerator().compareTo(Expression.ZERO) >= 0) {
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

PolynomialRoot.SimpleInterval = SimpleInterval;

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
  var n = this.polynomial.getDegree();
  var newPolynomial = this.polynomial.map((coefficient, degree) => k._pow(n - degree).multiply(coefficient)).primitivePart();
  var newInterval = null;
  if (!(isRational(k))) { // TODO: remove
    var root = this;
    newPolynomial = toPolynomialWithIntegerCoefficients(newPolynomial);
    newInterval = calculateNewInterval(newPolynomial, function (precision) {
      //var e1 = toSimpleInterval(k, precision);
      //debugger;
      var tmp = ExpressionParser.parse(toDecimalStringInternal(k, {significantDigits: precision}));
      var e = new SimpleInterval(ExpressionParser.parse('1-5*10**-' + precision), ExpressionParser.parse('1+5*10**-' + precision));
      return root.toDecimal(precision).multiply(new SimpleInterval(tmp, tmp).multiply(e));
    });
  } else {
    newInterval = this.interval.scale(k);
  }
  return new PolynomialRoot(newPolynomial, newInterval);
};

//TODO: remove (?)
PolynomialRoot.prototype.translate = function (k) {
  //console.assert(k instanceof Expression.Integer || isRational(k));//TODO: ???
  // z = x + k, x = z - k
  var newPolynomial = this.polynomial.subs(x => x.subtract(k)).primitivePart();
  // to avoid intervals, which include zero
  var root = this;
  var newInterval = null;
  if (!(isRational(k))) { // TODO: remove
    newPolynomial = toPolynomialWithIntegerCoefficients(newPolynomial);
    newInterval = calculateNewInterval(newPolynomial, function (precision) {
      var tmp = ExpressionParser.parse(toDecimalStringInternal(k, {significantDigits: precision}));
      var e = new SimpleInterval(ExpressionParser.parse('1-5*10**-' + precision), ExpressionParser.parse('1+5*10**-' + precision));
      return root.toDecimal(precision).add(new SimpleInterval(tmp, tmp).multiply(e));
    });
  } else {
    // to avoid intervals, which include zero
    newInterval = calculateNewInterval(newPolynomial, function (precision) {
      return root.toDecimal(precision).add(new SimpleInterval(k, k));
    });
  }
  return new PolynomialRoot(newPolynomial, newInterval);
};

PolynomialRoot.prototype.multiply = function (other) {
  var x = this;
  var y = other;
  // z = x * y, y = z / x
  //TODO: variable names
  var second = Polynomial.toPolynomial(y.polynomial.calcAt(ExpressionParser.parse("z/x")).multiply(ExpressionParser.parse("x")._pow(y.polynomial.getDegree())), ExpressionParser.parse("x"));
  var newPolynomial = Polynomial.resultant(x.polynomial, second, "z").primitivePart();
  var newInterval = calculateNewInterval(newPolynomial, function (precision) {
    return x.toDecimal(precision).multiply(y.toDecimal(precision));
  });
  return new PolynomialRoot(newPolynomial, newInterval);
};

PolynomialRoot.prototype.add = function (other) {
  var x = this;
  var y = other;
  if (x.polynomial.equals(y.polynomial) && x.interval.a.equals(y.interval.b.negate()) && x.interval.b.equals(y.interval.a.negate())) {
    return new PolynomialRoot(Polynomial.of(Expression.ONE).shift(1), new SimpleInterval(Expression.ONE.negate(), Expression.ONE));
  }
  //TODO: variable names
  var second = Polynomial.toPolynomial(y.polynomial.calcAt(ExpressionParser.parse("z-x")), ExpressionParser.parse("x"));
  var newPolynomial = Polynomial.resultant(x.polynomial, second, "z").primitivePart();
  var newInterval = calculateNewInterval(newPolynomial, function (precision) {
    return x.toDecimal(precision).add(y.toDecimal(precision));
  });
  return new PolynomialRoot(newPolynomial, newInterval);
};

//TODO: remove (?)
PolynomialRoot.prototype.negate = function () {
  return new PolynomialRoot(this.polynomial.subs(x => x.negate()), this.interval.negate());
};

//TODO: remove (?)
PolynomialRoot.prototype.inverse = function () {
  // z = 1/y, y = 1/z
  var newPolynomial = this.polynomial.subs(x => x.inverse());
  console.assert(this.interval.a.getNumerator().sign() === this.interval.b.getNumerator().sign());
  var newInterval = new SimpleInterval(this.interval.b.inverse(), this.interval.a.inverse());
  return new PolynomialRoot(newPolynomial, newInterval);
};

PolynomialRoot.prototype.sign = function () {
  if (this.interval.a.getNumerator().compareTo(Expression.ZERO) > 0) {
    return +1;
  }
  if (this.interval.b.getNumerator().compareTo(Expression.ZERO) < 0) {
    return -1;
  }
  throw new TypeError("should not happen");
};

PolynomialRoot.prototype._pow = function (n) {
  var gcd = function (a, b) { return b === 0 ? a : gcd(b, a % b); };
  var pow = function (x, count, accumulator) {
    if (!(count >= 0)) {
      throw new RangeError();
    }
    if (count > Number.MAX_SAFE_INTEGER) {
      throw new RangeError("NotSupportedError");
    }
    return (count < 1 ? accumulator : (2 * Math.floor(count / 2) !== count ? pow(x, count - 1, accumulator.multiply(x)) : pow(x._pow(2), Math.floor(count / 2), accumulator)));
  };
  var g = gcd(n, this.polynomial.getGCDOfTermDegrees());
  if (g === 1) {
    if (n === 0) {
      return new PolynomialRoot(Polynomial.of(Expression.ONE.negate(), Expression.ONE), new SimpleInterval(Expression.ZERO, Expression.TWO)); // x-1=0
    }
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
  return new PolynomialRoot(this.polynomial.subs(x => x._nthRoot(n)), new SimpleInterval(newInterval.a, newInterval.b));
};

PolynomialRoot.prototype._nthRoot = function (n) {
  var newPolynomial = this.polynomial.subs(x => x._pow(n));
  var root = this;
  var newInterval = calculateNewInterval(newPolynomial, function (precision) {
    //TODO: 
    //return root.toDecimal(precision).nthRoot(n);
    const value = ExpressionParser.parse(toDecimalStringInternal(Expression.NthRoot.makeRoot(new Expression.ExpressionPolynomialRoot(new LazyPolynomialRoot(new Expression.Symbol('α'), root)), n), {significantDigits: precision}));
    var a = value.multiply(ExpressionParser.parse('1-2**-' + precision));
    var b = value.multiply(ExpressionParser.parse('1+2**-' + precision));
    return new SimpleInterval(a, b);
  });
  return new PolynomialRoot(newPolynomial, newInterval);
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
  var intersection = function (a, b) {
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
  if (intersection(this.interval, other.interval) == null) {
    return false;
  }
  //TODO: ?
  //return this.polynomial.equals(other.polynomial) && intersection(this.interval, other.interval) != null && this.add(other.negate()).equals(Expression.ZERO);
  var interval = this.add(other.negate()).interval;
  return interval.a.getNumerator().compareTo(Expression.ZERO) <= 0 && interval.b.getNumerator().compareTo(Expression.ZERO) >= 0;
};

PolynomialRoot._calculateNewInterval = calculateNewInterval;//TODO: remove
LazyPolynomialRoot._calculateNewInterval = calculateNewInterval;//TODO: remove

function LazyPolynomialRoot(e, _root) {
  console.assert(e instanceof Expression);
  console.assert(_root instanceof PolynomialRoot);
  //TODO:
  //console.assert(Expression.isConstant(e) && !Expression.has(e, Expression.Complex));
  this.e = e; // internal symbolic expression with a "root" as a symbol
  this._root = _root;
}


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
    return makeExpressionWithPolynomialRoot(e.getNumerator().multiply(p.modularInverse(root.polynomial).calcAt(variable)), root, variable);
  }
}

  var c = function (x) {
    //!new 2020-08-27
    //TODO: remove
    if (true && !(x instanceof Expression.Integer) && !(x instanceof Expression.Multiplication && x.a === Expression.I && x.b === v)) {
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
            if (pn.calcAt(q).equals(Expression.ZERO)) {
              if (pn0.numberOfRoots({a: q.subtract(ExpressionParser.parse('0.5').divide(lc)), b: q.add(ExpressionParser.parse('0.5').divide(lc))}) === 1) {
                console.debug(q.toString());
                return q;
              }
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

LazyPolynomialRoot._makeExpressionWithPolynomialRoot = makeExpressionWithPolynomialRoot;

function simplifyExpressionWithPolynomialRoot(e, root) {
  
  //var variable = new Expression.Symbol('α');
  //TODO: ... see makeExpressionWithPolynomialRoot
  //var p = Polynomial.toPolynomial(e, variable);
  //p = p.divideAndRemainder(root.polynomial).remainder;
  //e = p.calcAt(variable);
  //return new LazyPolynomialRoot(e, root);
  return new LazyPolynomialRoot(makeExpressionWithPolynomialRoot(e, root, new Expression.Symbol('α')), root);
}
LazyPolynomialRoot.SimpleInterval = SimpleInterval;
LazyPolynomialRoot.create = function (polynomial, interval, options) {
  return new LazyPolynomialRoot(new Expression.Symbol('α'), new PolynomialRoot(polynomial, interval, options));
};
function fromRoot(root) {
  return new LazyPolynomialRoot(new Expression.Symbol('α'), root);
}
LazyPolynomialRoot.prototype.scale = function (k) {
  console.assert(Expression.isConstant(k));//TODO: ??? only some constants
  return simplifyExpressionWithPolynomialRoot(this.e.multiply(k), this._root);
};
LazyPolynomialRoot.prototype.translate = function (k) {
  console.assert(Expression.isConstant(k));//TODO: ??? only some constants
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
  //TODO: modInverse or upgrade numerator and denominator separately - ?
  var p = Polynomial.toPolynomial(e.getNumerator(), variable);
  var p2 = Polynomial.toPolynomial(e.getDenominator(), variable);
  var scale = Expression.ONE;
  if (p2.getDegree() === 0 && p2.hasIntegerCoefficients()) {
    scale = p2.getLeadingCoefficient();
    p2 = Polynomial.of(Expression.ONE);
  }
  var polynomial = p.subtract(Polynomial.of(ExpressionParser.parse('β')).multiply(p2));
  polynomial = toPolynomialWithIntegerCoefficients(polynomial);//TODO: ???
  var newPolynomial = Polynomial.resultant(polynomial, root.polynomial, 'β').primitivePart();
  if (scale !== Expression.ONE) {
    // "unscale"
    newPolynomial = newPolynomial.map((coefficient, degree) => scale._pow(degree).multiply(coefficient)).primitivePart();
  }
  var newInterval = PolynomialRoot._calculateNewInterval(newPolynomial, function (precision) {
    var tmp = ExpressionParser.parse(toDecimalStringInternal(new Expression.ExpressionPolynomialRoot(lazyRoot), {significantDigits: precision}));
    var epsilonInterval = new PolynomialRoot.SimpleInterval(ExpressionParser.parse('1-5*10**-' + precision), ExpressionParser.parse('1+5*10**-' + precision));
    return epsilonInterval.scale(tmp);
  });
  return new PolynomialRoot(newPolynomial, newInterval);
}
LazyPolynomialRoot.prototype.add = function (other) {
  if (this._root.equals(other._root)) {
    return simplifyExpressionWithPolynomialRoot(this.e.add(other.e), this._root);
  }
  var root = upgrade(this).add(upgrade(other));
  return fromRoot(root);
};
LazyPolynomialRoot.prototype.multiply = function (other) {
  if (this._root.equals(other._root)) {
    return simplifyExpressionWithPolynomialRoot(this.e.multiply(other.e), this._root);
  }
  var root = upgrade(this).multiply(upgrade(other));
  return fromRoot(root);
};
LazyPolynomialRoot.prototype._pow = function (n) {
  //TODO: modular exponentiation (?)
  return simplifyExpressionWithPolynomialRoot(this.e.getNumerator()._pow(n), this._root).multiply(simplifyExpressionWithPolynomialRoot(this.e.getDenominator()._pow(n), this._root).inverse());
};
LazyPolynomialRoot.prototype.inverse = function () {
  return simplifyExpressionWithPolynomialRoot(this.e.inverse(), this._root);
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
LazyPolynomialRoot.prototype.upgrade = function () {
  return fromRoot(upgrade(this));
};
LazyPolynomialRoot.prototype.sign = function () {
  if (this.e.equals(Expression.ZERO)) {
    return 0;
  }
  //return this.e;
  //?
  //TODO: ???
  var s = toDecimalStringInternal(new Expression.ExpressionPolynomialRoot(this), {significantDigits: 1});
  return s.startsWith('-') ? -1 : +1;
};

LazyPolynomialRoot.prototype.negate = function () {
  return new LazyPolynomialRoot(this.e.negate(), this._root);
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
  return simplifyExpressionWithPolynomialRoot(this.e._nthRoot(n), this._root);
};

LazyPolynomialRoot.PolynomialRoot = PolynomialRoot;//TODO: REMOVE!!!

export default LazyPolynomialRoot;
