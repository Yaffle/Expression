import Expression from './Expression.js';
import Polynomial from './Polynomial.js';
import toDecimalStringInternal from './toDecimalString.js';

function PolynomialRoot(polynomial, interval) {
  Expression.Symbol.call(this, "@");
  this.polynomial = polynomial;
  //TODO: https://www.wolframalpha.com/input/?i=x**5%2B7x**3%2Bx**2%2Bx%2B1%3D0
  this.interval = interval;
}
PolynomialRoot.prototype = Object.create(Expression.Symbol.prototype);

PolynomialRoot.prototype.toDecimal = function (precision) {
  return this.polynomial.getZero(this.interval, precision);
};

Expression.PolynomialRoot = PolynomialRoot;





function ExpressionWithPolynomialRoot(e, root) {
  this.e = e; // internal symbolic expression with a "root" as a symbol
  this.root = root;
}

function makeExpressionWithPolynomialRoot(e, root) {
  //TODO: use cases - ?
  if (e instanceof Expression.Integer) {
    return e;
  }
  if (e instanceof Expression.Division && e.a instanceof Expression.Integer && e.b instanceof Expression.Integer) {
    return e;
  }
  //!

  var en = e.getNumerator();
  if (en.equals(Expression.ZERO)) {
    return Expression.ZERO;
  }
  var v = root;

  //TODO: use polynomial from the start - ?
  if (Polynomial.toPolynomial(en, v).hasRoot(v)) {
    return Expression.ZERO;
  }
  var c = function (x) {
    return Polynomial.toPolynomial(x, v).divideAndRemainder(v.polynomial).remainder.calcAt(v);
  };
  var oldE = e;
  e = c(e.getNumerator()).divide(c(e.getDenominator()));
  if (e instanceof Expression.Integer) {
    return e;
  }
  if (e instanceof Expression.Division && e.a instanceof Expression.Integer && e.b instanceof Expression.Integer) {
    return e;
  }

  return new ExpressionWithPolynomialRoot(e, root);
}

ExpressionWithPolynomialRoot.prototype = Object.create(Expression.Symbol.prototype);

//Expression.prototype.isExact = function () {
//  return true;
//};
ExpressionWithPolynomialRoot.prototype.isExact = function () {
  //TODO: fix - ?
  return false;
};

ExpressionWithPolynomialRoot.prototype.negate = function () {
  return makeExpressionWithPolynomialRoot(this.e.negate(), this.root);
};
ExpressionWithPolynomialRoot.prototype.equals = function (other) {
  // optimization
  var s = other instanceof Expression.Integer && other.equals(Expression.ZERO) ? this : this.subtract(other);
  return s instanceof ExpressionWithPolynomialRoot ? false : s.equals(Expression.ZERO);
};
ExpressionWithPolynomialRoot.prototype.simplifyExpression = function () {
  return this;
};

ExpressionWithPolynomialRoot.prototype.toString = function (options) {
  options = options || {};
  //TODO:
  if (this.equals(Expression.ZERO)) {
    return Expression.ZERO.toString(options);
  }
  //return this.e.toString(options);
  if (options.fractionDigits !== -1 && options.fractionDigits != null) {
    console.debug('options.fractionDigits is deprecated, please use options.rounding');
  }
  var rounding = options.rounding != null ? options.rounding : (options.fractionDigits !== -1 && options.fractionDigits != null ? {fractionDigits: options.fractionDigits} : {fractionDigits: 3});
  var tmp = toDecimalStringInternal(this.e, rounding, undefined, undefined);
  return tmp;
};

ExpressionWithPolynomialRoot.prototype.toMathML = function (options) {
  options = options || {};
  //TODO:
  if (this.equals(Expression.ZERO)) {
    return Expression.ZERO.toMathML(options);
  }
  //return this.e.toMathML(options);
  if (options.fractionDigits !== -1 && options.fractionDigits != null) {
    console.debug('options.fractionDigits is deprecated, please use options.rounding');
  }
  var rounding = options.rounding != null ? options.rounding : (options.fractionDigits !== -1 && options.fractionDigits != null ? {fractionDigits: options.fractionDigits} : {fractionDigits: 3});
  var tmp = toDecimalStringInternal(this.e, rounding, Expression._decimalToMathML, Expression._complexToMathML);
  return tmp;
};

ExpressionWithPolynomialRoot.prototype.multiply = function (other) {
  if (other instanceof ExpressionWithPolynomialRoot) {
    if (this.root !== other.root) {
      throw new RangeError("NotSupportedError");
    }
    return this.multiply(other.e);
  }
  return makeExpressionWithPolynomialRoot(this.e.multiply(other), this.root);
};
ExpressionWithPolynomialRoot.prototype.divide = function (other) {
  if (other instanceof ExpressionWithPolynomialRoot) {
    if (this.root !== other.root) {
      throw new RangeError("NotSupportedError");
    }
    return this.divide(other.e);
  }
  return makeExpressionWithPolynomialRoot(this.e.divide(other), this.root);
};
ExpressionWithPolynomialRoot.prototype.add = function (other) {
  if (other instanceof ExpressionWithPolynomialRoot) {
    if (this.root !== other.root) {
      throw new RangeError("NotSupportedError");
    }
    return this.add(other.e);
  }
  return makeExpressionWithPolynomialRoot(this.e.add(other), this.root);
};

ExpressionWithPolynomialRoot.prototype.divideExpression = function (other) {
  return makeExpressionWithPolynomialRoot(other.divide(this.e), this.root);
};
ExpressionWithPolynomialRoot.prototype.multiplyExpression = function (other) {
  return makeExpressionWithPolynomialRoot(other.multiply(this.e), this.root);
};
ExpressionWithPolynomialRoot.prototype.addExpression = function (other) {
  return makeExpressionWithPolynomialRoot(other.add(this.e), this.root);
};

ExpressionWithPolynomialRoot.prototype.getPrecedence = function () {
  return 1000;
};
ExpressionWithPolynomialRoot.prototype.isRightToLeftAssociative = function () {
  return true;
};
ExpressionWithPolynomialRoot.prototype.isUnaryPlusMinus = function () {
  return true;
};

ExpressionWithPolynomialRoot.prototype._nthRoot = function (n) {//?
  return makeExpressionWithPolynomialRoot(this.e._nthRoot(n), this.root);
};

  // https://math.stackexchange.com/questions/309178/polynomial-root-finding
  function SturmSequence(f) {
    var d = f.derive();
    this.f = f;
    this.s = [];
    var fp = f;
    var fc = d;
    this.s.push(fp);
    this.s.push(fc);
    while (fc.getDegree() > 0) {
      var fn = fp.divideAndRemainder(fc).remainder.negate();
      if (fn.getDegree() >= 0) {
        var y = fn.getContent().inverse();
        if (y.isNegative()) {
          y = y.negate();
        }
        this.s.push(fn.scale(y));
        //?
        //this.s.push(fn);
      }
      fp = fc;
      fc = fn;
    }
  }

  SturmSequence.prototype.signChanges = function (x) {
    var result = 0;
    var sign = 0;
    for (var i = 0; i < this.s.length; i += 1) {
      var p = this.s[i];
      //var v = p.calcAt(x);
      //! 2018-10-15
      var v = Expression.ZERO;
      if (p.getDegree() >= 0) {
        var n = p.getDegree();
        var e = x.getDenominator();
        v = p.map(function (coefficient, degree) {
          return coefficient.multiply(Expression.pow(e, n - degree));
        }).calcAt(x.getNumerator());
      }
      //!
      //!?2020-01-17
      if (!(v instanceof Expression.Integer)) {
        //TODO: move to Expression._isPositive - ?
        if (v instanceof Expression.Addition && v.a instanceof Expression.Multiplication && v.b instanceof Expression.Integer) {
          if (v.a.isNegative() !== v.b.isNegative()) {
            var c = Expression.getConjugate(v);
            if (c.isNegative()) {
              c = c.negate();
            }
            v = v.multiply(c);//?
          } else {
            if (!(v.multiply(Expression.getConjugate(v)) instanceof Expression.Integer)) {
              throw new TypeError();
            }
            v = v.a.isNegative() ? Expression.ONE.negate() : Expression.ONE;//TODO: FIX
          }
        }
        if (!(v instanceof Expression.Integer)) {
          var isPositive = Expression._isPositive(v);
          if (isPositive === true) {
            v = Expression.ONE;
          }
          if (isPositive === false) {
            v = Expression.ONE.negate();
          }
        }
      }
      //!?
      var c = v.getNumerator().compareTo(Expression.ZERO);
      if (c !== 0) {
        if (sign === 0) {
          sign = c;
        } else {
          if (sign !== c) {
            sign = c;
            result += 1;
          }
        }
      }
    }
    return result;
  };

  SturmSequence.prototype.numberOfRoots = function (interval) {
    if (interval.a.equals(interval.b)) {
      throw new TypeError();
    }
    return this.signChanges(interval.a) - this.signChanges(interval.b);
  };

  SturmSequence.prototype.getRootIntervals = function () {
    var that = this;
    var interval = this.f.getRootsInterval();
    var getIntervals = function (interval) {
      var n = that.numberOfRoots(interval);
      if (n === 1) {
        return [interval];
      }
      if (n > 1) {
        var middle = interval.a.add(interval.b).divide(Expression.TWO);
        var a = getIntervals({a: interval.a, b: middle});
        var b = getIntervals({a: middle, b: interval.b});
        return a.concat(b);
      }
      return [];
    };
    return getIntervals(interval);
  };

  var abs = function (i) {
    return i.compareTo(Expression.ZERO) < 0 ? i.negate() : i;
  };

  Polynomial.prototype.getRootsInterval = function () {
    //TODO: only integer coefficients (?)
    // https://en.wikipedia.org/wiki/Sturm%27s_theorem#Number_of_real_roots
    // https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#Lagrange's_and_Cauchy's_bounds
    var max = null;
    //TODO: fix the iteration
    for (var i = 0; i < this.getDegree(); i += 1) {
      var c = abs(this.getCoefficient(i));
      if (max == null || max.compareTo(c) < 0) {
        max = c;
      }
    }
    var M = Expression.ONE.add(max.divide(abs(this.getLeadingCoefficient())));
    return {a: M.negate(), b: M};
  };

  //TODO: BigDecimal - ?, rounding - ?
  Polynomial.prototype.getZero = function (interval, precision) {
    var e = Expression.pow(Expression.TEN, precision); // epsilon^-1
    if (!(e instanceof Expression.Integer)) {
      throw new RangeError("epsilon^-1 is not an integer");
    }
    var a = interval.a;
    var b = interval.b;
    if (b.subtract(a).multiply(e).subtract(Expression.ONE).getNumerator().compareTo(Expression.ZERO) > 0) {
      //TODO:
      var tmp = true && precision >= 16 ? this.getZero(interval, Math.floor(precision / 4)) : interval;
      a = tmp.a;
      b = tmp.b;

      var n = this.getDegree();
      var p = this.map(function (coefficient, degree) {
        return coefficient.multiply(Expression.pow(e, n - degree));
      });
      var sa = a.getNumerator().multiply(e).truncatingDivide(a.getDenominator()).add(Expression.ONE);//?
      var sb = b.getNumerator().multiply(e).truncatingDivide(b.getDenominator());//?
      if (p.calcAt(sa).getNumerator().compareTo(Expression.ZERO) !== this.calcAt(a).getNumerator().compareTo(Expression.ZERO)) {
        return {a: a, b: sa.divide(e)};
      }
      if (p.calcAt(sb).getNumerator().compareTo(Expression.ZERO) !== this.calcAt(b).getNumerator().compareTo(Expression.ZERO)) {
        return {a: sb.divide(e), b: b};
      }
      a = sa;
      b = sb;
      // bisection method
      var pa = p.calcAt(a).getNumerator();
      var pb = p.calcAt(b).getNumerator();
      if (pa.compareTo(Expression.ZERO) === 0) {
        return {a: pa, b: pa};
      }
      if (pb.compareTo(Expression.ZERO) === 0) {
        return {a: pb, b: pb};
      }
      if (pa.compareTo(Expression.ZERO) === pb.compareTo(Expression.ZERO)) {
        throw new RangeError();//?
      }
      var cc = 0;
      var d = p.derive();
      while (b.subtract(a).compareTo(Expression.ONE) > 0) {// b - a > 1
        var middle = a.add(b).truncatingDivide(Expression.TWO);
        //console.log(eval(a.divide(e).toString()) + ' - ' + eval(b.divide(e).toString()));
        //?
        if (cc % 3 !== 2 && b.subtract(a).compareTo(e) < 0) {
          // TODO: better guesses
          // Newton's method
          var x = cc % 3 === 1 ? a : b;
          var c = d.calcAt(x);
          if (!c.equals(Expression.ZERO)) {
            x = x.multiply(c).subtract(p.calcAt(x)).truncatingDivide(c);
            if (x.compareTo(a) <= 0) {
              x = a.add(Expression.ONE);
            }
            if (x.compareTo(b) >= 0) {
              x = b.subtract(Expression.ONE);
            }
            //console.log("N: " + a + "  - " + x);
            middle = x;
          }
        }
        cc += 1;
        //?
        var v = p.calcAt(middle).getNumerator();
        if (v.compareTo(Expression.ZERO) === pb.compareTo(Expression.ZERO)) {
          b = middle;
        } else if (v.compareTo(Expression.ZERO) === pa.compareTo(Expression.ZERO)) {
          a = middle;
        } else {
          a = middle;
          b = middle;
        }
      }
      //console.log(cc);
      a = a.divide(e);
      b = b.divide(e);
    }
    return {a: a, b: b};
  };

  //var gcd = function (a, b) {
  //  return b.getDegree() === -1 ? a : gcd(b, a.divideAndRemainder(b).remainder);
  //};

  var gcd = function (a, b) {
    return b.getDegree() === -1 ? a : Polynomial.polynomialGCD(a, b);
  };

  Polynomial.prototype.hasRoot = function (polynomialRoot) {
    var f = this;
    if (f.getDegree() === -1) {
      return false;
    }
    var p = polynomialRoot.polynomial;
    var g = Polynomial.polynomialGCD(f, p);
    if (g.getDegree() < 1) {
      return false;
    }
    var i = polynomialRoot.interval;
    var sturmSequence = new SturmSequence(g);
    return sturmSequence.numberOfRoots(i) === 1;
  };

  Polynomial.prototype.hasIntegerCoefficients = function () {
    //TODO: optimize
    for (var i = 0; i <= this.getDegree(); i += 1) {
      if (!(this.getCoefficient(i) instanceof Expression.Integer)) {
        return false;
      }
    }
    return true;
  };

  //TODO: optimize - ?
  /*
  var _countMultiplicities = function (multiplicities, rootIntervals, f) {
    var d = f.derive();
    var g = gcd(f, d);
    var p = f.divideAndRemainder(g).quotient;
    //!
    p = p.scale(p.getContent().inverse());
    //!
    var sturmSequence = new SturmSequence(p);
    for (var i = 0; i < rootIntervals.length; i += 1) {
      if (sturmSequence.numberOfRoots(rootIntervals[i]) === 1) {
        multiplicities[i] += 1;
      }
    }
    if (g.getDegree() > 0) {
      _countMultiplicities(multiplicities, rootIntervals, g);
    }
  };
  */

  // Polynomial.toPolynomial(RPN("x^3-8x^2+21x-18"), RPN("x")).getZeros(3).result.toString()
  Polynomial.prototype.getZeros = function (precision, complex) {
    precision = precision || 1;
    complex = complex || false;
    if (this.getCoefficient(0).equals(Expression.ZERO)) {
      if (this.getLeadingCoefficient().equals(Expression.ZERO)) {
        throw new TypeError();
      }
      var i = 0;
      while (this.getCoefficient(i).equals(Expression.ZERO)) {
        i += 1;
      }
      var tmp = this.divideAndRemainder(Polynomial.of(Expression.ONE).shift(i)).quotient.getZeros(precision, complex);
      return {
        result: tmp.result.concat([Expression.ZERO]),
        multiplicities: tmp.multiplicities.concat([i])
      };
    }
    //TODO: test
    var content = this.getContent();
    var f = this.scale(content.getDenominator()).divideAndRemainder(Polynomial.of(content.getNumerator()), "throw").quotient;

    // https://en.wikipedia.org/wiki/Square-free_polynomial
    var d = f.derive();
    var a0 = gcd(f, d);
    var p = f.divideAndRemainder(a0).quotient;

    if (a0.getDegree() !== 0) {
      p = p.divideAndRemainder(gcd(p, a0)).quotient; // roots with multiplicity = 1 (?)
      var tmp1 = p.getZeros(precision, complex);
      var tmp2 = a0.getZeros(precision, complex);
      return {
        result: tmp1.result.concat(tmp2.result),
        multiplicities: tmp1.multiplicities.concat(tmp2.multiplicities.map(function (m) {
          return m + 1;
        }))
      };
    }

    if (p.getDegree() === 0) {
      return {
        result: [],
        multiplicities: []
      };
    }

    //!
    p = p.scale(p.getContent().inverse());
    //!

    if (!f.hasIntegerCoefficients()) {
      //?new
      var variable = new Expression.Symbol('$$')
      var e = f.calcAt(variable);
      var c = Expression.getConjugate(e);
      if (c != null) {
        var result = [];
        var tmp = Polynomial.toPolynomial(c.multiply(e), variable).getZeros(precision, complex);
        for (var i = 0; i < tmp.result.length; i += 1) {
          if (tmp.result[i] instanceof ExpressionWithPolynomialRoot ? this.hasRoot(tmp.result[i].root) : this.calcAt(tmp.result[i]).equals(Expression.ZERO)) {
            result.push(tmp.result[i]);
          } else {
            //TODO:?
            console.debug(tmp.result[i].root);
          }
        }
        return {result: result, multiplicities: new Array(result.length).fill(1)};
      }
      //?
      //return [];
      return {result: [], multiplicities: []};
    }

    //!new
    if (p.getDegree() === 3) {
      //?
    }
    //!

    // https://en.wikipedia.org/wiki/Sturm%27s_theorem
    var sturmSequence = new SturmSequence(p);
    var intervals = sturmSequence.getRootIntervals();

    // https://math.stackexchange.com/questions/309178/polynomial-root-finding
    // "it is guaranteed that there is a sign change inside every interval (because there are no repeated zeroes)"
    var result = new Array(intervals.length);
    for (var i = 0; i < intervals.length; i += 1) {
      var zero = p.getZero(intervals[i], precision);
      if (zero.a.equals(zero.b)) {
        result[i] = zero.a;//TODO: fix
      } else {
        //! p, not f, as f may have roots with multiplicity > 1
        var root = new PolynomialRoot(p, zero);
        result[i] = new ExpressionWithPolynomialRoot(root, root);
      }
    }
    //return result;

    //!new
    //var p = np;
    if (intervals.length !== p.getDegree() && true && complex) {
      //var p = stringToPolynomial("x^5+2*x^2+2*x+3");

      var e = p.calcAt(RPN("a+b*i"));
      var ce = Expression.getComplexConjugate(e);
      var pa = ce.add(e);//TODO: ?
      var pb = ce.subtract(e).multiply(Expression.I).divide(RPN('b'));
      pa = Polynomial.toPolynomial(pa, RPN('a'));
      pb = Polynomial.toPolynomial(pb, RPN('a'));

      while (pb.getDegree() !== 0) {
        //console.log(pa + '', pb + '');
        var tmp = pa.divideAndRemainder(pb).remainder;
        pa = pb;
        pb = tmp;
        if (false && pa.getDegree() < 3 && (pa.getDegree() === 1 || pa.getCoefficient(1).equals(Expression.ZERO))) {
        var candidates = Polynomial.toPolynomial(pb.getLeadingCoefficient().getNumerator(), RPN('b')).getZeros(undefined, false).result;
        if (candidates.length > 1) {//?
          /*
          var pb2 = pb.subtract(Polynomial.of(pb.getLeadingCoefficient()).shift(pb.getDegree()));
          var pa2 = pa;
          while (pb2.getDegree() > 0) {
            var tmp = pa2.divideAndRemainder(pb2).remainder;
            pa2 = pb2;
            pb2 = tmp;
          }
          */
          //if (pb2.getDegree() > 0) {
            //var bs = Polynomial.toPolynomial(pb2.calcAt(RPN('0')).getNumerator(), RPN('b')).getZeros(undefined, false).result;
            var pa2 = pa;
            pa2 = pa2.scale(pa2.getContent().inverse());
            var bs = candidates;
            for (var i = 0; i < bs.length; i += 1) {
              var b = bs[i];
              //var e = pa2.getCoefficient(0).negate().divide(pa2.getCoefficient(1));
              //var a = Polynomial.toPolynomial(e.getNumerator(), RPN('b')).calcAt(b).divide(Polynomial.toPolynomial(e.getDenominator(), RPN('b')).calcAt(b));
              var pp = pa2.map(function (coefficient) { return Polynomial.toPolynomial(coefficient, RPN('b')).calcAt(b); });
              //var pp = pa2;
              var roots = pp.getDegree() === 1 ? [pp.getCoefficient(0).negate().divide(pp.getCoefficient(1))] : pp.getroots();
              for (var j = 0; j < roots.length; j += 1) {
                var a = roots[j];
                //a = Polynomial.toPolynomial(a.getNumerator(), RPN('b')).calcAt(b).divide(Polynomial.toPolynomial(a.getDenominator(), RPN('b')).calcAt(b));
                var root = a.add(b.multiply(Expression.I));
                if (p.calcAt(root).equals(Expression.ZERO)) {
                  //console.log(a.toString({fractionDigits: 10}), b.toString({fractionDigits: 10}));
                  //console.log(root + '', a);
                  result.push(root);
                }
              }
            }
          //}
        }
        }
      }
      //console.log(pa + '', pb + '');

      pa = pa.scale(pa.getContent().inverse());
      var bs = Polynomial.toPolynomial(pb.calcAt(RPN('0')).getNumerator(), RPN('b')).getZeros(undefined, false).result;
      for (var i = 0; i < bs.length; i += 1) {
        var b = bs[i];
        var e = pa.getCoefficient(0).negate().divide(pa.getCoefficient(1));
        var a = Polynomial.toPolynomial(e.getNumerator(), RPN('b')).calcAt(b).divide(Polynomial.toPolynomial(e.getDenominator(), RPN('b')).calcAt(b));
        var root = a.add(b.multiply(Expression.I));
        if (p.calcAt(root).equals(Expression.ZERO)) {
          //console.log(a.toString({fractionDigits: 10}), b.toString({fractionDigits: 10}));
          //console.log(root + '', a);
          result.push(root);
        }
      }
    }
    //!

    if (false) {
    result = result.filter(function (x, index) {
      for (var j = index + 1; j < result.length; j += 1) {
        if (result[j].toString() === x.toString()) {
          return false;
        }
      }
      return true;
    });
    }

    //?
    var multiplicities = new Array(result.length).fill(1);
    //_countMultiplicities(multiplicities, intervals, a0);

    return {result: result, multiplicities: multiplicities};
  };
