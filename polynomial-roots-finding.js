import Expression from './Expression.js';
import Polynomial from './Polynomial.js';
import toDecimalStringInternal from './toDecimalString.js';

function PolynomialRoot(polynomial, interval) {
  Expression.Symbol.call(this, "@");
  this.polynomial = polynomial;
  this.interval = interval;
}
PolynomialRoot.prototype = Object.create(Expression.Symbol.prototype);

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
  var f = Polynomial.toPolynomial(en, v);
  if (f.hasRoot(v)) {
    return Expression.ZERO;
  }
  var c = function (x) {
    return Polynomial.toPolynomial(x, v).divideAndRemainder(v.polynomial).remainder.calcAt(v);
  };
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
  var tmp = toDecimalStringInternal(this.e, options.fractionDigits !== -1 && options.fractionDigits != null ? options.fractionDigits : 3, undefined, undefined);
  return tmp;
};

ExpressionWithPolynomialRoot.prototype.toMathML = function (options) {
  //TODO: remove
  var decimalToMathML = function (sign, number) {
    return (sign < 0 ? "<mrow>" : "") + (sign < 0 ? "<mo>&minus;</mo>" : "") + "<mn>" + number + "</mn>" + (sign < 0 ? "</mrow>" : "");
  };
  var complexToMathML = function (real, imaginary, imaginarySign) {
    return real + (imaginarySign >= 0 ? "<mo>+</mo>" : "") + imaginary + "<mo>&#x2062;</mo><mi>i</mi>";
  };

  //TODO:
  if (this.equals(Expression.ZERO)) {
    return Expression.ZERO.toMathML(options);
  }
  //return this.e.toMathML(options);
  var tmp = toDecimalStringInternal(this.e, options.fractionDigits !== -1 && options.fractionDigits != null ? options.fractionDigits : 3, decimalToMathML, complexToMathML);
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
      throw new Error();
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

  var gcd = function (a, b) {
    return b.getDegree() === -1 ? a : gcd(b, a.divideAndRemainder(b).remainder);
  };

  Polynomial.prototype.hasRoot = function (polynomialRoot) {
    var f = this;
    var p = polynomialRoot.polynomial;
    var g = gcd(f, p);
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

  // Polynomial.toPolynomial("x^3-8x^2+21x-18").getZeros();
  Polynomial.prototype.getZeros = function (precision) {
    //TODO: test
    var content = this.getContent();
    var f = this.scale(content.getDenominator()).divideAndRemainder(Polynomial.of(content.getNumerator()), "throw").quotient;

    if (!f.hasIntegerCoefficients()) {
      //return [];
      return {result: [], multiplicities: []};
    }

    // https://en.wikipedia.org/wiki/Square-free_polynomial
    var d = f.derive();
    var a0 = gcd(f, d);
    var p = f.divideAndRemainder(a0).quotient;
    
    //!
    p = p.scale(p.getContent().inverse());
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

    //?
    var multiplicities = new Array(intervals.length);
    for (var i = 0; i < intervals.length; i += 1) {
      multiplicities[i] = 1;
    }
    _countMultiplicities(multiplicities, intervals, a0);

    return {result: result, multiplicities: multiplicities};
  };
