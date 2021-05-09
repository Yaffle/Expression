import Expression from './Expression.js';
import Matrix from './Matrix.js';
import Polynomial from './Polynomial.js';
import toDecimalStringInternal from './toDecimalString.js';
import primeFactor from './primeFactor.js';
import './polynomialFactorization.js';
import nthRoot from './nthRoot.js';
import ExpressionParser from './ExpressionParser.js';
import ExpressionWithPolynomialRoot from './ExpressionWithPolynomialRoot.js';
import LazyPolynomialRoot from './PolynomialRoot.js';


Expression.ExpressionWithPolynomialRoot = ExpressionWithPolynomialRoot;



function ExpressionPolynomialRoot(root) {
  const polynomial = root instanceof LazyPolynomialRoot ? root._root.polynomial : root.polynomial;
  const interval = root instanceof LazyPolynomialRoot ? root._root.interval : root.interval;
  if (root.e != null && Expression.isConstant(root.e)) {
    return root.e;
  }
  if (polynomial.getDegree() === 1 || polynomial.getDegree() === 2 || (polynomial.getDegree() === 4 && false) || polynomial.getDegree() === polynomial.getGCDOfTermDegrees()) {//TODO: other - ? like biqudratic - ?
    var roots = polynomial.getroots();
    for (var rr of roots) {
      if (!Expression.has(rr, ExpressionPolynomialRoot)) {//?
        if (Expression._isPositive(rr.subtract(interval.a)) && Expression._isPositive(rr.subtract(interval.b).negate()) || rr.equals(interval.b)) {
          if (root.e == null || root.e instanceof Expression.Symbol && root.e.symbol === 'α') {//TODO: ???
            return rr;//TODO: MOVE!
          } else {
            var n = Polynomial.toPolynomial(root.e.getNumerator(), new Expression.Symbol('α')).calcAt(rr);
            var d = Polynomial.toPolynomial(root.e.getDenominator(), new Expression.Symbol('α')).calcAt(rr);
            return n.divide(d);
          }
        }
      }
    }
  }
  if (polynomial.getDegree() < 3) {
    throw new TypeError();
  }
  if (polynomial.getDegree() > 64 * 2) {
    throw new Error();//TODO: too long
  }
  Expression.Symbol.call(this, root.toString());
  this.root = root;
}
ExpressionPolynomialRoot.prototype = Object.create(Expression.Symbol.prototype);

ExpressionPolynomialRoot.SimpleInterval = LazyPolynomialRoot.SimpleInterval;

ExpressionPolynomialRoot.prototype.getAlpha = function () {
  return this.root instanceof LazyPolynomialRoot ? this.root._root : this.root;
};
ExpressionPolynomialRoot.prototype.getAlphaExpression = function () {
  return this.root instanceof LazyPolynomialRoot ? this.root.e : new Expression.Symbol('α');
};

ExpressionPolynomialRoot._create = function (polynomial, interval) {
  return new ExpressionPolynomialRoot(new LazyPolynomialRoot.PolynomialRoot(polynomial, new LazyPolynomialRoot.SimpleInterval(interval.a, interval.b)));
};

ExpressionPolynomialRoot.prototype.multiplyInteger = function (x) {
  return x.multiplyPolynomialRoot(this);
};
ExpressionPolynomialRoot.prototype.multiply = function (e) {
  return e.multiplyPolynomialRoot(this);
};
ExpressionPolynomialRoot.prototype.multiplyExpression = function (e) {
  if (e.equals(Expression.ONE)) {
    return this;
  }
  //?
  //TODO: fix
  if (Expression.isConstant(e) && !Expression.has(e, Expression.Complex)) {
    return this.multiply(e);
  }
  return Expression.Symbol.prototype.multiplyExpression.call(this, e);
};
ExpressionPolynomialRoot.prototype.multiplyAddition = function (e) { // for performance (?) when `e` is a constant
  if (Expression.isConstant(e) && !Expression.has(e, Expression.Complex)) {
    return this.multiplyExpression(e);
  }
  return Expression.Symbol.prototype.multiplyAddition.call(this, e);
};
Expression.prototype.multiplyPolynomialRoot = function (root) {
  if (Expression.isConstant(this) && !Expression.has(this, Expression.ExpressionWithPolynomialRoot) && !(Expression.has(this, Expression.Complex)) && !(Expression.has(this, Expression.Exponentiation))) {
    const k = this;
    if (k.equals(Expression.ZERO)) {
      return k;
    }
    if (k.equals(Expression.ONE)) {
      return root;
    }
    return new ExpressionPolynomialRoot(root.root.scale(k));
  }
  if (this instanceof Expression.Complex && !this.imaginary.equals(Expression.ONE)) {
    return this.divide(this.imaginary).multiply(this.imaginary.multiply(root));
  }
  //TODO: ?
  //throw new Error();
  return this.multiplyExpression(root);
};
ExpressionPolynomialRoot.prototype._pow = function (n) {
  return new ExpressionPolynomialRoot(this.root._pow(n));
};
ExpressionPolynomialRoot.prototype.pow = function (e) {
  if (e instanceof Expression.Integer) {
    return this._pow(e.toNumber());
  }
  //TODO: ?
  if (e instanceof Expression.Division && e.getDenominator() instanceof Expression.Integer) {
    //TODO: verify
    return this._nthRoot(e.getDenominator().toNumber()).pow(e.getNumerator());
  }
  return Expression.Symbol.prototype.pow.call(this, e);
};
ExpressionPolynomialRoot.prototype.multiplyPolynomialRoot = function (x) {
  var y = this;
  return new ExpressionPolynomialRoot(x.root.multiply(y.root));
};
ExpressionPolynomialRoot.prototype.add = function (e) {
  return e.addPolynomialRoot(this);
};
ExpressionPolynomialRoot.prototype.addPolynomialRoot = function (x) {
  var y = this;
  return new ExpressionPolynomialRoot(x.root.add(y.root));
};
Expression.prototype.addPolynomialRoot = function (root) {
  if (Expression.isConstant(this) && !Expression.has(this, Expression.ExpressionWithPolynomialRoot) && !(Expression.has(this, Expression.Complex)) && !(Expression.has(this, Expression.Exponentiation))) {
    var k = this;
    if (k.equals(Expression.ZERO)) { // for performance
      return root;
    }
    return new ExpressionPolynomialRoot(root.root.translate(k));
  }
  //throw new Error();
  return this.addExpression(root);
};
ExpressionPolynomialRoot.prototype.addExpression = function (e) {
  return this.add(e);//!?
};
ExpressionPolynomialRoot.prototype.divide = function (e) {
  //if (e.equals(Expression.ONE)) {
  //  return this;
  //}
  //if (!(e instanceof ExpressionPolynomialRoot) && !Expression.isConstant(e) || Expression.has(e, Expression.Matrix) || Expression.has(e, Expression.MatrixSymbol)) {
    //TODO: why - ?
  //  throw new Error();
  //}
  return this.multiply(e.inverse());
};
ExpressionPolynomialRoot.prototype.divideExpression = function (x) {
  return x.multiply(this.inverse());
};
ExpressionPolynomialRoot.prototype.inverse = function () {
  return new ExpressionPolynomialRoot(this.root.inverse());
};
ExpressionPolynomialRoot.prototype.sign = function () {
  return this.root.sign();
};
var toRadicalExpression = function (polynomialRoot) {
  if (polynomialRoot.polynomial.getDegree() === 1) {
    //TODO: ???
    return polynomialRoot.polynomial.getroots()[0];
  }
  //TODO: ?
  var g = polynomialRoot.polynomial.getGCDOfTermDegrees();
  if (g > 1) {
    var v = toRadicalExpression(polynomialRoot._pow(g));
    if (v != null) {
      if (g % 2 === 1 || polynomialRoot.sign() > 0) {
        return Expression.NthRoot.makeRoot(v, g);
      } else {
        return new Expression.Negation(Expression.NthRoot.makeRoot(v, g));
      }
    }
  }
  // convert to depressed:
  var h = polynomialRoot.polynomial._getShiftToDepressed();
  if (!h.equals(Expression.ZERO)) {
    var tmp = toRadicalExpression(polynomialRoot.translate(h));
    if (tmp != null) {
      return new Expression.Addition(tmp, h.negate());
    }
  }
  return null;
};
ExpressionPolynomialRoot.prototype.toString = function (options) {
  //return new ExpressionWithPolynomialRoot(this, this).toString(options);
  options = options || {};
  if (options.fractionDigits !== -1 && options.fractionDigits != null) {
    console.debug('options.fractionDigits is deprecated, please use options.rounding');
  }
  var rounding = options.rounding != null ? options.rounding : (options.fractionDigits !== -1 && options.fractionDigits != null ? {fractionDigits: options.fractionDigits} : {fractionDigits: 3});
  if (options.rounding == null) {
    var p = (this.root instanceof LazyPolynomialRoot ? this.root._root : this.root).polynomial;
    if (p.getDegree() / p.getGCDOfTermDegrees() < 10) { //TODO: REMOVE !!!
      var re = toRadicalExpression(this.root.toPolynomialRoot());
      if (re != null) {
        return re.toString(options);
      }
    }
  }
  return toDecimalStringInternal(this, rounding);
};
ExpressionPolynomialRoot.prototype.equals = function (other) {
  if (other instanceof ExpressionPolynomialRoot) {
    return this.root.equals(other.root);
  }
  // optimization
  if (other instanceof Expression.Integer) {
    return false;
  }
  /*if (Expression.isConstant(other)) {
    if (!this.root.polynomial.calcAt(other).equals(Expression.ZERO)) {
      return false;
    }
    var withinInterval = function (x, interval) {
      return Expression._isPositive(x.subtract(interval.a)) && Expression._isPositive(x.subtract(interval.b).negate());
    };
    return withinInterval(other, this.interval);
  }*/
  //TODO: optimize
  return this.subtract(other).equals(Expression.ZERO);
};
ExpressionPolynomialRoot.prototype.compare4MultiplicationComplex = function (x) {
  return -1;
  //return +1;
};
ExpressionPolynomialRoot.prototype.compare4MultiplicationNthRoot = function (x) {
  return 0;
};
ExpressionPolynomialRoot.prototype.compare4Multiplication = function (y) {
  if (y instanceof Expression.Complex) {
    return +1;
    //return -1;
  }
  if (y instanceof Expression.Integer) {
    return +1;
  }
  if (y instanceof ExpressionPolynomialRoot) {
    return 0;
  }
  if (y instanceof Expression.NthRoot) {
    return 0;//?
  }
  if (y instanceof Expression.Symbol) {
    return -1;
  }
  return Expression.Symbol.prototype.compare4Multiplication.call(this, y);
};
ExpressionPolynomialRoot.prototype.compare4MultiplicationSymbol = function (x) {
  return +1;
};
ExpressionPolynomialRoot.prototype.compare4Addition = function (y) {
  if (y instanceof ExpressionPolynomialRoot) {
    return 0;//?
  }
  if (y instanceof Expression.Symbol) {
    return +1;
  }
  if (y instanceof Expression.NthRoot) {
    return 0;//?
  }
  return Expression.Symbol.prototype.compare4Addition.call(this, y);
};
ExpressionPolynomialRoot.prototype.compare4AdditionSymbol = function (x) {
  return -1;
};

ExpressionPolynomialRoot.prototype._nthRoot = function (n) {//?
  if (this.root.sign() < 0) {
    //TODO: check
    return Expression.I.multiply(this.negate()._nthRoot(n));
  }
  return new ExpressionPolynomialRoot(this.root._nthRoot(n));
};

ExpressionPolynomialRoot.prototype.upgrade = function () {
  return new ExpressionPolynomialRoot(this.root.upgrade());//TODO: ?
};

ExpressionPolynomialRoot.prototype.isNegative = function () {
  return this.root.sign() < 0;
};

Expression.prototype.upgrade = function () { //TODO: remove !!!
  return this;
};

ExpressionPolynomialRoot.prototype.isExact = function () {
  //TODO: fix - ?
  return false;
};

ExpressionPolynomialRoot.prototype.negate = function () {
  return new ExpressionPolynomialRoot(this.root.negate()); // for performance
};

ExpressionPolynomialRoot.prototype.simplifyExpression = function () {//TODO: remove - ?
  return this;
};

ExpressionPolynomialRoot.prototype.toMathML = function (options) {
  options = options || {};
  if (options.fractionDigits !== -1 && options.fractionDigits != null) {
    console.debug('options.fractionDigits is deprecated, please use options.rounding');
  }
  var rounding = options.rounding != null ? options.rounding : (options.fractionDigits !== -1 && options.fractionDigits != null ? {fractionDigits: options.fractionDigits} : {fractionDigits: 3});
  if (options.rounding == null) {
    var p = (this.root instanceof LazyPolynomialRoot ? this.root._root : this.root).polynomial;
    if (p.getDegree() / p.getGCDOfTermDegrees() < 10) { //TODO: REMOVE !!!
      var re = toRadicalExpression(this.root.toPolynomialRoot());
      if (re != null) {
        return re.toMathML(options);
      }
    }
  }
  var tmp = toDecimalStringInternal(this, rounding, Expression._decimalToMathML, Expression._complexToMathML);
  return tmp;
};

//TODO: ?????
ExpressionPolynomialRoot.prototype.getPrecedence = function () {
  return 1000;
};
ExpressionPolynomialRoot.prototype.isRightToLeftAssociative = function () {
  return true;
};
ExpressionPolynomialRoot.prototype.isUnaryPlusMinus = function () {
  return true;//TODO: !?
};

//ExpressionPolynomialRoot.prototype.complexConjugate = function () {//TODO: test
//  return this;
//};

Expression.ExpressionPolynomialRoot = ExpressionPolynomialRoot;

Expression.toPolynomialRoot = function (e) {
  var x = e instanceof Expression.NthRoot ? e.a : e;//TODO: remove
  var n = e instanceof Expression.NthRoot ? e.n : 1;//TODO: remove
  var symbol = new Expression.Symbol('x');
  if (!(x.getDenominator() instanceof Expression.Integer)) {
    throw new TypeError();
  }
  var p = Polynomial.toPolynomial(Expression.getConjugateExpression(symbol._pow(n).subtract(x).getNumerator()), symbol);
  //TODO: remove:
  if (p.getDegree() <= 8 && (true || isSmall(p))) {//TODO: ?
    var factor = p.factorize();
    if (factor != null && factor.getDegree() < p.getDegree() && factor.getDegree() === 4) {//?
      var roots = Polynomial.polynomialGCD(factor, Polynomial.toPolynomial(symbol._pow(n).subtract(x), symbol)).getroots();
      for (const root of roots) {
        if (root._pow(n).equals(x)) {
          //debugger;
          return Expression._isPositive(root) || n % 2 !== 0 ? root : root.negate();
        }
      }
    }
  }
  if (true) {
    //var root = Expression.toPolynomialRoot(x)._nthRoot(n);
    
    //TODO: Expression#toPolynomialRoot() - ?
    //TODO: move up (!?)
    p = p.squareFreeFactors().a1;//TODO: which one (?)
    const isComplex = n === 2 && Expression.has(e.radicand, Expression.Complex);
    var zeros = p.getZeros(undefined, isComplex);
    if (zeros.length === 2) {
      //TODO: remove
      if (Expression._isPositive(zeros[1]) && Expression._isPositive(x)) {
        return zeros[1];
      }
    }
    //TODO: find zero only on interval
    for (var zero of zeros) {
      if (zero.root != null && Expression._isPositive(zero) || isComplex && Expression._isPositive(Expression.getComplexNumberParts(zero).real)) {
        if (zero._pow(n).equals(x)) {
          return zero;
        }
      }
    }
    //TODO: ?
  }
  console.error(e.toString());
  return undefined;
};



globalThis.SturmSequence = SturmSequence;//TODO: ???

  // https://math.stackexchange.com/questions/309178/polynomial-root-finding
  function SturmSequence(f) {
    //f = f.scale(f.getContent().inverse().abs());//?
    var d = f.derive();
    d = d.scale(d.getContent().inverse().abs());//?
    let s = [];
    var fp = f;
    var fc = d;
    s.push(fp);
    s.push(fc);
    while (fc.getDegree() > 0) {
      //var fn = fp.divideAndRemainder(fc).remainder.negate();
      //TODO: ?
      // https://en.wikipedia.org/wiki/Sturm%27s_theorem#Use_of_pseudo-remainder_sequences
      var fn = Polynomial.pseudoRemainder(fp, fc.getLeadingCoefficient().abs().equals(fc.getLeadingCoefficient()) ? fc : fc.negate()).negate();
      if (fn.getDegree() >= 0) {
        var y = fn.getContent().inverse();
        if (y.isNegative()) {
          y = y.negate();
        }
        fn = fn.scale(y);//!
        //s.push(fn.scale(y));
        //?
        s.push(fn);
      }
      fp = fc;
      fc = fn;
    }
    this.s = s;
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
        v = (e.equals(Expression.ONE) ? p : p.map(function (coefficient, degree) {
          return coefficient.multiply(Expression.pow(e, n - degree));
        })).calcAt(x.getNumerator());
      }
      //!
      var c = v.compareTo(Expression.ZERO);
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

  // interval - the half-open interval (a, b] (see Wikipedia's article)
  SturmSequence.prototype.numberOfRoots = function (interval) {
    if (interval.a.equals(interval.b)) {
      throw new TypeError();
    }
    return this.signChanges(interval.a) - this.signChanges(interval.b);
  };

  Polynomial.prototype.getRootIntervals = function () {
    var sturmSequence = new SturmSequence(this);
    var interval = {a: this.subs(x => x.negate()).getPositiveRealRootsBound().negate(), b: this.getPositiveRealRootsBound()};
    var getIntervals = function (interval, rootsAfterA, rootsAfterB) {
      //var n = sturmSequence.numberOfRoots(interval);
      var n = rootsAfterA - rootsAfterB;
      if (n === 1) {
        return [interval];
      }
      if (n > 1) {
        var middle = interval.a.add(interval.b).divide(Expression.TWO);
        var rootsAfterM = sturmSequence.signChanges(middle);
        var a = getIntervals({a: interval.a, b: middle}, rootsAfterA, rootsAfterM);
        var b = getIntervals({a: middle, b: interval.b}, rootsAfterM, rootsAfterB);
        return a.concat(b);
      }
      return [];
    };
    var negative = getIntervals({a: interval.a, b: Expression.ZERO}, sturmSequence.signChanges(interval.a), sturmSequence.signChanges(Expression.ZERO));
    var positive = getIntervals({a: Expression.ZERO, b: interval.b}, sturmSequence.signChanges(Expression.ZERO), sturmSequence.signChanges(interval.b));
    return negative.concat(positive);
  };

  Polynomial.prototype.getPositiveRealRootsBound = function () {
    //TODO: only integer coefficients (?)
    // https://en.wikipedia.org/wiki/Sturm%27s_theorem#Number_of_real_roots
    // https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#Bounds_of_positive_real_roots
    let M = null;
    //TODO: fix the iteration
    const n = this.getDegree();
    const an = this.getLeadingCoefficient();
    for (let i = 0; i <= this.getDegree() - 1; i += 1) {
      const v = this.getCoefficient(i).negate().truncatingDivide(an);
      if (v.sign() >= 0) {
        const c = Expression.TWO.multiply(Expression.Integer.fromBigInt(nthRoot(v.toBigInt(), n - i)).add(Expression.ONE));
        if (M == null || M.compareTo(c) < 0) {
          M = c;
        }
      }
    }
    if (M == null) {
      return Expression.ZERO;
    }
    //!2020-12-19
    // round to a power of two:
    M = Expression.TWO._pow(primeFactor._bitLength(M.getNumerator().toBigInt()) - primeFactor._bitLength(M.getDenominator().toBigInt()) + 1);
    //!
    return M;
  };

  //TODO: BigDecimal - ?, rounding - ?
  Polynomial.prototype.getZero = function (interval, precision) {
    var roundFloor = function (point, e) {
      var n = point.getNumerator().multiply(e);
      var d = point.getDenominator();
      console.assert(d.compareTo(Expression.ZERO) > 0);
      var q = n.truncatingDivide(d);
      var r = n.subtract(q.multiply(d));
      return r.compareTo(Expression.ZERO) >= 0 ? q : q.subtract(Expression.ONE);
    };
    const sign = function (v) {
      return Math.sign(v.getNumerator().compareTo(Expression.ZERO));
    };
    //const BASE = Expression.TEN;
    const BASE = Expression.TWO;
    const e = Expression.pow(BASE, precision); // epsilon^-1
    if (!(e instanceof Expression.Integer)) {
      throw new RangeError("epsilon^-1 is not an integer");
    }
    var a = interval.a;
    var b = interval.b;
    // (b - a) * 10**precision > min(abs(a), abs(b))
    // (b - a) * 10**fractionDigits > 1
    //TODO: fix to use precision, not fractionDigits
    // e * (b - a) > 0:
    if (e.multiply(a.getDenominator().multiply(b.getNumerator()).subtract(b.getDenominator().multiply(a.getNumerator()))).compareTo(b.getDenominator().multiply(a.getDenominator())) > 0) {
      //TODO:
      var tmp = true && precision >= 16 / Math.log10(BASE.toNumber()) ? this.getZero(interval, Math.floor(precision / 4)) : interval;
      a = tmp.a;
      b = tmp.b;

      var n = this.getDegree();
      var p = this.map(function (coefficient, degree) {
        return coefficient.multiply(Expression.pow(e, n - degree));
      });
      //p = p.primitivePart();//?
      const sa = roundFloor(a, e).add(Expression.ONE); // a.getNumerator().multiply(e).truncatingDivide(a.getDenominator()).add(Expression.ONE);//?
      const sb = roundFloor(b, e); // b.getNumerator().multiply(e).truncatingDivide(b.getDenominator());//?
      console.assert(sa.multiply(a.getDenominator()).subtract(a.getNumerator().multiply(e)).compareTo(Expression.ZERO) >= 0); // sa/e >= a
      console.assert(sb.multiply(b.getDenominator()).subtract(b.getNumerator().multiply(e)).compareTo(Expression.ZERO) <= 0); // sb/e <= b
      //TODO: bigdecimal - ?
      // remember values at boundaries to reuse in the loop:
      let pa = p.calcAt(sa);
      let pb = p.calcAt(sb);
      const spb = sign(pb);
      const spa = sign(pa);
      if (spa === 0) {
        return {a: sa.divide(e), b: sa.divide(e)};
      }
      if (spb === 0) {
        return {a: sb.divide(e), b: sb.divide(e)};
      }
      if (spa === spb) {
        if (spa !== (sign(this.calcAt(a)) || sign(this.calcAt(b).negate()) || spa)) {
          return {a: a, b: sa.divide(e)};
        }
        if (spb !== sign(this.calcAt(b))) {
          return {a: sb.divide(e), b: b};
        }
        throw new RangeError();//?
      }
      a = sa;
      b = sb;
      // bisection method
      var cc = 0;
      var d = p.derive();
      var width = b.subtract(a);
      while (width.compareTo(Expression.ONE) > 0) {// b - a > 1
        var middle = a.add(width.truncatingDivide(Expression.TWO));
        //console.log(eval(a.divide(e).toString()) + ' - ' + eval(b.divide(e).toString()));
        //?
        if (cc % 3 !== 2 && width.compareTo(a.abs()) < 0) {// TODO: test for the case when a < 0
          // TODO: better guesses
          // Newton's method
          var x = cc % 3 === 1 ? a : b;
          var px = x === a ? pa : (x === b ? pb : undefined);
          var c = d.calcAt(x);
          if (!c.equals(Expression.ZERO)) {
            x = x.subtract(px.truncatingDivide(c));
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
        var v = p.calcAt(middle);
        var sv = sign(v);
        if (sv === spb) {
          b = middle;
          pb = v;
        } else if (sv === spa) {
          a = middle;
          pa = v;
        } else {
          a = middle;
          b = middle;
          pa = v;
          pb = v;
        }
        width = b.subtract(a);
      }
      //console.debug(cc);
      a = a.divide(e);
      b = b.divide(e);
    }
    return {a: a, b: b};
  };

  Polynomial.prototype.hasRoot = function (polynomialRoot) {
    var f = this;
    if (f.equals(Polynomial.ZERO)) {
      return true;
    }
    //!new 2021-02-20 (TODO: CHECK)
    if (!f.hasIntegerCoefficients() && f.hasComplexCoefficients()) {
      return f.map(c => c instanceof Expression.Integer ? c : c.real).hasRoot(polynomialRoot) && f.map(c => c instanceof Expression.Integer ? Expression.ZERO : c.imaginary).hasRoot(polynomialRoot);
    }
    //!
    var p = polynomialRoot.polynomial;
    var g = null;
    //!
    if (!f.hasIntegerCoefficients()) {
      var variable = new Expression.Symbol('~');
      var ff = f.calcAt(variable);
      var tmp = Expression.getMultivariatePolynomial(ff);
      if (tmp != null && !tmp.v.equals(variable) && tmp.v instanceof Expression.Symbol) {
        g = Polynomial.polynomialGCD(Polynomial.toPolynomial(tmp.p.getContent(), variable), p);
      }
    }
    if (g == null) {
      g = Polynomial.polynomialGCD(f, p);
    }
    //!
    if (g.getDegree() < 1) {
      return false;
    }
    var i = polynomialRoot.interval;

    if (!g.hasIntegerCoefficients()) {
      //TODO: BUG?
      //?new
      var variable = new Expression.Symbol('$$');
      var e = g.calcAt(variable);
      var c = Expression.getComplexConjugate(e);
      if (c != null) {
        g = Polynomial.toPolynomial(c.multiply(e), variable).getSquareFreePolynomial();
      }
      //?
    }

    return g.numberOfRoots(i) === 1;
  };

  // interval - the half-open interval (a, b] (see Wikipedia's article)
  Polynomial.prototype.numberOfRoots = function (interval = null) {
    if (interval == null) {
      interval = {a: this.subs(x => x.negate()).getPositiveRealRootsBound().negate(), b: this.getPositiveRealRootsBound()};//TODO: use (-1/0; +1/0)
    }
    var sturmSequence = new SturmSequence(this);
    return sturmSequence.numberOfRoots(interval);
  };


  // Polynomial.toPolynomial(ExpressionParser.parse("x^3-8x^2+21x-18"), ExpressionParser.parse("x")).getZeros(3).toString()
  Polynomial.prototype.getZeros = function (precision = 0, complex = false) {
    if (this.getCoefficient(0).equals(Expression.ZERO)) {
      if (this.getLeadingCoefficient().equals(Expression.ZERO)) {
        throw new TypeError();
      }
      var i = 0;
      while (this.getCoefficient(i).equals(Expression.ZERO)) {
        i += 1;
      }
      var tmp = this.divideAndRemainder(Polynomial.of(Expression.ONE).shift(i)).quotient.getZeros(precision, complex);
      return tmp.concat(new Array(i).fill(Expression.ZERO));
    }
    //TODO: test
    var content = this.getContent();
    var f = this.scale(content.getDenominator()).divideAndRemainder(Polynomial.of(content.getNumerator()), "throw").quotient;

    // https://en.wikipedia.org/wiki/Square-free_polynomial
    var tmp = f.squareFreeFactors();
    var a0 = tmp.a0;
    var a1 = tmp.a1;

    if (a0.getDegree() !== 0) {
      var tmp1 = a1.getZeros(precision, complex); // roots with multiplicity = 1 (?)
      var tmp2 = a0.getZeros(precision, complex);
      var result = [];
      var previous = undefined;
      for (var i = 0; i < tmp2.length; i += 1) {
        var zero = tmp2[i];
        if (zero !== previous) {
          result.push(zero);
          previous = zero;
        }
        result.push(zero);
      }
      return tmp1.concat(result);
    }

    var p = f;
    if (p.getDegree() === 0) {
      return [];
    }

    //!
    p = p.scale(p.getContent().inverse());
    //!

    if (!f.hasIntegerCoefficients()) {
      //?new
      var variable = new Expression.Symbol('$$')
      var e = f.calcAt(variable);
      var c = Expression.getConjugateExpression(e);
      if (c != null && !e.equals(c)) {
        var result = [];
        var tmp = Polynomial.toPolynomial(c, variable).getZeros(precision, complex);
        console.time('checking roots');
        for (var i = 0; i < tmp.length; i += 1) {
          var zero = tmp[i];
          if (zero instanceof ExpressionPolynomialRoot && zero.root.e.equals(new Expression.Symbol('α')) ? f.hasRoot(zero.root._root) :
              zero instanceof ExpressionWithPolynomialRoot && zero.e === zero.root ? f.hasRoot(zero.root) :
              f.calcAt(zero).equals(Expression.ZERO)) {
            result.push(zero);
          } else {
            //TODO:?
            console.debug(zero.root);
          }
        }
        console.timeEnd('checking roots');
        return result;
      }
      //!new
      // u * x + v = t
      // u**n*x**n = a_n*x**n, u = a_n**(1/n)
      // u**(n-1)*x**(n-1)*v*n+u**(n-1)*x**(n-1) = a_(n-1)*x**(n-1)
      var u = Polynomial.of(Expression.ONE).shift(f.getDegree()).subtract(Polynomial.of(f.getLeadingCoefficient())).getroots();
      if (u.length !== 0) {
        u = u[0];
        var v = f.getCoefficient(f.getDegree() - 1).divide(u._pow(f.getDegree() - 1)).subtract(Expression.ZERO).divide(Expression.Integer.fromNumber(f.getDegree()));
        var x = new Expression.Symbol('$$');//?
        var pt = Polynomial.toPolynomial(p.calcAt(x.subtract(v).divide(u)).getNumerator(), x);
        if (pt.hasIntegerCoefficients()) {//TODO: ?
          return pt.getZeros(precision, complex).map(function (zero) {
            return zero.subtract(v).divide(u);
          });
        }
      }
      //!
      //?
      return [];
    }

    //!new
    if (p.getDegree() === 3) {
      //?
    }
    //!

    // https://en.wikipedia.org/wiki/Sturm%27s_theorem
    var intervals = p.getRootIntervals();

    // https://math.stackexchange.com/questions/309178/polynomial-root-finding
    // "it is guaranteed that there is a sign change inside every interval (because there are no repeated zeroes)"
    var result = new Array(intervals.length);
    var enableNewClass = false;
    for (var i = 0; i < intervals.length; i += 1) {
      var zero = p.getZero(intervals[i], precision);
      if (zero.a.equals(zero.b)) {
        result[i] = zero.a;//TODO: fix
      } else {
        //! p, not f, as f may have roots with multiplicity > 1
        if (!enableNewClass) {
          var root = new Expression.PolynomialRootSymbol(p, zero);
          result[i] = new ExpressionWithPolynomialRoot(root, root);
        } else {
          var root = new ExpressionPolynomialRoot(LazyPolynomialRoot.create(p, new LazyPolynomialRoot.SimpleInterval(zero.a, zero.b), {skipFactorization: true}));
          //result[i] = new ExpressionWithPolynomialRoot(new Expression.Symbol('$α'), root);
          result[i] = root;
        }
      }
    }
    //return result;

    //!new
    //var p = np;
    if (intervals.length !== p.getDegree() && true && complex) {
      //!new
      if (p.getDegree() > 4) {//?
        var factor = p.factorize();
        if (factor != null) {
          //TODO: remove double work
          return factor.getZeros(precision, complex).concat(p.divideAndRemainder(factor, "throw").quotient.getZeros(precision, complex));
        }
      }
      //!
      if (p.isEven()) {
        //debugger;
        const zeros = p.subs(x => x.squareRoot()).getZeros(precision, complex);
        for (var zero of zeros) {
          //var z = zero.squareRoot();
          // https://en.wikipedia.org/wiki/Complex_number#Square_root
          var squareRoot = function (z) {
            var tmp = Expression.getComplexNumberParts(z);
            var a = tmp.real;
            var b = tmp.imaginary;
            var aapbb = a._pow(2).add(b._pow(2)).squareRoot();
            var γ = a.add(aapbb).divide(Expression.TWO).squareRoot();
            var sign = (b.compareTo(Expression.ZERO) > 0 ? Expression.ONE : Expression.ONE.negate());
            var tmp = a.negate().add(aapbb).divide(Expression.TWO);
            //debugger;
            var δ = sign.multiply(tmp.squareRoot());
            return γ.add(δ.multiply(Expression.I));
          };
          //zero = zero instanceof ExpressionWithPolynomialRoot ? zero.upgrade() : zero;
          if (!Expression._isPositive(zero)) {
            var z = squareRoot(zero.e != null ? zero.upgrade() : zero);
            result.push(z);
            result.push(z.negate());
          }
        }
        return result;
      }
      //var p = stringToPolynomial("x^5+2*x^2+2*x+3");

      var e = p.calcAt(ExpressionParser.parse("a+b*i"));
      var ce = Expression.getComplexConjugate(e);
      var pa = ce.add(e);//TODO: ?
      var pb = ce.subtract(e).multiply(Expression.I).divide(ExpressionParser.parse('b'));
      const cpa = pa;
      const cpb = pb;
    if (true) {
      pa = Polynomial.toPolynomial(pa, ExpressionParser.parse('a'));
      pb = Polynomial.toPolynomial(pb, ExpressionParser.parse('a'));
      while (pa.getCoefficient(0).equals(Expression.ZERO)) {
        pa = pa.divideAndRemainder(Polynomial.of(Expression.ONE).shift(1)).quotient;//TODO: simplify
      }
      while (pb.getCoefficient(0).equals(Expression.ZERO)) {
        // a = 0, p(b*i) = 0
        var candidates = Polynomial.toPolynomial(pa.calcAt(Expression.ZERO), ExpressionParser.parse('b')).getZeros(undefined, false);
        for (var c of candidates) {
          var root = c.multiply(Expression.I);
          if (p.calcAt(root).equals(Expression.ZERO)) {
            result.push(root);
          }
        }
        pb = pb.divideAndRemainder(Polynomial.of(Expression.ONE).shift(1)).quotient;//TODO: simplify
      }
      //TODO: verify that no roots are lost
      console.assert(pa.getContent() instanceof Expression.Integer);
      pa = pa.divideAndRemainder(Polynomial.of(pa.getContent()), "throw").quotient;
      //TODO: verify that no roots are lost
      console.assert(pb.getContent() instanceof Expression.Integer);
      pb = pb.divideAndRemainder(Polynomial.of(pb.getContent()), "throw").quotient;

      var walk = function (p1, p2, condition) {
        if (p2.getDegree() > 0) {
          var lc = p2.getLeadingCoefficient();
          var c1 = condition.andZero(lc);
          var c2 = condition.andNotZero(lc);
          var simplifyCoefficients = function (p, condition) {
            if (condition.array.length === 1 && condition.array[0].operator === " == 0") {
              //console.assert(condition.array.length === 1 && condition.array[0].operator === " == 0");
              var zero = Polynomial.toPolynomial(condition.array[0].expression, ExpressionParser.parse('b'));
              return p.map(function (coefficient) {
                var n = Polynomial.toPolynomial(coefficient.getNumerator(), ExpressionParser.parse('b')).divideAndRemainder(zero).remainder.calcAt(ExpressionParser.parse('b'));
                var d = Polynomial.toPolynomial(coefficient.getDenominator(), ExpressionParser.parse('b')).divideAndRemainder(zero).remainder.calcAt(ExpressionParser.parse('b'));
                return n.divide(d);
              });
            }
            return p;
          };
          if (!c1.isFalse()) {
            walk(simplifyCoefficients(p1, c1), simplifyCoefficients(p2, c1), c1);
          }
          if (!c2.isFalse()) {
            console.assert(result.length < p.getDegree());
            var newp2 = p2.scale(p2.getContent().inverse());
            p1 = p1.scale(p1.map(c => c.getDenominator().inverse()).getContent().inverse());
            var r = Polynomial.pseudoRemainder(p1, newp2);
            walk(newp2, simplifyCoefficients(r, c2), c2);
          }
        } else {
          //TODO: ?
          condition = condition.andZero(p2.getLeadingCoefficient());
          if (!condition.isFalse()) {
            console.assert(condition.array.length === 1 && condition.array[0].operator === " == 0");
            const bPolynomial = Polynomial.toPolynomial(condition.array[0].expression, ExpressionParser.parse('b'));
            const getZeros1 = function (p) {//TODO: !? use everywhere (?)
              var factor = p.factorize();
              if (factor != null) {
                //TODO: remove double work
                return factor.getZeros(undefined, false).concat(p.divideAndRemainder(factor, "throw").quotient.getZeros(undefined, false));
              }
              return p.getZeros(undefined, false);
            };
            //TODO: fix for higher degrees (?)
            let candidates = bPolynomial.getDegree() < 3 ? bPolynomial.getroots() : getZeros1(bPolynomial);
            candidates = candidates.filter(c => Expression._isPositive(c));//!?
            for (const b of candidates) {
              const pp = p1.map(function (coefficient) { return Polynomial.toPolynomial(coefficient.getNumerator(), ExpressionParser.parse('b')).calcAt(b).divide(Polynomial.toPolynomial(coefficient.getDenominator(), ExpressionParser.parse('b')).calcAt(b)); });
              if (pp.getDegree() === 1 || pp.getCoefficient(1).equals(Expression.ZERO) && pp.getDegree() < 3) {
                const roots = pp.getroots();
                for (const a of roots) {
                  if (!Expression.has(a, Expression.Complex)) {
                    result.push(a.add(b.multiply(Expression.I)));
                    result.push(a.add(b.negate().multiply(Expression.I)));
                  }
                }
              }
            }
          }
        }
      };
      walk(pa, pb, Condition.TRUE);
      //console.log(pa + '', pb + '');
    }
      //!TODO: 
      //!new 2021-01-03
      if (result.length < p.getDegree()) {
        console.count('yyy');
        //debugger;
        var resultant = function (v1, v2) {
          var A = Polynomial.toPolynomial(cpa, ExpressionParser.parse(v1));
          var B = Polynomial.toPolynomial(cpb, ExpressionParser.parse(v1));
          return Polynomial.resultant(A, B, v2).primitivePart();
        };
        var bCandidates = resultant('a', 'b').getZeros();
        var aCandidates = resultant('b', 'a').getZeros();
        var unique = function (array) {
          //return Array.from(new Set(array));
          var result = [];
          for (const element of array) {
            if (result.indexOf(element) === -1) {
              result.push(element);
            }
          }
          return result;
        };
        bCandidates = unique(bCandidates);
        aCandidates = unique(aCandidates);
        bCandidates = bCandidates.filter(c => Expression._isPositive(c));//!?
        //console.log(bCandidates.map(x =>  typeof x.upgrade === 'function' ?  x.upgrade() : x).toString());
        //console.log(aCandidates.map(x =>  typeof x.upgrade === 'function' ?  x.upgrade() : x).toString());
        // https://en.wikipedia.org/wiki/Resultant#Application_to_polynomial_systems
        //debugger;
        for (var a of aCandidates) {
          for (var b of bCandidates) {
            var candidate = a.add(b.multiply(Expression.I));
            if (p.calcAt(candidate).equals(Expression.ZERO)) {
              result.push(candidate);
              result.push(a.add(b.negate().multiply(Expression.I)));
            }
          }
        }
      }
    }
    //!

    if (intervals.length !== p.getDegree() && true && complex) {
      //TODO: FIX!!!
      var strings = result.map(x => x.toMathML({rounding: {fractionDigits: 3}}));
      result = result.filter(function (x, index) {
        for (var j = index - 1; j >= 0; j -= 1) {
          if (strings[j] === strings[index]) {
            if (result[j].equals(result[index])) {
              return false;
            }
          }
        }
        return true;
      });
    }

    return result;
  };

Polynomial.getSylvesterMatrix = function (p, q) {
  var m = p.getDegree();
  var n = q.getDegree();
  return Matrix.Zero(n + m, n + m).map(function (element, i, j) {
    var index1 = m - (j - i);
    var index2 = n - (j - (i - n));
    return i < n ? (index1 < 0 || index1 > m ? Expression.ZERO : p.getCoefficient(index1)) : (index2 < 0 || index2 > n ? Expression.ZERO : q.getCoefficient(index2));
  });
};

//Polynomial.resultant = function (p, q, v2) {
//  //return Polynomial.toPolynomial(Polynomial.getSylvesterMatrix(p, q).determinant(), ExpressionParser.parse(v2));
//  return Polynomial.getSylvesterMatrix(p, q).map(e => new Expression.Polynomial(Polynomial.toPolynomial(e, ExpressionParser.parse(v2)))).determinant().polynomial;
//};

function subresultantPseudoRemainderSequence(A, B) {
  let first = true;
  let phi = Expression.ONE;
  var iterator = {
    next: function () {
      console.assert(A.getDegree() >= B.getDegree());
      // For the explanation and proof see Donald E. Knuth The Art of computer programming Third Edition, Volume 2 (Seminumerical algorithms), page 428.
      // https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Subresultant_pseudo-remainder_sequence
      if (!B.equals(Polynomial.ZERO)) {
        const d = A.getDegree() - B.getDegree();
        const scale = B.getLeadingCoefficient()._pow(d + 1);
        const α = first ? Expression.ONE : A.getLeadingCoefficient().multiply(phi._pow(d));
        const tmp = A.scale(scale).divideAndRemainder(B, "throw");
        const q = tmp.quotient;
        const R = tmp.remainder.divideAndRemainder(Polynomial.of(α), "throw").quotient;
        first = false;
        phi = d === 0 ? phi : phi.inverse()._pow(d - 1).multiply(B.getLeadingCoefficient()._pow(d));
        const value = {A: A, B: B, q: q, R: R, α: α};
        A = B;
        B = R;
        return {value: value, done: false};
      }
      return {value: undefined, done: true};
    }
  };
  iterator[globalThis.Symbol.iterator] = function () {
    return this;
  };
  return iterator;
}

Polynomial._subresultantPseudoRemainderSequence = subresultantPseudoRemainderSequence;

function primitivePseudoRemainderSequence(A, B) {
  var iterator = {
    next: function () {
      console.assert(A.getDegree() >= B.getDegree());
      // https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Primitive_pseudo-remainder_sequence
      if (!B.equals(Polynomial.ZERO)) {
        const tmp = A.scale(B.getLeadingCoefficient()._pow(A.getDegree() - B.getDegree() + 1)).divideAndRemainder(B, "throw");
        const q = tmp.quotient;
        const r = tmp.remainder;
        const α = r.getContent();
        const R = r.divideAndRemainder(Polynomial.of(α), "throw").quotient;
        const value = {A: A, B: B, q: q, R: R, α: α};
        A = B;
        B = R;
        return {value: value, done: false};
      }
      return {value: undefined, done: true};
    }
  };
  iterator[globalThis.Symbol.iterator] = function () {
    return this;
  };
  return iterator;
}

Polynomial._primitivePseudoRemainderSequence = primitivePseudoRemainderSequence;

Polynomial.resultant = function (A, B, v2) {
  A = A.map(c => new Expression.Polynomial(Polynomial.toPolynomial(c, ExpressionParser.parse(v2))));
  B = B.map(c => new Expression.Polynomial(Polynomial.toPolynomial(c, ExpressionParser.parse(v2))));
  if (A.getDegree() < B.getDegree()) {
    const tmp = A;
    A = B;
    B = tmp;
    //TODO: change the sign
  }
  const resultant2 = [];
  resultant2.push({base: Expression.ONE, exponent: -2 * B.getDegree()});
  var isPseudoRemainderSequence = true;
  for (var tmp of Polynomial._subresultantPseudoRemainderSequence(A, B)) {
    var A = tmp.A;
    var B = tmp.B;
    var R = tmp.R;
    var α = tmp.α;
    // https://en.wikipedia.org/wiki/Resultant#Properties
    // b_0**(deg(A - Q * B) - deg(A)) * res(A, B) = res(B, A - Q * B)
    //resultant = resultant.multiply(B.getLeadingCoefficient()._pow(A.getDegree() - Math.max(R.getDegree(), 0)));
    //resultant = resultant.multiply(α._pow(B.getDegree()));
    //resultant = resultant.divide(scale._pow(B.getDegree()));
    
    const previous = resultant2.pop();
    console.assert(previous.exponent === -2 * B.getDegree());
    resultant2.push({base: α.divide(previous.base._pow(2)), exponent: B.getDegree()});
    resultant2.push({base: B.getLeadingCoefficient(), exponent: Math.max(R.getDegree(), 0) + (isPseudoRemainderSequence ? (1 - B.getDegree()) * (A.getDegree() - B.getDegree()) : A.getDegree())});
    resultant2.push({base: B.getLeadingCoefficient(), exponent: 0 - 2 * Math.max(R.getDegree(), 0)});
  }
  let resultant = Expression.ONE;
  for (const x of resultant2) {
    if (x.exponent < 0) {
      resultant = resultant.divide(x.base._pow(-x.exponent));
    } else if (x.exponent > 0) {
      resultant = resultant.multiply(x.base._pow(x.exponent));
    }
  }
  resultant = resultant.polynomial;
  //resultant = Polynomial.toPolynomial(resultant, ExpressionParser.parse(v2));
  return resultant;
};

Polynomial.prototype.subs = function (variableMapFunction) {
  var variable = new Expression.Symbol('$x');//TODO:
  return Polynomial.toPolynomial(this.calcAt(variableMapFunction(variable)).getNumerator(), variable);
};


function GramSchmidt(rowVectorsMatrix) {
  if (false) {
    const V = rowVectorsMatrix;
    const n = V.cols();
    const k = V.rows();
    const U = new Array(k).fill(null).map(x => new Matrix.Vector(new Array(n).fill(Expression.ZERO)));;
    U[0] = V.row(0);
    for (let i = 1; i < k; i += 1) {
        U[i] = V.row(i);
        for (let j = 0; j < i; j += 1) {
            U[i] = U[i].subtract(U[j].scale(U[j].dot(U[i]).divide(U[j].dot(U[j]))));
        }
    }
    return Matrix.Zero(k, n).map((e, i, j) => U[i].e(j));
  }
  // https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process#Via_Gaussian_elimination
  /*
  function matrixWithVectorRows(vectors) {
    return Matrix.Zero(vectors.length, vectors[0].length).map(function (e, i, j) {
      console.assert(vectors[i].length === vectors[0].length);
      return vectors[i][j];
    });
  }
  var A = matrixWithVectorRows(vectors);
  */
  var A = rowVectorsMatrix;
  var matrix = A.multiply(A.transpose()).augment(A).toRowEchelon(Matrix.Gauss, "inverse").matrix;
  return matrix.slice(0, matrix.rows(), A.rows(), matrix.cols());
}

globalThis.GramSchmidtOrthogonalization = function (vectors) {
  var tmp = GramSchmidt(Matrix.Zero(vectors.length, Math.max.apply(null, vectors.map(vector => vector.dimensions()))).map(function (e, i, j) {
    return vectors[i].e(j);
  }));
  var result = [];
  for (var i = 0; i < vectors.length; i += 1) {
    result.push(tmp.row(i));
  }
  return result;
};
Matrix.prototype.orthogonalizeColumnVectors = function () {
  return GramSchmidt(this);
};

Expression.Complex.prototype.abs = function () {
  // https://en.wikipedia.org/wiki/Absolute_value#Complex_numbers
  return this.multiply(this.conjugate()).squareRoot();
};
Expression.Division.prototype.abs = function () {
  return this.getNumerator().abs().divide(this.getDenominator().abs());
};
Expression.prototype.abs = function () {//TODO: remove - ?
  if (this.compareTo(Expression.ZERO) < 0) {
    return this.negate();
  }
  return this;
};
Expression.prototype.compareTo = function (other) {//TODO: remove - ?
  if (other.equals(Expression.ZERO)) {
    if (Expression._isPositive(this)) {
      return +1;
    }
    if (Expression._isPositive(this.negate())) {
      return -1;
    }
    throw new TypeError(this.toString());
  }
  return this.subtract(other).getNumerator().compareTo(Expression.ZERO);
};
Expression.prototype.round = function () {//TODO: remove - ?
  //TODO: half away from zero - ?
  //console.log(this.getNumerator(), this.getDenominator());
  //return this.getNumerator().add(this.getDenominator().truncatingDivide(Expression.TWO)).truncatingDivide(this.getDenominator());
  return ExpressionParser.parse(toDecimalStringInternal(this, {fractionDigits: 0}));
};

//console.assert(GramSchmidt(new Matrix([[Expression.Integer.fromNumber(3), Expression.Integer.fromNumber(1)], [Expression.Integer.fromNumber(2), Expression.Integer.fromNumber(2)]])).toString() === '{{3,1},{-2/5,6/5}}');
//console.assert(GramSchmidt(new Matrix([[Expression.Integer.fromNumber(3), Expression.Integer.fromNumber(1)], [Expression.Integer.fromNumber(2), Expression.Integer.fromNumber(2)], [new Expression.Integer(0), new Expression.Integer(0)]])).toString() === '{{3,1},{-2/5,6/5},{0,0}}');
//throw new Error();


// Math.log2(Math.hypot.apply(null, coefficients))
Polynomial.prototype._log2hypot = function () {
  const polynomial = this;
  let max = Expression.ZERO;
  for (var i = 0; i < polynomial.a.size; i += 1) {
    var c = polynomial.a.coefficient(i).abs();
    if (c.compareTo(max) > 0) {
      max = c;
    }
  }
  //const maxBitLength = Math.max.apply(null, coefficients.map(c => c.equals(Expression.ZERO) ? 0 : primeFactor._bitLength(c.abs().toBigInt())));
  const maxBitLength = max.toNumber() < 1 / 0 ? 0 : primeFactor._bitLength(max.toBigInt());
  const k = maxBitLength < 1024 ? maxBitLength : Math.min(Math.floor(Math.log2(Number.MAX_SAFE_INTEGER + 1)), maxBitLength);
  const scale = Expression.TWO._pow(maxBitLength - k);
  const p2k = 2**k;
  const coefficients = new Array(polynomial.a.size);
  for (var i = 0; i < polynomial.a.size; i += 1) {
    coefficients[i] = polynomial.a.coefficient(i).truncatingDivide(scale).toNumber() / p2k;
  }
  const hypot = Math.hypot.apply(null, coefficients);
  const log2hypot = maxBitLength + Math.log2(hypot);
  return log2hypot;
};

Polynomial.prototype._log2OfBoundForCoefficientsOfFactor = function (factorDegreeBound, factorLeadingCoefficientBound) {
  // https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#:~:text=This%20bound%20is%20also%20useful%20to%20bound%20the%20coefficients%20of%20a%20divisor%20of%20a%20polynomial%20with%20integer%20coefficients:
  // see also
  // The art of computer programming. Vol.2: Seminumerical algorithms
  // exersize 20, page 458
  // which gives better result (~2 times smaller)
  if (factorDegreeBound == undefined) {
    factorDegreeBound = Math.floor(this.getDegree() / 2);
  }
  if (factorLeadingCoefficientBound == undefined) {
    factorLeadingCoefficientBound = this.getLeadingCoefficient().abs();
  }
  var log2 = function (integer) {
    var e = primeFactor._bitLength(integer.toBigInt());
    if (e <= 53) {
      return Math.log2(integer.toNumber());
    }
    return (e - 53) + Math.log2(integer.truncatingDivide(Expression.TWO._pow(e - 53)).toNumber());
  };
  var m = factorDegreeBound;
  return (m - Math.log2(Math.sqrt(Math.PI * Math.ceil(m / 2)))) + (log2(factorLeadingCoefficientBound) - log2(this.getLeadingCoefficient().abs())) + this._log2hypot();
};


Polynomial.prototype.isDivisibleBy = function (guess) {
  var tmp = this.divideAndRemainder(guess, "undefined");
  return tmp != null && tmp.remainder.equals(Polynomial.ZERO);
};