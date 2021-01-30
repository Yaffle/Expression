import Expression from './Expression.js';
import Matrix from './Matrix.js';
import Polynomial from './Polynomial.js';
import toDecimalStringInternal from './toDecimalString.js';
import bitLength from './bitLength.js';
import {BigFloat} from './BigDecimal/BigDecimal.js';

  var isRational = function (e) {
    return e.getNumerator() instanceof Expression.Integer && e.getDenominator() instanceof Expression.Integer;
  };

function SimpleInterval(a, b) {
  console.assert(isRational(a));
  console.assert(isRational(b));
  console.assert(a.getDenominator().multiply(b.getNumerator()).compareTo(b.getDenominator().multiply(a.getNumerator())) >= 0);
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
SimpleInterval.prototype.toString = function () {
  return '[' + this.a + ';' + this.b + ']';
};

function PolynomialRoot(polynomial, interval) {
  Expression.Symbol.call(this, "[root of " + polynomial + " near " + interval.a.add(interval.b).divide(Expression.TWO).toString() + "]");
  this.polynomial = polynomial;
  //TODO: https://www.wolframalpha.com/input/?i=x**5%2B7x**3%2Bx**2%2Bx%2B1%3D0
  this.interval = interval;
}
PolynomialRoot.prototype = Object.create(Expression.Symbol.prototype);

PolynomialRoot.prototype.toDecimal = function (precision) {
  var tmp = this.polynomial.getZero(this.interval, precision);
  return new SimpleInterval(tmp.a, tmp.b);
};

Expression.PolynomialRoot = PolynomialRoot;

function isSameRoot(x, y) {
  return x instanceof PolynomialRoot && y instanceof PolynomialRoot && x.polynomial.equals(y.polynomial) && x.interval.a.equals(y.interval.a) && x.interval.b.equals(y.interval.b);
}

function PolynomialRoot1(polynomial, interval) {
  PolynomialRoot.call(this, polynomial, interval);
  if (!polynomial.hasIntegerCoefficients()) {
    //?new
    var variable = new Expression.Symbol('$$');
    var e = polynomial.calcAt(variable);
    var c = Expression.getConjugateExpression(e);
    if (c != null && !c.equals(e)) {
      //TODO: what if multiple (?) - ?
      return new PolynomialRoot1(Polynomial.toPolynomial(c, variable), interval);
    }
    //?
  }
  if (polynomial.getLeadingCoefficient().compareTo(Expression.ZERO) < 0) {
    return new PolynomialRoot1(polynomial.negate(), interval);
  }
  if (!polynomial.hasIntegerCoefficients()) {
    throw new TypeError();
  }
  /*
  var rr = null;
  do {
    var rr = polynomial.doRationalRootTest();
    if (rr != null) {
      //TODO: >= 0 && <= 0
      if (Expression._isPositive(rr.subtract(interval.a)) && Expression._isPositive(rr.subtract(interval.b).negate())) {
        return rr;//TODO: MOVE!
      }
      polynomial = polynomial.divideAndRemainder(Polynomial.of(rr.getNumerator().negate(), rr.getDenominator()), "throw").quotient;//?
    }
  } while (rr != null);
  */
  if (polynomial.getDegree() === 1 || polynomial.getDegree() === 2 || polynomial.getDegree() === 4) {//TODO: other - ? like biqudratic - ?
    var roots = polynomial.getroots();
    for (var rr of roots) {
      if (!Expression.has(rr, ExpressionWithPolynomialRoot)) {//?
        if (Expression._isPositive(rr.subtract(interval.a)) && Expression._isPositive(rr.subtract(interval.b).negate()) || rr.equals(interval.b)) {
          return rr;//TODO: MOVE!
        }
      }
    }
  }
  var content = polynomial.getContent();
  if (!content.equals(Expression.ONE)) {
    return new PolynomialRoot1(polynomial.scale(content.inverse()), interval);
  }
  var factor = polynomial.factorize();
  // var factor = polynomial._factorByLLL({zero: new PolynomialRoot(polynomial, interval)});
  if (factor != null && !factor.equals(polynomial)) {
    if (factor.numberOfRoots(interval) !== 0) {
      return new PolynomialRoot1(factor, interval);
    } else {
      var otherFactor = polynomial.divideAndRemainder(factor, "throw").quotient;
      return new PolynomialRoot1(otherFactor, interval);
    }
  }
  if (!isRational(interval.a) || !isRational(interval.b)) {
    throw new TypeError();
  }
  if (interval.a.subtract(interval.b).getNumerator().compareTo(Expression.ZERO) >= 0) {
    throw new TypeError();
  }
  if ((interval.a.getNumerator().sign() || interval.b.getNumerator().sign()) !== (interval.b.getNumerator().sign() || interval.a.getNumerator().sign())) {
    throw new TypeError();
  }
  if (polynomial.numberOfRoots(interval) !== 1) {
    throw new TypeError();
  }
  if (!polynomial.getContent().equals(Expression.ONE)) {
    throw new TypeError();
  }
  if (polynomial.getDegree() < 3) {
    throw new TypeError();
  }
  if (polynomial.getDegree() > 64) {
    throw new Error();//TODO: too long
  }
}
PolynomialRoot1.prototype = Object.create(PolynomialRoot.prototype);

PolynomialRoot1.SimpleInterval = SimpleInterval;

  var calculateNewInterval = function (newPolynomial, zeroFunction) {
    var precision = 1;
    var guess = zeroFunction(precision);
    while ((guess.a.getNumerator().sign() !== guess.b.getNumerator().sign() && !newPolynomial.calcAt(Expression.ZERO).equals(Expression.ZERO)) || newPolynomial.numberOfRoots(guess) > 1) {
      precision *= 2;
      guess = zeroFunction(precision);
      if (precision > 1024) throw new Error();//TODO: ?
    }
    return guess;
  };
PolynomialRoot1.prototype.multiplyInteger = function (x) {
  return x.multiplyPolynomialRoot1(this);
};
PolynomialRoot1.prototype.multiply = function (e) {
  return e.multiplyPolynomialRoot1(this);
};
PolynomialRoot1.prototype.multiplyExpression = function (e) {
  //?
  //TODO: fix
  if (Expression.isConstant(e) && !Expression.has(e, Expression.Complex)) {
    return this.multiply(e);
  }
  return Expression.Symbol.prototype.multiplyExpression.call(this, e);
};
PolynomialRoot1.prototype.multiplyAddition = function (e) { // for performance (?) when `e` is a constant
  if (Expression.isConstant(e) && !Expression.has(e, Expression.Complex)) {
    return this.multiplyExpression(e);
  }
  return Expression.Symbol.prototype.multiplyAddition.call(this, e);
};
Expression.prototype.multiplyPolynomialRoot1 = function (root) {
  if (Expression.isConstant(this) && !(Expression.has(this, Expression.Complex))) {
    var k = this;
    if (k.equals(Expression.ZERO)) {
      return k;
    }
    if (k.equals(Expression.ONE)) {
      return root;
    }
    // z = k * x, x = z / k
    //TODO: fix: use division instead of multiplication by inverse (?)
    var newPolynomial = root.polynomial.subs(x => x.divide(k));
    newPolynomial = newPolynomial.scale(newPolynomial.getContent().inverse());
    //TODO: FIX!!!
    k = Expression.has(k, Expression.NthRoot) ? RPN(toDecimalStringInternal(k, {significantDigits: 10})) : k;
    var newInterval = Expression._isPositive(k) ? {a: root.interval.a.multiply(k), b: root.interval.b.multiply(k)} : {a: root.interval.b.multiply(k), b: root.interval.a.multiply(k)};
    return new PolynomialRoot1(newPolynomial, newInterval);
  }
  if (this instanceof Expression.Complex && !this.imaginary.equals(Expression.ONE)) {
    return this.divide(this.imaginary).multiply(this.imaginary.multiply(root));
  }
  //TODO: ?
  //throw new Error();
  return this.multiplyExpression(root);
};
PolynomialRoot1.prototype._pow = function (count) {
  //if (e instanceof Expression.Integer) {
    var degree = this.polynomial.getDegree();
    for (var i = 0; i <= degree; i += 1) {
      if (!this.polynomial.getCoefficient(i).equals(Expression.ZERO)) {
        if (i % count !== 0) {
          return Expression.prototype._pow.call(this, count);
        }
      }
    }
    //TODO: faster method
    var newInterval = {a: this.interval.a._pow(count), b: this.interval.b._pow(count)};
    if (count % 2 === 0 && this.interval.b.compareTo(Expression.ZERO) < 0) {
      //throw new Error();//TODO:
      newInterval = {a: this.interval.b._pow(count), b: this.interval.a._pow(count)};
    }
    return new PolynomialRoot1(this.polynomial.subs(x => x._nthRoot(count)), newInterval);
    //throw new Error("!");
  //}
  //TODO: ?
  //throw new Error();
};
PolynomialRoot1.prototype.pow = function (e) {
  if (e instanceof Expression.Integer) {
    return this._pow(e.toNumber());
  }
  //TODO: ?
  throw new Error();
};
PolynomialRoot1.prototype.multiplyPolynomialRoot1 = function (x) {
  //var x = this;
  var y = this;
  //TODO: variable names
  var newPolynomial = Polynomial.toPolynomial(Polynomial.resultant(x.polynomial, Polynomial.toPolynomial(y.polynomial.calcAt(RPN("z/x")).multiply(RPN("x")._pow(y.polynomial.getDegree())), RPN("x"))), RPN("z"));
  //var tmp = newPolynomial.squareFreeFactors();
  //newPolynomial = newPolynomial.getSquareFreePolynomial();//TODO: ?
  newPolynomial = newPolynomial.scale(newPolynomial.getContent().inverse());
  //var s = toDecimalStringInternal(new Expression.Multiplication(new PolynomialRoot(x.polynomial, x.interval), new PolynomialRoot(y.polynomial, y.interval)), {significantDigits: 10});
  //if (s == undefined) {
    //debugger;
  //  throw new Error();
  //}
  //TODO consider sign (?)
  //var interval = {a: RPN(s + '*(1-10**-9)'), b: RPN(s + '*(1+10**-9)')};
  //if (s.startsWith('-')) {
  //  interval = {a: interval.b, b: interval.a};
  //}
  var interval = calculateNewInterval(newPolynomial, function (precision) {
    return x.toDecimal(precision).multiply(y.toDecimal(precision));
  });
  //TODO: validate interval
  //var tmp = newPolynomial.squareFreeFactors();
  //while (tmp.a1.numberOfRoots(interval) === 0) {
  //  tmp = tmp.a0;
  //  tmp = tmp.squareFreeFactors();
  //}
  //newPolynomial = tmp.a1;
  var root = new PolynomialRoot1(newPolynomial, interval);
  return root;
};
PolynomialRoot1.prototype.add = function (e) {
  return e.addPolynomialRoot1(this);
};
PolynomialRoot1.prototype.addPolynomialRoot1 = function (x) {
  //var x = this;
  var y = this;
  //TODO: variable names
  var newPolynomial = Polynomial.toPolynomial(Polynomial.resultant(x.polynomial, Polynomial.toPolynomial(y.polynomial.calcAt(RPN("z-x")), RPN("x"))), RPN("z"));
  //newPolynomial = newPolynomial.getSquareFreePolynomial(); //TODO: fix
  //var zero = {a: x.interval.a.add(y.interval.a), b: x.interval.b.add(y.interval.b)};
  var zero = calculateNewInterval(newPolynomial, function (precision) {
    return x.toDecimal(precision).add(y.toDecimal(precision));
    //return toDecimalStringInternal(new Expression.Addition(x, y), {significantDigits: precision});
  });
  //TODO: validate interval
  var root = new PolynomialRoot1(newPolynomial, zero);
  //return new ExpressionWithPolynomialRoot(root, root);
  return root;
};
Expression.prototype.addPolynomialRoot1 = function (e) {
  if (Expression.isConstant(this) && !(Expression.has(this, Expression.Complex))) {
    if (this.equals(Expression.ZERO)) { // for performance
      return e;
    }
    // z = x + k, x = z - k
    var k = this;
    var root = e;
    var newPolynomial = root.polynomial.subs(x => x.subtract(k));
    newPolynomial = newPolynomial.scale(newPolynomial.getContent().inverse());
    var newInterval = null;
    if (!(isRational(k))) {
      newInterval = calculateNewInterval(newPolynomial, function (precision) {
        var tmp = RPN(toDecimalStringInternal(k, {significantDigits: precision}));
        var e = new SimpleInterval(RPN('1-5*10**-' + precision), RPN('1+5*10**-' + precision));
        return root.toDecimal(precision).add(new SimpleInterval(tmp, tmp).multiply(e));
      });
    } else {
      // to avoid intervals, which include zero
      newInterval = calculateNewInterval(newPolynomial, function (precision) {
        return root.toDecimal(precision).add(new SimpleInterval(k, k));
      });
    }
    return new PolynomialRoot1(newPolynomial, newInterval);
  }
  //throw new Error();
  return this.addExpression(e);
};
PolynomialRoot1.prototype.addExpression = function (e) {
  return this.add(e);//!?
};
PolynomialRoot1.prototype.divide = function (e) {
  if (e instanceof ExpressionWithPolynomialRoot) {
    return this.divide(e.upgrade());
  }
  if (!(e instanceof PolynomialRoot1) && !Expression.isConstant(e) || Expression.has(e, Expression.Matrix) || Expression.has(e, Expression.MatrixSymbol)) {
    throw new Error();
  }
  return this.multiply(e.inverse());
};
PolynomialRoot1.prototype.divideExpression = function (x) {
  return x.multiply(this.inverse());
};
PolynomialRoot1.prototype.inverse = function () {
  //?
  // z = 1/y, y = 1/z
  var root = this;
  var newPolynomial = root.polynomial.subs(x => x.inverse());
  newPolynomial = newPolynomial.scale(newPolynomial.getContent().inverse());
  var precision = 1;
  var interval = root.interval
  while (interval.a.getNumerator().compareTo(Expression.ZERO) <= 0 && interval.b.getNumerator().compareTo(Expression.ZERO) >= 0) {
    interval = root.polynomial.getZero(interval, precision);
    precision *= 2;
  }
  return new PolynomialRoot1(newPolynomial, {a: interval.b.inverse(), b: interval.a.inverse()});
};
PolynomialRoot1.prototype.toString = function (options) {
  //return new ExpressionWithPolynomialRoot(this, this).toString(options);
  options = options || {};
  var rounding = options.rounding != null ? options.rounding : {fractionDigits: 3};
  return toDecimalStringInternal(this, rounding);
};
PolynomialRoot1.prototype.equals = function (other) {
  if (other instanceof PolynomialRoot) {
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
    return this.polynomial.equals(other.polynomial) && intersection(this.interval, other.interval) != null && this.subtract(other).equals(Expression.ZERO);
  }
  if (Expression.isConstant(other)) {
    if (!this.polynomial.calcAt(other).equals(Expression.ZERO)) {
      return false;
    }
    var withinInterval = function (x, interval) {
      return Expression._isPositive(x.subtract(interval.a)) && Expression._isPositive(x.subtract(interval.b).negate());
    };
    return withinInterval(other, this.interval);
  }
  //TODO: optimize
  return this.subtract(other).equals(Expression.ZERO);
};
PolynomialRoot1.prototype.compare4MultiplicationComplex = function (x) {
  return -1;
};
PolynomialRoot1.prototype.compare4Multiplication = function (y) {
  if (y instanceof Expression.Complex) {
    return +1;
  }
  if (y instanceof Expression.Integer) {
    return +1;
  }
  if (y instanceof PolynomialRoot) {
    return 0;
  }
  return Expression.Symbol.prototype.compare4Multiplication.call(this, y);
};
PolynomialRoot1.prototype.compare4Addition = function (y) {
  if (y instanceof PolynomialRoot1) {
    return 0;//?
  }
  if (y instanceof Expression.Symbol) {
    return +1;
  }
  return Expression.Symbol.prototype.compare4Addition.call(this, y);
};
PolynomialRoot1.prototype.compare4AdditionSymbol = function (x) {
  return -1;
};

PolynomialRoot1.prototype._nthRoot = function (n) {//?
  // PolynomialRoot#_nthRoot - ?
  if (this.interval.a.getNumerator().compareTo(Expression.ZERO) < 0) {
    //TODO: check
    var newRoot = new PolynomialRoot(this.polynomial.subs(x => x.negate()), {a: this.interval.b.negate(), b: this.interval.a.negate()});
    return Expression.I.multiply(newRoot._nthRoot(n));
  }
  var newPolynomial = this.polynomial.subs(function (x) { return x._pow(n); });
  var root = this;
  var newInterval = calculateNewInterval(newPolynomial, function (precision) {
    //TODO: 
    //return root.toDecimal(precision).nthRoot(n);
    var a = RPN(toDecimalStringInternal(new Expression.NthRoot('', root, n), {significantDigits: precision}));
    var b = RPN(toDecimalStringInternal(new Expression.NthRoot('', root, n), {significantDigits: precision}));
    a = a.multiply(RPN('1-2**-' + precision));
    b = b.multiply(RPN('1+2**-' + precision));
    return new SimpleInterval(a, b);
  });
  var newRoot = new PolynomialRoot1(newPolynomial, newInterval);
  return newRoot;
};

Expression.PolynomialRoot1 = PolynomialRoot1;

Expression.toPolynomialRoot = function (e) {
  var x = e instanceof Expression.NthRoot ? e.a : e;//TODO: remove
  var n = e instanceof Expression.NthRoot ? e.n : 1;//TODO: remove
  var symbol = new Expression.Symbol('x');
  var p = Polynomial.toPolynomial(Expression.getConjugateExpression(symbol._pow(n).subtract(x)), symbol);
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
    var zeros = p.getZeros();
    //TODO: find zero only on interval
    for (var zero of zeros) {
      if (Expression._isPositive(zero.e) && zero._pow(n).equals(x)) {
        return zero;
      }
    }
    //TODO: ?
  }
};


function ExpressionWithPolynomialRoot(e, root) {
  this.e = e; // internal symbolic expression with a "root" as a symbol
  this.root = root;
}

Expression.ExpressionWithPolynomialRoot = ExpressionWithPolynomialRoot;

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
    //!new 2020-08-27
    //TODO: remove
    if (true && !(x instanceof Expression.Integer) && !(x instanceof Expression.Multiplication && x.a === Expression.I && x.b === v)) {
      var test = v.polynomial.divideAndRemainder(Polynomial.toPolynomial(x.subtract(new Expression.Symbol('$n')), v)).remainder;
      if (test.getDegree() === 0) {
        //(x**2-2)(x**2+x-1) = 0
        var pn = Polynomial.toPolynomial(Expression.getConjugateExpression(test.calcAt(Expression.ZERO).getNumerator()), new Expression.Symbol('$n'));
        //pn = pn.scale(pn.getLeadingCoefficient().inverse());
        pn = pn.scale(pn.getContent().inverse());
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
              if (pn.numberOfRoots({a: q.subtract(RPN('0.5').divide(lc)), b: q.add(RPN('0.5').divide(lc))}) === 1) {
                return q;
              }
            }
          }
        }
      }
    }
    var px = Polynomial.toPolynomial(x, v);
    if (false) {
      if (v.polynomial.getDegree() === 6 && v.polynomial.getCoefficient(1).equals(Expression.ZERO) && v.polynomial.getCoefficient(2).equals(Expression.ZERO) && v.polynomial.getCoefficient(4).equals(Expression.ZERO) && v.polynomial.getCoefficient(5).equals(Expression.ZERO)) {
        if (px.getDegree() >= 3) {
          var alpha = c(v._pow(3));
          var dv = Polynomial.of(alpha.negate(), Expression.ZERO, Expression.ZERO, Expression.ONE);
          return px.divideAndRemainder(dv).remainder.calcAt(v);
        }
      }
    }
    //!
    return px.divideAndRemainder(v.polynomial).remainder.calcAt(v);
  };
  var oldE = e;
  e = c(e.getNumerator()).divide(c(e.getDenominator()));
  if (!oldE.equals(e)) {
    var oldE1 = e;
    e = c(e.getNumerator()).divide(c(e.getDenominator())); // something may change after the previous step
  }

  if (!Expression.has(e, Expression.PolynomialRoot)) {
    return e;
  }

  return new ExpressionWithPolynomialRoot(e, root);
}

ExpressionWithPolynomialRoot.prototype = Object.create(Expression.Symbol.prototype);

/*
ExpressionWithPolynomialRoot.prototype.compare4Multiplication = function (y) {
  return y.compare4MultiplicationExpressionWithPolynomialRoot(this);
};
ExpressionWithPolynomialRoot.prototype.compare4MultiplicationExpressionWithPolynomialRoot = function (x) {
  return 0;//?
};
Expression.prototype.compare4MultiplicationExpressionWithPolynomialRoot = function (x) {
  return this.compare4MultiplicationSymbol(x);//?
};
Expression.Symbol.prototype.compare4MultiplicationExpressionWithPolynomialRoot = function (x) {
  return -1;
};
ExpressionWithPolynomialRoot.prototype.compare4MultiplicationSymbol = function (x) {
  return +1;
};
*/

//Expression.prototype.isExact = function () {
//  return true;
//};
ExpressionWithPolynomialRoot.prototype.isExact = function () {
  //TODO: fix - ?
  return false;
};

PolynomialRoot.prototype.isExact = function () {
  //TODO: fix - ?
  return false;
};

ExpressionWithPolynomialRoot.prototype.negate = function () {
  //return makeExpressionWithPolynomialRoot(this.e.negate(), this.root);
  return new ExpressionWithPolynomialRoot(this.e.negate(), this.root); // for performance
};
ExpressionWithPolynomialRoot.prototype.equals = function (other) {
  //!TODO: remove (hack to avoid error)
  if (this instanceof ExpressionWithPolynomialRoot && other instanceof ExpressionWithPolynomialRoot) {
    if (!isSameRoot(this.root, other.root)) {
      //var s1 = toDecimalStringInternal(new Expression.Addition(this.e, other.e.negate()), {fractionDigits: 3});
      //if (s1 != undefined && !s1.endsWith('000')) {
      //  return false;
      //}
      if (this.toMathML() !== other.toMathML()) {
        return false;//?
      }
      if (true) {
        return upgrade(this.e).subtract(upgrade(other.e)).equals(Expression.ZERO);
      }
      var s = toDecimalStringInternal(new Expression.Addition(this.e, other.e.negate()), {significantDigits: 1});
      //TODO: will it hang for zero?
      return s === '0';
    }
  }
  //!
  //if (Expression.has(other, Expression.NthRoot)) {
  //  return this.upgrade().equals(Expression.toPolynomialRoot(other));//!?
  //}
  // optimization
  var s = other instanceof Expression.Integer && other.equals(Expression.ZERO) ? this : this.subtract(other);
  return s instanceof ExpressionWithPolynomialRoot ? false : s.equals(Expression.ZERO);
};
ExpressionWithPolynomialRoot.prototype.simplifyExpression = function () {
  return this;
};

ExpressionWithPolynomialRoot.prototype.toString = function (options) {
  options = options || {};
  //TODO: return 'polynomial-root of x**2+2x+1 on [a; b]';
  //TODO:
  if (this.equals(Expression.ZERO)) {
    return Expression.ZERO.toString(options);
  }
  //return this.e.toString(options);
  if (options.fractionDigits !== -1 && options.fractionDigits != null) {
    console.debug('options.fractionDigits is deprecated, please use options.rounding');
  }
  var rounding = options.rounding != null ? options.rounding : (options.fractionDigits !== -1 && options.fractionDigits != null ? {fractionDigits: options.fractionDigits} : {fractionDigits: 3});
  //if (true) {
  //  return Expression.toDecimalString(this.e, Object.assign({}, options, {rounding: rounding}));
  //}
  var tmp = toDecimalStringInternal(this.e, rounding, undefined, undefined);
  return tmp;
};

PolynomialRoot.prototype.toMathML = function (options) {
  options = options || {};
  if (options.fractionDigits !== -1 && options.fractionDigits != null) {
    console.debug('options.fractionDigits is deprecated, please use options.rounding');
  }
  var rounding = options.rounding != null ? options.rounding : (options.fractionDigits !== -1 && options.fractionDigits != null ? {fractionDigits: options.fractionDigits} : {fractionDigits: 3});
  var tmp = toDecimalStringInternal(this, rounding, Expression._decimalToMathML, Expression._complexToMathML);
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
  if (true) {
    //return Expression.toDecimalString(this.e, );
    return this.e.toMathML(Object.assign({}, options, {rounding: rounding}));
  }
  var tmp = toDecimalStringInternal(this.e, rounding, Expression._decimalToMathML, Expression._complexToMathML);
  return tmp;
};

function upgrade(e) {
  var cache = null;//TODO: ?
  var root = null;
  return Expression._map(function (x) {
    return x instanceof PolynomialRoot && !(x instanceof PolynomialRoot1) ? (x === cache ? root : (cache = x, root = new PolynomialRoot1(x.polynomial, x.interval))) : x;
  }, e);
}

ExpressionWithPolynomialRoot.prototype.multiply = function (other) {
  if (other instanceof ExpressionWithPolynomialRoot) {
    if (this.root !== other.root && !isSameRoot(this.root, other.root)) {
      return upgrade(this.e).multiply(upgrade(other.e));
    }
    return this.multiply(other.e);
  }
  if (Expression.has(other, PolynomialRoot1)) {//TODO: ?
    return this.upgrade().multiply(other);
  }
  return makeExpressionWithPolynomialRoot(this.e.multiply(other), this.root);
};
ExpressionWithPolynomialRoot.prototype.divide = function (other) {
  if (other.equals(Expression.ONE)) {
    return this;
  }
  if (other instanceof ExpressionWithPolynomialRoot) {
    if (this.root !== other.root && !isSameRoot(this.root, other.root)) {
      return upgrade(this.e).divide(upgrade(other.e));
    }
    return this.divide(other.e);
  }
  if (Expression.has(other, PolynomialRoot1)) {//TODO: ?
    return this.upgrade().divide(other);
  }
  return makeExpressionWithPolynomialRoot(this.e.divide(other), this.root);
};
ExpressionWithPolynomialRoot.prototype.add = function (other) {
  if (other instanceof ExpressionWithPolynomialRoot) {
    if (this.root !== other.root && !isSameRoot(this.root, other.root)) {
      return upgrade(this.e).add(upgrade(other.e));
    }
    return this.add(other.e);
  }
  if (Expression.has(other, PolynomialRoot1)) {//TODO: ?
    return this.upgrade().add(other);
  }
  return makeExpressionWithPolynomialRoot(this.e.add(other), this.root);
};

ExpressionWithPolynomialRoot.prototype.divideExpression = function (other) {
  if (Expression.has(other, PolynomialRoot1)) {//TODO: ?
    return other.divide(this.upgrade());
  }
  return makeExpressionWithPolynomialRoot(other.divide(this.e), this.root);
};
ExpressionWithPolynomialRoot.prototype.multiplyExpression = function (other) {
  if (other.equals(Expression.ONE)) {
    return this;
  }
  if (Expression.has(other, PolynomialRoot1)) {//TODO: ?
    return other.multiply(this.upgrade());
  }
  return makeExpressionWithPolynomialRoot(other.multiply(this.e), this.root);
};
ExpressionWithPolynomialRoot.prototype.addExpression = function (other) {
  if (Expression.has(other, PolynomialRoot1)) {//TODO: ?
    return other.add(this.upgrade());
  }
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
  if (this.e === this.root) {//TODO: ?
    // PolynomialRoot#_nthRoot - ?
    if (this.root.interval.a.getNumerator().compareTo(Expression.ZERO) < 0) {
      //TODO: check
      var newRoot = new PolynomialRoot(this.root.polynomial.subs(x => x.negate()), {a: this.root.interval.b.negate(), b: this.root.interval.a.negate()});
      return Expression.I.multiply(new ExpressionWithPolynomialRoot(newRoot, newRoot)._nthRoot(n));
    }
    var precision = 3;
    var a = null;
    var b = null;
    do {
      a = RPN(toDecimalStringInternal(this.root.interval.a._nthRoot(n), {significantDigits: precision}));
      b = RPN(toDecimalStringInternal(this.root.interval.b._nthRoot(n), {significantDigits: precision}));
      precision *= 2;
      if (precision > 128) {
        debugger;
        throw new TypeError();
      }
    //TODO: fix !!!
    } while (a.equals(b));
    var newRoot = new PolynomialRoot(this.root.polynomial.subs(function (x) { return x._pow(n); }), {a: a, b: b});
    if (newRoot.polynomial.numberOfRoots(newRoot.interval) === 1) {
      return new ExpressionWithPolynomialRoot(newRoot, newRoot);
    } else {
      console.assert(false);
      debugger;
    }
  }
  if (true && n === 2) {
    return this.upgrade()._nthRoot(n);
  }
  return makeExpressionWithPolynomialRoot(this.e._nthRoot(n), this.root);
};
ExpressionWithPolynomialRoot.prototype.pow = function (count) {
  if (count instanceof Expression.Integer) {
    return makeExpressionWithPolynomialRoot(this.e._pow(count.value), this.root);
  } else {
    if (count instanceof Expression.Division && count.getDenominator() instanceof Expression.Integer) {
      return makeExpressionWithPolynomialRoot(this.e.pow(count.getNumerator()), this.root)._nthRoot(count.getDenominator().value);
    }
    //TODO: upgrade (?)
    return makeExpressionWithPolynomialRoot(this.e.pow(count), this.root);
  }
};
ExpressionWithPolynomialRoot.prototype._pow = function (count) {
  return makeExpressionWithPolynomialRoot(this.e.getNumerator()._pow(count), this.root).divide(makeExpressionWithPolynomialRoot(this.e.getDenominator()._pow(count), this.root));
};

//TODO: remove
ExpressionWithPolynomialRoot.prototype.upgrade = function () {
  return upgrade(this.e);
};

ExpressionWithPolynomialRoot.prototype.complexConjugate = function () {
  return makeExpressionWithPolynomialRoot(this.e.complexConjugate(), this.root);
};

  //TODO: remove
  Polynomial.pseudoRemainder2 = function (x, y) {
    var lcg = y.getLeadingCoefficient().abs();
    var n = x.getDegree() - y.getDegree();
    // assertion
    if (n < 0) {
      throw new RangeError();
    }
    //var sx = x.multiply(Polynomial.of(Expression.pow(lcg, n).multiply(lcg)));
    var sx = x.scale(lcg._pow(n + 1));
    return sx.divideAndRemainder(y, "throw").remainder;
  };

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
      var fn = Polynomial.pseudoRemainder2(fp, fc).negate();
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
      //!?2020-01-17
      if (!(v instanceof Expression.Integer)) {
        if (Expression._isPositive(v)) {
          v = Expression.ONE;
        } else if (Expression._isPositive(v.negate())) {
          v = Expression.ONE.negate();
        } else {
          throw new TypeError(v.toString());
        }
      }
      //!?
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
    var interval = this.getRootsInterval();
    var getIntervals = function (interval) {
      var n = sturmSequence.numberOfRoots(interval);
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
    //!2020-12-19
    // rount to a power of two:
    M = Expression.TWO._pow(bitLength(M.getNumerator().toBigInt()) - bitLength(M.getDenominator().toBigInt()) + 1);
    //!
    return {a: M.negate(), b: M};
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
    //const BASE = Expression.TEN;
    const BASE = Expression.TWO;
    var e = Expression.pow(BASE, precision); // epsilon^-1
    if (!(e instanceof Expression.Integer)) {
      throw new RangeError("epsilon^-1 is not an integer");
    }
    var a = interval.a;
    var b = interval.b;
    // (b - a) * 10**p > 1
    //TODO: fix to use precision, not fractionDigits
    if (e.multiply(a.getDenominator().multiply(b.getNumerator()).subtract(b.getDenominator().multiply(a.getNumerator()))).compareTo(b.getDenominator().multiply(a.getDenominator())) > 0) {
      //TODO:
      var tmp = true && precision >= 16 / Math.log10(BASE.toNumber()) ? this.getZero(interval, Math.floor(precision / 4)) : interval;
      a = tmp.a;
      b = tmp.b;

      var n = this.getDegree();
      var p = this.map(function (coefficient, degree) {
        return coefficient.multiply(Expression.pow(e, n - degree));
      });
      var sa = roundFloor(a, e).add(Expression.ONE); // a.getNumerator().multiply(e).truncatingDivide(a.getDenominator()).add(Expression.ONE);//?
      var sb = roundFloor(b, e); // b.getNumerator().multiply(e).truncatingDivide(b.getDenominator());//?
      console.assert(sa.divide(e).subtract(a).getNumerator().compareTo(Expression.ZERO) >= 0);
      console.assert(sb.divide(e).subtract(b).getNumerator().compareTo(Expression.ZERO) <= 0);
      //TODO: bigdecimal - ?
      if (p.calcAt(sa).getNumerator().compareTo(Expression.ZERO) !== (this.calcAt(a).getNumerator().compareTo(Expression.ZERO) || this.calcAt(b).negate().getNumerator().compareTo(Expression.ZERO) || p.calcAt(sa).getNumerator().compareTo(Expression.ZERO))) {
        return {a: a, b: sa.divide(e)};
      }
      if (p.calcAt(sb).getNumerator().compareTo(Expression.ZERO) !== this.calcAt(b).getNumerator().compareTo(Expression.ZERO)) {
        if (p.calcAt(sb).getNumerator().compareTo(Expression.ZERO) === 0) {
          return {a: sb.divide(e), b: sb.divide(e)};
        }
        return {a: sb.divide(e), b: b};
      }
      a = sa;
      b = sb;
      // bisection method
      var pa = Math.sign(p.calcAt(a).getNumerator().compareTo(Expression.ZERO));
      var pb = Math.sign(p.calcAt(b).getNumerator().compareTo(Expression.ZERO));
      if (pa === 0) {
        return {a: a.divide(e), b: a.divide(e)};
      }
      if (pb === 0) {
        return {a: b.divide(e), b: b.divide(e)};
      }
      if (pa === pb) {
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
        var v = Math.sign(p.calcAt(middle).getNumerator().compareTo(Expression.ZERO));
        if (v === pb) {
          b = middle;
        } else if (v === pa) {
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

  Polynomial.prototype.numberOfRoots = function (interval = null) {
    if (interval == null) {
      interval = this.getRootsInterval();//TODO: use (-1/0; +1/0)
    }
    var sturmSequence = new SturmSequence(this);
    return sturmSequence.numberOfRoots(interval);
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

  // Polynomial.toPolynomial(RPN("x^3-8x^2+21x-18"), RPN("x")).getZeros(3).toString()
  Polynomial.prototype.getZeros = function (precision = 1, complex = false) {
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
        for (var i = 0; i < tmp.length; i += 1) {
          var zero = tmp[i];
          if (zero instanceof ExpressionWithPolynomialRoot && zero.e === zero.root ? f.hasRoot(zero.root) : f.calcAt(zero).equals(Expression.ZERO)) {
            result.push(zero);
          } else {
            //TODO:?
            console.debug(zero.root);
          }
        }
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
          if (!Expression._isPositive(zero instanceof ExpressionWithPolynomialRoot ? zero.e : zero)) {
            var z = squareRoot(zero.e != null ? zero.upgrade() : zero);
            result.push(z);
            result.push(z.negate());
          }
        }
        return result;
      }
      //var p = stringToPolynomial("x^5+2*x^2+2*x+3");

      var e = p.calcAt(RPN("a+b*i"));
      var ce = Expression.getComplexConjugate(e);
      var pa = ce.add(e);//TODO: ?
      var pb = ce.subtract(e).multiply(Expression.I).divide(RPN('b'));
      const cpa = pa;
      const cpb = pb;
    if (true) {
      pa = Polynomial.toPolynomial(pa, RPN('a'));
      pb = Polynomial.toPolynomial(pb, RPN('a'));
      while (pa.getCoefficient(0).equals(Expression.ZERO)) {
        pa = pa.divideAndRemainder(Polynomial.of(Expression.ONE).shift(1)).quotient;//TODO: simplify
      }
      while (pb.getCoefficient(0).equals(Expression.ZERO)) {
        // a = 0, p(b*i) = 0
        var candidates = Polynomial.toPolynomial(pa.calcAt(Expression.ZERO), RPN('b')).getZeros(undefined, false);
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
              var zero = Polynomial.toPolynomial(condition.array[0].expression, RPN('b'));
              return p.map(function (coefficient) {
                var n = Polynomial.toPolynomial(coefficient.getNumerator(), RPN('b')).divideAndRemainder(zero).remainder.calcAt(RPN('b'));
                var d = Polynomial.toPolynomial(coefficient.getDenominator(), RPN('b')).divideAndRemainder(zero).remainder.calcAt(RPN('b'));
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
            var newp2 = p2.scale(lc.inverse());
            var r = p1.divideAndRemainder(newp2).remainder;
            walk(newp2, simplifyCoefficients(r, c2), c2);
          }
        } else {
          //TODO: ?
          condition = condition.andZero(p2.getLeadingCoefficient());
          if (!condition.isFalse()) {
            console.assert(condition.array.length === 1 && condition.array[0].operator === " == 0");
            const bPolynomial = Polynomial.toPolynomial(condition.array[0].expression, RPN('b'));
            //TODO: fix for higher degrees (?)
            let candidates = bPolynomial.getDegree() < 3 ? bPolynomial.getroots() : bPolynomial.getZeros(undefined, false);
            candidates = candidates.filter(c => Expression._isPositive(c instanceof ExpressionWithPolynomialRoot ? c.e : c));//!?
            for (const b of candidates) {
              const pp = p1.map(function (coefficient) { return Polynomial.toPolynomial(coefficient.getNumerator(), RPN('b')).calcAt(b).divide(Polynomial.toPolynomial(coefficient.getDenominator(), RPN('b')).calcAt(b)); });
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
        debugger;
        var bCandidates = Polynomial.toPolynomial(Polynomial.resultant(Polynomial.toPolynomial(cpa, RPN('a')), Polynomial.toPolynomial(cpb, RPN('a'))), RPN('b')).getZeros();
        var aCandidates = Polynomial.toPolynomial(Polynomial.resultant(Polynomial.toPolynomial(cpa, RPN('b')), Polynomial.toPolynomial(cpb, RPN('b'))), RPN('a')).getZeros();
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
        bCandidates = bCandidates.filter(c => Expression._isPositive(c instanceof ExpressionWithPolynomialRoot ? c.e : c));//!?
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
      result.sort((a, b) => (a instanceof ExpressionWithPolynomialRoot ? 1 : 0) - (b instanceof ExpressionWithPolynomialRoot ? 1 : 0));
      var strings = result.map(x => x.toMathML());
    result = result.filter(function (x, index) {
      for (var j = index - 1; j >= 0; j -= 1) {
        if (strings[j] === strings[index]) {
          if (result[j] instanceof ExpressionWithPolynomialRoot) {
            if (result[j].equals(result[index])) {
              return false;
            }
          } else if (result[index] instanceof ExpressionWithPolynomialRoot) {
            if (result[index].equals(result[j])) {
              return false;
            }
          } else if (result[j].equals(result[index])) {
            return false;
          }
        }
      }
      return true;
    });
    }

    //?
    //var multiplicities = new Array(result.length).fill(1);
    //_countMultiplicities(multiplicities, intervals, a0);

    return result;
  };

Polynomial.getSylvesterMatrix = function (A, B) {
  var d = A.getDegree();
  var e = B.getDegree();
  return Matrix.Zero(d + e, d + e).map(function (element, i, j) {
    return j < e ? (i - j < 0 || i - j > d ? Expression.ZERO : A.getCoefficient(i - j)) : (i + e - j < 0 || i + e - j > e ? Expression.ZERO : B.getCoefficient(i + e - j));
  });
};

Polynomial.resultant = function (A, B) {
  return Polynomial.getSylvesterMatrix(A, B).determinant();
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
    throw new RangeError();
  }
  return this.subtract(other).getNumerator().compareTo(Expression.ZERO);
};
Expression.prototype.round = function () {//TODO: remove - ?
  //TODO: half away from zero - ?
  //console.log(this.getNumerator(), this.getDenominator());
  //return this.getNumerator().add(this.getDenominator().truncatingDivide(Expression.TWO)).truncatingDivide(this.getDenominator());
  return RPN(toDecimalStringInternal(this, {fractionDigits: 0}));
};

//console.assert(GramSchmidt(new Matrix([[Expression.Integer.fromNumber(3), Expression.Integer.fromNumber(1)], [Expression.Integer.fromNumber(2), Expression.Integer.fromNumber(2)]])).toString() === '{{3,1},{-2/5,6/5}}');
//console.assert(GramSchmidt(new Matrix([[Expression.Integer.fromNumber(3), Expression.Integer.fromNumber(1)], [Expression.Integer.fromNumber(2), Expression.Integer.fromNumber(2)], [new Expression.Integer(0), new Expression.Integer(0)]])).toString() === '{{3,1},{-2/5,6/5},{0,0}}');
//throw new Error();

function BD(value, rounding) {
  if (!(value instanceof BigFloat)) {
    throw new RangeError();
  }
  this.value = value;
  this.rounding = rounding;
}
BD.prototype.add = function (other) {
  return new BD(BigFloat.add(this.value, other.value, this.rounding), this.rounding);
};
BD.prototype.subtract = function (other) {
  return new BD(BigFloat.subtract(this.value, other.value, this.rounding), this.rounding);
};
BD.prototype.multiply = function (other) {
  return new BD(BigFloat.multiply(this.value, other.value, this.rounding), this.rounding);
};
BD.prototype.divide = function (other) {
  return new BD(BigFloat.divide(this.value, other.value, this.rounding), this.rounding);
};
BD.prototype.simplifyExpression = function () {
  return this;
};
BD.prototype.equals = function (other) {
  if (other === Expression.ZERO) {
    return BigFloat.equal(this.value, BigFloat.round(BigFloat.BigFloat(0), this.rounding));
  }
  throw new TypeError(other);
};
BD.prototype.getNumerator = function () {
  return this;//?
};
BD.prototype.abs = function () {
  return new BD(BigFloat.lessThan(this.value, BigFloat.round(BigFloat.BigFloat(0), this.rounding)) ? BigFloat.unaryMinus(this.value) : this.value, this.rounding);
};
BD.prototype.compareTo = function (other) {
  return BigFloat.lessThan(this.value, other.value) ? -1 : (BigFloat.lessThan(other.value, this.value) ? +1 : 0);
};
BD.prototype.round = function () {
  return new BD(BigFloat.round(this.value, {maximumFractionDigits: 0, roundingMode: 'half-even'}), this.rounding);
};
BD.prototype._pow = function (n) {
  if (n !== 2) {
    throw new RangeError(n);
  }
  return this.multiply(this);
};

function LLL(latticeBasis) {
  //const δ = RPN('3/4');
  const δ = latticeBasis.e(0, 0).rounding == null ? RPN('3/4') : new BD(BigFloat.divide(BigFloat.BigFloat(3), BigFloat.BigFloat(4), latticeBasis.e(0, 0).rounding), latticeBasis.e(0, 0).rounding);
  const ONE_HALF = latticeBasis.e(0, 0).rounding == null ? RPN('1/2') : new BD(BigFloat.divide(BigFloat.BigFloat(1), BigFloat.BigFloat(2), latticeBasis.e(0, 0).rounding), latticeBasis.e(0, 0).rounding);
  // https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm#LLL_algorithm_pseudocode
  let B = latticeBasis;
  const n = B.rows() - 1;
  let Bstar = GramSchmidt(B); // B⃰ - not supported in IE 11
  let μ = (i, j) => B.row(i).dot(Bstar.row(j)).divide(Bstar.row(j).dot(Bstar.row(j)));
  let k = 1;
  while (k <= n) {
    for (let j = k - 1; j >= 0; j -= 1) {
      if (μ(k, j).abs().compareTo(ONE_HALF) > 0) {
        //B[k] = B[k].subtract(B[j].scale(Math.floor(μ(k, j))));
        var a = μ(k, j).round();
        //console.log(a);
        B = B.map(function (e, row, column) {
          return row === k ? B.e(k, column).subtract(a.multiply(B.e(j, column))) : e;
        });
        Bstar = GramSchmidt(B);
      }
    }
    if (Bstar.row(k).dot(Bstar.row(k)).compareTo((δ.subtract(μ(k, k - 1).abs()._pow(2))).multiply(Bstar.row(k - 1).dot(Bstar.row(k - 1)))) >= 0) {
      k = k + 1;
    } else {
      B = B.swapRows(k, k - 1);
      Bstar = GramSchmidt(B);
      k = Math.max(k - 1, 1);
    }
  }
  return B;
}
globalThis.LLL = LLL;//!
// https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm#Example_from_'%22%60UNIQ--postMath-00000048-QINU%60%22'
//console.assert(LLL(RPN('{{1,-1,3},{1,0,5},{1,2,6}}').matrix.transpose()).transpose().toString() === '{{0,1,-1},{1,0,0},{0,1,2}}');

// https://en.wikipedia.org/wiki/Lenstra–Lenstra–Lovász_lattice_basis_reduction_algorithm#Example_from_'"%60UNIQ--postMath-0000004C-QINU%60"'

//var M = '{{-2+2i,7+3i,7+3i,-5+4i},{3+3i,-2+4i,6+2i,-1+4i},{2+2i,-8+0i,-9+1i,-7+5i},{8+2i,-9+0i,6+3i,-4+4i}}'
//LLL(RPN(M).matrix.transpose()).transpose().toString()

// https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm#Applications

//console.assert(LLL(RPN('{{1,0,0,10000*(1.618034)**2},{0,1,0,10000*(1.618034)},{0,0,1,10000}}').matrix.map(e => e.round())).slice(0, 1, 0, 3).toString() === '{{-1,1,1}}'); // or '{{1,-1,-1}}'

Polynomial.prototype._log2hypot = function () {
  const polynomial = this;
  const coefficients = new Array(polynomial.getDegree() + 1).fill(0).map((x, k) => polynomial.getCoefficient(k));
  const maxBitLength = Math.max.apply(null, coefficients.map(c => c.equals(Expression.ZERO) ? 0 : bitLength(c.abs().toBigInt())));
  const scale = Expression.TWO._pow(maxBitLength);
  const hypot = Math.hypot.apply(null, coefficients.map(c => Number(toDecimalStringInternal(c.divide(scale), {fractionDigits: 15}))));
  const log2hypot = maxBitLength + Math.log2(hypot);
  return log2hypot;
};

globalThis.precision1 = 1;
function factorByLLL(polynomial, options) {
  options = options || {};
  // https://en.wikipedia.org/wiki/Factorization_of_polynomials#:~:text=LLL%20factorization
  var zeros = polynomial.getZeros();
  if (zeros.length === 0) {
    zeros = polynomial.getZeros(undefined, true);
  }
  if (zeros.length !== 0) {
    //var n = Math.ceil(polynomial.getDegree() / 2);//!!! should be polynomial.getDegree() - ?
    //n = Math.min(n, 8);//TODO: REMOVE!!!
    var log2hypot = polynomial._log2hypot();
    for (var n = (polynomial.isEven() && options.evenFlag ? polynomial.getDegree() / 2 : 2); n <= (polynomial.isEven() && options.evenFlag ? polynomial.getDegree() / 2 : polynomial.getDegree() - 2) && n <= 8; n += 1) {
      for (var index = zeros.length - 1; index >= Math.max(zeros.length - n, 0) && (!polynomial.isEven() || !options.evenFlag || index === zeros.length - 1); index -= 1) {
        var zero = zeros[index];
        if (zero instanceof ExpressionWithPolynomialRoot) {//TODO: remove
          zero = zero.e;
        }
        //TODO: what is the precision - ?
        // https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#:~:text=This%20bound%20is%20also%20useful%20to%20bound%20the%20coefficients%20of%20a%20divisor%20of%20a%20polynomial%20with%20integer%20coefficients:
        // https://oeis.org/search?q=1+2+6+20+70+252&sort=&language=&go=Search
        const precision = Math.ceil(((n * Math.log2(Math.sqrt(2 + Math.sqrt(3))) + log2hypot / (polynomial.isEven() && options.evenFlag ? 2 : 1) * n)) / Math.log2(2) * globalThis.precision1);
        //const precision = globalThis.precision1;
        // 52
        // 244
        console.log(n, n * Math.log2(Math.sqrt(2 + Math.sqrt(3))), log2hypot, precision);
        const rounding = {maximumFractionDigits: Math.ceil(precision * 2), roundingMode: 'floor'};//?
        const ZERO = new BD(BigFloat.round(BigFloat.BigFloat(0), rounding), rounding);
        var matrix = Matrix.Zero(n + 1, n + 2).map(function (e, i, j) {
          if (i === j) {
            return Expression.ONE;
          }
          if (j !== n + 1) {
            return Expression.ZERO;
          }
          return Expression.pow(Expression.TWO, precision).multiply(zero._pow(n - i)).round();
        }).map(function (e, i, j) {
          var tmp = e;
          if (!(tmp instanceof Expression.Integer)) {
            throw new TypeError();
          }
          if (tmp.value === Expression.ZERO) {
            return ZERO;
          }
          return new BD(BigFloat.round(BigFloat.BigFloat(tmp.value), rounding), rounding);
          //return tmp;
        });
        var reduced = LLL(matrix);
        reduced = reduced.map(e => RPN(e.value.toFixed(0)));
        var guess = Polynomial.of();
        for (var j = 0; j <= n; j += 1) {
          guess = guess.add(Polynomial.of(reduced.e(0, j)).shift(n - j));
        }
        while (guess.getCoefficient(0).equals(Expression.ZERO)) {
          guess = guess.divideAndRemainder(Polynomial.of(Expression.ZERO, Expression.ONE)).quotient;
        }
        if (guess.getDegree() > 0) {
          console.debug(guess.toString(), reduced.e(0, n + 1), reduced.e(0, n + 2), log2hypot, guess._log2hypot());
          if (polynomial.divideAndRemainder(guess).remainder.equals(Polynomial.ZERO)) {
            if (guess.getDegree() > polynomial.getDegree() / 2) {//? wrong choose of zero
              guess = polynomial.divideAndRemainder(guess).quotient;
            }
            return guess.getLeadingCoefficient().compareTo(Expression.ZERO) < 0 ? guess.negate() : guess;
          }
        }
      }
    }
  }
  return undefined;
}

Polynomial.prototype._factorByLLL = function (options) {
  return factorByLLL(this, options);
};
