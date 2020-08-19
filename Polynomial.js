/*global hit*/

  import Expression from './Expression.js';
  import Matrix from './Matrix.js';

  function PolynomialTerm(degree, coefficient) {
    if (degree < 0) {
      throw new RangeError();
    }
    if (degree > 9007199254740991) {
      throw new RangeError("NotSupportedError");
    }
    if (coefficient.equals(Expression.ZERO)) {
      throw new RangeError();
    }
    this.degree = degree;
    this.coefficient = coefficient;
  }
  function PolynomialData(length) {
    this.array = new Array(length);
    this.i = 0;
    this.length = length; // public
  }
  PolynomialData.prototype.add = function (degree, coefficient) {
    var k = this.i;
    if (k >= this.length) {
      throw new RangeError();
    }
    this.array[k] = new PolynomialTerm(degree, coefficient);
    this.i = k + 1;
  };
  PolynomialData.prototype.degree = function (i) {
    return this.array[i].degree;
  };
  PolynomialData.prototype.coefficient = function (i) {
    return this.array[i].coefficient;
  };
  PolynomialData.prototype.trim = function () {
    if (this.i !== this.length) {
      var array = new Array(this.i);
      for (var i = 0; i < this.i; i += 1) {
        array[i] = this.array[i];
      }
      this.array = array;
      this.length = this.i;
    }
    return this;
  };
  PolynomialData.prototype.getCoefficient = function (degree) {
    var from = 0;
    var to = this.array.length;
    while (from < to) {
      var middle = from + Math.floor((to - from) / 2);
      var y = this.array[middle].degree;
      if (y < degree) {
        from = middle + 1;
      } else if (y > degree) {
        to = middle;
      } else {
        return this.array[middle].coefficient;
      }
    }
    return Expression.ZERO;
  };

  // Polynomial([a0, a1, a2, ...., an]);
  // an*x^n+ an-1 x ^n-1 +... + a0
  function Polynomial(a) {
    this.a = a;
  }

  Polynomial.of = function () {
    var newData = new PolynomialData(arguments.length);
    for (var i = 0; i < arguments.length; i += 1) {
      var a = arguments[i];
      if (!a.equals(Expression.ZERO)) {
        newData.add(i, a);
      }
    }
    return new Polynomial(newData.trim());
  };
  Polynomial.from = function (array) {
    var newData = new PolynomialData(array.length);
    for (var i = 0; i < array.length; i += 1) {
      var a = array[i];
      if (!a.equals(Expression.ZERO)) {
        newData.add(i, a);
      }
    }
    return new Polynomial(newData.trim());
  };

  Polynomial.ZERO = Polynomial.of();

  Polynomial.prototype.getDegree = function () {
    return this.a.length === 0 ? -1 : this.a.degree(this.a.length - 1);
  };
  Polynomial.prototype.getCoefficient = function (degree) {
    if (degree > this.getDegree()) {
      throw new RangeError();
    }
    return this.a.getCoefficient(degree);
  };
  Polynomial.prototype.getLeadingCoefficient = function () {
    return this.a.length === 0 ? Expression.ZERO : this.a.coefficient(this.a.length - 1);
  };

  Polynomial.prototype.map = function (mapFunction) {//?
    var newData = new PolynomialData(this.a.length);
    for (var i = 0; i < this.a.length; i += 1) {
      var coefficient = this.a.coefficient(i);
      var degree = this.a.degree(i);
      var c = mapFunction(coefficient, degree);
      if (!c.equals(Expression.ZERO)) {
        newData.add(degree, c);
      }
    }
    return new Polynomial(newData.trim());
  };

  Polynomial.prototype.equals = function (p) {
    var i = this.a.length;
    if (i !== p.a.length) {
      return false;
    }
    while (--i >= 0) {
      if (this.a.degree(i) !== p.a.degree(i) || !this.a.coefficient(i).equals(p.a.coefficient(i))) {
        return false;
      }
    }
    return true;
  };

  Polynomial.prototype.add = function (p) {
    var i = 0;
    var j = 0;
    var newData = new PolynomialData(this.a.length + p.a.length);
    while (i < this.a.length && j < p.a.length) {
      var x = this.a.degree(i);
      var y = p.a.degree(j);
      if (x < y) {
        newData.add(x, this.a.coefficient(i));
        i += 1;
      } else if (x > y) {
        newData.add(y, p.a.coefficient(j));
        j += 1;
      } else {
        var c = this.a.coefficient(i).add(p.a.coefficient(j));
        if (!c.equals(Expression.ZERO)) {
          newData.add(x, c);
        }
        i += 1;
        j += 1;
      }
    }
    while (i < this.a.length) {
      newData.add(this.a.degree(i), this.a.coefficient(i));
      i += 1;
    }
    while (j < p.a.length) {
      newData.add(p.a.degree(j), p.a.coefficient(j));
      j += 1;
    }
    return new Polynomial(newData.trim());
  };

  Polynomial.prototype.multiply = function (p) {
    if (this.a.length === 0 || p.a.length === 0) {
      return Polynomial.ZERO;
    }
    var result = Polynomial.ZERO;
    for (var i = 0; i < this.a.length; i += 1) {
      var xd = this.a.degree(i);
      var xc = this.a.coefficient(i);
      var newData = new PolynomialData(p.a.length);
      for (var j = 0; j < p.a.length; j += 1) {
        var yd = p.a.degree(j);
        var yc = p.a.coefficient(j);
        newData.add(xd + yd, xc.multiply(yc));
      }
      result = result.add(new Polynomial(newData.trim()));
    }
    return result;
  };

  Polynomial.prototype.shift = function (n) { // *= x**n, n >= 0
    if (n < 0) {
      throw new TypeError();
    }
    var newData = new PolynomialData(this.a.length);
    for (var i = 0; i < this.a.length; i += 1) {
      newData.add(this.a.degree(i) + n, this.a.coefficient(i));
    }
    return new Polynomial(newData.trim());
  };

  Polynomial.prototype.divideAndRemainder = function (p, w) {
    w = w || undefined;
    if (p.equals(Polynomial.ZERO)) {
      throw new TypeError("ArithmeticException");
    }
    var quotient = Polynomial.ZERO;
    var remainder = this;
    while (remainder.getDegree() >= p.getDegree()) {
      var n = remainder.getDegree() - p.getDegree();
      var lcr = remainder.getLeadingCoefficient();
      var lcp = p.getLeadingCoefficient();
      var q = lcr.divide(lcp);
      if (w != undefined) {
        if (q instanceof Expression.Division) {
          if (w === "throw") {
            throw new RangeError(); // AssertionError
          }
          if (w === "undefined") {
            return undefined;
          }
          throw new RangeError();
        }
      }
      var pq = Polynomial.of(q);
      quotient = quotient.add(pq.shift(n));
      //TODO: optimize - ?
      remainder = remainder.subtract(pq.multiply(p).shift(n));
      if (remainder.getDegree() - p.getDegree() === n) {
        // to avoid the infite loop
        throw new TypeError("there is a some problem with the expression evaluation");//!
      }
    }
    return {quotient: quotient, remainder: remainder};
  };

  Polynomial.pseudoRemainder = function (x, y) {
    var lcg = y.getLeadingCoefficient();
    var n = x.getDegree() - y.getDegree();
    // assertion
    if (n < 0) {
      throw new RangeError();
    }
    var sx = x.multiply(Polynomial.of(Expression.pow(lcg, n).multiply(lcg)));
    return sx.divideAndRemainder(y, "throw").remainder;
  };

  Polynomial.polynomialGCD = function (a, b) {
    //TODO: fix (place condition for degrees earlier - ?)
    if (a.getDegree() < b.getDegree()) {
      //!!!
      var tmp = a;
      a = b;
      b = tmp;
    }

    var contentA = a.getContent();
    var contentB = b.getContent();
    var ppA = a.divideAndRemainder(Polynomial.of(contentA), "throw").quotient;
    var ppB = b.divideAndRemainder(Polynomial.of(contentB), "throw").quotient;
    var A = ppA;
    var B = ppB;
    while (!B.equals(Polynomial.ZERO)) {
      var r = Polynomial.pseudoRemainder(A, B);
      //! 2018-05-12
      if (!r.equals(Polynomial.ZERO)) {
        // https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Primitive_pseudo-remainder_sequence
        // ExpressionParser.parse("(x^8+x^6-3x^4-3x^3+8x^2+2x-5)/(3x^6+5x^4-4x^2-9x+21)")
        r = r.divideAndRemainder(Polynomial.of(r.getContent()), "throw").quotient;
        //console.log(r.toString());
      }
      //!
      A = B;
      B = r;
    }
    var c = contentA.gcd(contentB);
    return Polynomial.of(c).multiply(A.divideAndRemainder(Polynomial.of(A.getContent()), "throw").quotient);
  };

  Polynomial.prototype.calcAt = function (point) {//!!!
    var n = Expression.ZERO;
    var lastDegree = -1;
    var i = this.a.length;
    while (--i >= -1) {
      var degree = i === -1 ? 0 : this.a.degree(i);
      var coefficient = i === -1 ? Expression.ZERO : this.a.coefficient(i);
      if (i !== this.a.length - 1) {
        n = n.multiply(Expression.pow(point, lastDegree - degree));
      }
      if (i !== -1) {
        n = n.add(coefficient);
      }
      lastDegree = degree;
    }
    return n;
  };

  Polynomial.prototype.getContent = function () {
    if (this.a.length === 0) {
      throw new RangeError();
    }
    var x = this.a.coefficient(this.a.length - 1);
    var denominator = x.getDenominator();
    var numerator = x.getNumerator();
    var i = this.a.length - 1;
    while (--i >= 0) {
      var y = this.a.coefficient(i);
      denominator = denominator.lcm(y.getDenominator());
      numerator = numerator.gcd(y.getNumerator());
    }
    var c = numerator.divide(denominator);
    return x.isNegative() && !numerator.isNegative() ? c.negate() : c;
  };

  // add, multiply, divideAndRemainder

  Polynomial.prototype.negate = function () {
    //TODO: fix
    return this.map(function (coefficient, degree) {
      return coefficient.negate();
    });
  };

  Polynomial.prototype.subtract = function (l) {
    return this.add(l.negate());
  };

  Polynomial.prototype.scale = function (x) {
    return this.map(function (coefficient, degree) {
      return coefficient.multiply(x);
    });
  };

  Polynomial.toPolynomial = function (e, v) {
    if (e instanceof Expression.Division) {
      throw new RangeError();
    }
    var coefficients = Expression.getCoefficients(e, v);
    var newData = new PolynomialData(coefficients.length);
    for (var i = 0; i < coefficients.length; i += 1) {
      var x = coefficients[i];
      var d = x.degree.toNumber();
      var c = x.coefficient;
      newData.add(d, c);
    }
    return new Polynomial(newData.trim());
  };

  Polynomial.prototype.toExpression = function (variableSymbol) {
    var i = this.a.length;
    var result = undefined;
    while (--i >= 0) {
      var degree = this.a.degree(i);
      var coefficient = this.a.coefficient(i);
      var v = degree === 0 ? undefined : (degree === 1 ? variableSymbol : new Expression.Exponentiation(variableSymbol, Expression.Integer.fromNumber(degree)));
      var current = v == undefined ? coefficient : (coefficient.equals(Expression.ONE) ? v : new Expression.Multiplication(coefficient, v));
      result = result == undefined ? current : new Expression.Addition(result, current);
    }
    return result == undefined ? Expression.ZERO : result;
  };

  // return a first founded root to simplify and as the next call may be called with reduced coefficients
  Polynomial.prototype.doRationalRootTest = function () {
    var np = this;
    var an = np.getLeadingCoefficient();
    var a0 = np.getCoefficient(0);
    a0 = Expression._expandTrigonometry(a0);//!

    //TODO: http://en.wikipedia.org/wiki/Polynomial_remainder_theorem

    // http://scask.ru/g_book_mav.php?id=26
    var hasIntegerCoefficients = np._hasIntegerCoefficients();
    // f(k) = -q(a)(a - k)
    var filter = [];
    if (hasIntegerCoefficients) {
      if (np.calcAt(Expression.ONE).equals(Expression.ZERO)) {
        return Expression.ONE;
      }
      filter.push({k: Expression.ONE, fk: np.calcAt(Expression.ONE)});
      if (np.calcAt(Expression.ONE.negate()).equals(Expression.ZERO)) {
        return Expression.ONE.negate();
      }
      filter.push({k: Expression.ONE.negate(), fk: np.calcAt(Expression.ONE.negate())});
    }

    //!new 2020-01-13
    if (hasIntegerCoefficients) {
      var zeros = np.getZeros(1).result;
      for (var i = 0; i < zeros.length; i += 1) {
        var zero = zeros[i];
        var candidate = Expression.Integer.fromString(zero.multiply(an).toString({rounding: {fractionDigits: 0}})).divide(an);
        if (np.calcAt(candidate).equals(Expression.ZERO)) {
          return candidate;
        }
      }
      return null;
    }

    /*
    TODO:
    k = k
    f(k) = fk
    t = x + k
    x = t - k
    f(t) = f(x + k) = an * x**n + ... + f(k)
    */


    var result = null;
    // p/q
    //TODO: forEach -> some ?
    Expression.everyDivisor(a0, function (p) {
      return Expression.everyDivisor(an, function (q) {
        var sign = -3;
        while ((sign += 2) < 3) {
          var sp = sign === -1 ? p.negate() : p;

          var i = 0;
          while (i < filter.length && !filter[i].k.multiply(sp).subtract(q).equals(Expression.ZERO) && filter[i].fk.remainder(filter[i].k.multiply(sp).subtract(q)).equals(Expression.ZERO)) {
            i += 1;
          }
          var filteredOut = i < filter.length;

          if (//sp.gcd(q).equals(Expression.ONE) &&
              !filteredOut) {//?
            var x = Polynomial.of(sp.negate(), q);
            var z = np.divideAndRemainder(x, "undefined");
            var r = z == undefined ? undefined : z.remainder.map(function (x) {
              return x.simplifyExpression();
            });
            if (r != undefined && r.equals(Polynomial.ZERO)) {
              result = sp.divide(q);
              return false;
            }
          }
        }
        return true;
      });
    });

    return result;
  };

  Polynomial.prototype._hasIntegerCoefficients = function () {
    for (var i = 0; i <= this.getDegree(); i += 1) {
      if (!(this.getCoefficient(i) instanceof Expression.Integer)) {
        return false;
      }
    }
    return true;
  };

  Polynomial.prototype._hasIntegerLikeCoefficients = function () {
    var isIntegerLike = function (c) {
      if (c instanceof Expression.Integer) {
        return true;
      }
      if (c instanceof Expression.Symbol) {
        return true;
      }
      if (c instanceof Expression.Addition) {
        return isIntegerLike(c.a) && isIntegerLike(c.b);
      }
      if (c instanceof Expression.Multiplication) {
        return isIntegerLike(c.a) && isIntegerLike(c.b);
      }
      if (c instanceof Expression.Exponentiation) {
        return isIntegerLike(c.a) && c.b instanceof Expression.Integer;
      }
      return false;
    };
    var integerLikeCoefficients = true;
    for (var i = 0; i <= this.getDegree(); i += 1) {
      integerLikeCoefficients = integerLikeCoefficients && isIntegerLike(this.getCoefficient(i));
    }
    return integerLikeCoefficients;
  };

  Polynomial.prototype._getFactorByKroneckersMethod = function () {
    // https://ru.wikipedia.org/wiki/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%9A%D1%80%D0%BE%D0%BD%D0%B5%D0%BA%D0%B5%D1%80%D0%B0
    // https://ru.wikipedia.org/wiki/%D0%98%D0%BD%D1%82%D0%B5%D1%80%D0%BF%D0%BE%D0%BB%D1%8F%D1%86%D0%B8%D0%BE%D0%BD%D0%BD%D1%8B%D0%B9_%D0%BC%D0%BD%D0%BE%D0%B3%D0%BE%D1%87%D0%BB%D0%B5%D0%BD_%D0%9B%D0%B0%D0%B3%D1%80%D0%B0%D0%BD%D0%B6%D0%B0
    // https://en.wikipedia.org/wiki/Vandermonde_matrix
    var np = this;
    if (!np._hasIntegerLikeCoefficients()) {
      return undefined;
    }
    var n = np.getDegree();
    var ys = new Array(n);
    var total = 1;
    for (var i = 0; i <= Math.floor(n / 2); i += 1) {
      var bi = Expression.Integer.fromNumber(i);
      var y = np.calcAt(bi);
      if (y.equals(Expression.ZERO)) {
        return Polynomial.of(bi.negate(), Expression.ONE);
      }
      var attachNegative = function (array) {
        var result = new Array(array.length * 2);
        for (var i = 0; i < array.length; i += 1) {
          result[i * 2] = array[i];
          result[i * 2 + 1] = array[i].negate();
        }
        return result;
      };
      var divisors = [];
      Expression.everyDivisor(y, function (d) {
        divisors.push(d);
        return true;
      });
      // let the first be positive as two sets with different signs give polynomials that differ only in sign of coefficients
      ys[i] = i === 0 ? divisors : attachNegative(divisors);
      total *= ys[i].length;
      var V = Matrix.Zero(i + 1, i + 1).map(function (e, i, j) {
        return Expression.pow(Expression.Integer.fromNumber(i), j);
      });
      var inv = V.inverse();
      //scale?
      inv = inv.scale(V.determinant());
      //?
      var u = new Array(i + 1);
      for (var j = 0; j < i + 1; j += 1) {
        u[j] = 0;
      }
      u[0] = -1;
      for (var j = 0; j < total; j += 1) {
        var k = 0;
        u[k] += 1;
        while (u[k] === ys[k].length) {
          u[k] = 0;
          k += 1;
          u[k] += 1;
        }
        var y = Matrix.Zero(i + 1, 1).map(function (e, i, j) {
          return ys[i][u[i]];
        });
        var s = inv.multiply(y);
        var polynomialFromVector = function (s) {
          var c = new Array(s.rows());
          for (var j = 0; j < s.rows(); j += 1) {
            c[j] = s.e(j, 0);
          }
          return Polynomial.from(c);
        };
        var g = polynomialFromVector(s);
        //if (g.getDegree() > 0 && np.divideAndRemainder(g).remainder.equals(Polynomial.ZERO)) {
        //  return g;
        //}
        if (g.getDegree() > 0) {
          var gc = g.getContent();
          g = g.scale(gc.getDenominator()).divideAndRemainder(Polynomial.of(gc.getNumerator()), "throw").quotient;
          var t = np.divideAndRemainder(g, "undefined");
          if (t != undefined && t.remainder.equals(Polynomial.ZERO)) {
            return g;
          }
        }
      }
    }
    return undefined;
  };

  var counter = 0;//TODO: remove

  Polynomial.prototype.getroots = function (callback) {
    //TODO: merge hit and callback
    callback = callback || undefined;
    var np = this;

    var roots = [];


    //!new 2018-12-24
    //TODO: fix (?Polynomial#getContent()?)
      var ct = Expression.ONE;
      var t = Expression.ZERO;
      while (t != null) {
        var t = Expression.getConjugate(np.getCoefficient(np.getDegree()));
        if (t != undefined) {
          np = np.scale(t);
          ct = ct.multiply(t);
        }
      }
    //!

    //!new 2020-07-11
    np = np.map(function (x) {
      return x.simplifyExpression();
    });
    //!

    var content = np.getContent();
    if (!content.equals(Expression.ONE)) {
      np = np.scale(content.getDenominator()).divideAndRemainder(Polynomial.of(content.getNumerator()), "throw").quotient;
      //np = np.divideAndRemainder(Polynomial.of(content), "throw").quotient;
    }
    if (!ct.equals(Expression.ONE)) {
      content = content.divide(ct);
    }



    // x = 0
    while (np.getCoefficient(0).equals(Expression.ZERO)) {
      np = np.divideAndRemainder(Polynomial.of(Expression.ZERO, Expression.ONE), "throw").quotient;
      roots.push(Expression.ZERO);
    }
    if ((!content.equals(Expression.ONE) && !content.equals(Expression.ONE.negate())) || roots.length > 0) {
      if (roots.length > 0) {
        if (typeof hit === "function") {
          hit({getroots: {special: "0"}});
        }
      }
      if (callback != undefined) {
        callback({content: content, roots: roots, newPolynomial: np, type: "factorOutTheGreatestCommonFactor"});
      }
    }

    if (np.getDegree() === 1) {
      roots.push(np.getCoefficient(0).negate().divide(np.getCoefficient(1)));
      np = Polynomial.of(np.getLeadingCoefficient());
      if (typeof hit === "function") {
        hit({getroots: {linear: ""}});
      }
      if (callback != undefined) {
        callback({content: content, roots: roots, newPolynomial: np, type: "solveLinearEquation"});
      }
      return roots;
    }

    var nthRootInternal = function (n, x) {
      if (x instanceof Expression.Division) {
        var sa1 = nthRootInternal(n, x.a);
        var sb1 = nthRootInternal(n, x.b);
        return sa1 == undefined || sb1 == undefined ? undefined : sa1.divide(sb1);
      }
      if (x instanceof Expression.Exponentiation) {
        var N = Expression.Integer.fromNumber(n);
        if (x.b instanceof Expression.Integer) {
          if (x.b.remainder(N).equals(Expression.ZERO)) {
            return x.a.pow(x.b.divide(N));
          }
          //return undefined;
        }
        if (x.a instanceof Expression.Integer || x.a === Expression.E) {//?
          return x.a.pow(x.b.divide(N));
        }
        if (x.b instanceof Expression.Division && x.b.a instanceof Expression.Integer && x.b.a.remainder(N).equals(Expression.ZERO)) {//TODO:
          return x.a.pow(x.b.divide(N));
        }
      }
      if (x instanceof Expression.Multiplication) {
        var sa = nthRootInternal(n, x.a);
        var sb = nthRootInternal(n, x.b);
        return sa == undefined || sb == undefined ? undefined : sa.multiply(sb);
      }
      if (x instanceof Expression.Complex || (Expression.isConstant(x) && Expression.has(x, Expression.Complex))) {
        //TODO: - ?
        //var real = x.real;
        //var imaginary = x.imaginary;
        var c = Expression.getComplexConjugate(x);
        var real = x.add(c).divide(Expression.TWO);
        var imaginary = x.subtract(c).multiply(Expression.I.negate()).divide(Expression.TWO);
        if (n === 2) {
          var m = real.multiply(real).add(imaginary.multiply(imaginary)).squareRoot();
          var a = nthRootInternal(2, real.add(m).divide(Expression.TWO));
          if (a != undefined) {
            var b = imaginary.divide(Expression.TWO.multiply(a));
            var result = a.add(b.multiply(Expression.I));
            return result;
          }
        }
        if (real.equals(Expression.ZERO) && n % 2 === 0) {
          var c = nthRootInternal(Math.floor(n / 2), x);
          if (c != undefined) {
            return nthRootInternal(2, c);
          }
        }
        if (real.equals(Expression.ZERO) && n % 2 === 1) {
          //?
          var c = nthRootInternal(n, imaginary);
          if (c != undefined) {
            return c.multiply((n % 4 === 1 ? Expression.I : Expression.I.negate()));
          }
        }
        //!new 2020-07-24
        if (x instanceof Expression.Complex && !imaginary.equals(Expression.ZERO)) {//?TODO: ?
          // https://en.wikipedia.org/wiki/Complex_number#Modulus_and_argument
          var rho = real._pow(2).add(imaginary._pow(2)).squareRoot();
          try {
            var phi = Expression.TWO.multiply(imaginary.divide(rho.add(real)).arctan());
            return rho._nthRoot(n).multiply(Expression.I.multiply(phi.divide(Expression.Integer.fromNumber(n))).exp());
          } catch (error) {
            //TODO: ?
            console.debug(error);
          }
        }
      }
      if (x instanceof Expression.Addition) {
        var lastFactor = undefined;
        var e = 0;
        var result = Expression.ONE;
        var rest = Expression.ONE;
        var t = x;
        while (!t.equals(Expression.ONE) && !t.equals(Expression.ONE.negate())) {
          var f = Expression.simpleDivisor(t);
          if (e === 0) {
            lastFactor = f;
            e += 1;
          } else if (f.equals(lastFactor)) {
            e += 1;
            if (e === n) {
              e = 0;
              result = result.multiply(lastFactor);
            }
          } else if (e !== 0) {
            rest = rest.multiply(Expression.pow(lastFactor, e));
            lastFactor = f;
            e = 1;
          }
          t = t.divide(f);
        }
        if (result !== Expression.ONE) {
          if (e !== 0) {
            rest = rest.multiply(Expression.pow(lastFactor, e));
          }
          if (t.equals(Expression.ONE.negate())) {
            rest = rest.multiply(t);
          }
          var rn = nthRootInternal(n, rest);
          if (rn != undefined) {
            return result.multiply(rn);
          }
        }
      }
      if (x instanceof Expression.Exponentiation && x.a instanceof Expression.Symbol) {
        var b = x.b.divide(Expression.Integer.fromNumber(n));
        return b.equals(Expression.ONE) ? x.a : new Expression.Exponentiation(x.a, b);
      }
      if (!Expression.isConstant(x) && x.isNegative() && (n === 2 || n % 2 !== 0)) {
        x = x.negate();
        var c = nthRootInternal(n, x);
        return c == null ? null : Expression.ONE.negate()._nthRoot(n).multiply(c);
      }
      if ((x instanceof Expression.Integer || x instanceof Expression.Complex) && x.isNegative() && n % 2 === 0) {//?
        var c = x instanceof Expression.Integer ? x._nthRoot(2) : nthRootInternal(2, x);
        return c == null ? null : nthRootInternal(n / 2, c);
      }
      if (Expression.has(x, Expression.Sin) || Expression.has(x, Expression.Cos)) {
        //?
        var tmp = nthRootInternal(2, Expression._replaceSinCos(x));
        if (tmp != null) {
          return Expression._replaceBySinCos(tmp).simplifyExpression();//?
        }
      }
      var y = undefined;
      try {
        y = x._nthRoot(n);
      } catch (error) {
        //TODO:
        console.error(error);
      }
      if (y == undefined) {//?
        var a = x;
        var ac = Expression.getConjugate(a.getNumerator());
        if ((n === 3 || n === 2) && ac != undefined && ac.divide(a.getDenominator()).multiply(a) instanceof Expression.Integer) {//TODO: ?
          //TODO: a > 0 - ?
          var a = x;
          var tmp = new Expression.Symbol('x')._pow(n).subtract(a);
          var polynomial = Polynomial.toPolynomial(Expression.getConjugate(tmp.getNumerator()).divide(tmp.getDenominator()).multiply(tmp), new Expression.Symbol('x'));
          var tmp2 = polynomial.getZeros(1);//TODO: ?
          return tmp2.result[1];//TODO: which one - ?
        }
      }
      return y;
    };

    var nthRoot = function (n, x, np) {
      if (n === 1) {
        return x;
      }
      var y = nthRootInternal(n, x);
      if (y == undefined) {
        if (!(x instanceof Expression.Integer)) {
          if (typeof hit === "function") {
            hit({nthRoot: (n === 2 ? "squareRoot" : (n === 3 ? "cubeRoot" : n + "-root")) + ":" + x.toString() + ":" + np.toString()});
          }
        }
      }
      return y;
    };

    var continueWithNewPolynomial = function (roots, np) {
      var rs = np.getroots(callback != undefined ? function (info) {
        var xxx = Object.assign({}, info, {content: content.multiply(info.content), roots: roots.concat(info.roots)});
        callback(xxx);
      } : undefined);
      for (var i = 0; i < rs.length; i += 1) {
        roots.push(rs[i]);
      }
    };

    if (np.getDegree() >= 2) {
      var gcd = function (a, b) {
        return b === 0 ? a : gcd(b, a % b);
      };
      var g = np.getDegree();
      for (var i = 1; i <= np.getDegree(); i += 1) {
        if (!np.getCoefficient(i).equals(Expression.ZERO)) {
          g = gcd(g, i);
        }
      }
      if (g >= 2) {
        var allZeros = g === np.getDegree();
        if (typeof hit === "function") {
          if (g === np.getDegree()) {
            hit({getroots: {allZeros: ""}});
          } else {
            hit({getroots: g % 2 === 0 ? (np.getDegree() === 4 ? {biquadratic: ""} : {even: ""}) : {xyz: np.toString()}});
          }
        }
        // t = x^g
        var newData = new Array(Math.floor((np.getDegree() + g) / g));
        var k = 0;
        for (var i = 0; i <= np.getDegree(); i += g) {
          newData[k] = np.getCoefficient(i);
          k += 1;
        }
        var q = Polynomial.from(newData);
        var qRoots = q.getroots();
        var n = np.getDegree();//TODO: 2018-02-04
        //var ok = false;//?
        for (var k = 0; k < qRoots.length; k += 1) {
          var qRoot = qRoots[k];
          var s = nthRoot(g, qRoot, np);
          if (s != undefined) {
            var d = null;
            if ((!allZeros || g > 4) && g <= 24 && ((17896830 >> g) & 1) === 1) {
              d = Polynomial.of(Expression.ONE).shift(g).add(Polynomial.of(qRoot.negate()));
              var c = Expression.E.pow(Expression.I.multiply(Expression.TWO.multiply(Expression.PI)).divide(Expression.Integer.fromNumber(g)));
              for (var i = 0; i < g; i += 1) {
                var root = Expression.pow(c, i).multiply(s);
                roots.push(root);
              }
            } else {
              roots.push(s);
              d = Polynomial.of(s.negate(), Expression.ONE);
              if (g % 2 === 0) {
                roots.push(s.negate());
                d = Polynomial.of(s.multiply(s).negate(), Expression.ZERO, Expression.ONE);
              }
            }
            var quotient = np.divideAndRemainder(d).quotient;
            console.assert(np.subtract(quotient.multiply(d)).getDegree() < 0);
            np = quotient;
          }
          //ok = ok || Expression.has(qRoot, Expression.Complex);//?
        }
        if (callback != undefined) {
          var type = allZeros ? (g === 2 ? "applyDifferenceOfSquaresRule" : (g === 3 ? "applyDifferenceOfCubesRule" : "applyDifferenceOfNthPowersRule")) : "t = x^g";
          callback({content: content, roots: roots, newPolynomial: np, type: type, g: g});//TODO: ?
        }
        var ok = true;
        if (n !== np.getDegree() && ok && np.getDegree() > 0) {
          continueWithNewPolynomial(roots, np);
        }
        return roots;
      }
    }

    //! new: solution of quadratic equations
    if (np.getDegree() === 2) {
      var a = np.getCoefficient(2);
      var b = np.getCoefficient(1);
      var c = np.getCoefficient(0);

      var D = b.multiply(b).subtract(Expression.TWO.multiply(Expression.TWO).multiply(a).multiply(c));
      D = D.simplifyExpression();
      var sD = nthRoot(2, D, np);
      if (typeof hit === "function") {
        hit({getroots: {quadratic: (sD == undefined ? (D instanceof Expression.Integer ? D.compareTo(Expression.ZERO) : "?" + D.toString()) : "OK")}});
      }
      if (sD != undefined) {
        var x1 = b.negate().subtract(sD).divide(Expression.TWO.multiply(a));
        var x2 = b.negate().add(sD).divide(Expression.TWO.multiply(a));
        roots.push(x1);
        roots.push(x2);
        np = Polynomial.of(np.getLeadingCoefficient());
        if (callback != undefined) {
          callback({content: content, roots: roots, newPolynomial: np, type: "solveQuadraticEquation"});
        }
        return roots;
      }
    }

    //TODO: odd degrees ?
    if (np.getDegree() >= 4 && np.getDegree() % 2 === 0) {
      var middle = Math.floor(np.getDegree() / 2);
      var j = 1;
      while (j < middle + 1 && np.getCoefficient(middle + j).equals(Expression.ZERO) && np.getCoefficient(middle - j).equals(Expression.ZERO)) {
        j += 1;
      }
      if (j < middle + 1 && !np.getCoefficient(middle + j).equals(Expression.ZERO) && !np.getCoefficient(middle - j).equals(Expression.ZERO)) {
        var jj = Expression.Integer.fromNumber(j);
        var mj = np.getCoefficient(middle + j).divide(np.getCoefficient(middle - j));
        var isQuasiPalindromic = true;
        for (var i = 2; i < middle + 1; i += 1) {
          isQuasiPalindromic = isQuasiPalindromic && np.getCoefficient(middle + i).pow(jj).subtract(np.getCoefficient(middle - i).pow(jj).multiply(mj.pow(Expression.Integer.fromNumber(i)))).equals(Expression.ZERO);
        }
        if (isQuasiPalindromic) {
          //TODO: fix
          if (typeof hit === "function") {
            hit({getroots: {quasiPalindromic: np.getDegree()}});
          }
        }
        if (isQuasiPalindromic && np.getDegree() <= 53) { // log2(9007199254740991 + 1)
          var substitute = function (m, np) { // t = mx + 1 / x
            // https://stackoverflow.com/a/15302448/839199
            var choose = function (n, k) {
              return k === 0 ? 1 : Math.floor((n * choose(n - 1, k - 1)) / k);
            };
            var p = function (n, i, mpi) {
              return n - 2 * i >= 0 ? p(n - 2 * i, 1, m).scale(Expression.Integer.fromNumber(choose(n, i)).multiply(mpi).negate()).add(p(n, i + 1, mpi.multiply(m))) : Polynomial.of(Expression.ONE).shift(n);
            };
            var f = function (n, i) {
              return i <= n ? p(n - i, 1, m).scale(np.getCoefficient(i)).add(f(n, i + 1)) : Polynomial.ZERO;
            };
            return f(Math.floor(np.getDegree() / 2), 0);
          };
          var m = j === 1 ? mj : nthRoot(j, mj, np); // TODO: check the result of nthRoot - ?
          // t = mx + 1 / x
          var pt = substitute(m, np);
          //var pt = Polynomial.of(np.getCoefficient(2).subtract(Expression.ONE.add(Expression.ONE).multiply(m).multiply(np.getCoefficient(0))), np.getCoefficient(1), np.getCoefficient(0));
          var ptRoots = pt.getroots();
          for (var i = 0; i < ptRoots.length; i += 1) {
            var ptRoot = ptRoots[i];
            // mx^2 - tx + 1 = 0
            var u = Polynomial.of(Expression.ONE, ptRoot.negate(), m);
            var uRoots = u.getroots();
            for (var j = 0; j < uRoots.length; j += 1) {
              var root = uRoots[j];
              np = np.divideAndRemainder(Polynomial.of(root.negate(), Expression.ONE)).quotient;
              roots.push(root);
            }
          }
          if (callback != undefined) {
            callback({content: content, roots: roots, newPolynomial: np, type: "solvePalindromicEquaion"});
          }
          return roots;
        }
      }
    }

    if (np.getDegree() >= 2) {//?
      // (ax+b)**n = a**n*x**n + n*a**(n-1)*x**(n-1)*b + ...
      //?
      // a**n
      // n*a**(n-1)*b
      var n = np.getDegree();
      var hasZeroCoefficients = function (np) {
        for (var i = 0; i <= np.getDegree(); i += 1) {
          if (np.getCoefficient(i).equals(Expression.ZERO)) {
            return true;
          }
        }
        return false;
      };
      if (!hasZeroCoefficients(np)) {
        var g = np.getCoefficient(n - 1).divide(np.getCoefficient(n)).divide(Expression.Integer.fromNumber(n));
        var ok = true;
        for (var k = np.getDegree() - 1; k >= 1 && ok; k -= 1) {
          ok = g.equals(np.getCoefficient(k - 1).divide(np.getCoefficient(k)).multiply(Expression.Integer.fromNumber(n - k + 1)).divide(Expression.Integer.fromNumber(k)));
        }
        if (ok) {
          var root = g.negate();
          for (var k = 0; k < n; k += 1) {
            roots.push(root);
          }
          np = Polynomial.of(np.getLeadingCoefficient());
          if (callback != undefined) {
            callback({content: content, roots: roots, newPolynomial: np, type: "(ax+b)**n"});//TODO:
          }
          return roots;
        }
      }
    }

    if (np.getDegree() >= 2) {
      var root = np.doRationalRootTest();
      if (root != null) {
        //np = np.divideAndRemainder(Polynomial.of(root.getNumerator().negate(), root.getDenominator())).quotient;
        np = np.divideAndRemainder(Polynomial.of(root.negate(), Expression.ONE)).quotient;
        roots.push(root);
        if (typeof hit === "function") {
          hit({getroots: {rational: ""}});
        }
        if (callback != undefined) {
          callback({content: content, roots: roots, newPolynomial: np, type: "useTheRationalRootTest"});
        }
        if (np.getDegree() > 0) {
          continueWithNewPolynomial(roots, np);
        }
        return roots;
      }
    }

    if (np.getDegree() === 3) {
      // https://en.wikipedia.org/wiki/Cubic_function#Algebraic_solution
      var a = np.getCoefficient(3);
      var b = np.getCoefficient(2);
      var c = np.getCoefficient(1);
      var d = np.getCoefficient(0);
      var TWO = Expression.TWO;
      var THREE = TWO.add(Expression.ONE);
      var FOUR = TWO.add(TWO);
      var NINE = THREE.multiply(THREE);
      var EIGHTEEN = NINE.multiply(TWO);
      var TWENTY_SEVEN = NINE.multiply(THREE);
      // 18*a*b*c*d-4*b^3*d+b^2*c^2-4*a*c^3-27*a^2*d^2
      var delta = EIGHTEEN.multiply(a).multiply(b).multiply(c).multiply(d)
                  .subtract(FOUR.multiply(b.pow(THREE)).multiply(d))
                  .add(b.pow(TWO).multiply(c.pow(TWO)))
                  .subtract(FOUR.multiply(a).multiply(c.pow(THREE)))
                  .subtract(TWENTY_SEVEN.multiply(a.pow(TWO)).multiply(d.pow(TWO)));
      // b^2-3*a*c
      var delta0 = b.pow(TWO).subtract(THREE.multiply(a).multiply(c));
      if (typeof hit === "function") {
        hit({getroots: {cubic: (delta instanceof Expression.Integer ? delta.compareTo(Expression.ZERO) : "?") + "-" + (delta0 instanceof Expression.Integer ? delta0.compareTo(Expression.ZERO) : "?")}});
      }
      if (delta.equals(Expression.ZERO)) {
        if (delta0.equals(Expression.ZERO)) {
          //TODO: link to a^3+3a^2b+3ab^2+b^3=0 - ?
          // -b/(3*a)
          var root = b.negate().divide(THREE.multiply(a));
          roots.push(root);
          roots.push(root);
          roots.push(root);
          np = Polynomial.of(np.getLeadingCoefficient());
          if (callback != undefined) {
            callback({content: content, roots: roots, newPolynomial: np, type: "solveCubicEquation"});
          }
          return roots;
        } else {
          // https://github.com/nicolewhite/algebra.js/blob/master/src/equations.js
          // https://en.wikipedia.org/wiki/Cubic_function#Multiple_roots.2C_.CE.94_.3D_0
          // a double root
          // 9*a*d-b*c
          var root = NINE.multiply(a).multiply(d).subtract(b.multiply(c)).divide(TWO.multiply(delta0));
          roots.push(root);
          roots.push(root);
          var p = Polynomial.of(root.negate(), Expression.ONE);
          np = np.divideAndRemainder(p.multiply(p)).quotient;
          roots.push(np.getCoefficient(0).negate().divide(np.getCoefficient(1)));
          np = Polynomial.of(np.getLeadingCoefficient());
          if (callback != undefined) {
            callback({content: content, roots: roots, newPolynomial: np, type: "solveCubicEquation"});
          }
          return roots;
        }
      } else {
        // 2*b^3-9*a*b*c+27*a^2*d
        var delta1 = TWO.multiply(b.pow(THREE)).subtract(NINE.multiply(a).multiply(b).multiply(c)).add(TWENTY_SEVEN.multiply(a.pow(TWO)).multiply(d));
        var C = undefined;
        if (delta0.equals(Expression.ZERO)) {
          C = nthRoot(3, delta1, np);
        } else {
          var tmp = nthRoot(2, delta1.pow(TWO).subtract(FOUR.multiply(delta0.pow(THREE))), np);
          if (tmp != undefined) {
            C = nthRoot(3, delta1.add(tmp).divide(TWO), np);
          }
        }
        if (C != undefined) {
          var root = b.add(C).add(delta0.divide(C)).negate().divide(THREE.multiply(a));
          roots.push(root);

          if (true) {
            var C1 = C.multiply(Expression.ONE.negate().add(THREE.squareRoot().multiply(Expression.I)).divide(TWO)); // C*(-1+sqrt(3)*i)/2
            var root1 = b.add(C1).add(delta0.divide(C1)).negate().divide(THREE.multiply(a));
            roots.push(root1);

            var C2 = C.multiply(Expression.ONE.negate().subtract(THREE.squareRoot().multiply(Expression.I)).divide(TWO)); // C*(-1-sqrt(3)*i)/2
            var root2 = b.add(C2).add(delta0.divide(C2)).negate().divide(THREE.multiply(a));
            roots.push(root2);

            np = Polynomial.of(np.getLeadingCoefficient());
            if (callback != undefined) {
              callback({content: content, roots: roots, newPolynomial: np, type: "solveCubicEquation"});
            }
          } else {
            np = np.divideAndRemainder(Polynomial.of(root.negate(), Expression.ONE)).quotient;
            if (callback != undefined) {
              callback({content: content, roots: roots, newPolynomial: np, type: "solveCubicEquation"});
            }
            continueWithNewPolynomial(roots, np);
          }

          return roots;
        }
      }
    }

    if (np.getDegree() >= 4) {
      if (true) {
        //TODO: performance
        var isSmall = function () {
          var c = 1;
          for (var i = 0; i <= np.getDegree(); i += 1) {
            var k = np.getCoefficient(i);
            //TODO: ilog2(k)
            if (!k.equals(Expression.ZERO)) {
              c *= k instanceof Expression.Integer ? Math.floor(Math.log(Math.abs(k.toNumber()) + 0.5) / Math.log(2)) + 1 : 1;
            }
          }
          return c <= 4 * 1024;
        };
        //console.time('Kronecker\'s method');
        var g = isSmall() ? np._getFactorByKroneckersMethod() : undefined; // too slow
        //console.timeEnd('Kronecker\'s method');
        if (g != undefined) {
          var h = np.divideAndRemainder(g).quotient;
          var gNew = null;
          var gRoots = g.getroots(function (x) {
            gNew = x.newPolynomial;
          });
          for (var i = 0; i < gRoots.length; i += 1) {
            roots.push(gRoots[i]);
            //np = np.divideAndRemainder(Polynomial.of(gRoots[i].negate(), Expression.ONE)).quotient;//TODO: optimize somehow - ?
          }
          if (gRoots.length > 0) {
            np = np.divideAndRemainder(g.divideAndRemainder(gNew).quotient).quotient;
          }
          var hNew = null;
          var hRoots = h.getroots(function (x) {
            hNew = x.newPolynomial;
          });
          for (var i = 0; i < hRoots.length; i += 1) {
            roots.push(hRoots[i]);
            //np = np.divideAndRemainder(Polynomial.of(hRoots[i].negate(), Expression.ONE)).quotient;//TODO: optimize somehow - ?
          }
          if (hRoots.length > 0) {
            np = np.divideAndRemainder(h.divideAndRemainder(hNew).quotient).quotient;
          }
          if (hRoots.length > 0 || gRoots.length > 0) {
            if (typeof hit === "function") {
              hit({getroots: {methodOfKronecker: np.toString()}});
            }
            if (callback != undefined) {
              callback({content: content, roots: roots, newPolynomial: np, type: "methodOfKronecker"});//?
            }
          }
          return roots;
        }
      }
    }

    //!2018-12-23
    if (np.getDegree() > 2) {
      // https://en.wikipedia.org/wiki/Square-free_polynomial
      //var gcd = function (a, b) {
      //  return b.getDegree() === -1 ? a : gcd(b, a.divideAndRemainder(b).remainder);
      //};
      var f = np;
      var d = f.derive();
      var a0 = Polynomial.polynomialGCD(f, d);
      if (a0.getDegree() > 0) {
        //TODO: merge with a code for Kronecker's method
        a0 = a0.scale(a0.getContent().inverse());//!?
        var q = f.divideAndRemainder(a0).quotient;
        var a0r = a0.getroots();
        for (var i = 0; i < a0r.length; i += 1) {
          var root = a0r[i];
          while (np.calcAt(root).equals(Expression.ZERO)) {
            roots.push(root);
            np = np.divideAndRemainder(Polynomial.of(root.negate(), Expression.ONE)).quotient;//TODO: optimize somehow - ?
          }
          q = q.divideAndRemainder(Polynomial.of(root.negate(), Expression.ONE)).quotient;//TODO: optimize somehow - ?
        }
        if (a0r.length > 0) {
          if (typeof hit === "function") {
            hit({getroots: {squareFreeFactorization: np.toString()}});
          }
          if (callback != undefined) {
            //TODO: better details, t = sqrt(3), show the polynomial, ...
            callback({content: content, roots: roots, newPolynomial: np, type: "squareFreeFactorization"});//?
          }
          continueWithNewPolynomial(roots, np);
          return roots;
        }
      }
    }
    //!

    //TODO: ???
    //TODO: move up
    if (np.getDegree() >= 3) {
      for (var i = 0; i <= np.getDegree(); i += 1) {
        if (Expression.has(np.getCoefficient(i), Expression.SquareRoot)) {
          var c = null;
          Expression._map(function (x) {
            if (c == null) {//TODO: fix - ?
              if (x instanceof Expression.SquareRoot && x.a instanceof Expression.Integer) {
                c = x;
              }
            }
            return x;
          }, np.getCoefficient(i));
          var tmp = new Expression.Symbol('_t');//?
          // substitute
          var p = np.map(function (coefficient) {
            var s1 = Expression.ZERO;
            // Expression._map does not work here as it goes into Expression.Exponentiation: x**2 -> x**(t**2). It throws an exception.
            for (var s of coefficient.summands()) {
              var t = Expression.ONE;
              for (var x of s.factors()) {
                if (x.equals(c)) {
                  t = t.multiply(tmp);
                } else if (x instanceof Expression.Integer) {
                  var exp = 0;
                  while (x.gcd(c.a).equals(c.a)) {
                    exp += 2;
                    x = x.divide(c.a);
                  }
                  t = t.multiply(x.multiply(tmp._pow(exp)));
                  //?
                  //var q = x.truncatingDivide(c.a);
                  //var r = x.subtract(q.multiply(c.a));
                  //return q.multiply(tmp.multiply(tmp)).add(r);
                } else {
                  t = t.multiply(x);
                }
              }
              s1 = s1.add(t);
            }
            return s1;
          });
          var a = "_x" + (++counter);
          var newp = Polynomial.toPolynomial(p.calcAt(new Expression.Symbol(a)), tmp);
          if (newp.getDegree() > 1) {
            var roots1 = newp.getroots();
            var rootsCount = roots.length;
            for (var j = 0; j < roots1.length; j += 1) {
              //TODO: check
              var roots2 = Polynomial.toPolynomial(roots1[j].subtract(c).getNumerator(), new Expression.Symbol(a)).getroots();
              for (var k = 0; k < roots2.length; k += 1) {
                var root = roots2[k];
                roots.push(root);
                np = np.divideAndRemainder(Polynomial.of(root.negate(), Expression.ONE)).quotient;//TODO: optimize somehow - ?
              }
            }
            if (roots.length > rootsCount) {
              if (typeof hit === "function") {
                hit({getroots: {methodOfIntroducingANewVariable: ""}});
              }
              if (callback != undefined) {
                //TODO: better details, t = sqrt(3), show the polynomial, ...
                callback({content: content, roots: roots, newPolynomial: np, type: "methodOfIntroducingANewVariable", t: c});//?
              }
              continueWithNewPolynomial(roots, np);
              return roots;//TODO: when it is not all roots
            }
          }
        }
      }
    }

    if (np.getDegree() >= 4) {
      //TODO: fix
      if (typeof hit === "function") {
        hit({getroots: {other: np.getDegree()}});
      }
    }

    return roots;
  };

  Polynomial.prototype.derive = function () {
    var newData = new PolynomialData(this.a.length - 1);
    for (var i = 1; i < this.a.length; i += 1) {
      var n = this.a.degree(i);
      var c = this.a.coefficient(i);
      newData.add(n - 1, c.multiply(Expression.Integer.fromNumber(n)));
    }
    return new Polynomial(newData);
  };

  Polynomial.prototype.getSquareFreePolynomial = function () {
    // https://en.wikipedia.org/wiki/Square-free_polynomial
    var p = this;
    var zero = 0;
    while (p.getCoefficient(0).equals(Expression.ZERO)) {
      p = p.divideAndRemainder(Polynomial.of(Expression.ONE).shift(1)).quotient;
      zero += 1;
    }
    var f = p;
    var d = f.derive();
    var g = d.getDegree() === -1 ? Polynomial.of(Expression.ONE) : Polynomial.polynomialGCD(f, d);
    return f.divideAndRemainder(g).quotient.multiply(Polynomial.of(Expression.ONE).shift(zero > 0 ? 1 : 0));
  };

  export default Polynomial;
