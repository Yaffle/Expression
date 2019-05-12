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
      throw new RangeError();
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
      throw new RangeError("ArithmeticException");
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
      remainder = remainder.subtract(p.multiply(pq).shift(n));
    }
    return {quotient: quotient, remainder: remainder};
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
      var d = Number.parseInt(x.degree.toString(), 10);
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

  Polynomial.prototype.doRationalRootTest = function (callback) {
    var np = this;
    var an = np.getLeadingCoefficient();
    var a0 = np.getCoefficient(0);

    //TODO: http://en.wikipedia.org/wiki/Polynomial_remainder_theorem
    //var fp1 = getBigInteger(np.calcAt(1));
    //var fm1 = getBigInteger(np.calcAt(-1));

    // p/q
    //TODO: forEach -> some ?
    Expression.everyDivisor(a0, function (p) {
      return Expression.everyDivisor(an, function (q) {
        var sign = -3;
        while ((sign += 2) < 3) {
          var sp = sign === -1 ? p.negate() : p;

          //if (// fp1.remainder(sp.subtract(q)).equals(ZERO) &&
                // fm1.remainder(sp.add(q)).equals(ZERO) &&
                // sp.gcd(q).equals(ONE)) {//?
            var x = Polynomial.of(sp.negate(), q);
            var z = np.divideAndRemainder(x, "undefined");
            while (z != undefined && z.remainder.equals(Polynomial.ZERO)) {
              np = z.quotient;
              np = np.scale(q);
              if (!callback(sp, q)) {
                np = undefined;
                return false;
              }
              if (np.getDegree() === 1) {
                sp = np.getCoefficient(0).negate();
                q = np.getCoefficient(1);
                np = Polynomial.of(q);
                if (!callback(sp, q)) {
                  np = undefined;
                  return false;
                }
                return false;// or divide
              }
              if (np.getDegree() < 1) {
                return false;
              }
              //TODO: !!!
              //an = np.getLeadingCoefficient();//!
              //a0 = np.getCoefficient(0);//!

              // fp1 = fp1.divide(q.subtract(sp));
              // fm1 = fm1.divide(q.negate().subtract(sp));
              z = np.divideAndRemainder(x, "undefined");
            }
          //}
        }
        return true;
      });
    });

    return np;
  };

  Polynomial.prototype._getFactorByKroneckersMethod = function () {
    // https://ru.wikipedia.org/wiki/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%9A%D1%80%D0%BE%D0%BD%D0%B5%D0%BA%D0%B5%D1%80%D0%B0
    // https://ru.wikipedia.org/wiki/%D0%98%D0%BD%D1%82%D0%B5%D1%80%D0%BF%D0%BE%D0%BB%D1%8F%D1%86%D0%B8%D0%BE%D0%BD%D0%BD%D1%8B%D0%B9_%D0%BC%D0%BD%D0%BE%D0%B3%D0%BE%D1%87%D0%BB%D0%B5%D0%BD_%D0%9B%D0%B0%D0%B3%D1%80%D0%B0%D0%BD%D0%B6%D0%B0
    // https://en.wikipedia.org/wiki/Vandermonde_matrix
    var np = this;
    var isInteger = function (c) {
      if (c instanceof Expression.Integer) {
        return true;
      }
      if (c instanceof Expression.Symbol) {
        return true;
      }
      if (c instanceof Expression.Addition) {
        return isInteger(c.a) && isInteger(c.b);
      }
      if (c instanceof Expression.Multiplication) {
        return isInteger(c.a) && isInteger(c.b);
      }
      if (c instanceof Expression.Exponentiation) {
        return isInteger(c.a) && c.b instanceof Expression.Integer;
      }
      return false;
    };
    var integerCoefficients = true;
    for (var i = 0; i <= np.getDegree(); i += 1) {
      integerCoefficients = integerCoefficients && isInteger(np.getCoefficient(i));
    }
    if (!integerCoefficients) {
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

  Polynomial.prototype.getroots = function (callback) {
    //TODO: merge hit and callback
    callback = callback || undefined;
    var np = this;

    var roots = [];


    //!new 2018-12-24
    //TODO: fix (?Polynomial#getContent()?)
      var t = Expression.ZERO;
      while (t != null) {
        var t = Expression.getConjugate(np.getCoefficient(np.getDegree()));
        if (t != undefined) {
          np = np.scale(t);
        }
      }
    //!

    var content = np.getContent();
    if (!content.equals(Expression.ONE)) {
      np = np.scale(content.getDenominator()).divideAndRemainder(Polynomial.of(content.getNumerator()), "throw").quotient;
      //np = np.divideAndRemainder(Polynomial.of(content), "throw").quotient;
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
      np = Polynomial.of(np.getCoefficient(1));
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
        return x.b.remainder(N).equals(Expression.ZERO) ? x.a.pow(x.b.divide(N)) : undefined;
      }
      if (x instanceof Expression.Multiplication) {
        var sa = nthRootInternal(n, x.a);
        var sb = nthRootInternal(n, x.b);
        return sa == undefined || sb == undefined ? undefined : sa.multiply(sb);
      }
      if (x instanceof Expression.Complex) {
        //TODO: - ?
        if (n === 2) {
          var m = x.real.multiply(x.real).add(x.imaginary.multiply(x.imaginary)).squareRoot();
          var g = nthRootInternal(2, x.real.add(m).divide(Expression.TWO));
          var d = nthRootInternal(2, x.real.negate().add(m).divide(Expression.TWO));
          if (g != undefined && d != undefined) {
            var result = g.add((x.imaginary.compareTo(Expression.ZERO) < 0 ? d.negate() : d).multiply(Expression.I));
            return result;
          }
        }
        if (x.real.equals(Expression.ZERO) && n % 2 === 0) {
          var c = nthRootInternal(Math.floor(n / 2), x);
          if (c != undefined) {
            return nthRootInternal(2, c);
          }
        }
        if (x.real.equals(Expression.ZERO) && n % 2 === 1) {
          //?
          var c = nthRootInternal(n, x.imaginary);
          if (c != undefined) {
            return c.multiply((n % 4 === 1 ? Expression.I : Expression.I.negate()));
          }
        }
      }
      if (x instanceof Expression.Addition) {
        var lastFactor = undefined;
        var e = 0;
        var result = Expression.ONE;
        var rest = Expression.ONE;
        Expression.everySimpleDivisor(x, function (f) {
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
          return result != undefined;
        });
        if (result !== Expression.ONE) {
          if (e !== 0) {
            rest = rest.multiply(Expression.pow(lastFactor, e));
          }
          var rn = nthRootInternal(n, rest);
          if (rn != undefined) {
            return result.multiply(rn);
          }
        }
      }
      var y = undefined;
      try {
        y = x._nthRoot(n);
      } catch (error) {
        //TODO:
      }
      return y;
    };

    var nthRoot = function (n, x, np) {
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
        np = np.divideAndRemainder(Polynomial.of(rs[i].negate(), Expression.ONE)).quotient;
      }
      return np;
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
        var ok = false;//?
        for (var k = 0; k < qRoots.length; k += 1) {
          var qRoot = qRoots[k];
          var s = nthRoot(g, qRoot, np);
          if (s != undefined) {
            roots.push(s);
            np = np.divideAndRemainder(Polynomial.of(s.negate(), Expression.ONE)).quotient;
            if (g % 2 === 0) {
              roots.push(s.negate());
              np = np.divideAndRemainder(Polynomial.of(s, Expression.ONE)).quotient;
            }
          }
          ok = ok || Expression.has(qRoot, Expression.Complex);//?
        }
        if (callback != undefined) {
          callback({content: content, roots: roots, newPolynomial: np, type: "t = x^g", g: g, allZeros: allZeros});//TODO: ?
        }
        if (n !== np.getDegree() && ok) {
          np = continueWithNewPolynomial(roots, np);
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
      var sD = nthRoot(2, D, np);
      if (typeof hit === "function") {
        hit({getroots: {quadratic: (sD == undefined ? (D instanceof Expression.Integer ? D.compareTo(Expression.ZERO) : "?" + D.toString()) : "OK")}});
      }
      if (sD != undefined) {
        var x1 = b.negate().subtract(sD).divide(Expression.TWO.multiply(a));
        var x2 = b.negate().add(sD).divide(Expression.TWO.multiply(a));
        roots.push(x1);
        roots.push(x2);
        np = Polynomial.of(a);
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
        if (isQuasiPalindromic && np.getDegree() <= 53) { // Math.log2(9007199254740991 + 1)
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

    if (np.getDegree() > 2) {
      var n = roots.length;
      np = np.doRationalRootTest(function (sp, q) {
        roots.push(sp.divide(q));
        return true;
      });

      if (n !== roots.length) {
        if (typeof hit === "function") {
          hit({getroots: {rational: ""}});
        }
        if (callback != undefined) {
          callback({content: content, roots: roots, newPolynomial: np, type: "useTheRationalRootTest"});
        }
        if (np.getDegree() > 0) {
          np = continueWithNewPolynomial(roots, np);
        }
        return roots;
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
          np = Polynomial.of(Expression.pow(root.negate().getDenominator(), n));
          if (callback != undefined) {
            callback({content: content, roots: roots, newPolynomial: np, type: "(ax+b)**n"});//TODO:
          }
          return roots;
        }
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
      if (delta0.equals(Expression.ZERO)) {
        //TODO: link to a^3+3a^2b+3ab^2+b^3=0 - ?
        if (delta.equals(Expression.ZERO)) {
          // -b/(3*a)
          var root = b.negate().divide(THREE.multiply(a));
          roots.push(root);
          roots.push(root);
          roots.push(root);
          var p = Polynomial.of(root.negate(), Expression.ONE);
          np = np.divideAndRemainder(p.multiply(p).multiply(p)).quotient;
          if (callback != undefined) {
            callback({content: content, roots: roots, newPolynomial: np, type: "solveCubicEquation"});
          }
          return roots;
        } else {
          // 2*b^3-9*a*b*c+27*a^2*d
          var delta1 = TWO.multiply(b.pow(THREE)).subtract(NINE.multiply(a).multiply(b).multiply(c)).add(TWENTY_SEVEN.multiply(a.pow(TWO)).multiply(d));
          var C = nthRoot(3, delta1, np);
          if (C != undefined) {
            var root = np.getCoefficient(2).add(C).add(delta0.divide(C)).negate().divide(THREE.multiply(np.getCoefficient(3)));
            roots.push(root);
            np = np.divideAndRemainder(Polynomial.of(root.negate(), Expression.ONE)).quotient;
            if (callback != undefined) {
              callback({content: content, roots: roots, newPolynomial: np, type: "solveCubicEquation"});
            }
            np = continueWithNewPolynomial(roots, np);
            return roots;
          }
        }
      } else {
        // https://github.com/nicolewhite/algebra.js/blob/master/src/equations.js
        // https://en.wikipedia.org/wiki/Cubic_function#Multiple_roots.2C_.CE.94_.3D_0
        if (delta.equals(Expression.ZERO)) {
          // a double root
          // 9*a*d-b*c
          var root = NINE.multiply(a).multiply(d).subtract(b.multiply(c)).divide(TWO.multiply(delta0));
          roots.push(root);
          roots.push(root);
          var p = Polynomial.of(root.negate(), Expression.ONE);
          np = np.divideAndRemainder(p.multiply(p)).quotient;
          roots.push(np.getCoefficient(0).negate().divide(np.getCoefficient(1)));
          np = Polynomial.of(np.getCoefficient(1));
          if (callback != undefined) {
            callback({content: content, roots: roots, newPolynomial: np, type: "solveCubicEquation"});
          }
          return roots;
        } else {
          //?
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
              c *= k instanceof Expression.Integer ? Math.floor(Math.log(Math.abs(Number.parseInt(k.toString(), 10)) + 0.5) / Math.log(2)) + 1 : 1;
            }
          }
          return c <= 4 * 1024;
        };
        //console.time('Kronecker\'s method');
        var g = isSmall() ? np._getFactorByKroneckersMethod() : undefined; // too slow
        //console.timeEnd('Kronecker\'s method');
        if (g != undefined) {
          var h = np.divideAndRemainder(g).quotient;
          var gr = g.getroots();
          for (var i = 0; i < gr.length; i += 1) {
            roots.push(gr[i]);
            np = np.divideAndRemainder(Polynomial.of(gr[i].negate(), Expression.ONE)).quotient;//TODO: optimize somehow - ?
          }
          var hr = h.getroots();
          for (var i = 0; i < hr.length; i += 1) {
            roots.push(hr[i]);
            np = np.divideAndRemainder(Polynomial.of(hr[i].negate(), Expression.ONE)).quotient;//TODO: optimize somehow - ?
          }
          if (hr.length > 0 || gr.length > 0) {
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
      var gcd = function (a, b) {
        return b.getDegree() === -1 ? a : gcd(b, a.divideAndRemainder(b).remainder);
      };
      var f = np;
      var d = f.derive();
      var a0 = gcd(f, d);
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
            hit({getroots: {squareFreeRoot: np.toString()}});
          }
          np = continueWithNewPolynomial(roots, np);
          return roots;
        }
      }
    }
    //!

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

  export default Polynomial;
