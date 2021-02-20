/*global hit*/

  import Expression from './Expression.js';
  import ExpressionParser from './ExpressionParser.js';
  import Matrix from './Matrix.js';
  import Heap from './Heap.js';

  function PolynomialTerm(degree, coefficient) {
    if (degree < 0 || degree > Number.MAX_SAFE_INTEGER || Math.floor(degree) !== degree) {
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
    if (k >= this.array.length || k < 0) {
      throw new RangeError();
    }
    if (k >= 1 && degree <= this.array[k - 1].degree) {
      throw new RangeError();
    }
    this.array[k] = Object.freeze(new PolynomialTerm(degree, coefficient));
    this.i = k + 1;
  };
  PolynomialData.prototype.degree = function (i) {
    return this.array[i].degree;
  };
  PolynomialData.prototype.coefficient = function (i) {
    return this.array[i].coefficient;
  };
  PolynomialData.prototype.trim = function () {
    if (this.i !== this.array.length) {
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
    this.a = a.trim();
  }

  Polynomial.of = function () {
    var newData = new PolynomialData(arguments.length);
    for (var i = 0; i < arguments.length; i += 1) {
      var a = arguments[i];
      if (!a.equals(Expression.ZERO)) {
        newData.add(i, a);
      }
    }
    return new Polynomial(newData);
  };
  Polynomial.from = function (array) {
    var newData = new PolynomialData(array.length);
    for (var i = 0; i < array.length; i += 1) {
      var a = array[i];
      if (!a.equals(Expression.ZERO)) {
        newData.add(i, a);
      }
    }
    return new Polynomial(newData);
  };

  Polynomial.ZERO = Polynomial.of();

  Polynomial.prototype.getDegree = function () {
    return this.a.length === 0 ? -1 : this.a.degree(this.a.length - 1);
  };
  Polynomial.prototype.getCoefficient = function (degree) {
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
    return new Polynomial(newData);
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
    if (p.a.length === 0) {
      return this;
    }
    if (this.a.length === 0) {
      return p;
    }
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
    return new Polynomial(newData);
  };

  Polynomial.prototype.multiply = function (p) {
    if (this.a.length === 0 || p.a.length === 0) {
      return Polynomial.ZERO;
    }
    var result = new FastAdditionPolynomial();
    if (this.a.length <= p.a.length) {
      for (var i = this.a.length - 1; i >= 0; i -= 1) {
        var d = this.a.degree(i);
        var c = this.a.coefficient(i);
        result.add(c, d, p, true);
      }
    } else {
      for (var i = p.a.length - 1; i >= 0; i -= 1) {
        var d = p.a.degree(i);
        var c = p.a.coefficient(i);
        result.add(c, d, this, false);
      }
    }
    return result.toPolynomial();
  };

  Polynomial.prototype.shift = function (n) { // *= x**n, n >= 0
    if (!(n >= 0)) {
      throw new TypeError();
    }
    var newData = new PolynomialData(this.a.length);
    for (var i = 0; i < this.a.length; i += 1) {
      newData.add(this.a.degree(i) + n, this.a.coefficient(i));
    }
    return new Polynomial(newData);
  };

  function FastAdditionPolynomial() {
    // see http://www.cecm.sfu.ca/~mmonagan/teaching/TopicsinCA11/johnson.pdf
    // see https://en.wikipedia.org/wiki/K-way_merge_algorithm#Heap
    this.maxHeap = new Heap(function (a, b) {
      return b.degree - a.degree;
    });
    this.degree = -1;
    this.leadingCoefficient = Expression.ZERO;
  }
  FastAdditionPolynomial.prototype.add = function (scale, shift, polynomial, fromLeft) {
    var iterator = {
      i: polynomial.a.length,
      next: function () {
        if (this.i === 0) {
          return undefined;
        }
        this.i -= 1;
        let c = polynomial.a.coefficient(this.i);
        return Object.freeze({
          degree: polynomial.a.degree(this.i) + shift,
          coefficient: fromLeft ? scale.multiply(c) : c.multiply(scale),
          iterator: iterator
        });
      }
    };
    var newEntry = iterator.next();
    if (newEntry.degree > this.degree && this.degree !== -1) {
      throw new RangeError();
    } else if (newEntry.degree < this.degree) {
      this.maxHeap.push(newEntry);
    } else {
      if (this.degree === -1) {
        this.leadingCoefficient = newEntry.coefficient;
        this.degree = newEntry.degree;
      } else {
        this.leadingCoefficient = this.leadingCoefficient.add(newEntry.coefficient);
      }
      newEntry = newEntry.iterator.next();
      if (newEntry != undefined) {
        this.maxHeap.push(newEntry);
      }
      // Computation of the leading coefficient:
      while (this.maxHeap.size() > 0 && this.leadingCoefficient.equals(Expression.ZERO)) {
        var tmp = this.maxHeap.peek();
        this.leadingCoefficient = tmp.coefficient;
        this.degree = tmp.degree;
        var next = tmp.iterator.next();
        if (next != undefined) {
          this.maxHeap.pop();
          this.maxHeap.push(next);
        } else {
          this.maxHeap.pop();
        }
        while (this.maxHeap.size() > 0 && this.degree === this.maxHeap.peek().degree) {
          var tmp = this.maxHeap.peek();
          this.leadingCoefficient = this.leadingCoefficient.add(tmp.coefficient);
          var next = tmp.iterator.next();
          if (next != undefined) {
            this.maxHeap.pop();
            this.maxHeap.push(next);
          } else {
            this.maxHeap.pop();
          }
        }
      }
      if (this.maxHeap.size() === 0 && this.leadingCoefficient.equals(Expression.ZERO)) {
        this.degree = -1;
      }
    }
  };
  FastAdditionPolynomial.prototype.getDegree = function () {
    return this.degree;
  };
  FastAdditionPolynomial.prototype.getLeadingCoefficient = function () {
    return this.leadingCoefficient;
  };
  FastAdditionPolynomial.prototype.toPolynomial = function () {
    var terms = [];
    var ONE = Polynomial.of(Expression.ONE);
    while (this.getDegree() >= 0) {
      var degree = this.getDegree();
      var coefficient = this.getLeadingCoefficient();
      this.add(coefficient.negate(), degree, ONE);
      terms.push({
        degree: degree,
        coefficient: coefficient
      });
    }
    return Polynomial.fromTerms(terms);
  };



  Polynomial.prototype.divideAndRemainder = function (p, w) {
    w = w || undefined;
    if (p.equals(Polynomial.ZERO)) {
      throw new TypeError("ArithmeticException");
    }
    if (this.equals(Polynomial.ZERO)) {
      return {quotient: this, remainder: this};
    }
    var quotient = [];
    var remainder = new FastAdditionPolynomial();
    remainder.add(Expression.ONE, 0, this, true);
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
      //TODO: optimize - ?
      quotient.push({coefficient: q, degree: n});
      //TODO: optimize - ?
      remainder.add(q.negate(), n, p, true);
      if (remainder.getDegree() - p.getDegree() === n) {
        // to avoid the infite loop
        throw new TypeError("there is a some problem with the expression evaluation");//!
      }
    }
    var q = new PolynomialData(quotient.length);
    for (var i = quotient.length - 1; i >= 0; i -= 1) {
      q.add(quotient[i].degree, quotient[i].coefficient);
    }
    return {quotient: new Polynomial(q), remainder: remainder.toPolynomial()};
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
    if (point instanceof Expression.Division &&
        point.getNumerator() instanceof Expression.Integer &&
        point.getDenominator() instanceof Expression.Integer &&
        this.hasIntegerCoefficients()) {
      var n = this.getDegree();
      var p = this.map(function (coefficient, degree) {
        return coefficient.multiply(Expression.pow(point.getDenominator(), n - degree));
      });
      return p.calcAt(point.getNumerator()).divide(point.getDenominator()._pow(n));
    }
    var n = Expression.ZERO;
    var lastDegree = -1;
    var i = this.a.length;
    while (--i >= -1) {
      var degree = i === -1 ? 0 : this.a.degree(i);
      var coefficient = i === -1 ? Expression.ZERO : this.a.coefficient(i);
      if (i !== this.a.length - 1) {
        n = Expression.pow(point, lastDegree - degree).multiply(n);
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
    if (x.equals(Expression.ONE)) {
      return this;
    }
    if (Expression.has(x, Expression.Matrix) || Expression.has(x, Expression.MatrixSymbol)) {
      throw new TypeError();
    }
    return this.map(function (coefficient, degree) {
      return coefficient.multiply(x);
    });
  };

  Polynomial.fromTerms = function (terms) {
    var data = new PolynomialData(terms.length);
    for (var i = terms.length - 1; i >= 0; i -= 1) {
      var term = terms[i];
      var d = term.degree;
      var c = term.coefficient;
      data.add(d, c);
    }
    return new Polynomial(data);
  };

  Polynomial.toPolynomial = function (e, v) {
    if (e instanceof Expression.Division) {
      throw new RangeError();
    }
    var terms = Expression.getCoefficients(e, v);
    return Polynomial.fromTerms(terms);
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
    if (np.getCoefficient(0).equals(Expression.ZERO)) {
      return Expression.ZERO;//!TODO: test
    }

    var an = np.getLeadingCoefficient();
    var a0 = np.getCoefficient(0);
    a0 = Expression._expandTrigonometry(a0);//!
    if (np.getDegree() === 1) {
      return a0.negate().divide(an);
    }

    //TODO: http://en.wikipedia.org/wiki/Polynomial_remainder_theorem

    // http://scask.ru/g_book_mav.php?id=26
    var hasIntegerCoefficients = np.hasIntegerCoefficients();
    // f(k) = q(k)(k - a)
    var fp1 = null;
    var fm1 = null;
    if (hasIntegerCoefficients) {
      fp1 = np.calcAt(Expression.ONE);
      if (fp1.equals(Expression.ZERO)) {
        return Expression.ONE;
      }
      fm1 = np.calcAt(Expression.ONE.negate());
      if (fm1.equals(Expression.ZERO)) {
        return Expression.ONE.negate();
      }
    }
    var filter = function (n, d) {
      if (fp1 != null) {
        if (d.subtract(n).equals(Expression.ZERO)) {
          return false;
        }
        if (!fp1.remainder(d.subtract(n)).equals(Expression.ZERO)) {
          return false;
        }
      }
      if (fm1 != null) {
        if (d.add(n).equals(Expression.ZERO)) {
          return false;
        }
        if (!fm1.remainder(d.add(n)).equals(Expression.ZERO)) {
          return false;
        }
      }
      return true;
    };

    //!new 2020-01-13
    if (hasIntegerCoefficients) {
      var tmp = np.squareFreeFactors();
      if (tmp.a0.getDegree() > 0) {
        return tmp.a0.doRationalRootTest() || tmp.a1.doRationalRootTest();
      }
      var zeros = np.getZeros();
      for (var i = 0; i < zeros.length; i += 1) {
        var zero = zeros[i];
        if (i === 0 || zero !== zeros[i - 1]) {
          var candidate = zero.root != null ? Expression.Integer.fromString(toDecimalStringInternal(an.equals(Expression.ONE) ? zero.root : new Expression.Multiplication(zero.root, an), {fractionDigits: 0})).divide(an) : zero;
          if (filter(candidate.getNumerator(), candidate.getDenominator()) && np.calcAt(candidate).equals(Expression.ZERO)) {
            return candidate;
          }
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


          if (//sp.gcd(q).equals(Expression.ONE) &&
              filter(sp, q)) {//?
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

  Polynomial.prototype.hasIntegerCoefficients = function () {
    for (var i = 0; i < this.a.length; i += 1) {
      if (!(this.a.coefficient(i) instanceof Expression.Integer)) {
        return false;
      }
    }
    return true;
  };
  Polynomial.prototype.hasComplexCoefficients = function () {
    for (var i = 0; i < this.a.length; i += 1) {
      if (!(this.a.coefficient(i) instanceof Expression.Complex) && !(this.a.coefficient(i) instanceof Expression.Integer)) {
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

  Polynomial.prototype._canBeFactored = function (depth) {
    // https://en.wikipedia.org/wiki/Eisenstein%27s_criterion
    if (!this.hasIntegerCoefficients()) {
      //throw new Error();
      return true;
    }
    if (this.getCoefficient(0).equals(Expression.ZERO)) {
      return true;
    }
    var content = this.getContent();
    if (!content.equals(Expression.ONE) && !content.equals(Expression.ONE.negate())) {
      return this.scale(content.inverse())._canBeFactored();
    }
    // p divides each a_i for 0 ≤ i < n
    var g = this.subtract(Polynomial.of(this.getLeadingCoefficient()).shift(this.getDegree())).getContent();
    // p does not divide a_n
    var x = null;
    do {
      x = g.gcd(this.getLeadingCoefficient());
      g = g.truncatingDivide(x);
    } while (!x.equals(Expression.ONE));
    var x = null;
    g = g.abs();//?
    while (!g.equals(Expression.ONE)) {
      var p = g.primeFactor();
      // p**2 does not divide a_0
      if (!this.getCoefficient(0).remainder(p._pow(2)).equals(Expression.ZERO)) {
        return false;
      }
      g = g.truncatingDivide(p);
    }
    if (depth == undefined) {
      // see https://en.wikipedia.org/wiki/Eisenstein%27s_criterion#Indirect_(after_transformation)
      //TODO: ?
      if (!this.subs(x => x.add(Expression.Integer.fromNumber(3)))._canBeFactored(1)) {
        return false;
      }
      //TODO: (too slow)
      //if (!this.subs(x => x.inverse())._canBeFactored(1)) {
      //  return false;
      //}
    }
    return true;
  };

  Polynomial.prototype.isEven = function () {
    for (var i = 1; i <= this.getDegree(); i += 2) {
      if (!this.getCoefficient(i).equals(Expression.ZERO)) {
        return false;
      }
    }
    return true;
  };

  Polynomial.prototype._getFactorByKroneckersMethod = function () {
    // https://ru.wikipedia.org/wiki/Метод_Кронекера
    // https://ru.wikipedia.org/wiki/Интерполяционный_многочлен_Лагранжа
    // https://en.wikipedia.org/wiki/Vandermonde_matrix
    var np = this;
    if (np.getDegree() < 2) {
      return undefined;
    }
    if (!np._hasIntegerLikeCoefficients()) {
      return undefined;
    }
    var n = np.getDegree();
    var ys = new Array(n);
    var total = 1;
    var x = function (i) {
      if (np.isEven()) {
        return Expression.Integer.fromNumber((i % 2 === 0 ? 0 - i / 2 : (i + 1) / 2));
      }
      return Expression.Integer.fromNumber(i);
    };
    var findSomeSmallIntegerPoints = function (polynomial, count) {
      let candidates = [];
      for (let j = 0; j < count * 2; j += 1) {
        const i = x(j);
        candidates.push({
          i: i,
          y: np.calcAt(i)
        });
      }
      if (polynomial.hasIntegerCoefficients()) {
        var d = function (n) {
          var count = 0;
          if (n.abs().toNumber() > Math.pow(Number.MAX_SAFE_INTEGER, 2)) {
            return 1/0;
          }
          Expression.everyDivisor(n.abs(), function (d) {
            count += 1;
            return true;
          });
          return count;
        };
        //var zeros = polynomial.multiply(polynomial).derive().getZeros();
        var zeros = polynomial.getZeros().concat(polynomial.derive().getZeros());
        for (const zero of zeros) {
          var i = Expression.Integer.fromString(toDecimalStringInternal(zero.root != null ? zero.root : zero, {fractionDigits: 0}));
          //TODO: siblings
          for (var v = -n; v <= n; v += 1) {
            var y = np.calcAt(i.add(Expression.Integer.fromNumber(v)));
            candidates.push({
              i: i.add(Expression.Integer.fromNumber(v)),
              y: y
            });
          }
        }
        candidates.sort((a, b) => a.i.compareTo(b.i));
        // remove duplicates:
        let unique = [];
        for (let i = 0; i < candidates.length; i += 1) {
          if (unique.length === 0 || !candidates[i].i.equals(unique[unique.length - 1].i)) {
            unique.push(candidates[i]);
          }
        }
        candidates = unique;
        //candidates.sort((a, b) => a.y._pow(2).compareTo(b.y._pow(2)));
        candidates = candidates.map(function (c) {
          return {i: c.i, ni: c.i.toNumber(), y: c.y.toNumber(), d: d(c.y)};
        });
        candidates.sort((a, b) => (a.d - b.d) || (Math.abs(a.y) - Math.abs(b.y)));
      }
      return candidates.map(c => c.i);
    };
    var integerPoints = findSomeSmallIntegerPoints(np, Math.floor(n / 2) + 1);
    for (var i = 0; i <= Math.floor(n / 2); i += 1) {
      const bi = integerPoints[i];
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
      if (Expression.isConstant(y) && y.abs().toNumber() > 2**106) {//!
        throw new Error("TOO SLOW");
      }
      var divisors = [];
      Expression.everyDivisor(y, function (d) {
        divisors.push(d);
        return true;
      });
      // let the first be positive as two sets with different signs give polynomials that differ only in sign of coefficients
      ys[i] = i === 0 ? divisors : attachNegative(divisors);
      total *= ys[i].length;
      var V = Matrix.Zero(i + 1, i + 1).map(function (e, i, j) {
        return integerPoints[i]._pow(j);
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
        if (j > 1e5) throw new Error("TOO SLOW");
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
          if (!gc.equals(Expression.ONE) && !gc.equals(Expression.ONE.negate())) {//?
          g = g.scale(gc.getDenominator()).divideAndRemainder(Polynomial.of(gc.getNumerator()), "throw").quotient;
          var t = np.divideAndRemainder(g, "undefined");
          if (t != undefined && t.remainder.equals(Polynomial.ZERO)) {
            return g;
          }
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
        var t = Expression.getConjugate(np.getLeadingCoefficient());
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
      if (x instanceof Expression.ExpressionWithPolynomialRoot) {
        return undefined;//?
      }
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
        //TODO: different cases (1+sqrt(2)+sqrt(3)) - (?)
        var ac = Expression.getConjugateExpression(a.getNumerator());
        if ((n === 3 || n === 2) && ac instanceof Expression.Integer) {//TODO: ?
          if (n === 2 && Expression._isPositive(x.negate())) {
            return Expression.I.multiply(nthRootInternal(2, x.negate()));
          }
          //TODO: a > 0 - ?
          var a = x;
          var tmp = new Expression.Symbol('x')._pow(n).subtract(a).getNumerator();
          var polynomial = Polynomial.toPolynomial(Expression.getConjugateExpression(tmp), new Expression.Symbol('x'));
          var tmp2 = polynomial.getZeros();//TODO: ?
          for (const zero of tmp2) {
            if (zero._pow(n).equals(x)) {
              return zero;
            }
          }
          //TODO: fix
          return null;
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

    var continueWithNewPolynomial = function (roots, np, newPolynomialVariable) {
      var rs = np.getroots(callback != undefined ? function (info) {
        var xxx = Object.assign({}, info, {content: content.multiply(info.content), roots: roots.concat(info.roots), newPolynomialVariable: newPolynomialVariable});
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
      for (var i = 1; i <= np.getDegree() && g >= 2; i += 1) {
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
        var qRoots = [];
        if (!allZeros) {
          if (callback != undefined) {
            callback({content: content, roots: roots, newPolynomial: q, type: "t = x^g", g: g, newPolynomialVariable: new Expression.Symbol('t')});//TODO: ?
          }
          continueWithNewPolynomial(qRoots, q, new Expression.Symbol('t'));
        } else {
          qRoots = q.getroots();
        }
        var n = np.getDegree();//TODO: 2018-02-04
        //var ok = false;//?
        for (var k = 0; k < qRoots.length; k += 1) {
          var qRoot = qRoots[k];
          var s = nthRoot(g, qRoot, np);
          if (s != undefined) {
            var d = null;
            if ((!allZeros || g > 4) && g <= 24 && ((17896830 >> g) % 2) === 1) {
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
        if (!allZeros) {
          if (callback != undefined) {
            callback({content: content, roots: roots, newPolynomial: np, type: "x = t^(1/g)", g: g});//TODO: ?
          }
        } else {
          var type = (g === 2 ? "applyDifferenceOfSquaresRule" : (g === 3 ? "applyDifferenceOfCubesRule" : "applyDifferenceOfNthPowersRule"));
          if (callback != undefined) {
            callback({content: content, roots: roots, newPolynomial: np, type: type, g: g});//TODO: ?
          }
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
        if (sD.equals(Expression.ZERO)) {
          var x12 = b.negate().divide(Expression.TWO.multiply(a));
          roots.push(x12);
          roots.push(x12);
          //TODO: different details (?)
        } else {
          var x1 = b.negate().subtract(sD).divide(Expression.TWO.multiply(a));
          var x2 = b.negate().add(sD).divide(Expression.TWO.multiply(a));
          roots.push(x1);
          roots.push(x2);
        }
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
        for (var i = 2; i < middle + 1 && isQuasiPalindromic; i += 1) {
          isQuasiPalindromic = isQuasiPalindromic && np.getCoefficient(middle + i).pow(jj).subtract(np.getCoefficient(middle - i).pow(jj).multiply(mj.pow(Expression.Integer.fromNumber(i)))).equals(Expression.ZERO);
        }
        if (isQuasiPalindromic) {
          //TODO: fix
          if (typeof hit === "function") {
            hit({getroots: {quasiPalindromic: np.getDegree()}});
          }
        }
        if (isQuasiPalindromic && np.getDegree() <= Math.log2(Number.MAX_SAFE_INTEGER + 1)) {
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
              //np = np.divideAndRemainder(Polynomial.of(root.negate(), Expression.ONE)).quotient;//TODO: optimize
              roots.push(root);
            }
            np = np.divideAndRemainder(u.scale(u.getLeadingCoefficient().inverse())).quotient;
            //TODO: multiply by "newU"
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
      var hasZeroCoefficient = function (np) {
        for (var i = 0; i <= np.getDegree(); i += 1) {
          if (np.getCoefficient(i).equals(Expression.ZERO)) {
            return true;
          }
        }
        return false;
      };
      if (!hasZeroCoefficient(np)) {
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

    //TODO: depressed for all degrees (?) in Polynomial#getZeros() - ?
    if (!np.hasIntegerCoefficients() && np.getDegree() === 3) {//?
      // convert to depressed (?)
      // x = t - b / (n * a)
      var n = np.getDegree();
      var a = np.getLeadingCoefficient();
      var b = np.getCoefficient(n - 1);
      if (!b.equals(Expression.ZERO)) {
        var f = x => x.subtract(b.divide(Expression.Integer.fromNumber(n).multiply(a)));
        var p = np.subs(f);
        if (p.hasIntegerCoefficients()) {//?
          var zeros = p.getroots();
          if (zeros.length === 0) {//?
            zeros = p.getZeros();
            debugger;
          }
          if (zeros.length === np.getDegree()) {//?
            for (const zero of zeros) {
              roots.push(f(zero));
            }
            np = Polynomial.of(np.getLeadingCoefficient());
            if (callback != undefined) {
              callback({content: content, roots: roots, newPolynomial: np, type: "solveCubicEquation"});
            }
            return roots;
          } else {
            debugger;
          }
        }
      }
    }
   if (np.hasIntegerCoefficients()) {//?
      // convert to depressed (?)
      // x = t - b / (n * a)
      var n = np.getDegree();
      var a = np.getLeadingCoefficient();
      var b = np.getCoefficient(n - 1);
      if (!b.equals(Expression.ZERO)) {
        var f = x => x.subtract(b.divide(Expression.Integer.fromNumber(n).multiply(a)));
        var depressed = np.subs(f);
        var depressedRoots = [];
        if (callback != undefined) {
          // https://en.wikipedia.org/wiki/Algebraic_equation#Elimination_of_the_sub-dominant_term
          callback({content: content, roots: roots, newPolynomial: depressed, type: "eliminationOfTheSubDominantTerm", b: b, n: n, a: a, newPolynomialVariable: new Expression.Symbol('t')});
        }
        continueWithNewPolynomial(depressedRoots, depressed, new Expression.Symbol('t'));//TODO: ?
          for (const depressedRoot of depressedRoots) {
            roots.push(f(depressedRoot));
          }
          if (depressedRoots.length = depressed.getDegree()) {
            np = Polynomial.of(np.getLeadingCoefficient());
          } else {
            //?
            throw new TypeError();
          }
          if (callback != undefined) {
            callback({content: content, roots: roots, newPolynomial: np, type: "t = x - b/(n*a)", g: "?"});//TODO: ?
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
      var THREE = Expression.Integer.fromNumber(3);
      var substitute = function (y) {
        return y.subtract(b.divide(THREE.multiply(a)))
      };
      var tmp = np.subs(substitute);
      var depressed = tmp.scale(tmp.getLeadingCoefficient().inverse());
      var p = depressed.getCoefficient(1);
      var q = depressed.getCoefficient(0);
      var discriminant = p.divide(THREE)._pow(3).add(q.divide(Expression.TWO)._pow(2));
      if (typeof hit === "function") {
        hit({getroots: {cubic: (discriminant instanceof Expression.Integer ? discriminant.compareTo(Expression.ZERO) : "?") + "-" + (p instanceof Expression.Integer ? p.compareTo(Expression.ZERO) : "?")}});
      }
      var minusOneOverTwo = Expression.ONE.negate().divide(Expression.TWO);
      var iSqrtOfThreeOverTwo = Expression.I.multiply(THREE.squareRoot()).divide(Expression.TWO);
      var cbrtOfMinusOne1 = minusOneOverTwo.subtract(iSqrtOfThreeOverTwo); // (-1-sqrt(3)*i)/2
      var cbrtOfMinusOne2 = minusOneOverTwo.add(iSqrtOfThreeOverTwo); // (-1+sqrt(3)*i)/2
      if (q.equals(Expression.ZERO) && p.equals(Expression.ZERO)) {
        //TODO: link to a^3+3a^2b+3ab^2+b^3=0 - ?
        // -b/(3*a)
        var root = substitute(Expression.ZERO);
        roots.push(root);
        roots.push(root);
        roots.push(root);
        np = Polynomial.of(np.getLeadingCoefficient());
        if (callback != undefined) {
          callback({content: content, roots: roots, newPolynomial: np, type: "solveCubicEquation"});
        }
        return roots;
      } else if (q.equals(Expression.ZERO)) {
        roots.push(substitute(Expression.ZERO));
        var tmp = nthRoot(2, p.negate(), np);
        roots.push(substitute(tmp));
        roots.push(substitute(tmp.negate()));
        np = Polynomial.of(np.getLeadingCoefficient());
        if (callback != undefined) {
          callback({content: content, roots: roots, newPolynomial: np, type: "solveCubicEquation"});
        }
        return roots;
      } else if (p.equals(Expression.ZERO)) {
        //TODO: should not reach this point (?) - should be solved by previous methods
        var tmp = nthRoot(3, q.negate(), np);
        roots.push(substitute(tmp));
        roots.push(substitute(tmp.multiply(cbrtOfMinusOne1)));
        roots.push(substitute(tmp.multiply(cbrtOfMinusOne2)));
        np = Polynomial.of(np.getLeadingCoefficient());
        if (callback != undefined) {
          callback({content: content, roots: roots, newPolynomial: np, type: "solveCubicEquation"});
        }
        return roots;
      } else if (discriminant.equals(Expression.ZERO)) {
        // https://github.com/nicolewhite/algebra.js/blob/master/src/equations.js
        // https://en.wikipedia.org/wiki/Cubic_equation#Multiple_root
        // a double root
        var t23 = THREE.multiply(q).divide(Expression.TWO.multiply(p)).negate();
        var root23 = substitute(t23);
        roots.push(root23);
        roots.push(root23);
        var t1 = t23.multiply(Expression.TWO).negate();
        var root1 = substitute(t1);
        roots.push(root1);
        np = Polynomial.of(np.getLeadingCoefficient());
        if (callback != undefined) {
          callback({content: content, roots: roots, newPolynomial: np, type: "solveCubicEquation"});
        }
        return roots;
      } else {
        // https://en.wikipedia.org/wiki/Cubic_equation#Cardano's_formula
        // 2*b^3-9*a*b*c+27*a^2*d
        var tmp = nthRoot(2, discriminant, np);
        if (tmp != undefined) {
          var C = nthRoot(3, q.negate().divide(Expression.TWO).add(tmp), np);
          if (C != undefined) {
            var rootFromC = function (C) {
              return substitute(C.subtract(p.divide(THREE.multiply(C))));
            };
            roots.push(rootFromC(C));
            roots.push(rootFromC(C.multiply(cbrtOfMinusOne1)));
            roots.push(rootFromC(C.multiply(cbrtOfMinusOne2)));
            np = Polynomial.of(np.getLeadingCoefficient());
            if (callback != undefined) {
              callback({content: content, roots: roots, newPolynomial: np, type: "solveCubicEquation"});
            }
            return roots;
          }
        }
      }
    }

    //!2018-12-23
    if (np.getDegree() > 2) {
      // https://en.wikipedia.org/wiki/Square-free_polynomial
      var tmp = np.squareFreeFactors();
      var a0 = tmp.a0;
      var a1 = tmp.a1;
      if (a0.getDegree() > 0) {
        //TODO: merge with a code for Kronecker's method
        //TODO: factorization - ?
        var newA0 = a0;
        var a0r = a0.getroots(function (x) {
          newA0 = x.newPolynomial;
        });
        var previousRoot = null;
        for (var i = 0; i < a0r.length; i += 1) {
          var root = a0r[i];
            roots.push(root);
            if (previousRoot == null || !previousRoot.equals(root)) {
              roots.push(root);
            }
            previousRoot = root;
        }
        // find roots of a1 at first (for better performance):
        var newA1 = a1;
        var a1Roots = a1.getroots(function (x) {
          newA1 = x.newPolynomial;
        });
        for (var i = 0; i < a1Roots.length; i += 1) {
          roots.push(a1Roots[i]);
        }
        if (newA0 != null) {
          //TODO: test
          np = newA1.multiply(newA0).multiply(Polynomial.polynomialGCD(newA0, np));
        }
        if (a0r.length > 0 || a1Roots.length > 0) {
          if (typeof hit === "function") {
            hit({getroots: {squareFreeFactorization: np.toString()}});
          }
          if (callback != undefined) {
            //TODO: better details, t = sqrt(3), show the polynomial, ...
            callback({content: content, roots: roots, newPolynomial: np, type: "squareFreeFactorization"});//?
          }
          //continueWithNewPolynomial(roots, np);
          return roots;
        }
      }
    }
    //!

    if (np.getDegree() >= 4) {
      if (true) {
        //TODO: !!! show correct method name in details
        var g = np.factorize();
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
              //TODO: better details
              callback({content: content, roots: roots, newPolynomial: np, type: "methodOfKronecker"});//?
            }
          }
          return roots;
        }
      }
    }

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
        if (c != null) {
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
          //!2020-12-11
          var gcdOfDegrees = function () {
            var gcd = function (a, b) {
              return b === 0 ? a : gcd(b, a % b);
            };
            var g = newp.getDegree();
            for (var i = 1; i <= newp.getDegree() && g >= 2; i += 1) {
              if (!newp.getCoefficient(i).equals(Expression.ZERO)) {
                g = gcd(g, i);
              }
            }
            return g;
          };
          var g = gcdOfDegrees();
          if (g >= 2) {
            var newData = new Array(Math.floor((newp.getDegree() + g) / g));
            var k = 0;
            for (var i1 = 0; i1 <= newp.getDegree(); i1 += g) {
              newData[k] = newp.getCoefficient(i1).multiply(c._pow(i1 - k));
              k += 1;
            }
            newp = Polynomial.from(newData);
          }
          //!
          if (newp.getDegree() > 1) {
            var roots1 = [];
            var rr = null;
            //TODO: ?
            while ((rr = newp.doRationalRootTest()) != null) {
              roots1.push(rr);
              newp = newp.divideAndRemainder(Polynomial.of(rr.negate(), Expression.ONE), "throw").quotient;
            }
            //TODO: details - ?
            //var roots1 = newp.getroots();
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
    }

    if (np.getDegree() === 4) {
      // https://en.wikipedia.org/wiki/Quartic_function
      // 1. coverting to a depressed quartic:
      var a_4 = np.getCoefficient(4);
      var p = np.scale(a_4.inverse());
      var b = p.getCoefficient(3);
      var y = new Expression.Symbol('$y');
      var substitute = function (y) {
        return y.subtract(b.divide(Expression.Integer.fromNumber(4)))
      };
      var e = p.calcAt(substitute(y));
      var sp = Polynomial.toPolynomial(e.getNumerator(), y);//TODO: ?
      var p = sp.getCoefficient(2);
      var q = sp.getCoefficient(1);
      var r = sp.getCoefficient(0);
      //var Q = nthRootInternal(3, );

      //var root = y0.subtract(b.divide(Expression.Integer.fromNumber(4)));
      
      // https://en.wikipedia.org/wiki/Quartic_function#Descartes'_solution
      const pU = Polynomial.of(q._pow(2).negate(), p._pow(2).subtract(Expression.TWO.add(Expression.TWO).multiply(r)), Expression.TWO.multiply(p), Expression.ONE);
      //TODO: when depressed cubic equation (?)
      //TODO: only one cubic equation root (?)
      const pURoots = pU.getCoefficient(1).equals(Expression.ZERO) && pU.getCoefficient(2).equals(Expression.ZERO) || pU.getCoefficient(0).equals(Expression.ZERO) ? pU.getroots() : [pU.doRationalRootTest() || Expression.ZERO];
      for (const U of pURoots) {
        const u = nthRootInternal(2, U);
        if (!U.equals(Expression.ZERO) && u != null) {
          const s = u.negate();
          const t = p.add(u._pow(2)).add(q.divide(u)).divide(Expression.TWO);
          const v = p.add(u._pow(2)).subtract(q.divide(u)).divide(Expression.TWO);
          //TODO: details (factorization of two quadratic)
          const p1 = Polynomial.of(t, s, Expression.ONE);
          const p2 = Polynomial.of(v, u, Expression.ONE);
          var p12 = Polynomial.of(Expression.ONE);
          for (const root of p1.getroots()) {
            const s = substitute(root);
            roots.push(s);
            p12 = p12.multiply(Polynomial.of(s.negate(), Expression.ONE));
          }
          for (const root of p2.getroots()) {
            const s = substitute(root);
            roots.push(s);
            p12 = p12.multiply(Polynomial.of(s.negate(), Expression.ONE));
          }
          if (p12.getDegree() > 0) {
            np = np.divideAndRemainder(p12).quotient;
            if (callback != undefined) {
              callback({content: content, roots: roots, newPolynomial: np, type: "solveQuarticEcuation"});//?
            }
            continueWithNewPolynomial(roots, np);
            return roots;
          }
        }
      }
      //debugger;
      //console.log(pURoots);
    }

    if (!np.hasIntegerCoefficients() && np.getDegree() > 0) {
      //?new
      var variable = new Expression.Symbol('$$');
      var e = np.calcAt(variable);
      var c = Expression.getConjugate(e);
      if (c != null) {
        var result = [];
        var newP = null;
        var ceRoots = Polynomial.toPolynomial(c.multiply(e), variable).getSquareFreePolynomial().getroots(function (x) {
          newP = x.newPolynomial;
        });//TODO: details - ?
        for (var i = 0; i < ceRoots.length; i += 1) {
          var root = ceRoots[i];
          if (np.calcAt(root).equals(Expression.ZERO)) {
            roots.push(root);
            //TODO: this also filters out duplicate roots (?)
            //TODO: factorization may be not good, so this also needed to calc polynomial to continue the factorization
            np = np.divideAndRemainder(Polynomial.of(root.negate(), Expression.ONE)).quotient;//TODO: optimize somehow - ?
          } else {
            //TODO:?
            console.debug(root);
          }
        }
        if (ceRoots.length > 0) {
          if (callback != undefined) {
            //TODO: better details, t = sqrt(3), show the polynomial, ...
            callback({content: content, roots: roots, newPolynomial: np, type: "multiplyByConjugates", t: c});//?
          }
          continueWithNewPolynomial(roots, np);
        }
        return roots;
      }
    }

    if (np.getDegree() >= 0) {
      //TODO: fix
      if (typeof hit === "function") {
        hit({getroots: {other: np.getDegree()}});
      }
    }

    return roots;
  };

  Polynomial.prototype.derive = function () {
    var newData = new PolynomialData(this.a.length);
    for (var i = 0; i < this.a.length; i += 1) {
      var n = this.a.degree(i);
      if (n >= 1) {
        var c = this.a.coefficient(i);
        newData.add(n - 1, Expression.Integer.fromNumber(n).multiply(c));
      }
    }
    return new Polynomial(newData);
  };

  Polynomial.prototype.getSquareFreePolynomial = function () {
    //TODO: remove (it is not good for performance of the factoring)
    return this.divideAndRemainder(this.squareFreeFactors().a0).quotient;
  };

  // f = a_1 * a_2**2 * a_3**3 * ... * a_n**n
  // a1 = a_1
  // a0 = a_2**1 * a_3**2 * ... * a_n**(n-1)
  // returns factor a1 - square free factor, a0 - a factor where the degree of every coprime factor is less by one
  Polynomial.prototype.squareFreeFactors = function () {
    // https://en.wikipedia.org/wiki/Square-free_polynomial
    var p = this;
    var zero = 0;
    while (p.getCoefficient(0).equals(Expression.ZERO)) {
      p = p.divideAndRemainder(Polynomial.of(Expression.ONE).shift(1)).quotient;
      zero += 1;
    }
    if (p.getDegree() !== 0) {
      var f = p;
      var d = f.derive();
      var a0 = Polynomial.polynomialGCD(f, d);
      if (a0.getDegree() !== 0) {
        if (Expression.isConstant(a0.getLeadingCoefficient())) {//TODO: ?
          a0 = a0.scale(a0.getLeadingCoefficient().inverse());//TODO: ?
        }
        a0 = a0.scale(a0.getContent().inverse());//?
        var b1 = f.divideAndRemainder(a0).quotient;
        var g1 = Polynomial.polynomialGCD(b1, a0);
        var a1 = b1.divideAndRemainder(g1).quotient;
        return {a1: a1.shift(zero === 1 ? 1 : 0), a0: a0.shift(zero > 1 ? zero - 1 : 0)};
      }
    }
    return {a1: p.shift(zero === 1 ? 1 : 0), a0: Polynomial.of(Expression.ONE).shift(zero > 1 ? zero - 1 : 0)};
  };

  export default Polynomial;

Polynomial.prototype.factorize = function () {
  if (this.getDegree() !== 3 && !this._canBeFactored()) {//TODO: ?
    //TODO: details - ?
    return undefined;
  }
  if (this.getCoefficient(0).equals(Expression.ZERO)) {//?
    return Polynomial.of(Expression.ZERO, Expression.ONE);
  }
  var content = this.getContent();//?TODO: ?
  if (!content.equals(Expression.ONE) && !content.equals(Expression.ONE.negate())) {
    throw new RangeError();
  }
  var tmp = this.squareFreeFactors();
  if (tmp.a0.getDegree() > 0) {
    //TODO: ?
    if (tmp.a1.getDegree() > 0) {
      return tmp.a1;
    }
    return tmp.a0;
  }
  if (this.getDegree() < 2) {
    return null;
  }
  //!new
  if (this.isEven()) {
    var f = this.subs(x => x.squareRoot()).factorize();
    if (f != undefined) {
      return f.subs(x => x._pow(2));
    }
    if (this.hasIntegerCoefficients() && !(this.getCoefficient(0).abs().squareRoot() instanceof Expression.Integer)) {
      return null;
    }
    //TODO: it should be factored into a product of (ax**(n/2)+...+bx+c) and (ax**(n/2)+...-bx+c)
    //if (this._getFactorByKroneckersMethod() != undefined) {
    //  debugger;
    //}
  }
if (this.getDegree() < 4 || !this.hasIntegerCoefficients()) {//? avoid Polynomial#getZeros(...)
  var root = this.doRationalRootTest();
  if (root != null) {
    return Polynomial.of(root.getNumerator(), root.getDenominator().negate());
  }
  if (this.getDegree() < 4) {
    // https://math.stackexchange.com/questions/1138416/how-do-i-show-a-cubic-polynomial-does-not-factorise#answer-1138428
    return null;
  }
}
  //!
  //TODO: performance
  var isSmall = function (np) {
    if (np.getDegree() >= 12) {
      return false;
    }
    var c = 1;
    for (var i = 0; i <= np.getDegree(); i += 1) {
      var k = np.getCoefficient(i);
      //TODO: ilog2(k)
      if (!k.equals(Expression.ZERO)) {
        c *= k instanceof Expression.Integer ? Math.floor(Math.log2(Math.abs(k.toNumber()) + 0.5)) + 1 : 1;
      }
    }
    return c <= 32 * 4 * 1024;
  };
  var np = this;
  //console.time('Kronecker\'s method');
  var g = undefined;
  //TODO: ?
  if (!np.hasIntegerCoefficients()) {
    //TODO: ?
    return np._getFactorByKroneckersMethod();
  }
  var realRootsCount = np.numberOfRoots();
  console.assert(realRootsCount <= np.getDegree());
  //TODO: try "np._factorByUsingZeros(2)" - ?
    //g = np._factorByUsingZeros(4, false);//?
  if (realRootsCount > 0) {
    var start1 = Date.now();
    var maxDegree = np.getDegree() >= 32 ? 3 : (np.getDegree() >= 24 ? 5 : (np.getDegree() >= 16 ? 7 : 8));
    var g2 = np._factorByUsingZeros(maxDegree, false, np.isEven() && np.hasIntegerCoefficients() ? np.getDegree() / 2 : 1);//?
    var end1 = Date.now();
    g = g2;
    if (end1 - start1 > 100) {
      //debugger;
    }
    if (g != undefined) {
      return g;
    }
    if (g == undefined && realRootsCount >= np.getDegree() - 2 && realRootsCount <= maxDegree) {
      //debugger;
      return g;//?
    }
    if (g == undefined && realRootsCount === np.getDegree() && maxDegree >= Math.floor(np.getDegree() / 2)) {
      //debugger;
      return g;//?
    }
  }
  if (isSmall(np) || realRootsCount === 0) {
    try {
      g = np._getFactorByKroneckersMethod(); // too slow
    } catch (error) {
      if (error.message !== "TOO SLOW") {
        throw error;
      }
      g = np._factorByLLL();//TODO: ?
    }
  } else if (np.getDegree() <= 6 && realRootsCount >= 4 || np.getDegree() <= 4 && realRootsCount >= 2) {
    g = np._factorByUsingZeros();//?
  } else {
    g = np._factorByLLL({evenFlag: true});
  }
  //console.timeEnd('Kronecker\'s method');
  return g;
};

//TODO: remove - ?
Polynomial.prototype._factorByUsingZeros = function (maxDegree = 1/0, tryComplex = true, minDegree = 1) {
  // https://lowrey.me/es6-javascript-combination-generator/
  const combinations = function(elements, length) {
    let i = 0;
    let remaining = null;
    const tmp = {
      next: function () {
        for (; i < elements.length;) {
          if (length === 1) {
            const result = [elements[i]];
            i += 1;
            return {value: result, done: false};
          } else {
            if (remaining == null) {
              remaining = combinations(elements.slice(i + 1, elements.length), length - 1)[globalThis.Symbol.iterator]();
            }
            const next = remaining.next();
            if (next.done) {
              remaining = null;
              i += 1;
            } else {
              const result = [elements[i]].concat(next.value);
              return {value: result, done: false};
            }
          }
        };
        return {value: undefined, done: true};
      }
    };
    tmp[globalThis.Symbol.iterator] = function () {
      return this;
    };
    return tmp;
  };
  //var maxDegree = 4;
  //var zeros = "-5.274,-2.730,-1.676,-0.868,0.868,1.676,2.730,5.274".split(',');
  var zeros = this.getZeros(Math.ceil(2 * this._log2hypot()));//TODO: complex - ?
  if (zeros.length < Math.min(maxDegree, this.getDegree() - 2) && tryComplex) {//TODO: ???
      zeros = this.getZeros(undefined, true);
  }
  var c1 = 0;
  var c2 = 0;
  for (var degree = minDegree; degree <= maxDegree && degree <= this.getDegree() / 2; degree += 1) {
    for (var combination of combinations(zeros, degree)) {
      var a = this.getLeadingCoefficient();
      // guess = a(x-z_0)(x-z_1)(x-z_2)(x-z_3) = ... with rounded to integer coefficients
      var product = a;
      for (var zero of combination) {
        //TODO: fix -> instanceOf - ?
        product = product.multiply(zero.e != null ? zero.e : zero);
      }
      var test = toDecimalStringInternal(product, {fractionDigits: 3});
      //debugger;
      c1 += 1;
      if (!/\.0*[1-9]/.test(test)) { // product is an integer
        c2 += 1;
        var guess = Polynomial.of(a);
        for (var zero of combination) {
          //guess = guess.multiply(Polynomial.of(ExpressionParser.parse(zero.toString()).negate(), Expression.ONE));
          guess = guess.multiply(Polynomial.of(zero.e != null ? zero.e.negate() : zero.negate(), Expression.ONE));
        }
        var tmp2 = guess;
        guess = guess.map(function (coefficient) {
          return ExpressionParser.parse(toDecimalStringInternal(coefficient, {fractionDigits: 0}));
        });
        guess = guess.scale(guess.getContent().inverse());
        var tmp = this.divideAndRemainder(guess, "undefined");
        if (tmp != undefined && tmp.remainder.equals(Polynomial.ZERO)) {
          console.info(guess + '', c1, c2, this + '');
          return guess;
        }
      }
    }
  }
  console.log(c1, c2);
  return null;
};
