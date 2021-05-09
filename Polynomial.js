/*global hit*/

  import Expression from './Expression.js';
  import Matrix from './Matrix.js';
  import Heap from './Heap.js';
  import primeFactor from './primeFactor.js';
  //import bitLength from './bitLength.js';

  const nextPrime = primeFactor._nextPrime;

  //const freeze = Object.freeze; - too slow
  const freeze = x => x;

  function PolynomialData(length) {
    if (length !== Math.floor(length)) {
      throw new RangeError();
    }
    this.degrees = new Array(length);
    this.coefficients = new Array(length);
    this.size = 0; // public
  }
  PolynomialData.prototype.add = function (degree, coefficient) {
    if (degree < 0 || degree > Number.MAX_SAFE_INTEGER || Math.floor(degree) !== degree) {
      throw new RangeError("NotSupportedError");
    }
    const k = this.size;
    if (k > 0 && degree >= this.degrees[k - 1]) {
      throw new RangeError();
    }
    if (k >= this.degrees.length || k >= this.coefficients.length || k < 0) {
      throw new RangeError();
    }
    this.degrees[k] = degree;
    this.coefficients[k] = coefficient;
    this.size = k + 1;
  };
  PolynomialData.prototype.degree = function (i) {
    return this.degrees[i];
  };
  PolynomialData.prototype.coefficient = function (i) {
    return this.coefficients[i];
  };
  PolynomialData.prototype.trim = function () {
    if (this.size !== this.degrees.length || this.size !== this.coefficients.length) {
      this.degrees = this.degrees.slice(0, this.size);
      this.coefficients = this.coefficients.slice(0, this.size);
    }
    return this;
  };
  PolynomialData.prototype.getCoefficient = function (degree) {
    var from = 0;
    var to = this.size;
    while (from < to) {
      var middle = from + Math.floor((to - from) / 2);
      var y = this.degrees[middle];
      if (y > degree) {
        from = middle + 1;
      } else if (y < degree) {
        to = middle;
      } else {
        return this.coefficients[middle];
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
    for (var i = arguments.length - 1; i >= 0; i -= 1) {
      var a = arguments[i];
      if (!a.equals(Expression.ZERO)) {
        newData.add(i, a);
      }
    }
    return new Polynomial(newData);
  };
  Polynomial.from = function (array) {
    var newData = new PolynomialData(array.length);
    for (var i = array.length - 1; i >= 0; i -= 1) {
      var degree = i;
      var coefficient = array[i];
      if (!coefficient.equals(Expression.ZERO)) {
        newData.add(degree, coefficient);
      }
    }
    return new Polynomial(newData);
  };

  Polynomial.ZERO = Polynomial.of();

  Polynomial.prototype.getDegree = function () {
    return this.a.size === 0 ? -1 : this.a.degree(0);
  };
  Polynomial.prototype.getCoefficient = function (degree) {
    return this.a.getCoefficient(degree);
  };
  Polynomial.prototype.getLeadingCoefficient = function () {
    return this.a.size === 0 ? Expression.ZERO : this.a.coefficient(0);
  };

  Polynomial.prototype.map = function (mapFunction) {//?
    const length = mapFunction.length;
    const newData = new PolynomialData(this.a.size);
    for (let i = 0; i < this.a.size; i += 1) {
      const degree = this.a.degree(i);
      const coefficient = this.a.coefficient(i);
      const c = length === 1 ? mapFunction(coefficient) : mapFunction(coefficient, degree);
      if (!c.equals(Expression.ZERO)) {
        newData.add(degree, c);
      }
    }
    return new Polynomial(newData);
  };

  Polynomial.prototype.equals = function (p) {
    if (p === Expression.ZERO) {
      return this.a.size === 0;
    }
    if (this.a.size !== p.a.size) {
      return false;
    }
    var i = 0;
    var j = 0;
    while (i < this.a.size && j < p.a.size) {
      if (this.a.degree(i) !== p.a.degree(j) || !this.a.coefficient(i).equals(p.a.coefficient(j))) {
        return false;
      }
      i += 1;
      j += 1;
    }
    return true;
  };

  Polynomial.prototype.add = function (p) {
    if (p.a.size === 0) {
      return this;
    }
    if (this.a.size === 0) {
      return p;
    }
    var newData = new PolynomialData(this.a.size + p.a.size);
    let i = -1;
    let x = -2;
    let xc = Expression.ZERO;
    let j = -1;
    let y = -2;
    let yc = Expression.ZERO;
    while (x !== -1 || y !== -1) {
      const d = Math.max(x, y);
      const flag = x !== y;
      const c = x > y ? xc : (y > x ? yc : (d === -2 ? Expression.ZERO : xc.add(yc)));
      if (x === d) {
        i += 1;
        if (i < this.a.size) {
          x = this.a.degree(i);
          xc = this.a.coefficient(i);
        } else {
          x = -1;
          xc = Expression.ZERO;
        }
      }
      if (y === d) {
        j += 1;
        if (j < p.a.size) {
          y = p.a.degree(j);
          yc = p.a.coefficient(j);
        } else {
          y = -1;
          yc = Expression.ZERO;
        }
      }
      if (flag || (d !== -2 && !c.equals(Expression.ZERO))) {
        newData.add(d, c);
      }
    }
    return new Polynomial(newData);
  };

  function multiplyInternal(A, B, fromLeft) {
    var npmp1 = A.getDegree() + B.getDegree() + 1;
    if ((A.a.size + B.a.size + Math.min(A.a.size, B.a.size)) * 4 >= npmp1) { // "dense"
      var result = new Array(npmp1).fill(Expression.ZERO);
      for (var i = 0; i < A.a.size; i += 1) {
        var d = A.a.degree(i);
        var c = A.a.coefficient(i);
        for (var j = 0; j < B.a.size; j += 1) {
          var bd = B.a.degree(j);
          var bj = B.a.coefficient(j);
          var degree = d + bd;
          var coefficient = fromLeft ? c.multiply(bj) : bj.multiply(c);
          result[degree] = result[degree].add(coefficient);
        }
      }
      return Polynomial.from(result);
    } else {
      var result = new FastAdditionPolynomial();
      for (var i = 0; i < A.a.size; i += 1) {
        var d = A.a.degree(i);
        var c = A.a.coefficient(i);
        result.add(c, d, B, fromLeft);
      }
      return result.toPolynomial();
    }
  }
  Polynomial.prototype.multiply = function (p) {
    if (this.a.size === 0 || p.a.size === 0) {
      return Polynomial.ZERO;
    }
    if (p.a.size === 1) {
      return this.shift(p.a.degree(0)).scale(p.a.coefficient(0));
    }
    if (this.a.size <= p.a.size) {
      return multiplyInternal(this, p, true);
    } else {
      return multiplyInternal(p, this, false);
    }
  };

  Polynomial.prototype.shift = function (n) { // *= x**n, n >= 0
    if (!(n >= 0)) {
      throw new TypeError();
    }
    if (n === 0) {
      return this;
    }
    var newData = new PolynomialData(this.a.size);
    for (var i = 0; i < this.a.size; i += 1) {
      newData.add(this.a.degree(i) + n, this.a.coefficient(i));
    }
    return new Polynomial(newData);
  };

  //note: no Map is needed
  function FastAdditionPolynomial() {
    // see http://www.cecm.sfu.ca/~mmonagan/teaching/TopicsinCA11/johnson.pdf
    // see https://en.wikipedia.org/wiki/K-way_merge_algorithm#Heap
    this.maxHeap = new Heap((a, b) => b.degree - a.degree);
    this.degree = -1;
    this.leadingCoefficient = Expression.ZERO;
  }
  FastAdditionPolynomial.prototype.add = function (scale, shift, polynomial, fromLeft) {
    let i = 0;
    var iterator = {
      next: function () {
        if (i === polynomial.a.size) {
          return undefined;
        }
        const d = polynomial.a.degree(i);
        const c = polynomial.a.coefficient(i);
        i += 1;
        return freeze({
          degree: d + shift,
          coefficient: (c === Expression.ONE ? scale : (fromLeft ? scale.multiply(c) : c.multiply(scale))),
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
      while (this.maxHeap.size() > 0 && (this.degree === this.maxHeap.peek().degree || this.leadingCoefficient.equals(Expression.ZERO))) {
        var tmp = this.maxHeap.peek();
        this.leadingCoefficient = this.leadingCoefficient.add(tmp.coefficient);
        this.degree = tmp.degree;
        var next = tmp.iterator.next();
        if (next != undefined) {
          //this.maxHeap.pop();
          //this.maxHeap.push(next);
          this.maxHeap.replace(next);
        } else {
          this.maxHeap.pop();
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
      this.add(coefficient.negate(), degree, ONE, true);
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
    if (this.getDegree() < p.getDegree()) {
      return {quotient: Polynomial.ZERO, remainder: this};
    }
    var quotient = [];
    //var remainder = this;
    var remainder = new FastAdditionPolynomial();
    remainder.add(Expression.ONE, 0, this, true);
    const minusP = p.negate();
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
      quotient.push({degree: n, coefficient: q});
      //TODO: optimize - ?
      remainder.add(q, n, minusP, true);
      //if (!(q instanceof Expression.Integer)) {
      //  remainder = remainder.add(Polynomial.of(q).shift(n).multiply(minusP));
      //} else {
      //  remainder = remainder.add(minusP.shift(n).scale(q));
      //}
      if (remainder.getDegree() - p.getDegree() === n) {
        // to avoid the infite loop
        throw new TypeError("there is a some problem with the expression evaluation");//!
      }
    }
    remainder = remainder.toPolynomial();
    return {quotient: Polynomial.fromTerms(quotient), remainder: remainder};
  };

  Polynomial.pseudoRemainder = function (A, B) {
    var d = A.getDegree() - B.getDegree();
    // assertion
    if (d < 0) {
      throw new RangeError();
    }
    return A.scale(B.getLeadingCoefficient()._pow(d + 1)).divideAndRemainder(B, "throw").remainder;
  };

  Polynomial.polynomialGCD = function (a, b) {
    //!optimization 2021-03-15
    if (b.getDegree() > 0 && a.getDegree() > 0) {
      var i = 0;
      while (b.getCoefficient(i).equals(Expression.ZERO)) {
        i += 1;
      }
      var j = 0;
      while (a.getCoefficient(j).equals(Expression.ZERO)) {
        j += 1;
      }
      if (i !== 0 || j !== 0) {
        //TODO: optimize (?)
        return Polynomial.polynomialGCD(a.divideAndRemainder(Polynomial.of(Expression.ONE).shift(j), "throw").quotient, b.divideAndRemainder(Polynomial.of(Expression.ONE).shift(i), "throw").quotient).shift(Math.min(i, j));
      }
    }
    //!
    //TODO: fix (place condition for degrees earlier - ?)
    if (a.getDegree() < b.getDegree()) {
      //!!!
      var tmp = a;
      a = b;
      b = tmp;
    }
    if (b.equals(Polynomial.ZERO)) {
      return a;
    }

    var contentA = a.getContent();
    var contentB = b.getContent();
    var ppA = a.divideAndRemainder(Polynomial.of(contentA), "throw").quotient;
    var ppB = b.divideAndRemainder(Polynomial.of(contentB), "throw").quotient;
    var ppGCD = gcdOfPrimitivePolynomials(ppA, ppB);
    var contentGCD = contentA.gcd(contentB);
    return Polynomial.of(contentGCD).multiply(ppGCD);
  };

  function gcdOfPrimitivePolynomials(A, B) {
    console.assert(A.getDegree() >= B.getDegree());
    //if (B.getDegree() === 1) {
    //  return A.divideAndRemainder(B).remainder.equals(Polynomial.ZERO) ? B : Polynomial.of(Expression.ONE);
    //}
    if (A.hasIntegerCoefficients() && B.hasIntegerCoefficients() && B.getDegree() > 2) {
      if (A.isDivisibleBy(B)) {
        return B;
      }
      return gcdByModularAlgorithm(A, B);
    }
    //TODO: 
    if (A._hasIntegerLikeCoefficients() && B._hasIntegerLikeCoefficients() && B.getDegree() > 2) {//TODO: why isn't it working or slow for other cases - ?
      return gcdUsingRemainderSequence(A, B, Polynomial._subresultantPseudoRemainderSequence).primitivePart();
    }
    return gcdUsingRemainderSequence(A, B, Polynomial._primitivePseudoRemainderSequence);
  }

  Polynomial.debug = function (R) {
  };

  function gcdUsingRemainderSequence(A, B, generator) {
    var lastNonzeroRemainder = null;
    for (var tmp of generator(A, B)) {
      lastNonzeroRemainder = tmp.B;
      Polynomial.debug(tmp.R);
    }
    return lastNonzeroRemainder
  }

  Polynomial._gcdUsingRemainderSequence = gcdUsingRemainderSequence;

  function ChineseRemainderTheorem(m1, m2) {
    // https://en.wikipedia.org/wiki/Chinese_remainder_theorem#Case_of_two_moduli
    // x = r1 (mod m1)
    // x = r2 (mod m2)
    //TODO: fix
    const c = m1.modInverse(m2);
    //const m1m2 = m1.multiply(m2);
    return function (r1, r2) {
      return r1.add((r2.subtract(r1)).multiply(c).modulo(m2).multiply(m1));
    };
  }

  /*
   const i = Expression.Integer.fromNumber;
   console.assert(ChineseRemainderTheorem(i(36), i(1), i(119), i(23)).equals(i(2416)));
  */

  function ChineseRemainderTheoremForPolynomialCoefficients(p1, p2, m1, m2) {
    var n = Math.max(p1.getDegree(), p2.getDegree());
    const solution = ChineseRemainderTheorem(m1, m2);
    var c = new Array(n);
    for (var i = 0; i <= n; i += 1) {
      var r1 = p1.getCoefficient(i);
      var r2 = p2.getCoefficient(i);
      c[i] = solution(r1, r2);
    }
    return Polynomial.from(c);
  }

  function gcdByModularAlgorithm(a, b) {
    // https://www3.risc.jku.at/education/courses/ws2011/ca/3-gcd.pdf
    console.assert(a.hasIntegerCoefficients() && a.getContent().equals(Expression.ONE) && a.getDegree() > 0);
    console.assert(b.hasIntegerCoefficients() && b.getContent().equals(Expression.ONE) && b.getDegree() > 0);
    if (!a.getCoefficient(0).equals(Expression.ZERO) &&
        !b.getCoefficient(0).equals(Expression.ZERO) &&
        a.getLeadingCoefficient().gcd(b.getLeadingCoefficient()).abs().compareTo(a.getCoefficient(0).gcd(b.getCoefficient(0)).abs()) > 0) {//?
      var i = function (polynomial) {
        var coefficients = new Array(polynomial.a.size);
        for (var j = 0; j < polynomial.a.size; j += 1) {
          coefficients[j] = {degree: polynomial.getDegree() - polynomial.a.degree(j), coefficient: polynomial.a.coefficient(j)};
        }
        coefficients.reverse();
        return Polynomial.fromTerms(coefficients).primitivePart();
      };
      return i(gcdByModularAlgorithm(i(a), i(b)));
    }
    const d = a.getLeadingCoefficient().gcd(b.getLeadingCoefficient());
    const maxGCDDegree = Math.min(a.getDegree(), b.getDegree());
    const bound = function (a, b, maxGCDDegree) {
      //TODO: not necessary to use power of two
      const maxGCDLeadingCoefficient = a.getLeadingCoefficient().gcd(b.getLeadingCoefficient()).abs();
      const logarithmOfCoefficientBound = Math.min(a._log2OfBoundForCoefficientsOfFactor(maxGCDDegree, maxGCDLeadingCoefficient), b._log2OfBoundForCoefficientsOfFactor(maxGCDDegree, maxGCDLeadingCoefficient));
      //TODO: multiply by gcd of leading coefficients (?)
      return Expression.TWO._pow(1 + Math.ceil(logarithmOfCoefficientBound));
    };
    let M = bound(a, b, maxGCDDegree);
    //if (M === 1 / 0) {
    //  return gcdUsingPrimitivePseudoRemainderSequence(a, b);//TODO: !!!
    //}
    let p = Expression.Integer.fromBigInt(Math.min(M.toNumber(), Math.sqrt((Number.MAX_SAFE_INTEGER + 1) / 2))); // TODO: should we divide on n - ?
    //let p = Expression.Integer.fromNumber(3);//TODO: remove
    let counter = 0;
    let g = null;
    while (g == null) {
      let P = Expression.Integer.fromNumber(0);
      while (P.compareTo(M) < 0) {
        counter += 1;
        if (counter > 50 && p._pow(counter - 50).compareTo(M) >= 0) {
          throw new TypeError("!!!");
        }
        do {
          p = Expression.Integer.fromBigInt(nextPrime(p.toBigInt()));
        } while (d.remainder(p).equals(Expression.ZERO));
        const cp = Polynomial._gcdOfPolynomialsOverFiniteField(a, b, p);
        console.assert(cp.getLeadingCoefficient().equals(Expression.ONE));
        const gp = cp.scale(d.remainder(p)).mod(p); //TODO: should it do .mod(p) ? as the book does not tell this
        if (gp.getDegree() < 1) {
          return Polynomial.of(Expression.ONE);
        }
        if (gp.getDegree() > (2 / 3) * a.getDegree()) { // see Donald Knuth's book
          //TODO: ?
          //return gcdOfPrimitivePolynomials(a, b);
        }
        if (g == null || gp.getDegree() < g.getDegree()) {//TODO: why ?
          if (g != null) {
            //debugger;
          }
          g = gp.mod2(p);
          P = p;
          M = bound(a, b, gp.getDegree()); //!?
        } else if (gp.getDegree() === g.getDegree()) {//TODO: why ?
          var oldG = g;
          g = ChineseRemainderTheoremForPolynomialCoefficients(g, gp, P, p);
          P = P.multiply(p);
          g = g.mod2(P);
          //does not work (?), the number of iteration is small anyway (?)
          if (P.compareTo(M) < 0 && oldG.equals(g)) {
            // optimization from the book:
            // "Whenever g remains unchanged for a series of iterations through the while–loop, we might apply the test in step (5) and exit if the outcome is positive."
            g = g.primitivePart();
            // TODO: scale a and b instead (?) for the test
            if (a.isDivisibleBy(g) && b.isDivisibleBy(g)) {
              //console.log('counter:', counter);
              return g;
            }
            g = oldG; //TODO: ???
          } else {
            //debugger;
          }
        }
        //console.debug(counter, g.primitivePart().toString(), P.toString());
      }
      g = g.primitivePart();
      if (!a.isDivisibleBy(g) || !b.isDivisibleBy(g)) {
        g = null;
      }
    }
    //console.debug('counter:', counter);
    return g;
  }

  Polynomial._gcdByModularAlgorithm = gcdByModularAlgorithm;

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
    for (var i = 0; i < this.a.size; i += 1) {
      var degree = this.a.degree(i);
      var coefficient = this.a.coefficient(i);
      if (!n.equals(Expression.ZERO)) {
        n = Expression.pow(point, lastDegree - degree).multiply(n).add(coefficient);
      } else {
        n = coefficient;
      }
      lastDegree = degree;
    }
    if (!n.equals(Expression.ZERO)) {
      n = Expression.pow(point, lastDegree - 0).multiply(n);
    }
    return n;
  };

  Polynomial.prototype.getContent = function () {
    if (this.a.size === 0) {
      return Expression.ONE;
    }
    var denominator = Expression.ONE;
    var numerator = Expression.ZERO;
    for (var i = 0; i < this.a.size; i += 1) {
      var y = this.a.coefficient(i);
      denominator = denominator.lcm(y.getDenominator());
      numerator = numerator.gcd(y.getNumerator());
    }
    var c = numerator.divide(denominator);
    var x = this.a.coefficient(0);
    return x.isNegative() && !numerator.isNegative() || numerator.isNegative() && !x.isNegative() ? c.negate() : c;
  };

  // add, multiply, divideAndRemainder

  Polynomial.prototype.negate = function () {
    //TODO: fix
    return this.map(coefficient => coefficient.negate());
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
    return this.map(coefficient => coefficient.multiply(x));
  };

  Polynomial.fromTerms = function (terms) {
    var newData = new PolynomialData(terms.length);
    for (var i = 0; i < terms.length; i += 1) {
      var term = terms[i];
      newData.add(term.degree, term.coefficient);
    }
    return new Polynomial(newData);
  };

  Polynomial.toPolynomial = function (e, v) {
    if (e instanceof Expression.Division) {
      throw new RangeError();
    }
    var terms = Expression.getCoefficients(e, v);
    return Polynomial.fromTerms(terms);
  };

  Polynomial.prototype.toExpression = function (variableSymbol) {
    var result = undefined;
    for (var i = 0; i < this.a.size; i += 1) {
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
      const toInteger = function (zero, scale) {
        //TODO: ?
        const fraction = zero instanceof Expression.ExpressionPolynomialRoot ? zero.root._root.toDecimal(0).a : zero.toDecimal().a;
        return scale.multiply(fraction.getNumerator()).truncatingDivide(fraction.getDenominator());
      };
      var zeros = np.getZeros(primeFactor._bitLength(an.abs().toBigInt())); // note: for fractionDigits (!)
      for (var i = 0; i < zeros.length; i += 1) {
        var zero = zeros[i];
        if (i === 0 || zero !== zeros[i - 1]) {
          var candidate = zero.root != null ? toInteger(zero, an).divide(an) : (zero.getNumerator() instanceof Expression.Integer && zero.getDenominator() instanceof Expression.Integer ? zero : Expression.ZERO);
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

  Polynomial.prototype._testCoefficients = function (f) {
    for (var i = 0; i < this.a.size; i += 1) {
      if (!f(this.a.coefficient(i))) {
        return false;
      }
    }
    return true;
  };
  Polynomial.prototype.hasIntegerCoefficients = function () {
    return this._testCoefficients(c => c instanceof Expression.Integer);
  };
  Polynomial.prototype.hasComplexCoefficients = function () {
    return this._testCoefficients(c => (c instanceof Expression.Complex) || (c instanceof Expression.Integer));
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
    return this._testCoefficients(c => isIntegerLike(c));
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
    return this.getGCDOfTermDegrees() % 2 === 0;
  };

  Polynomial.prototype._factorByKroneckersMethod = function () {
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
          if (n.abs().toNumber() > 3317044064679887385961981) {
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
            g = g.primitivePart();
            if (np.isDivisibleBy(g)) {
              return g;
            }
          }
        }
      }
    }
    return undefined;
  };

  var counter = 0;//TODO: remove

  Polynomial.prototype.getGCDOfTermDegrees = function () {
    var gcd = function (a, b) {
      return b === 0 ? a : gcd(b, a % b);
    };
    var g = this.getDegree();
    for (var i = 1; i <= this.getDegree() && g >= 2; i += 1) {
      if (!this.getCoefficient(i).equals(Expression.ZERO)) {
        g = gcd(g, i);
      }
    }
    return g;
  };

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
        //TODO: !?
        return undefined;//?
        //return x._nthRoot(n);
      }
      if (x instanceof Expression.ExpressionPolynomialRoot) {
        //return undefined;//?
        return x._nthRoot(n);
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
      var g = np.getGCDOfTermDegrees();
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
              // https://en.wikipedia.org/wiki/Root_of_unity
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
            console.assert(np.subtract(quotient.multiply(d)).map(c => c.simplifyExpression()).getDegree() < 0);
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
      var h = np._getShiftToDepressed();
      if (!h.equals(Expression.ZERO)) {
        var f = x => x.subtract(h);
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

    if (!np.hasIntegerCoefficients() && np.getDegree() > 2) {
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
    var newData = new PolynomialData(this.a.size);
    for (var i = 0; i < this.a.size; i += 1) {
      var n = this.a.degree(i);
      var c = this.a.coefficient(i);
      if (n >= 1) {
        var nc = Expression.Integer.fromNumber(n).multiply(c);
        if (!nc.equals(Expression.ZERO)) {
          newData.add(n - 1, nc);
        }
      }
    }
    return new Polynomial(newData);
  };

  Polynomial.prototype.getSquareFreePolynomial = function () {
    //TODO: remove (it is not good for performance of the factoring)
    return this.divideAndRemainder(this.squareFreeFactors().a0).quotient;
  };
  Polynomial.prototype.isSquareFreePolynomial = function () {
    return this.squareFreeFactors().a0.equals(Polynomial.of(Expression.ONE));
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
      p = p.divideAndRemainder(Polynomial.of(Expression.ONE).shift(1), "throw").quotient;
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
        if (f.hasIntegerCoefficients()) {//?
          a0 = a0.scale(a0.getContent().inverse());//?
        }
        var b1 = f.divideAndRemainder(a0, "throw").quotient;
        var g1 = Polynomial.polynomialGCD(b1, a0);
        var a1 = b1.divideAndRemainder(g1, "throw").quotient;
        return {a1: a1.shift(zero === 1 ? 1 : 0), a0: a0.shift(zero > 1 ? zero - 1 : 0)};
      }
    }
    return {a1: p.shift(zero === 1 ? 1 : 0), a0: Polynomial.of(Expression.ONE).shift(zero > 1 ? zero - 1 : 0)};
  };

  export default Polynomial;
  
Polynomial.prototype._findGoodSubstitution = function () {
  var g = this.getCoefficient(0).abs();
  for (var i = 0; i < this.getDegree(); i += 1) {
    g = g.gcd(this.getCoefficient(i).abs());
  }
  if (g.toNumber() > 2**53) {
    return null;//!!!
  }
  if (g.toNumber() === 1) {
    return null;
  }
  var p = g.primeFactor();
  var before = this._log2hypot();
  var after = this.subs(x => x.multiply(p));
  after = after.scale(after.getContent().inverse());
  after = after._log2hypot();
  if (after < before) {
    return p;
  }

  /*
  var g = this.getCoefficient(0).abs();
  for (var i = 1; i <= this.getDegree(); i += 1) {
    g = g.gcd(this.getCoefficient(i).abs());
  }
  g = g.abs();
  //while (!g.equals(Expression.ONE)) {
    //var p = g.primeFactor();
    //g = g.truncatingDivide(p);
  //}
  var g = Expression.getConstant(this.getCoefficient(0).abs().gcd(g._pow(this.getDegree()))._nthRoot(this.getDegree()));
  for (var i = 1; i < this.getDegree(); i += 1) {
    var c = this.getCoefficient(i).gcd(g._pow(this.getDegree() - i));
    if (!c.equals(Expression.ZERO)) {
      g = g.gcd(Expression.getConstant(c.abs()._nthRoot(this.getDegree() - i)));
    }
  }
  if (!g.equals(Expression.ONE)) {
    //debugger;
    return g;
  }
  //TODO: denominator
  return null;*/
};

Polynomial.prototype._getShiftToDepressed = function () { // for testing (?)
  var n = this.getDegree();
  var a = this.getLeadingCoefficient();
  var b = this.getCoefficient(n - 1);
  var h = b.divide(Expression.Integer.fromNumber(n).multiply(a));
  //var f = x => x.subtract(h);
  //var depressed = this.subs(f); 
  //console.log(depressed);
  return h;
};

Polynomial.prototype.factorize = function () {
  //if (this.getDegree() !== 3 && !this._canBeFactored()) {//TODO: ?
    //TODO: details - ?
  //  return undefined;
  //}
  if (this.getCoefficient(0).equals(Expression.ZERO)) {//?
    return Polynomial.of(Expression.ZERO, Expression.ONE);
  }
  if (this.getDegree() === 3) {
    //console.log(this.toString());
  }
  var content = this.getContent();//?TODO: ?
  if (!content.equals(Expression.ONE) && !content.equals(Expression.ONE.negate())) {
    throw new RangeError();
  }
  if (this.getDegree() === 1) {
    return null;
  }
  if (this.getDegree() === 2) {
    const discriminant = function (p) {
      console.assert(p.getDegree() === 2);
      var a = p.getCoefficient(2);
      var b = p.getCoefficient(1);
      var c = p.getCoefficient(0);
      return b.multiply(b).subtract(Expression.Integer.fromNumber(4).multiply(a).multiply(c));
    };
    var d = discriminant(this);
    //TODO: not integer discriminants (?)
    if (d instanceof Expression.Integer && !d._isPerfectSquare()) {
      return null;
    }
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
  if (this.hasIntegerCoefficients() && false) {
    var k = this._findGoodSubstitution();
    if (k != null) {
      var newp = this.subs(x => x.multiply(k));
      newp = newp.scale(newp.getContent().inverse());
      var f = newp.factorize();
      if (f != undefined) {
        return f.subs(x => x.divide(k));
      }
      return null;//!?
    }
  }
  var n = this.getDegree();
  var b = this.getCoefficient(n - 1);
  if (!b.equals(Expression.ZERO)) {
    var a = this.getLeadingCoefficient();
    var f = x => x.subtract(b.divide(Expression.Integer.fromNumber(n).multiply(a)));
    var depressed = this.subs(f);
    if (!this.hasIntegerCoefficients() || depressed.isEven()) {//TODO: ?
      depressed = depressed.scale(depressed.getContent().inverse());//TODO: ?
      //debugger;
      var factor = depressed.factorize();
      if (factor != null) {
        return factor.subs(x => x.add(b.divide(Expression.Integer.fromNumber(n).multiply(a)))).primitivePart();
      }
      return null;
    }
  }
  //!
  //!new
  if (this.isEven()) {
    if (this.hasIntegerCoefficients() && (!this.getCoefficient(0).abs()._isPerfectSquare() || !this.getCoefficient(this.getDegree()).abs()._isPerfectSquare())) {
      var f = this.subs(x => x.squareRoot()).factorize();
      if (f != undefined) {
        return f.subs(x => x._pow(2));
      }
      return null;
    }
    //TODO: it should be factored into a product of (ax**(n/2)+...+bx+c) and (ax**(n/2)+...-bx+c)
    //if (this._factorByKroneckersMethod() != undefined) {
    //  debugger;
    //}
  }
  var g = this.getGCDOfTermDegrees();
  if (g > 2) {
    var f = this.subs(x => x._nthRoot(g)).factorize();
    if (f != undefined) {
      return f.subs(x => x._pow(g));
    }
    //?
  }
  //!
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
  var np = this;
  //console.time('Kronecker\'s method');
  //TODO: ?
  if (!np.hasIntegerCoefficients()) {
    return np._factorByKroneckersMethod();
  }
  return np._factorizeOverTheIntegers();
  //console.timeEnd('Kronecker\'s method');
};


  // https://lowrey.me/es6-javascript-combination-generator/
  Polynomial._combinations = function(elements, length) {
    let i = 0;
    let remaining = null;
    const tmp = {
      next: function () {
        for (; i < elements.length;) {
          if (length === 1) {
            const result = [];
            result.unshift(elements[i]);
            i += 1;
            return {value: result, done: false};
          } else {
            if (remaining == null) {
              remaining = Polynomial._combinations(elements.slice(i + 1, elements.length), length - 1)[globalThis.Symbol.iterator]();
            }
            const next = remaining.next();
            if (next.done) {
              remaining = null;
              i += 1;
            } else {
              const result = next.value;
              result.unshift(elements[i]);
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



Polynomial.prototype.primitivePart = function () {
  //return this.divideAndRemainder(Polynomial.of(this.getContent()), "throw").quotient;
  //TODO: test if this is faster:
  //var content = this.getContent();
  //return this.scale(content.getDenominator()).divideAndRemainder(Polynomial.of(content.getNumerator()), "throw").quotient;
  return this.scale(this.getContent().inverse());
};



Polynomial.prototype.modularInverse = function (m) {
  if (this.hasIntegerCoefficients() && m.hasIntegerCoefficients()) {
    //TODO: ???
  }
  //TODO: ?
  // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Polynomial_extended_Euclidean_algorithm
  var a = this;
  var oldR = m;
  var r = a.getDegree() >= m.getDegree() ? Polynomial.pseudoRemainder(a, m).primitivePart() : a; //? faster then swap (?)
  //var r = a;
  /*if (oldR.getDegree() < r.getDegree()) {
    var tmp = oldR;
    oldR = r;
    r = tmp;
  }*/
  var oldT = Polynomial.of(Expression.ZERO);
  var t = Polynomial.of(Expression.ONE);
  //var oldS = Polynomial.of(Expression.ONE);
  //var s = Polynomial.of(Expression.ZERO);
  //TODO: fix for different sequences (?)
  for (var tmp of Polynomial._subresultantPseudoRemainderSequence(oldR, r)) {
    //?
    var scale = r.getLeadingCoefficient()._pow(oldR.getDegree() - r.getDegree() + 1);
    var quotient = tmp.q;
    var α = tmp.α;
    var newT = oldT.scale(scale).subtract(quotient.multiply(t)).divideAndRemainder(Polynomial.of(α), "throw").quotient;
    [oldT, t] = [t, newT];
    //var newS = oldS.scale(scale).subtract(quotient.multiply(s)).divideAndRemainder(Polynomial.of(α), "throw").quotient;
    //[oldS, s] = [s, newS];
    var newR = tmp.R;
    [oldR, r] = [r, newR];
  }
  var gcd = oldR;
  //TODO: FIX
  if (gcd.getDegree() !== 0) {
    throw new TypeError();//?
  }
  oldT = oldT.primitivePart().scale(m.getLeadingCoefficient());
  var r = oldT.multiply(this).divideAndRemainder(m, "throw").remainder;
  if (r.getDegree() !== 0) {
    throw new TypeError();//?
  }
  oldT = oldT.scale(r.getCoefficient(0).inverse());
  //oldT = oldT.scale(gcd.getLeadingCoefficient().inverse());
  //oldS = oldS.primitivePart().scale(m.getLeadingCoefficient());
  //var r = oldS.multiply(this).divideAndRemainder(m, "throw").remainder;
  //if (r.getDegree() !== 0) {
  //  throw new TypeError();//?
  //}
  //oldS = oldS.scale(r.getCoefficient(0).inverse());
  //oldS = oldS.scale(gcd.getLeadingCoefficient().inverse());
  return oldT;
};
