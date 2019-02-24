/*jslint plusplus: true, vars: true, indent: 2 */

// Thanks to Eduardo Cavazos
// see also https://github.com/dharmatech/Symbolism/blob/master/Symbolism/Symbolism.cs

// public API: 
// Expression.prototype.add
// Expression.prototype.subtract
// Expression.prototype.multiply
// ...
// protected API:
// Expression.prototype.addExpression
// Expression.prototype.addInteger

  import BigInteger from './BigInteger.js';
  import Polynomial from './Polynomial.js';
  import Matrix from './Matrix.js';
  import primeFactor from './primeFactor.js';

  var pow = function (x, count, accumulator) {
    if (count < 0) {
      throw new RangeError();
    }
    if (count > 9007199254740991) {
      throw new RangeError("NotSupportedError");
    }
    return (count < 1 ? accumulator : (2 * Math.floor(count / 2) !== count ? pow(x, count - 1, accumulator.multiply(x)) : pow(x.multiply(x), Math.floor(count / 2), accumulator)));
  };

  var matrixInN = function (matrix, n) {
    // https://stackoverflow.com/a/15302448/839199
    var binomialCoefficient = function (n, k) { // binomail coefficient
      return k === 0 ? Expression.ONE : n.multiply(binomialCoefficient(n.subtract(Expression.ONE), k - 1)).divide(Integer.fromNumber(k));
    };
    //Note: experimental
    // {{1,0,0},{0,1,1},{0,0,1}}^n === {{1,0,0},{0,1,n},{0,0,1}}
    // A
    var a = matrix;
    // A^(n-1)
    var symbolName = "aa";
    var anm1 = matrix.map(function (e, i, j) {
      return new Expression.Symbol(symbolName + "_(" + i + "," + j + ")");
    });
    var anm1previous = anm1.map(function (e, i, j) {
      return Expression.ZERO;
    });
    var an = undefined;
    while (!anm1.eql(anm1previous)) {
      anm1previous = anm1;
      // A^(n) = A^(n-1) * A;
      an = anm1.multiply(a);
      anm1 = an.map(function (e, i, j) {
        // an: {{1,0,0},{0,1,1+aa_23},{0,0,1}}
        // a_n = n + a_(n-1)
        // a_n = k * a_(n-1) + k**(n-(m+1)) * choose(n-1, m)
        // =>
        // a_n = k**(n-(m+1)) * choose(n-1, m)
        if (!(e instanceof Integer)) {
          var m = Polynomial.toPolynomial(e.getNumerator(), n).getDegree();
          var previous = anm1.e(i, j);
          var p = Polynomial.toPolynomial(e.getNumerator(), previous);
          var k = p.getLeadingCoefficient().divide(e.getDenominator());
          //TODO: remove `k instanceof Integer`
          if (p.getDegree() === 1 &&
              a.e(i, j).equals(Expression.ZERO) && //TODO: remove
              k instanceof Integer &&
              e.equals(k.multiply(previous).add(k.pow(n).divide(k.pow(Integer.fromNumber(m + 1))).multiply(binomialCoefficient(n.subtract(Expression.ONE), m))))) {
            console.log("!", e.toString());
            // a.e(i, j).add()
            return k.pow(n).divide(k.pow(Integer.fromNumber(m + 2))).multiply(binomialCoefficient(n.subtract(Expression.ONE), m + 1));
          }
        }
        // a_n = a_(n-1)
        if (e.equals(anm1.e(i, j))) {
          return a.e(i, j);
        }
        // a_n = k * a_(n-1) => a_n = k**(n - 1) * a_1
        if (anm1.e(i, j) instanceof Expression.Symbol && !e.equals(Expression.ZERO)) {
          var k = e.divide(anm1.e(i, j));
          if (k instanceof Integer ||
              k instanceof Expression.Complex ||
              (k.getNumerator() instanceof Integer && e.getDenominator() instanceof Integer)) {
            return k.pow(n.subtract(Expression.TWO)).multiply(a.e(i, j));
          }
        }
        // a_n = a_(n-1) + b => a_n = a_1 + b*(n-1)
        var sub = e.subtract(anm1.e(i, j));
        if (sub instanceof Integer) {
          return a.e(i, j).add(sub.multiply(n.subtract(Expression.TWO)));
        }
        // a_n = d**(n-1) + d * a_(n-1)
        if (e instanceof Expression.Division && e.getDenominator() instanceof Integer) {
          var d = e.getDenominator();
          if (e.equals(d.pow(n.subtract(Expression.ONE)).add(d.multiply(anm1.e(i, j))))) {
            return d.pow(n.subtract(Expression.TWO)).multiply(n.subtract(Expression.TWO).add(a.e(i, j)));
          }
        }
        // a_n = (-1)**(n-1) - a_(n-1)
        var d = Expression.ONE.negate();
        if (e.equals(d.pow(n.subtract(Expression.ONE)).add(d.multiply(anm1.e(i, j))))) {
          return d.pow(n.subtract(Expression.TWO)).multiply(n.subtract(Expression.TWO).add(a.e(i, j)));
        }
        
        //?
        var d = Expression.I;
        if (e.equals(d.pow(n.subtract(Expression.ONE)).add(d.multiply(anm1.e(i, j))))) {
          return d.pow(n.subtract(Expression.TWO)).multiply(n.subtract(Expression.TWO).add(a.e(i, j)));
        }
        
        
        return anm1.e(i, j);
      });
    }
    var ok2 = true;
    for (var i = 0; i < anm1.rows(); i += 1) {
      for (var j = 0; j < anm1.cols(); j += 1) {
        var e = anm1.e(i, j);
        if (e instanceof Expression.Symbol && e.symbol.slice(0, symbolName.length) === symbolName) {
          ok2 = false;
        }
      }
    }
    return !ok2 ? undefined : an;
  };

  var enableEX = true;
  var enable2X = true;

  Expression.prototype.powExpression = function (x) {
    var y = this;

    if (x instanceof Expression.Matrix && y instanceof Expression.Symbol && (y.symbol === "t" || y.symbol === "T")) {
      return x.transpose();
    }

    //!
    if (y instanceof Division && y.a instanceof Integer && y.b instanceof Integer) {
      if (typeof hit === "function") {
        hit({powExpression: y.toString()});
      }
      var n = Number.parseInt(y.b.toString(), 10);
      //if (n >= 2 && n <= 9007199254740991) {//TODO:
        var q = y.a.truncatingDivide(y.b);
        var r = y.a.subtract(q.multiply(y.b));
        return x.pow(q).multiply(x.pow(r)._nthRoot(n));
      //}
    }
    //!

    //!new 2017-05-08
    if (enableEX) {
      if (x === Expression.E || (enable2X && x instanceof Integer && x.compareTo(Expression.ONE) > 0 && integerPrimeFactor(x).compareTo(x) === 0)) {
        var isValid = function (y) {
          if (y instanceof Addition) {
            return isValid(y.a) && isValid(y.b);
          }
          if (y instanceof Integer && x === Expression.E) {
            return true;
          }
          return (y instanceof Expression.Symbol || y instanceof Multiplication && y.a instanceof Integer && y.b instanceof Expression.Symbol);
        };
        if (isValid(y)) {
          if (y.isNegative()) {
            return Expression.ONE.divide(new Expression.Exponentiation(x, y.negate()));
          }
          return new Expression.Exponentiation(x, y);
        }
      }
      if (enable2X && x instanceof Integer && x.compareTo(Expression.ONE) > 0) {
        if (y instanceof Addition && (y.a instanceof Integer || y.b instanceof Integer)) {
          return x.pow(y.a).multiply(x.pow(y.b));
        }
        var xf = integerPrimeFactor(x);
        return xf.pow(y).multiply(x.divide(xf).pow(y));
      }
    }
    //!

    if (enableEX) {
      //TODO: - ?
      if (x instanceof Integer && x.equals(Expression.ONE)) {
        return Expression.ONE;
      }
      if (x instanceof Division) {
        return x.a.pow(y).divide(x.b.pow(y));
      }
      if (enable2X) {
        if (x instanceof Multiplication) {
          return x.a.pow(y).multiply(x.b.pow(y));
        }
      }
    }

    var yn = y.getNumerator();
    if (x === Expression.E && yn instanceof Multiplication && yn.a instanceof Expression.Complex && yn.a.real.equals(Expression.ZERO) && yn.b instanceof Expression.Symbol) {
      var t = y.multiply(Expression.I.negate());
      return t.cos().add(Expression.I.multiply(t.sin()));
    }

    //TODO: 
    if (x instanceof Expression.Matrix && x.matrix.isSquare() && y instanceof Expression.Symbol && (y.symbol === "n" || y.symbol === "k")) {
      var tmp = matrixInN(x.matrix, y);
      if (tmp != undefined) {
        return new Expression.Matrix(tmp);
      }

      //! 2018-08-26
      if (true) {
        var tmp = Expression.getEigenvalues(x.matrix);
        var eigenvalues = tmp.eigenvalues;
        var multiplicities = tmp.multiplicities;
        if (Expression.sum(multiplicities) === x.matrix.cols()) {
          var tmp2 = Expression.getEigenvectors(x.matrix, eigenvalues);
          var eigenvectors = tmp2.eigenvectors;
          if (eigenvectors.length === x.matrix.cols()) {
            var tmp = Expression.diagonalize(x.matrix, eigenvalues, multiplicities, eigenvectors);
            var L = tmp.L;
            var SL = matrixInN(L, y);
            if (SL != undefined) {
              if (Expression.callback != undefined) {
                Expression.callback(new Expression.Event("pow-using-diagonalization", x));
              }
              if (Expression.callback != undefined) {
                //TODO more details (A=P*D*P^-1 => A^n=P*D*P^-1 * ... * P*D*P^-1=P*D^n*P^1
                Expression.callback(new Expression.Event("diagonalize", x));
              }
              return new Expression.Matrix(tmp.T.multiply(SL).multiply(tmp.T_INVERSED));
            }
          } else {
            var tmp = Expression.getFormaDeJordan(x.matrix, eigenvalues, multiplicities);
            var JN = matrixInN(tmp.J, y);
            if (JN != undefined) {
              if (Expression.callback != undefined) {
                Expression.callback(new Expression.Event("pow-using-Jordan-normal-form", x));
              }
              if (Expression.callback != undefined) {
                //TODO more details (A=P*D*P^-1 => A^n=P*D*P^-1 * ... * P*D*P^-1=P*D^n*P^1
                Expression.callback(new Expression.Event("Jordan-decomposition", x));
              }
              //TODO: details !!!
              return new Expression.Matrix(tmp.P.multiply(JN).multiply(tmp.P_INVERSED));
            }
          }
        }
      }
      //!
    }
    
    if (Expression.ExponentiationOfMinusOne != null) {
      if (x instanceof Integer && x.compareTo(Expression.ZERO) < 0) {
        if (y instanceof Expression.Symbol && (y.symbol === "n" || y.symbol === "k")) {
          return new Expression.ExponentiationOfMinusOne(Expression.ONE.negate(), y).multiply(x.negate().pow(y));
        }
        if (y instanceof Addition && y.a instanceof Expression.Symbol && (y.a.symbol === "n" || y.a.symbol === "k") && y.b instanceof Integer) {
          return new Expression.ExponentiationOfMinusOne(Expression.ONE.negate(), y.a).multiply(Expression.ONE.negate().pow(y.b)).multiply(x.negate().pow(y));
        }
        if (y instanceof Multiplication) {
          return x.pow(y.a).pow(y.b);
        }
      }
    }

    if (Expression.ExponentiationOfImaginaryUnit != null) {
      if (x instanceof Expression.Complex && x.equals(Expression.I.negate())) {
        return Expression.ONE.negate().pow(y).multiply(x.negate().pow(y));
      }
      if (x instanceof Expression.Complex && x.equals(Expression.I)) {//TODO: -i, other complex numbers - ?
        if (y instanceof Expression.Symbol && (y.symbol === "n" || y.symbol === "k")) {
         return new Expression.ExponentiationOfImaginaryUnit(Expression.I, y);
        }
        if (y instanceof Addition && y.a instanceof Expression.Symbol && (y.a.symbol === "n" || y.a.symbol === "k") && y.b instanceof Integer) {
          var t = Expression.I.pow(y.b);
          return new Expression.ExponentiationOfImaginaryUnit(Expression.I, t instanceof Expression.Complex ? y.a.add(Expression.ONE) : y.a).multiply(t instanceof Expression.Complex ? t.divide(Expression.I) : t);
        }
        if (y instanceof Multiplication) {
          return x.pow(y.a).pow(y.b);
        }
      }
    }

    if (x === Expression.E && y instanceof Expression.Matrix && y.matrix.isSquare()) {
      // https://en.wikipedia.org/wiki/Matrix_exponential#Using_the_Jordan_canonical_form
      var tmp = Expression.getEigenvalues(y.matrix);
      var eigenvalues = tmp.eigenvalues;
      var multiplicities = tmp.multiplicities;
      if (Expression.sum(multiplicities) === y.matrix.cols()) {
        var tmp = Expression.getFormaDeJordan(y.matrix, eigenvalues, multiplicities);
        var D = tmp.J.map(function (e, i, j) {
          return i === j ? e : Expression.ZERO;
        });
        var N = tmp.J.map(function (e, i, j) {
          return i !== j ? e : Expression.ZERO;
        });
        var exp = function (N) {
          // https://en.wikipedia.org/wiki/Matrix_exponential#Nilpotent_case
          var z = Matrix.Zero(N.cols(), N.cols());
          var s = z;
          var p = Matrix.I(N.cols());
          var k = 0;
          var f = 1;
          while (!p.eql(z)) {
            s = s.add(p).scale(Expression.ONE.divide(Integer.fromNumber(f)));
            p = p.multiply(N);
            k += 1;
            f *= k;
          }
          return s;
        };
        if (Expression.callback != undefined) {
          Expression.callback(new Expression.Event("exponential-using-Jordan-canonical-form", y));
        }
        if (Expression.callback != undefined) {
          Expression.callback(new Expression.Event("Jordan-decomposition", y));
        }
        return new Expression.Matrix(tmp.P.multiply(D.map(function (e, i, j) {
          return i === j ? Expression.E.pow(e) : Expression.ZERO;
        }).multiply(exp(N))).multiply(tmp.P_INVERSED));
      }
    }
    
    throw new RangeError("NotSupportedError");
  };

  // compare two expression, which are factors (multiplicaiton operands) of terms (addition operands)
  Expression.prototype.compare4Addition = function (y) {
    var x = this;
    if (x instanceof Expression.Symbol && y instanceof Integer) {
      return +1;
    }
    if (x instanceof Integer && y instanceof Expression.Symbol) {
      return -1;
    }
    if (x instanceof Integer && y instanceof Integer) {
      return x.compareTo(y);
    }
    if (x instanceof Expression.Symbol && y instanceof Expression.Symbol) {
      return x.symbol < y.symbol ? -1 : (y.symbol < x.symbol ? +1 : 0);
    }
    //!
    if (x instanceof Expression.Matrix && y instanceof MatrixSymbol) {
      return +1;
    }
    if (x instanceof MatrixSymbol && y instanceof Expression.Matrix) {
      return -1;
    }
    if (x instanceof Expression.Matrix && y instanceof Expression.Matrix) {
      if (x.matrix.rows() === y.matrix.rows() &&
          x.matrix.cols() === y.matrix.cols()) {
        var rows = x.matrix.rows();
        var cols = x.matrix.cols();
        for (var i = 0; i < rows; i += 1) {
          for (var j = 0; j < cols; j += 1) {
            var c = x.matrix.e(i, j).compare4Addition(y.matrix.e(i, j));
            if (c !== 0) {
              return c;
            }
          }
        }
      }
    }
    
    //!new 2016-12-17
    //NOTE: the `x instanceof Addition || y instanceof Addition` should be used before `x instanceof Multiplication || y instanceof Multiplication`
    if (x instanceof Addition || y instanceof Addition) {
      return Addition.compare4Addition(x, y);
    }

    //!new 2016-10-09
    if (x instanceof Multiplication || y instanceof Multiplication) {
      return Multiplication.compare4Addition(x, y);
    }
    
    //!new 2017-02-10
    if (x instanceof Expression.Matrix || y instanceof Expression.Symbol) {
      return -1;
    }

    //!new 2018-10-11
    if (x instanceof Integer && y instanceof Expression.Function) {
      return -1;
    }

    //!new 2018-11-12
    if (x instanceof Division && y instanceof Division) {
      return x.a.compare4Addition(y.a) || x.b.compare4Addition(y.b);//?
    }
    if (x instanceof Expression && y instanceof Division) {
      return +1;//?
    }
    if (x instanceof Division && y instanceof Expression) {
      return -1;//?
    }

    //!2019-02-18
    if (x instanceof Integer && y instanceof Expression.Complex) {
      return -1;//?
    }
    if (x instanceof Expression.Complex && y instanceof Integer) {
      return +1;//?
    }
    //!

    //!
    throw new RangeError();
  };
  
  var compare = function (x, y) {
    return x.compare4Addition(y);
  };

  var compare4Multiplication = function (x, y) {
    //TODO: Exponentiation + Exponentiation, Exponentiation + Symbol, Symbol + Exponentiation
    return x.compare4Multiplication(y);
  };

  var getBase = function (x) {
    return x instanceof Exponentiation ? x.a : x;
  };
  var getExponent = function (x) {
    return x instanceof Exponentiation ? x.b : Expression.ONE;
  };

  var getConstant = function (e) {
    if (e instanceof Integer) {
      return e;
    } else if (e instanceof Expression.Complex) {
      return e;
    } else if (e instanceof Multiplication) {
      var c = undefined;
      var x = e;
      for (var multiplications = x.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
        var t = getConstant(y);
        c = c == undefined ? t : t.multiply(c);
      }
      if (c != undefined) {
        return c;
      }
    } else if (e instanceof Addition) { // -5*x+15
      var c = undefined;
      for (var additions = e.summands(), x = additions.next().value; x != null; x = additions.next().value) {
        var t = getConstant(x);
        //c = c == undefined ? t : integerGCD(t, c);
        c = c == undefined ? t : complexGCD(t, c);
      }
      if (c != undefined) {
        return c;
      }
    }
    return Expression.ONE;
  };
  var getTerm = function (x) {
  // TODO: fix performance ?
    if (x instanceof Integer) {
      return undefined;
    } else if (x instanceof Expression.Complex) {
      return undefined;
    } else if (x instanceof Multiplication) {
      var terms = [];
      for (var multiplications = x.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
        var t = getTerm(y);
        if (t != undefined) {
          terms.push(t);
        }
      }
      var result = undefined;
      for (var j = terms.length - 1; j >= 0; j -= 1) {
        result = result == undefined ? terms[j] : new Multiplication(result, terms[j]);
      }
      return result;
    }
    return x;
  };

  var multiplyByInteger = function (x, y) {
    if (x.compareTo(Expression.ZERO) === 0) {
      return x;
    }
    if (x.compareTo(Expression.ONE) === 0) {
      return y;
    }
    return new Multiplication(x, y);
  };
  
  Expression.prototype.multiplyExpression = function (x) {
    var y = this;

    if (x instanceof Expression && y instanceof Multiplication) {
      return x.multiply(y.a).multiply(y.b);
    }
    if (x instanceof Multiplication && y instanceof Expression) {
      var c = compare4Multiplication2(x.b, y);
      if (c === 0) {
        return x.a.multiply(x.b.multiply(y));
      }
      return c > 0 ? x.a.multiply(y).multiply(x.b) : new Multiplication(x, y);
    }

    //!
    if (x instanceof IdentityMatrix && y instanceof MatrixSymbol) {
      return y;
    }
    if (y instanceof IdentityMatrix && x instanceof MatrixSymbol) {
      return x;
    }
    //!
    // rest

    var c = 0;
    if (x instanceof Integer && y instanceof Expression.Symbol) {
      return multiplyByInteger(x, y);
    }
    if (x instanceof Expression.Symbol && y instanceof Integer) {
      return multiplyByInteger(y, x);
    }
    if (x instanceof Expression.Symbol && y instanceof Expression.Symbol) {
      c = compare4Multiplication(x, y);
      if (c === 0) {
        return x.pow(Expression.TWO);
      }
      return c > 0 ? new Multiplication(y, x) : new Multiplication(x, y);
    }
    if (x instanceof Integer && y instanceof Exponentiation) {
      return multiplyByInteger(x, y);
    }
    if (x instanceof Exponentiation && y instanceof Integer) {
      return multiplyByInteger(y, x);
    }
    if (x instanceof Exponentiation && y instanceof Expression.Symbol) {
      c = compare4Multiplication(getBase(x), y);
      if (c === 0) {
        return y.pow(getExponent(x).add(Expression.ONE));
      }
      return c > 0 ? new Multiplication(y, x) : new Multiplication(x, y);
    }
    if (x instanceof Expression.Symbol && y instanceof Exponentiation) {
      c = compare4Multiplication(x, getBase(y));
      if (c === 0) {
        return x.pow(getExponent(y).add(Expression.ONE));
      }
      return c > 0 ? new Multiplication(y, x) : new Multiplication(x, y);
    }
    if (x instanceof Exponentiation && y instanceof Exponentiation) {
      c = compare4Multiplication(getBase(x), getBase(y));
      if (c === 0) {
        return getBase(x).pow(getExponent(x).add(getExponent(y)));
      }
      return c > 0 ? new Multiplication(y, x) : new Multiplication(x, y);
    }

    if (x instanceof SquareRoot && y instanceof SquareRoot) {
      return x.a.multiply(y.a).squareRoot();
    }
    if (x instanceof CubeRoot && y instanceof CubeRoot) {
      return x.a.multiply(y.a).cubeRoot();
    }
    if (x instanceof NthRoot && y instanceof NthRoot) {
      if (x.n === y.n) {
        return x.a.multiply(y.a)._nthRoot(x.n);
      }
      return x.n < y.n ? new Multiplication(x, y) : new Multiplication(y, x);
    }

    //!
    if (x instanceof MatrixSymbol && y instanceof Expression.Matrix) {
      return new Multiplication(x, y);
    }
    if (x instanceof Expression.Matrix && y instanceof MatrixSymbol) {
      return new Multiplication(x, y);
    }
    if (has(x, MatrixSymbol) && y instanceof Expression.Matrix) { // X^2*{{1,2},{3,4}}
      return new Multiplication(x, y);
    }
    if (x instanceof Expression.Matrix && has(y, MatrixSymbol)) { // {{1,2},{3,4}}*X^2
      return new Multiplication(x, y);
    }
    
    //!
    //throw new RangeError();
    if (x instanceof Integer && y instanceof Expression) {
      if (x.compareTo(Expression.ZERO) === 0) {
        return x;
      }
      if (x.compareTo(Expression.ONE) === 0) {
        return y;
      }
    }
    if (x instanceof Expression && y instanceof Integer) {
      if (y.compareTo(Expression.ZERO) === 0) {
        return y;
      }
      if (y.compareTo(Expression.ONE) === 0) {
        return x;
      }
    }

    if (x instanceof Expression.Complex && y instanceof Expression.ExponentiationOfImaginaryUnit) {
      //!hack
      //TODO: remove
      if (x.real.equals(Expression.ZERO) && x !== Expression.I) {
        return x.imaginary.multiply(y.multiply(Expression.I));
      }
    }
    if (x instanceof Expression.ExponentiationOfImaginaryUnit && y instanceof Expression.Complex) {
      //!hack
      //TODO: remove
      if (y.real.equals(Expression.ZERO) && y !== Expression.I) {
        return y.imaginary.multiply(x.multiply(Expression.I));
      }
    }

    var cmp = compare4Multiplication(getBase(x), getBase(y));
    if (cmp === 0) {
      return getBase(x).pow(getExponent(x).add(getExponent(y)));
    }
    if (cmp < 0) {
      return new Multiplication(x, y);
    }
    if (cmp > 0) {
      return new Multiplication(y, x);
    }

  };

  function Iterator() {
  }
  if (typeof Symbol === "function") {
    Iterator.prototype[Symbol.iterator] = function () {
      return this;
    };
    Object.defineProperty(Iterator.prototype, "done", {
      get: function () {
        return this.value == null;
      }
    });
  }

  function TermFactorsIterator(e) {
    this.value = undefined;
    this.e = e;
  }
  TermFactorsIterator.prototype = Object.create(Iterator.prototype);
  TermFactorsIterator.prototype.next = function () {
    this.value = this.e instanceof Multiplication ? this.e.b : (this.e instanceof Integer || this.e instanceof Expression.Complex ? null : this.e);
    this.e = this.e instanceof Multiplication ? this.e.a : undefined;
    return this;
  };

  function termFactors(e) {
    return new TermFactorsIterator(e);
  }

  var compare4Addition = function (x, y) {
    // undefined | Symbol | Exponentiation | Multiplication
    var i = termFactors(x);
    var j = termFactors(y);
    var a = i.next().value;
    var b = j.next().value;
    while (a != null && b != null) {

      //!
      // x^3*y^2, x^2*y^3
      var cmp = 0 - compare(getBase(a), getBase(b));
      if (cmp === 0) {
        cmp = compare(getExponent(a), getExponent(b));
      }
      if (cmp !== 0) {
        return cmp;
      }
      a = i.next().value;
      b = j.next().value;
    }
    return a != null ? +1 : (b != null ? -1 : 0);
  };

  Expression.prototype.addExpression = function (x) {
    var y = this;
    if (x.equals(Expression.ZERO)) {
      return y;
    }
    if (y.equals(Expression.ZERO)) {
      return x;
    }

    //!2019-02-16
    if (x instanceof Multiplication && x.b instanceof IdentityMatrix) {
      var t = getIdentityMatrixCoefficient(y);
      if (t != null) {
        return x.a.add(t).multiply(x.b);
      }
    } else if (x instanceof IdentityMatrix) {
      var t = getIdentityMatrixCoefficient(y);
      if (t != null) {
        return Expression.ONE.add(t).multiply(x);
      }
    }
    if (y instanceof Multiplication && y.b instanceof IdentityMatrix) {
      var t = getIdentityMatrixCoefficient(x);
      if (t != null) {
        return t.add(y.a).multiply(y.b);
      }
    } else if (y instanceof IdentityMatrix) {
      var t = getIdentityMatrixCoefficient(x);
      if (t != null) {
        return t.add(Expression.ONE).multiply(y);
      }
    }
    //!2019-02-16

    // rest

    var i = x.summands();
    var j = y.summands();
    var a = i.next().value;
    var b = j.next().value;
    var s = [];
    //a + b, compare4Addition("a", "b") > 0
    while (a != null && b != null) {
      var c = compare4Addition(a, b);
      if (c < 0) {
        s.push(a);
        a = i.next().value;
      } else if (c > 0) {
        s.push(b);
        b = j.next().value;
      } else {
        var constant = getConstant(a).add(getConstant(b));
        var term = getTerm(a);
        var last = term == undefined ? constant : constant.multiply(term);
        if (!last.equals(Expression.ZERO)) {
          s.push(last);
        }
        a = i.next().value;
        b = j.next().value;
      }
    }
    while (a != null) {
      s.push(a);
      a = i.next().value;
    }
    while (b != null) {
      s.push(b);
      b = j.next().value;
    }
    if (s.length === 0) {
      return Expression.ZERO;
    }
    var accumulator = s[s.length - 1];
    for (var k = s.length - 2; k >= 0; k -= 1) {
      var currentValue = s[k];
      accumulator = new Addition(accumulator, currentValue);
    }
    return accumulator;
  };

  var pseudoRemainder = function (x, y) {
    var lcg = y.getLeadingCoefficient();
    var n = x.getDegree() - y.getDegree();
    // assertion
    if (n < 0) {
      throw new RangeError();
    }
    var sx = x.multiply(Polynomial.of(lcg.pow(Integer.fromNumber(n).add(Expression.ONE))));
    return sx.divideAndRemainder(y, "throw").remainder;
  };

  var divideByInteger = function (e, f) {
    if (f.equals(Expression.ZERO)) {
      throw new RangeError("ArithmeticException");
    }
    var result = Expression.ZERO;
    for (var additions = e.summands(), x = additions.next().value; x != null; x = additions.next().value) {
      var rest = Expression.ONE;
      var t = undefined;
      // TODO: check, fix?
      for (var multiplications = x.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
        var z = y;
        if (z instanceof Integer || z instanceof Expression.Complex) {
          if (t != undefined) {
            console.warn("!");
            t = t.multiply(z);
          } else {
            t = z;
          }
        } else {
          if (rest === Expression.ONE) {
            rest = z;
          } else {
            rest = z.multiply(rest);
          }
        }
      }
      if (!(t instanceof Expression.Complex)) {
      if (!(t instanceof Integer)) {
        throw new RangeError();
      }
      }
      //result = result.add(t.divide(f).multiply(rest));
      result = result.add(t.truncatingDivide(f).multiply(rest));
    }
    return result;
  };

  Expression.getCoefficients = function (e, v) {
    var result = [];
    for (var additions = e.summands(), x = additions.next().value; x != null; x = additions.next().value) {
      var d = Expression.ZERO;
      var c = Expression.ONE;
      for (var multiplications = x.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
        var t = y;
        for (var variables = getVariableInternal(t), ve = variables.next().value; ve != null; ve = variables.next().value) {
          if (ve.v.equals(v)) {
            d = d.add(ve.e);
          } else {
            c = c.multiply(ve.e === Expression.ONE ? ve.v : ve.v.pow(ve.e));
          }
        }
      }
      var tmp = {
        coefficient: c,
        degree: d
      };
      var k = result.length - 1;
      while (k >= 0 && tmp.degree.compareTo(result[k].degree) < 0) {
        k -= 1;
      }
      if (k >= 0 && tmp.degree.compareTo(result[k].degree) === 0) {
        result[k].coefficient = tmp.coefficient.add(result[k].coefficient);
      } else {
        result.push(tmp);
        var i = result.length - 1;
        while (i >= k + 2) {
          result[i] = result[i - 1];
          i -= 1;
        }
        result[k + 1] = tmp;
      }
    }
    return result;
  };

  //TODO: remove
  var getFirstAdditionOperand = function (x) {
    var result = x;
    while (result instanceof Addition) {
      result = result.a;
    }
    return result;
  };
  //TODO: remove
  var getLastMultiplicationOperand = function (x) {
    var result = x;
    while (result instanceof Multiplication) {
      result = result.b;
    }
    return result;
  };

  function VIterator(v) {
    if (v == undefined) {
      throw new Error();
    }
    this.value = undefined;
    this.v = v;
  }
  VIterator.prototype = Object.create(Iterator.prototype);
  VIterator.prototype.next = function () {
    this.value = this.v;
    this.v = undefined;
    return this;
  };

  function VariablesIterator(v, additions) {
    if (additions == undefined) {
      throw new Error();
    }
    this.value = undefined;
    this.v = v;
    this.additions = additions;
  }
  VariablesIterator.prototype = Object.create(Iterator.prototype);
  VariablesIterator.prototype.next = function () {
    var x = this.additions.next().value;
    var value = null;
    if (x == null) {
      value = null;
    } else if (x instanceof Expression.Symbol) {
      value = {v: new Exponentiation(this.v, x), e: Expression.ONE};
    } else if (x instanceof Multiplication && x.a instanceof Integer && x.b instanceof Expression.Symbol) {
      value = {v: new Exponentiation(this.v, x.b), e: x.a};
    } else if (x instanceof Integer) {
      value = {v: this.v, e: x};
    } else {
      throw new RangeError();
    }
    this.value = value;
    return this;
  };

  var getVariableInternal = function (t) {
    if (t instanceof Expression.ExponentiationOfMinusOne) {//TODO: ?
      return new VIterator({v: t, e: Expression.ONE});
    }
    if (t instanceof Expression.ExponentiationOfImaginaryUnit) {//TODO: ?
      return new VIterator({v: t, e: Expression.ONE});
    }
    var v = getBase(t);
    var e = getExponent(t);

    //!new 2017-05-08
    if (enableEX) {
      if (!(e instanceof Integer)) {
        var additions = e.summands();
        return new VariablesIterator(v, additions);
      }
    }
    //!
    return new VIterator({v: v, e: e});
  };

  var getVariable = function (e) {
    //? square roots at first
    for (var additions = e.summands(), x = additions.next().value; x != null; x = additions.next().value) {
      for (var multiplications = x.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
        if (y instanceof NthRoot) {
        //TODO: assert(y instanceof Integer)
          return y;
        }
      }
    }
    //?

    var result = getVariableInternal(getLastMultiplicationOperand(getFirstAdditionOperand(e))).next().value.v;
    //!?
    //if (result instanceof NthRoot) {
    //  return undefined;
    //}
    //
    if (result instanceof Expression.Complex) {
      return undefined;//?
    }
    if (result instanceof Integer) {
      return undefined;//?
    }
    return result;
  };

  Expression.getVariable = getVariable;

  var integerGCD = function (x, y) {
    if (x.compareTo(Expression.ZERO) < 0) {
      x = x.negate();
    }
    if (y.compareTo(Expression.ZERO) < 0) {
      y = y.negate();
    }
    var a = x.value;
    var b = y.value;
    while (!BigInteger.equal(b, BigInteger.BigInt(0))) {
      var t = BigInteger.remainder(a, b);
      a = b;
      b = t;
    }
    return new Integer(a);
  };

  var getIntegerContent = function (x) {
    return x instanceof Expression.Complex ? integerGCD(x.real, x.imaginary) : x;
  };

  var complexGCD = function (a, b) {//?
    return integerGCD(getIntegerContent(a), getIntegerContent(b));
  };

  // http://www-troja.fjfi.cvut.cz/~liska/ca/node33.html
  var gcd = function (a, b, v) {
    if (v == undefined) {
      return complexGCD(getConstant(a), getConstant(b));
    }
    return polynomialGCD(Polynomial.toPolynomial(a, v), Polynomial.toPolynomial(b, v)).calcAt(v);
  };

  var polynomialGCD = function (a, b) {
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
      var r = pseudoRemainder(A, B);
      //! 2018-05-12
      if (!r.equals(Polynomial.ZERO)) {
        // https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor#Primitive_pseudo-remainder_sequence
        // RPN("(x^8+x^6-3x^4-3x^3+8x^2+2x-5)/(3x^6+5x^4-4x^2-9x+21)")
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

  // ! new 21.12.2013 (square roots)

  var getConjugateFactor = function (e) {
    var r = 1 / 0;
    var p = undefined;
    for (var additions = e.summands(), x = additions.next().value; x != null; x = additions.next().value) {
      for (var multiplications = x.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
        if (y instanceof NthRoot) {
          var degree = y.getDegree();
          if (r > degree) {
            p = y.a;
            r = degree;
          } else if (r === degree) {
            //TODO: assert(y instanceof Integer)
            if (y.a instanceof Integer) {
              var z = integerGCD(p, y.a);
              if (z.compareTo(Expression.ONE) !== 0) {
                p = z;//!
              }
            } else {
              throw new TypeError();
            }
          }
        }
      }
    }
    return {p: p, r: r};
  };

  // TODO: test
  var getConjugate = function (a) {
    var e = undefined;
    e = Expression.getComplexConjugate(a);
    if (e != undefined) {
      return e;
    }
    e = Expression.getNthRootConjugate(a);
    if (e != undefined) {
      return e;
    }
    return undefined;
  };
  
  Expression.getConjugate = getConjugate;

  // https://en.wikipedia.org/wiki/Conjugate_(square_roots)
  Expression.getNthRootConjugate = function (e) {
  //TODO: fix
  //if (true) return undefined;
    var tmp = getConjugateFactor(e);
    var p = tmp.p;
    var r = tmp.r;
    if (p == undefined) {
      return undefined;
    }
    var polynomial = Polynomial.of();
    for (var additions = e.summands(), x = additions.next().value; x != null; x = additions.next().value) {
      var degree = 0;
      var coefficient = Expression.ONE;
      for (var multiplications = x.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
        if (y instanceof NthRoot && r === y.getDegree()) {
          if (y.a instanceof Integer) {
            var j = 0;
            var a = y.a;
            while (a.remainder(p).equals(Expression.ZERO)) {
              a = a.truncatingDivide(p);
              j += 1;
            }
            coefficient = coefficient.multiply(a._nthRoot(r));
            degree += j;
          } else {
            throw new TypeError();
          }
        } else {
          coefficient = coefficient.multiply(y);
        }
      }
      //TODO: efficiency ?
      polynomial = polynomial.add(Polynomial.of(coefficient).shift(degree));
    }
    var n = r;
    var x = p;
    var e = polynomial;
    var conjugate = Polynomial.of(n % 2 === 1 ? Expression.ONE : Expression.ONE.negate());
    while (e.getDegree() > 0) {
      var k = e.getDegree();
      var ak = Polynomial.of(e.getLeadingCoefficient());
      var m = Polynomial.of(Expression.ONE).shift(n - k);
      var p = e.multiply(m).subtract(ak.shift(n)).add(ak.multiply(Polynomial.of(x)));
      while (p.getDegree() >= k) {
        var d = p.getDegree();
        var lc = p.getLeadingCoefficient();
        m = m.multiply(ak);
        p = p.multiply(ak);
        var t = Polynomial.of(lc.negate()).shift(d - k);
        m = m.add(t);
        p = p.add(e.multiply(t));
      }
      e = p;
      conjugate = conjugate.multiply(m);
    }
    return conjugate.calcAt(x._nthRoot(n));
  };

  // without the checks
  Expression.collectLinearEquationVariables = function (e) {
    if (e instanceof Division) {
      throw new RangeError();
    }
    var list = [];
    for (var additions = e.summands(), x = additions.next().value; x != null; x = additions.next().value) {
      var v = undefined;
      var c = Expression.ONE;
      var NO_VARIABLE = "";
      for (var multiplications = x.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
        if (y instanceof Expression.Symbol && v == undefined) {
          v = y;
        } else {
          if (!(y instanceof Integer) && !(y instanceof NthRoot)) {
            if (v == undefined) {
              v = NO_VARIABLE;
            }
          }
          c = c.multiply(y);
        }
      }
      if (v == undefined) {
        v = NO_VARIABLE;
      }
      var variable = v === NO_VARIABLE ? "" : v.toString();
      list.push({c: c, v: variable});
    }
    return list;
  };

  var has = function (e, Class) {
    if (e instanceof Class) {
      return true;
    }
    if (e instanceof BinaryOperation) {
      if (has(e.b, Class)) {
        return true;
      }
      return has(e.a, Class);
    }
    if (e instanceof Negation) {
      return has(e.b, Class);
    }
    if (e instanceof Expression.Function) {
      return has(e.a, Class);
    }
    return false;//?
  };
  Expression.has = has;

  Expression.prototype.divideExpression = function (x) {
    var y = this;
    
    //if (Expression.getIdentityMatrixCoefficient(x) != undefined) {
    //  if (y instanceof Expression.Matrix) {
    //    return Expression.getIdentityMatrixCoefficient(x).divide(y);
    //  }
    //  return Expression.makeIdentityMatrixWithCoefficient(Expression.getIdentityMatrixCoefficient(x).divide(y));
    //}
    //if (Expression.getIdentityMatrixCoefficient(y) != undefined) {
    //  if (x instanceof Expression.Matrix) {
    //    return x.divide(Expression.getIdentityMatrixCoefficient(y));
    //  }
    //  return Expression.makeIdentityMatrixWithCoefficient(x.divide(Expression.getIdentityMatrixCoefficient(y)));
    //}

    //if (has(x, IdentityMatrix)) {//?
    //  throw new RangeError("NotSupportedError");
    //}
    //if (has(x, MatrixSymbol)) {
    //  throw new RangeError("NotSupportedError");
    //}

    if (x instanceof Multiplication && x.b instanceof IdentityMatrix) {
      return x.a.divide(y).multiply(x.b);
    } else if (x instanceof IdentityMatrix) {
      return Expression.ONE.divide(y).multiply(x);
    }
    if (y instanceof Multiplication && y.b instanceof IdentityMatrix) {
      return x.divide(y.a).multiply(y.b);
    } else if (y instanceof IdentityMatrix) {
      return x.multiply(y);
    }

    if (has(y, MatrixSymbol)) {
      throw new RangeError("NotSupportedError");
    }

    if (x instanceof Expression.Matrix && y instanceof Expression.Matrix) {
      // TODO: callback ???
      return new Expression.Matrix(x.matrix.multiply(y.matrix.inverse()));
    }
    if (x instanceof Expression.Matrix && y instanceof Expression) {
      //return new Expression.Matrix(x.matrix.scale(y.inverse()));
      return x.multiply(y.inverse());
    }
    if (x instanceof Expression && y instanceof Expression.Matrix) {
      if (Expression.callback != undefined) {
        Expression.callback(new Expression.Event(y.matrix.getDeterminantEventType("inverse").type, y));
      }
      //return new Expression.Matrix(y.matrix.inverse().scale(x));
      return new Expression.Matrix(y.matrix.inverse()).multiply(x);
    }

    if (y.equals(Expression.ZERO)) {
      //TODO: fix?
      throw new RangeError("ArithmeticException");
    }
    if (x.equals(Expression.ZERO)) {
      return Expression.ZERO;
    }
    if (y.equals(Expression.ONE)) {
      return x;
    }

    //!!! new (21.12.2013)
    if (true) { //TODO: remove hack!
      var e = getConjugate(y);
      if (e != undefined) {
        return x.multiply(e).divide(y.multiply(e));
      }
    }
    
    if (true) {//!new 2017-11-25
      var c = Expression.getComplexConjugate(x);
      if (c != undefined) {
        var a = x.add(c).divide(Expression.TWO);
        var b = x.subtract(c).multiply(Expression.I.negate()).divide(Expression.TWO);
        var g = (a.equals(Expression.ZERO) ? y : a.gcd(y)).gcd(b.equals(Expression.ZERO) ? y : b.gcd(y));
        if (!g.equals(Expression.ONE)) {
          x = a.divide(g).add(b.divide(g).multiply(Expression.I));
          y = y.divide(g);
        }
        if (y.isNegative()) {
          x = x.negate();
          y = y.negate();
        }
        return y.equals(Expression.ONE) ? x : new Division(x, y);
      }
    }//!

    // check if y is instance of Integer to avoid issues with nth-roots (?) - see a test
    //TODO: investigate
    var v = y instanceof Integer ? undefined : getVariable(x);//???
    //var v = getVariable(x);//???
    //TODO: move?

    // gcd
    if (v == undefined) {
      var g = complexGCD(getConstant(x), getConstant(y));
      if (!g.equals(Expression.ONE)) {
        x = divideByInteger(x, g);
        y = divideByInteger(y, g);
        return x.divide(y);
      } 
    } else {
      var px = Polynomial.toPolynomial(x, v);
      var py = Polynomial.toPolynomial(y, v);
      //!TODO: remove - performance optimization
      var t = px.divideAndRemainder(py, "undefined");
      if (t != undefined && t.remainder.equals(Polynomial.ZERO)) {
        return t.quotient.calcAt(v);
      }
      //!
      var g = polynomialGCD(px, py);
      if (g.getDegree() !== 0 || !g.getLeadingCoefficient().equals(Expression.ONE)) { // g !== 1
        var x2 = px.divideAndRemainder(g, "throw").quotient;
        var y2 = py.divideAndRemainder(g, "throw").quotient;
        return x2.calcAt(v).divide(y2.calcAt(v));
      }
    }

    //var lc = getConstant(getLeadingCoefficient(y, v));
    //var lc = getConstant(getLeadingCoefficient(y, getVariable(y)));
    var lc = getConstant(getFirstAdditionOperand(y));
    if (lc.compareTo(Expression.ZERO) < 0) {
      return x.negate().divide(y.negate());
    }
    return new Division(x, y);
  };

  function Expression() {
    throw new Error("Do not call for better performance");
  }

  Expression.callback = undefined;
  Expression.Event = function (type, data, second) {
    second = second == undefined ? undefined : second;
    this.type = type;
    this.data = data;
    this.second = second;
  };

  Expression.prototype.compare4Multiplication = function (y) {
    throw new Error(this.toString());
  };
  Expression.prototype.compare4MultiplicationInteger = function (x) {
    throw new Error();
  };
  Expression.prototype.compare4MultiplicationSymbol = function (x) {
    throw new Error();
  };
  Expression.prototype.compare4MultiplicationNthRoot = function (x) {
    throw new Error();
  };

  Expression.prototype.negate = function () {
    return Expression.ONE.negate().multiply(this);
  };
  Expression.prototype.add = function (y) {
    return y.addExpression(this);
  };
  Expression.prototype.subtract = function (y) {
    return this.add(y.negate());
  };
  Expression.prototype.divide = function (y) {
    return y.divideExpression(this);
  };
  Expression.prototype.multiply = function (y) {
    return y.multiplyExpression(this);
  };
  Expression.prototype.pow = function (y) {
    return y.powExpression(this);
  };
  Expression.prototype.getDenominator = function () {
    //TODO: FIX!!!!
    return this instanceof Division ? this.b : Expression.ONE;
  };
  Expression.prototype.getNumerator = function () {
    //TODO: FIX!!!!
    return this instanceof Division ? this.a : this;
  };
  Expression.prototype.inverse = function () {
    return Expression.ONE.divide(this);
  };

  //TODO: use in Expression#getCoefficients -?
  var variables = function (e) {
    var result = [];
    for (var additions = e.summands(), x = additions.next().value; x != null; x = additions.next().value) {
      for (var multiplications = x.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
        for (var variables = getVariableInternal(y), ve = variables.next().value; ve != null; ve = variables.next().value) {
          if (!(ve.v instanceof Integer)) {
            result.push(ve.v);
          }
        }
      }
    }
    return result;
  };

  //TODO: remove - performance optimization
  var getCommonVariable = function (x, y) {
    var a = variables(x);
    for (var i = 0; i < a.length; i += 1) {
      if (a[i] instanceof NthRoot) {
        return a[i];
      }
    }
    var b = variables(y);
    for (var i = 0; i < b.length; i += 1) {
      if (b[i] instanceof NthRoot) {
        return b[i];
      }
    }
    for (var i = 0; i < a.length; i += 1) {
      for (var j = 0; j < b.length; j += 1) {
        if (a[i].equals(b[j])) {
          return a[i];
        }
      }
    }
    return null;
  };

  // TODO: fix or remove ?
  Expression.prototype.gcd = function (x) {
    //return gcd(this, x, getVariable(this) || getVariable(x));
    return gcd(this, x, getCommonVariable(this, x));
  };
  Expression.prototype.lcm = function (x) {
    return this.divide(this.gcd(x)).multiply(x);
  };

  //TODO: merge with ExpressionParser.js ?!?
  var precedence = {
    binary: {
      ".^": 5,
      "^": 5,
      "*": 3,
      "/": 3,
      "+": 2,
      "-": 2
    },
    unary: {
      "-": 5//HACK
    }
  };

  Expression.Symbol = function (symbol) {
    //Expression.call(this);
    this.symbol = symbol;
  };

  Expression.Symbol.prototype = Object.create(Expression.prototype);

  Expression.Symbol.prototype.compare4Multiplication = function (y) {
    return y.compare4MultiplicationSymbol(this);
  };
  Expression.Symbol.prototype.compare4MultiplicationInteger = function (x) {
    return -1;
  };
  Expression.Symbol.prototype.compare4MultiplicationSymbol = function (x) {
    return x.symbol < this.symbol ? -1 : (this.symbol < x.symbol ? +1 : 0);
  };
  Expression.Symbol.prototype.compare4MultiplicationNthRoot = function (x) {
    return -1;
  };

  Expression.Symbol.prototype.toString = function (options) {
    return this.symbol;
  };


  Expression.prototype.addInteger = function (x) {
    return this.addExpression(x);
  };
  Expression.prototype.multiplyInteger = function (x) {
    return this.multiplyExpression(x);
  };
  Expression.prototype.divideInteger = function (x) {
    return this.divideExpression(x);
  };

  function Integer(value) {
    //Expression.call(this);
    this.value = value;
  }

  Integer.prototype = Object.create(Expression.prototype);
  
  Integer.prototype.powExpression = function (x) {
    var y = this;
    if (x instanceof IdentityMatrix) {
      return new IdentityMatrix(x.symbol);
    }
    if (x instanceof MatrixSymbol) {
      if (y.equals(Expression.ZERO)) {
        return Expression.ONE;
      }
      if (y.equals(Expression.ONE)) {
        return x;
      }
      return new Exponentiation(x, y);//?
    }
    if (y.compareTo(Expression.ZERO) < 0) {
      return Expression.ONE.divide(x.pow(y.negate()));
    }
    if (x instanceof Expression.Matrix) {
      if (y.compareTo(Expression.ONE) > 0) {
        if (Expression.callback != undefined) {
          Expression.callback(new Expression.Event("pow", x, new Expression.Matrix(Matrix.I(1).map(function () { return y; }))));
        }
      }
      return new Expression.Matrix(x.matrix.pow(Number.parseInt(y.toString(), 10)));
    }
    if (y.equals(Expression.ZERO)) {
      return Expression.ONE;
    }
    if (y.equals(Expression.ONE)) {
      return x;
    }

    if (x instanceof Expression.Symbol) {
      return new Exponentiation(x, y);
    }
    if (x instanceof Exponentiation) {
      return x.a.pow(x.b.multiply(y));
    }
    if (x instanceof Integer && (x.compareTo(Expression.ZERO) === 0 || x.compareTo(Expression.ONE) === 0 || x.compareTo(Expression.ONE.negate()) === 0)) {
      return y.remainder(Expression.TWO).compareTo(Expression.ZERO) === 0 ? x.multiply(x) : x;
    }
    // assert(x instanceof Operation || x instanceof Integer);
    return Expression.pow(x, Number.parseInt(y.toString(), 10));
  };

  Integer.prototype.compare4Multiplication = function (y) {
    return y.compare4MultiplicationInteger(this);
  };
  Integer.prototype.compare4MultiplicationInteger = function (x) {
    return x.compareTo(this);
    //return 0;
  };
  Integer.prototype.compare4MultiplicationSymbol = function (x) {
    return +1;
  };
  Integer.prototype.compare4MultiplicationNthRoot = function (x) {
    return +1;
  };

  Integer.prototype.negate = function () {
    return new Integer(BigInteger.unaryMinus(this.value));
  };
  Integer.prototype.compareTo = function (y) {
    return BigInteger.lessThan(this.value, y.value) ? -1 : (BigInteger.greaterThan(this.value, y.value) ? +1 : 0);
  };
  Integer.prototype.add = function (y) {
    return y.addInteger(this);
  };
  Integer.prototype.addInteger = function (x) {
    return new Integer(BigInteger.add(x.value, this.value));
  };
  Integer.prototype.multiply = function (y) {
    return y.multiplyInteger(this);
  };
  Integer.prototype.multiplyInteger = function (x) {
    return new Integer(BigInteger.multiply(x.value, this.value));
  };
  Integer.prototype.divide = function (y) {
    return y.divideInteger(this);
  };
  //! for performance only
  Integer.prototype.divideInteger = function (x) {
    var y = this;
    if (y.equals(Expression.ZERO)) {
      //TODO: fix?
      throw new RangeError("ArithmeticException");
    }
    var gInteger = integerGCD(x, y);
    if (y.compareTo(Expression.ZERO) < 0) {
      gInteger = gInteger.negate();
    }
    x = x.truncatingDivide(gInteger);
    y = y.truncatingDivide(gInteger);
    return y.compareTo(Expression.ONE) === 0 ? x : new Division(x, y);
  };
  Integer.prototype.truncatingDivide = function (y) {
    return new Integer(BigInteger.divide(this.value, y.value));
  };
  Integer.prototype.remainder = function (y) {
    return new Integer(BigInteger.remainder(this.value, y.value));
  };
  Integer.prototype.toString = function (options) {
    return this.value.toString();
  };

  Integer.fromNumber = function (n) {
    return new Integer(BigInteger.BigInt(n));
  };
  Integer.fromString = function (s) {
    return new Integer(BigInteger.BigInt(s));
  };

  Expression.ZERO = Integer.fromNumber(0);
  Expression.ONE = Integer.fromNumber(1);
  Expression.TWO = Integer.fromNumber(2);
  Expression.TEN = Integer.fromNumber(10);


  
  Expression.Matrix = function (matrix) {
    //Expression.call(this);
    this.matrix = matrix;
  };

  Expression.Matrix.prototype = Object.create(Expression.prototype);

  Expression.Matrix.prototype.equals = function (x) {
    if (!(x instanceof Expression.Matrix)) {
      return false;
    }
    return this.matrix.eql(x.matrix);
  };

  Expression.Matrix.prototype.compare4Multiplication = function (x) {
    if (x instanceof Expression.Matrix) {
      return 0;
    }
    return +1;
  };
  Expression.Matrix.prototype.compare4MultiplicationNthRoot = function (x) {
    return +1;
  };

  Expression.Matrix.prototype.multiply = function (y) {
    return y.multiplyMatrix(this);
  };
  Expression.prototype.multiplyMatrix = function (x) {
    var t = getIdentityMatrixCoefficient(this);
    if (t != undefined) {
      return new Expression.Matrix(x.matrix.scale(t));
    }
    return this.multiplyExpression(x);
  };
  Expression.Matrix.prototype.multiplyExpression = function (x) {
    var t = getIdentityMatrixCoefficient(x);
    if (t != undefined) {
      return new Expression.Matrix(this.matrix.scale(t));
    }
    return Expression.prototype.multiplyExpression.call(this, x);
  };
  Expression.Matrix.prototype.multiplyMatrix = function (x) {
    if (Expression.callback != undefined) {
      Expression.callback(new Expression.Event("multiply", x, this));
    }
    return new Expression.Matrix(x.matrix.multiply(this.matrix));
  };
  Expression.Matrix.prototype.compare4MultiplicationSymbol = function (x) {
    return +1;
  };
  Expression.Matrix.prototype.multiplyDivision = Expression.Matrix.prototype.multiplyExpression;
  Expression.Matrix.prototype.add = function (y) {
    return y.addMatrix(this);
  };
  Expression.Matrix.prototype.addMatrix = function (x) {
    return new Expression.Matrix(x.matrix.add(this.matrix));
  };

  var isScalar = function (x) {
    if (x instanceof Integer) {
      return true;
    }
    if (x instanceof Expression.Complex) {
      return true;
    }
    if (x instanceof MatrixSymbol) {
      return false;
    }
    if (x instanceof Expression.Symbol) {
      return true;
    }
    if (x instanceof BinaryOperation) {
      return isScalar(x.a) && isScalar(x.b);
    }
    if (x instanceof Negation) {
      return isScalar(x.b);
    }
    if (x instanceof Expression.Function) {
      return isScalar(x.a);
    }
    return false;//?
  };
  
  Expression.isScalar = isScalar;

  var getIdentityMatrixCoefficient = function (x) {
    var t = undefined;
    if (x instanceof Multiplication && x.b instanceof IdentityMatrix) {
      t = x.a;
    } else if (x instanceof IdentityMatrix) {
      t = Expression.ONE;
    } else if (isScalar(x)) {
      t = x;
    } else if (x instanceof Addition) {
      if (Expression.has(x, IdentityMatrix)) {//TODO: fix
        var ca = getIdentityMatrixCoefficient(x.a);
        var cb = getIdentityMatrixCoefficient(x.b);
        if (ca != undefined && cb != undefined) {
          t = ca.add(cb);
        }
      }
    }
    return t;
  };

  Expression.prototype.addMatrix = function (x) {
    var t = getIdentityMatrixCoefficient(this);
    if (t != undefined) {
      //?
      if (x.matrix.isSquare()) {
        return new Expression.Matrix(Matrix.I(x.matrix.rows()).scale(t)).add(x);
      } else {
        throw new RangeError("NonSquareMatrixException");
      }
    }
    return this.addExpression(x);
  };
  Expression.Matrix.prototype.addExpression = function (x) {
    var t = getIdentityMatrixCoefficient(x);
    if (t != undefined) {
      //?
      if (this.matrix.isSquare()) {
        return this.add(new Expression.Matrix(Matrix.I(this.matrix.rows()).scale(t)));
      } else {
        throw new RangeError("NonSquareMatrixException");
      }
    }
    return Expression.prototype.addExpression.call(this, x);
  };

  Expression.Matrix.prototype.toString = function (options) {
    return this.matrix.toString(setTopLevel(true, options));
  };

  Expression.Matrix.prototype.isExact = function () {
    return this.matrix.isExact();
  };

  function BinaryOperation(a, b) {
    //Expression.call(this);
    this.a = a;
    this.b = b;
  }

  BinaryOperation.prototype = Object.create(Expression.prototype);

  BinaryOperation.prototype.isNegation = function () {
    // TODO: What about NonSimplifiedExpression(s) ?
    //if (this instanceof Multiplication && this.a instanceof NonSimplifiedExpression && this.a.e instanceof Integer && this.a.e.equals(Expression.ONE.negate())) {
    //  return true;
    //}
    return (this instanceof Multiplication && this.a instanceof Integer && this.a.equals(Expression.ONE.negate()));
  };

  var setTopLevel = function (isTopLevel, options) {
    return options == undefined ? {isTopLevel: isTopLevel} : Object.assign({}, options, {isTopLevel: isTopLevel});
  };

  Expression.setTopLevel = setTopLevel;

  BinaryOperation.prototype.toString = function (options) {
    var a = this.a;
    var b = this.b;
    var isSubtraction = false;
    // TODO: check
    /*
    if (Expression.simplification && this instanceof Addition && a.isNegative()) {
      var tmp = b;
      b = a;
      a = tmp;
    }*/

    if (this instanceof Addition && b.isNegative()) {
      isSubtraction = true;
      b = b.negateCarefully();//?
    }
    var fa = a.getPrecedence() + (a.isRightToLeftAssociative() ? -1 : 0) < this.getPrecedence();
    var fb = this.getPrecedence() + (this.isRightToLeftAssociative() ? -1 : 0) >= b.getPrecedence();
    if (options != undefined && options.isTopLevel != undefined && options.isTopLevel === false) {
      fa = fa || a.isUnaryPlusMinus();
    }
    fb = fb || b.isUnaryPlusMinus(); // 1*-3 -> 1*(-3)
    fb = fb || (this.unwrap() instanceof Exponentiation && b.unwrap() instanceof Exponentiation); // 2^3^4
    fa = fa || (this.unwrap() instanceof Exponentiation && a.unwrap() instanceof Expression.Function); // cos(x)^(2+3)
    var s = isSubtraction ? "-" : this.getS();
    //TODO: fix spaces (matrix parsing)
    if (this.isNegation()) {
      // assert(fa === false);
      return "-" + (fb ? "(" : "") + b.toString(setTopLevel(fb, options)) + (fb ? ")" : "");
    }
    return (fa ? "(" : "") + a.toString(setTopLevel(fa || options == undefined || options.isTopLevel, options)) + (fa ? ")" : "") + s + (fb ? "(" : "") + b.toString(setTopLevel(fb, options)) + (fb ? ")" : "");
  };

  //?
  Expression.prototype.unwrap = function () {
    return this;
  };

  function Exponentiation(a, b) {
    BinaryOperation.call(this, a, b);
  }

  Exponentiation.prototype = Object.create(BinaryOperation.prototype);

  function Multiplication(a, b) {
    BinaryOperation.call(this, a, b);
  }

  Multiplication.prototype = Object.create(BinaryOperation.prototype);

  Multiplication.prototype.multiply = function (y) {
    return y.multiplyExpression(this);
  };
  //TODO:
  var compare4Multiplication2 = function (x, y) {//TODO: fix
    if (x instanceof Integer && y instanceof Exponentiation) {
      return -1;//?
    }
    if (x instanceof Exponentiation && y instanceof Integer) {
      return +1;//?
    }
    return compare4Multiplication(getBase(x), getBase(y));
  };

  function Negation(b) {
    //Expression.call(this);
    this.b = b;
  }

  Negation.prototype = Object.create(Expression.prototype);

  Expression.prototype.equalsNegation = function (x) {
    return false;
  };
  Negation.prototype.equalsNegation = function (b) {
    return this.b.equals(b.b);
  };
  Negation.prototype.equals = function (b) {
    return b.equalsNegation();
  };
  Negation.prototype.toString = function (options) {
    var b = this.b;
    var fb = this.getPrecedence() + (this.isRightToLeftAssociative() ? -1 : 0) >= b.getPrecedence();
    fb = fb || b.isUnaryPlusMinus();
    // assert(fa === false);
    return "-" + (fb ? "(" : "") + b.toString(setTopLevel(fb, options)) + (fb ? ")" : "");
  };

  function Subtraction(a, b) {
    BinaryOperation.call(this, a, b);
  }

  Subtraction.prototype = Object.create(BinaryOperation.prototype);

  Subtraction.prototype.getS = function () {
    return "-";
  };

  //

  function Addition(a, b) {
    BinaryOperation.call(this, a, b);
  }

  Addition.prototype = Object.create(BinaryOperation.prototype);
  Addition.prototype.multiply = function (y) {
    return y.multiplyAddition(this);
  };
  Expression.prototype.multiplyAddition = function (x) {
    return x.a.multiply(this).add(x.b.multiply(this));
  };
  Addition.prototype.multiplyExpression = function (x) {
    return x.multiply(this.a).add(x.multiply(this.b));
  };

  function Division(a, b) {
    BinaryOperation.call(this, a, b);
  }

  Division.prototype = Object.create(BinaryOperation.prototype);
  Division.prototype.multiply = function (y) {
    return y.multiplyDivision(this);
  };
  Expression.prototype.multiplyDivision = function (x) {
    return x.a.multiply(this).divide(x.b);
  };
  Division.prototype.multiplyExpression = function (x) {
    return x.multiply(this.a).divide(this.b);
  };
  Division.prototype.add = function (y) {
    return y.addDivision(this);
  };
  Expression.prototype.addDivision = function (x) {
    return x.a.add(this.multiply(x.b)).divide(x.b);
  };
  Division.prototype.addExpression = function (x) {
    return x.multiply(this.b).add(this.a).divide(this.b);
  };
  Division.prototype.divide = function (y) {
    return this.a.divide(this.b.multiply(y));
  };
  Division.prototype.divideExpression = function (x) {
    return x.multiply(this.b).divide(this.a);
  };
  //? not needed, but looks appropriate
  Division.prototype.multiplyAddition = function (x) {
    return x.multiply(this.a).divide(this.b);
  };

  // TODO: move
  Expression.prototype.equals = function (b) {
    throw new RangeError();//?
  };
  Expression.prototype.equalsInteger = function () {
    return false;
  };
  Integer.prototype.equals = function (y) {
    // TODO: fix
    if (y == undefined) {
      return false;
    }
    return y.equalsInteger(this);
  };
  Integer.prototype.equalsInteger = function (x) {
    return x.compareTo(this) === 0;
  };
  Expression.Symbol.prototype.equals = function (b) {
    return b instanceof Expression.Symbol && this.symbol === b.symbol;
  };
  BinaryOperation.prototype.equals = function (b) {
    return b instanceof BinaryOperation && this.getS() === b.getS() && this.a.equals(b.a) && this.b.equals(b.b);
  };

  function MatrixSymbol(symbol) {//TODO: only for square matrix !!!
    Expression.Symbol.call(this, symbol);
  }
  MatrixSymbol.prototype = Object.create(Expression.Symbol.prototype);

  Exponentiation.prototype.inverse = function () {
    return this.pow(Expression.ONE.negate());
  };
  MatrixSymbol.prototype.inverse = function () {//TODO: only for square matrix !!!
    return this.pow(Expression.ONE.negate());
  };
  MatrixSymbol.prototype.compare4Multiplication = function (y) {
    return y.compare4MultiplicationMatrixSymbol(this);
  };
  Expression.prototype.compare4MultiplicationMatrixSymbol = function (x) {
    return +1;
  };
  Expression.Matrix.prototype.compare4MultiplicationMatrixSymbol = function (x) {
    return x instanceof IdentityMatrix ? +1 : -1;//?
  };
  MatrixSymbol.prototype.compare4MultiplicationMatrixSymbol = function (x) {
    var c = Expression.Symbol.prototype.compare4MultiplicationSymbol.call(this, x);
    return c === +1 ? -1 : c;
  };
  MatrixSymbol.prototype.compare4MultiplicationSymbol = function (x) {
    return -1;
  };
  MatrixSymbol.prototype.equals = function (b) {
    return b instanceof MatrixSymbol && Expression.Symbol.prototype.equals.call(this, b);
  };
  //...

  Expression.MatrixSymbol = MatrixSymbol;

  function IdentityMatrix(symbol) {
    MatrixSymbol.call(this, symbol);
  }
  IdentityMatrix.prototype = Object.create(MatrixSymbol.prototype);
  IdentityMatrix.prototype.multiply = function (y) {
    return y.multiplyIdentityMatrix(this);
  };
  IdentityMatrix.prototype.multiplyAddition = function (x) {
    if (isScalar(x)) {
      return new Multiplication(x, this);
    }
    return Expression.prototype.multiplyAddition.call(this, x);
  };
  Expression.prototype.multiplyIdentityMatrix = function (x) {
    return this.multiplyExpression(x);
  };
  IdentityMatrix.prototype.multiplyIdentityMatrix = function (x) {
    return new IdentityMatrix(this.symbol);
  };
  IdentityMatrix.prototype.addMatrix = function (x) {
    return x.add(new Expression.Matrix(Matrix.I(x.matrix.rows())));
  };
  IdentityMatrix.prototype.add = function (y) {
    return y.addIdentityMatrix(this);
  };
  Expression.prototype.addIdentityMatrix = function (x) {
    return this.addExpression(x);//?
  };
  Expression.Matrix.prototype.addIdentityMatrix = function (x) {
    return new Expression.Matrix(Matrix.I(this.matrix.rows())).add(this);
  };
  IdentityMatrix.prototype.multiplyDivision = function (x) {
    if (isScalar(x)) {
      return new Multiplication(x, this);
    }
    return Expression.prototype.multiplyExpression.call(this, x);
  };

  Expression.IdentityMatrix = IdentityMatrix;

  BinaryOperation.prototype.getS = function () {
    throw new Error("abstract");
  };
  Exponentiation.prototype.getS = function () {
    return "^";
  };
  Multiplication.prototype.getS = function () {
    return "*";
  };
  Negation.prototype.getS = function () {
    return "-";
  };
  Addition.prototype.getS = function () {
    return "+";
  };
  Division.prototype.getS = function () {
    return "/";
  };

  Expression.Function = function (name, a) {
    //Expression.call(this);
    this.name = name;
    this.a = a;
  };
  Expression.Function.prototype = Object.create(Expression.prototype);
  Expression.Function.prototype.toString = function (options) {
    //?
    return this.name + "(" + this.a.toString(setTopLevel(true, options)) + ")";
  };
  Expression.Function.prototype.equals = function (b) {
    return b instanceof Expression.Function && this.name === b.name && this.a.equals(b.a);
  };

  Negation.prototype.isUnaryPlusMinus = function () {
    return true;
  };
  BinaryOperation.prototype.isUnaryPlusMinus = function () {
    return this.isNegation();
  };
  Expression.Function.prototype.isUnaryPlusMinus = function () {
    return false;//!
  };
  Expression.prototype.isUnaryPlusMinus = function () {
    return false;
  };
  Integer.prototype.isUnaryPlusMinus = function () {//?
    return this.compareTo(Expression.ZERO) < 0;
  };

  Negation.prototype.getPrecedence = function () {
    return precedence.unary["-"];
  };
  BinaryOperation.prototype.getPrecedence = function () {
    return this.isNegation() ? precedence.unary["-"] : precedence.binary[this.getS()];
  };
  Expression.Function.prototype.getPrecedence = function () {
    return precedence.unary["-"];
  };
  Expression.prototype.getPrecedence = function () {
    return 1000;
  };
  Integer.prototype.getPrecedence = function () {//?
    return this.compareTo(Expression.ZERO) < 0 ? precedence.unary["-"] : 1000;
  };

  Expression.prototype.isNegative = function () {
    var x = this;
    if (x instanceof Integer) {
      return x.compareTo(Expression.ZERO) < 0;
    }
    if (x instanceof Expression.Complex) {
      return x.real.compareTo(Expression.ZERO) < 0 || (x.real.compareTo(Expression.ZERO) === 0 && x.imaginary.compareTo(Expression.ZERO) < 0);
    }
    if (x instanceof Addition) {
      return x.a.isNegative();
      //return x.a.isNegative() && x.b.isNegative();
    }
    if (x instanceof Multiplication) {
      return x.a.isNegative() !== x.b.isNegative();
    }
    if (x instanceof Division) {
      return x.a.isNegative() !== x.b.isNegative();
    }
    if (x instanceof Negation) {
      //return !x.b.isNegative();
      return true;
    }
    return false;
  };

  //TODO: remove
  Expression.prototype.negateCarefully = function () {
    if (this instanceof Integer) {
      return this.negate();
    }
    if (this instanceof Addition) {
      return new Addition(this.a.negateCarefully(), this.b.negateCarefully());
    }
    if (this instanceof Multiplication) {
      return this.b.isNegative() ? new Multiplication(this.a, this.b.negateCarefully()) : (this.a.negateCarefully().equals(Expression.ONE) ? this.b : new Multiplication(this.a.negateCarefully(), this.b));
    }
    if (this instanceof Division) {
      return this.b.isNegative() ? new Division(this.a, this.b.negateCarefully()) : new Division(this.a.negateCarefully(), this.b);
    }
    if (this instanceof Negation) {
      return this.b;//!
    }
    return this.negate();
  };

  // https://en.wikipedia.org/wiki/Nth_root#Simplified_form_of_a_radical_expression
  // https://en.wikipedia.org/wiki/Factorization#Sum.2Fdifference_of_two_cubes

  function NthRoot(name, a, n) {
    Expression.Function.call(this, name, a);
    this.n = n;
  }

  NthRoot.prototype = Object.create(Expression.Function.prototype);

  NthRoot.prototype.compare4Multiplication = function (y) {
    return y.compare4MultiplicationNthRoot(this);
  };
  NthRoot.prototype.compare4MultiplicationInteger = function (x) {
    return -1;
  };
  NthRoot.prototype.compare4MultiplicationSymbol = function (x) {
    return +1;
  };
  NthRoot.prototype.compare4MultiplicationNthRoot = function (x) {
    return x.n < this.n ? -1 : (x.n > this.n ? + 1 : 0);
  };

  NthRoot.prototype.toString = function (options) {
    var fa = this.a.getPrecedence() < this.getPrecedence();
    return (fa ? "(" : "") + this.a.toString(setTopLevel(fa || options == undefined || options.isTopLevel, options)) + (fa ? ")" : "") + "^" + (this.n === 2 ? "0.5" : "(1/" + this.n + ")");
  };

  NthRoot.prototype.getDegree = function () {
    return this.n;
  };

  Expression.prototype._nthRoot = function (n) {
    if (n > 9007199254740991) {
      throw new RangeError("NotSupportedError");
    }
    var x = this;

    if (n === 2) {
      if (x instanceof Addition) {
        if ((x.a instanceof SquareRoot || x.a instanceof Multiplication && x.a.a instanceof Integer && x.a.b instanceof SquareRoot) && x.b instanceof Integer) {
          // (5^0.5+3)/2 or (-5^0.5+3)/2
          var u = x.a;
          var v = x.b;
          // (a+b)^2 = aa+2ab+bb = u+v
          // 2ab = u, b = u/(2a)
          // aa+bb = v, 4aaaa - 4vaa + uu = 0
          // t = sqrt(v*v-u*u);
          // a = sqrt(v+t)/sqrt(2)
          // b = sqrt(v-t)/sqrt(2)
          var tt = v.multiply(v).subtract(u.multiply(u));
          var t = tt.compareTo(Expression.ZERO) >= 0 ? tt.squareRoot() : undefined;
          if (t != undefined && (t instanceof Integer)) {//?
            var aa = v.add(t);
            var a = aa.compareTo(Expression.ZERO) >= 0 ? aa.squareRoot().divide(Expression.TWO.squareRoot()) : undefined;
            if (a != undefined) {
              var b = u.divide(Expression.TWO.multiply(a));
              return a.add(b);
            }
          }
        }
        //TODO: https://brownmath.com/alge/nestrad.htm  - √(√392 + √360)
        if ((x.a instanceof SquareRoot || x.a instanceof Multiplication && x.a.a instanceof Integer && x.a.b instanceof SquareRoot) && 
            (x.b instanceof SquareRoot || x.b instanceof Multiplication && x.b.a instanceof Integer && x.b.b instanceof SquareRoot)) {
          var a = x.a;
          var b = x.b;
          var aa = a.multiply(a);
          var bb = b.multiply(b);
          var g = aa.gcd(bb).squareRoot();
          if (!g.equals(Expression.ONE)) {
            var v = a.divide(g).add(b.divide(g)).squareRoot().multiply(g.squareRoot());
            if (typeof hit === "function") {
              hit({rootFromAddition: x.toString()});
            }
            return v;
          }
        }
      }
    }

    //?
    if (x instanceof NthRoot) {
      if (typeof hit === "function") {
        hit({rootFromRoot: ""});
      }
      return x.a._nthRoot(x.n * n);
    }
    if (x instanceof Division) {
      return x.a._nthRoot(n).divide(x.b._nthRoot(n));
    }
    if (x instanceof Multiplication) {
      return x.a._nthRoot(n).multiply(x.b._nthRoot(n));
    }
    if (x instanceof Integer) {
      if (x.compareTo(Expression.ZERO) < 0) {
        if (n % 2 === 0) {
          throw new RangeError("NotSupportedError");
        }
      }
      if (x.compareTo(Expression.ZERO) === 0) {
        return x;
      }
      var makeRoot = function (i, n) {
        return n === 1 ? i : (n === 2 ? new SquareRoot(i) : (n === 3 ? new CubeRoot(i) : new NthRoot("n" + "-root", i, n)));
      };
      var roots = {};
      var i = x.compareTo(Expression.ZERO) < 0 ? x.negate() : x;
      while (i.compareTo(Expression.ONE) > 0) {
        var d = integerPrimeFactor(i);
        var e = 0;
        while (i.remainder(d).compareTo(Expression.ZERO) === 0) {
          i = i.truncatingDivide(d);
          e += 1;
        }
        while (e !== 0) {
          var g = e;
          while (n % g !== 0) {
            g -= 1;
          }

          var e1 = Math.floor(e / g);
          var k = Math.floor(n / g);
          roots["~" + k] = (roots["~" + k] || Expression.ONE).multiply(Expression.pow(d, e1));

          e = e - g * e1; // e = e % g;
        }
      }
      var y = Expression.ONE;
      //for (var j = 1; j <= n; j += 1) {
      //}
      // `for-in` loop is used to have better performance for "a sparse array"
      for (var jj in roots) {//TODO: fix the iteration order
        if (Object.prototype.hasOwnProperty.call(roots, jj)) {
          var j = Number.parseInt(jj.slice("~".length), 10);
          var r = roots[jj];
          y = y.multiply(makeRoot(r, j));
        }
      }
      return x.compareTo(Expression.ZERO) < 0 ? y.negate() : y;
    }
    if (x instanceof Expression.Matrix) {
      if (typeof hit === "function") {
        hit(n === 2 ? {squareRoot: "matrix"} : (n === 3 ? {cubeRoot: "matrix"} : {nthRoot: "Matrix^(1/" + n + ")"}));
      }
      if (Expression.callback != undefined) {
        Expression.callback(new Expression.Event("diagonalize", x));
      }
      var tmp = Expression.getEigenvalues(x.matrix);
      var eigenvalues = tmp.eigenvalues;
      var multiplicities = tmp.multiplicities;
      if (Expression.sum(multiplicities) === x.matrix.cols()) {
        var tmp2 = Expression.getEigenvectors(x.matrix, eigenvalues);
        var eigenvectors = tmp2.eigenvectors;
        if (eigenvectors.length === x.matrix.cols()) {
          var tmp = Expression.diagonalize(x.matrix, eigenvalues, multiplicities, eigenvectors);
          var L = tmp.L;
          var SL = L.map(function (e, i, j) {
            return i === j ? e._nthRoot(n) : e;
          });
          return new Expression.Matrix(tmp.T.multiply(SL).multiply(tmp.T_INVERSED));
        }
      }
      //TODO: using Jordan normal form -?
    }
    throw new RangeError("NotSupportedError");
  };

  function SquareRoot(a) {
    NthRoot.call(this, "sqrt", a, 2);
  }

  SquareRoot.prototype = Object.create(NthRoot.prototype);

  Expression.prototype.squareRoot = function () {
    return this._nthRoot(2);
  };

  function CubeRoot(a) {
    NthRoot.call(this, "cbrt", a, 3);
  }

  CubeRoot.prototype = Object.create(NthRoot.prototype);

  Expression.prototype.cubeRoot = function () {
    return this._nthRoot(3);
  };

  Expression.Rank = function (matrix) {
    Expression.Function.call(this, "rank", matrix);
  };
  Expression.Rank.prototype = Object.create(Expression.Function.prototype);

  Expression.prototype.rank = function () {
    var x = this;
    if (!(x instanceof Expression.Matrix)) {
      throw new RangeError("NotSupportedError");//?
    }
    //!
    if (Expression.callback != undefined) {
      Expression.callback(new Expression.Event("rank", x));
    }
    return Integer.fromNumber(x.matrix.rank());
  };
  Expression.Determinant = function (matrix) {
    Expression.Function.call(this, "determinant", matrix);
  };
  Expression.Determinant.prototype = Object.create(Expression.Function.prototype);
  Expression.prototype.determinant = function () {
    var x = this;
    if (!(x instanceof Expression.Matrix)) {
      throw new RangeError("NotSupportedError");//?
    }
    //!
    if (Expression.callback != undefined) {
      Expression.callback(new Expression.Event(x.matrix.getDeterminantEventType("determinant").type, x));
    }
    return x.matrix.determinant();
  };
  Expression.RowReduce = function (matrix) {
    Expression.Function.call(this, "row-reduce", matrix);
  };
  Expression.RowReduce.prototype = Object.create(Expression.Function.prototype);
  Expression.prototype.rowReduce = function () {
    var x = this;
    if (!(x instanceof Expression.Matrix)) {
      throw new RangeError("NotSupportedError");//?
    }
    //!
    if (Expression.callback != undefined) {
      Expression.callback(new Expression.Event("row-reduce".type, x));
    }
    //TODO: Matrix.GaussMontante
    return new Expression.Matrix(x.matrix.toRowEchelon(Matrix.GaussJordan, "", null).matrix);
  };
  Expression.Transpose = function (matrix) {
    Expression.Function.call(this, "transpose", matrix);
  };
  Expression.Transpose.prototype = Object.create(Expression.Function.prototype);
  Expression.prototype.transpose = function () {
    var x = this;
    if (!(x instanceof Expression.Matrix)) {
      throw new RangeError("NotSupportedError");//?
    }
    return new Expression.Matrix(x.matrix.transpose());
  };
  Expression.Adjugate = function (matrix) {
    Expression.Function.call(this, "adjugate", matrix);
  };
  Expression.Adjugate.prototype = Object.create(Expression.Function.prototype);
  Expression.prototype.adjugate = function () {
    var x = this;
    if (!(x instanceof Expression.Matrix)) {
      throw new RangeError("NotSupportedError");//?
    }
    if (Expression.callback != undefined) {
      Expression.callback(new Expression.Event("adjugate", x));
    }
    if (x.matrix.rows() === 1 && x.matrix.cols() === 1) {
      return new Expression.Matrix(Matrix.I(1));
    }
    //TODO: optimize
    var C = x.matrix.map(function (element, i, j, matrix) {
      return ((i + j) - 2 * Math.floor((i + j) / 2) === 1 ? Expression.ONE.negate() : Expression.ONE).multiply(matrix.minorMatrix(i, j).determinant());
    });
    var CT = new Expression.Matrix(C.transpose());
    return CT;
    //return new Expression.Matrix(a.matrix.inverse().scale(a.matrix.determinant()));
  };

  Expression.NoAnswerExpression = function (matrix, name, second) {
    Expression.Function.call(this, name, matrix);
    this.second = second;
  };
  Expression.NoAnswerExpression.prototype = Object.create(Expression.Function.prototype);
  //TODO: remove secondArgument (?)
  Expression.prototype.transformNoAnswerExpression = function (name, second) {
    second = second == undefined ? undefined : second;
    if (!(this instanceof Expression.Matrix)) {
      throw new RangeError("NotSupportedError");//?
    }
    if (name === "solve") {
      if (Expression.callback != undefined) {
        Expression.callback(new Expression.Event("solve", this));
      }
    }
    return new Expression.NoAnswerExpression(this, name, second);
  };

  //Expression.NoAnswerExpression.prototype.multiplyExpression =
  //Expression.NoAnswerExpression.prototype.multiplyMatrix =
  //Expression.NoAnswerExpression.prototype.multiplySymbol =
  //Expression.NoAnswerExpression.prototype.multiplyInteger =
  Expression.NoAnswerExpression.prototype.multiply = function () {
    throw new RangeError("NotSupportedError");
  };
  
  //TODO: remove (only for second)
  Expression.NoAnswerExpression.prototype.toString = function (options) {
    if (this.second == undefined) {
      return Expression.Function.prototype.toString.call(this, options);
    }
    return this.a.toString(setTopLevel(true, options)) + " " + this.name + " " + this.second.toString(setTopLevel(true, options));
  };


  Expression.ElementWisePower = function (a, b) {
    BinaryOperation.call(this, a, b);
  };
  Expression.ElementWisePower.prototype = Object.create(BinaryOperation.prototype);
  Expression.ElementWisePower.prototype.getS = function () {
    return ".^";
  };
  Expression.prototype.elementWisePower = function (e) {
    if (!(this instanceof Expression.Matrix)) {
      throw new RangeError("NotSupportedError");//?
    }
    return new Expression.Matrix(this.matrix.map(function (element, i, j) {
      return element.pow(e);
    }));
  };

  Expression.prototype.isRightToLeftAssociative = function () {
    var x = this;
    if (x instanceof Integer) {
      return x.compareTo(Expression.ZERO) < 0;
    }
    if (x instanceof Negation) {
      return true;
    }
    if (x instanceof BinaryOperation) {
      if (x.isNegation()) {
        return true;
      }
      return x instanceof Exponentiation;
    }
    return false;
  };

  var integerPrimeFactor = function (n) {
    return new Integer(primeFactor(n.value));
  };

  //?
  var simpleDivisor = function (e) {
    if (e instanceof Division) {
      throw new RangeError();
    }
    if (e instanceof Expression.Matrix) {
      throw new RangeError();
    }
    if (e instanceof Expression.Symbol) {
      return e;
    }
    if (e instanceof Integer) {
      var x = e;
      var i = x.compareTo(Expression.ZERO) < 0 ? x.negate() : x;
      if (i.compareTo(Expression.ONE) > 0) {
        return integerPrimeFactor(i);
      }
      return null;
    }
    if (e instanceof Expression.Complex) {
      //TODO:
      var g = integerGCD(e.real, e.imaginary);
      var t = simpleDivisor(g);
      if (t != null) {
        return t;
      }
      if (typeof hit === "function") {
        hit({everySimpleDivisor: e.toString()});
      }
      return e;
    }
    var v = getVariable(e);
    if (v != undefined) {
      var np = Polynomial.toPolynomial(e, v);

      var content = np.getContent();
      var t = simpleDivisor(content);
      if (t != null) {
        return t;
      }

      //?
      if (np.getCoefficient(0).equals(Expression.ZERO)) {
        return v;
      }

      if (np.getDegree() >= 2) {
        var t = null;
        np = np.doRationalRootTest(function (sp, q) {
          t = v.multiply(q).subtract(sp);
          return false;
        });
        if (t == null) {
          t = np._getFactorByKroneckersMethod();
          if (t != null) {
            t = t.calcAt(v);
            np = null;
          }
        }
        if (np == undefined) {
          return t;
        }
      }

      e = np.calcAt(v);
      if (e.isNegative()) {//TODO: remove - ?
        e = e.negate();//!?
      }
      return e;
    }
    throw new RangeError();//?
  };

  Expression.everySimpleDivisor = function (e, callback) {//TODO: remove
    while (!e.equals(Expression.ONE) && !e.equals(Expression.ONE.negate())) {
      var t = simpleDivisor(e);
      if (!callback(t)) {
        return false;
      }
      e = e.divide(t);
    }
    return true;
  };

  Expression.everyDivisor = function (e, callback) {
    var divisors = [];
    if (!callback(Expression.ONE)) {
      return false;
    }
    return Expression.everySimpleDivisor(e, function (d) {
      var rec = function (start, n, s) {
        if (n >= 0) {
          var x = divisors[n];
          var d = x.d;
          var e = x.e;
          for (var i = start; i <= e; i += 1) {
            if (!rec(0, n - 1, s.multiply(Expression.pow(d, i)))) {
              return false;
            }
          }
        } else {
          if (!callback(s)) {
            return false;
          }
        }
        return true;
      };
      if (divisors.length === 0 || !divisors[divisors.length - 1].d.equals(d)) {
        divisors.push({
          d: d,
          e: 0
        });
      }
      divisors[divisors.length - 1].e += 1;
      return rec(divisors[divisors.length - 1].e, divisors.length - 1, Expression.ONE);
    });
  };

  Expression.Integer = Integer;
  Expression.NthRoot = NthRoot;
  Expression.SquareRoot = SquareRoot;
  Expression.CubeRoot = CubeRoot;
  Expression.Negation = Negation;
  Expression.Subtraction = Subtraction;
  Expression.BinaryOperation = BinaryOperation;
  Expression.Exponentiation = Exponentiation;
  Expression.Multiplication = Multiplication;
  Expression.Addition = Addition;
  Expression.Division = Division;
  Expression.pow = function (x, count) {
    return pow(x, count, Expression.ONE);
  };

  // --- 




  
  Expression.Equality = function (a, b) {
    BinaryOperation.call(this, a, b);
  };

  Expression.Equality.prototype = Object.create(BinaryOperation.prototype);
  Expression.Equality.prototype.getS = function () {
    return "=";
  };

  function AdditionIterator(e) {
    if (e == undefined) {
      throw new Error();
    }
    this.value = undefined;
    this.e = e;
  }
  AdditionIterator.prototype = Object.create(Iterator.prototype);
  AdditionIterator.prototype.next = function () {
    this.value = this.e instanceof Addition ? this.e.b : this.e;
    this.e = this.e instanceof Addition ? this.e.a : undefined;
    return this;
  };

  function MultiplicationIterator(e) {
    if (e == undefined) {
      throw new Error();
    }
    this.e = e;
  }
  MultiplicationIterator.prototype = Object.create(Iterator.prototype);
  MultiplicationIterator.prototype.next = function () {
    this.value = this.e instanceof Multiplication ? this.e.b : this.e;
    this.e = this.e instanceof Multiplication ? this.e.a : null;
    return this;
  };

  Expression.prototype.summands = function () {
    return new AdditionIterator(this);
  };

  Expression.prototype.factors = function () {
    return new MultiplicationIterator(this);
  };
  
  var splitX = function (e) {
    var scalar = undefined;
    var l = undefined;
    var r = undefined;
    for (var additions = e.summands(), x = additions.next().value; x != null; x = additions.next().value) {
      var state = 0;
      for (var multiplications = x.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
        var factor = y;
        if (!(factor instanceof Integer) && !(factor instanceof Expression.Symbol) && !(factor instanceof Expression.Matrix)) {
          throw new RangeError("NotSupportedError");
        }
        var s = factor instanceof Expression.Symbol ? factor.toString() : "";
        if (s === "X") {
          state = 1;
        } else {
          if (isScalar(factor)) {
            scalar = scalar == undefined ? factor: factor.multiply(scalar);
          } else {
            if (state === 0) {
              r = r == undefined ? factor : factor.multiply(r);
            }
            if (state === 1) {
              l = l == undefined ? factor : factor.multiply(l);
            }
          }
        }
      }
    }
    scalar = scalar == undefined ? Expression.ONE : scalar;
    return {s: scalar, l: l, r: r};
  };
  var groupX = function (a, b) {
    var tmp1 = splitX(a);
    var tmp2 = splitX(b);
    var s1 = tmp1.s;
    var l1 = tmp1.l;
    var r1 = tmp1.r;
    var s2 = tmp2.s;
    var l2 = tmp2.l;
    var r2 = tmp2.r;
    if (r1 == undefined && r2 == undefined) {
      l1 = l1 == undefined ? new IdentityMatrix("I") : l1;
      l2 = l2 == undefined ? new IdentityMatrix("I") : l2;
      return new Multiplication(s1.multiply(l1).add(s2.multiply(l2)), new Expression.Symbol("X"));
    }
    if (l1 == undefined && l2 == undefined) {
      r1 = r1 == undefined ? new IdentityMatrix("I") : r1;
      r2 = r2 == undefined ? new IdentityMatrix("I") : r2;
      return new Multiplication(new Expression.Symbol("X"), s1.multiply(r1).add(s2.multiply(r2)));
    }
    throw new RangeError("NotSupportedError");
  };

  //?
  var getExpressionWithX = function (e) {
    var withX = undefined;
    var withoutX = undefined;
    for (var additions = e.summands(), x = additions.next().value; x != null; x = additions.next().value) {
      var summand = x;
      var hasX = false;
      for (var multiplications = x.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
        var factor = y;
        var factorBase = getBase(factor);
        if (!(factorBase instanceof Integer) && !(factorBase instanceof Expression.Symbol)) {
          if (!(factorBase instanceof Expression.Matrix)) {//?
          if (!(factorBase instanceof NthRoot)) {//?TODO: remove - ?
            throw new RangeError("NotSupportedError");
          }
          }
        }
        if (factorBase instanceof Expression.Symbol) {
          var s = factorBase.toString();
          if (s === "X") {
            if (hasX) {
              throw new RangeError("NotSupportedError");
            }
            hasX = true;
          }
        }
      }
      if (hasX) {
        if (withX != undefined) {
          withX = groupX(withX, summand);
          //throw new RangeError("NotSupportedError");
        } else {
          withX = summand;
        }
      }
      if (!hasX) {
        withoutX = withoutX == undefined ? summand.negate() : withoutX.subtract(summand);
      }
    }
    return {withX: withX, withoutX: withoutX};
  };

  var isConstant = function (e) {
    if (e instanceof Expression.Integer) {
      return true;
    } else if (e instanceof Expression.Complex) {
      return true;
    } else if (e instanceof Expression.NthRoot) {
      return isConstant(e.a);
    } else if (e instanceof Expression.Multiplication) {
      return isConstant(e.a) && isConstant(e.b);
    } else if (e instanceof Expression.Addition) {
      return isConstant(e.a) && isConstant(e.b);
    }
    return false;
  };
  
  Expression.isConstant = isConstant;

  Expression.isSingleVariablePolynomial = function (e) {
    return Expression.getSingleVariablePolynomial(e) != undefined;
  };

  Expression.getSingleVariablePolynomial = function (e, allowMultivariate) {
    allowMultivariate = allowMultivariate == undefined ? false : allowMultivariate;
    if (e instanceof Expression.Division) {
      return undefined;
    }
    //var v = Expression.getVariable(e);
    // To avoid square roots / nth roots:
    var v = getVariableInternal(getLastMultiplicationOperand(getFirstAdditionOperand(e))).next().value.v;
    if (v instanceof NthRoot || v instanceof Integer || v instanceof Expression.Complex) {
      v = undefined;
    }
    if (v == undefined) {
      //throw new Error("undefined");
      return undefined;
    }
    var p = Polynomial.toPolynomial(e, v);
    //TODO: iteration by sparse coefficients
    for (var i = 0; i <= p.getDegree(); i += 1) {
      var c = p.getCoefficient(i);
      if (!isConstant(c)) {
        if (!allowMultivariate) {
          return undefined;
        }
        var pc = Expression.getMultivariatePolynomial(c);
        if (pc == undefined) {
          return undefined;
        }
      }
    }
    return {p: p, v: v};
  };

  //TODO: remove - ?
  Expression.getMultivariatePolynomial = function (e) {
    return Expression.getSingleVariablePolynomial(e, true);
  };
  Expression.isMultivariatePolynomial = function (e) {
    return Expression.getMultivariatePolynomial(e) != undefined;
  };

  // TODO: NotSupportedError
  Expression.prototype.transformEquality = function (b) {
    var e = this.subtract(b);
    var tmp = getExpressionWithX(e);
    var withX = tmp.withX;
    var withoutX = tmp.withoutX;
    if (withX == undefined) {
      if (e.getDenominator() instanceof Integer && !(e.getNumerator() instanceof Expression.Matrix)) {
        //TODO: tests
        var tmp = Expression.getSingleVariablePolynomial(e.getNumerator());
        if (tmp != undefined) {
          var p = tmp.p;
          var v = tmp.v;
          var m = Matrix.Zero(1, p.getDegree() + 1).map(function (e, i, j) {
            return p.getCoefficient(p.getDegree() - j);
          });
          return new Expression.NoAnswerExpression(new Expression.Matrix(m), "polynomial-roots", {variable: v});
        }
      }
      if (e instanceof Expression.Matrix && Expression.SystemOfEquations != null) {
        if (this instanceof Expression.Matrix && b instanceof Expression.Matrix) {
          return Expression.SystemOfEquations.from(this, b);
        }
        //TODO: other things - ?
      }
      throw new RangeError("NotSupportedError");
    }

    if (withoutX == undefined) {
      withoutX = Expression.ZERO;//?
    }
    //console.log(withX.toString() + "=" + withoutX.toString());

    var left = withX;
    var right = withoutX;

    var isToTheLeft = false;
    var x = withX;
    for (var multiplications = x.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
      var factor = y;
      var factorBase = getBase(factor);
      //if (!(factorBase instanceof Integer) && !(factorBase instanceof Expression.Symbol)) {
      //  if (!(factorBase instanceof Expression.Matrix)) {//?
      //    throw new RangeError("NotSupportedError");
      //  }
      //}
      var isX = false;
      if (factorBase instanceof Expression.Symbol) {
        var s = factorBase.toString();
        if (s === "X") {
          isX = true;
          isToTheLeft = true;
        }
      }
      if (!isX) {
        var f = factor.inverse();
        //  console.log(isToTheLeft, f.toString());
        if (isToTheLeft) {
          right = f.multiply(right);
          //left = f.multiply(left);
        } else {
          right = right.multiply(f);
          //left = left.multiply(f);
        }
      } else {
        left = factor;
      }
    }

    //console.log(left.toString() + "=" + right.toString());

    return new Expression.Equality(left, right);
  };

  Expression.simplifications = [];
  Expression.prototype.simplifyExpression = function () {
    var e = this;
    for (var i = 0; i < Expression.simplifications.length; i += 1) {
      e = Expression.simplifications[i](e);
    }
    return e;
  };

  Expression.prototype.isExact = function () {
    return true;
  };

  Expression.getComplexConjugate = function (e) {
    return undefined;
  };

  Expression.Complex = function () {
  };

  Expression.PI = new Expression.Symbol("\u03C0"); // PI
  Expression.E = new Expression.Symbol("e"); // EulerNumber
  Expression.I = new Expression.Symbol("i"); // ImaginaryUnit

  Expression.prototype.addPosition = function () {
    return this;
  };

  //! 2018-09-30
  Expression.SystemOfEquations = function (array) {
    this.array = array;
  };
  Expression.SystemOfEquations.from = function (s, b) {
    return new Expression.NoAnswerExpression({matrix: null}, "system-of-equations", {s: s, b: b});
  };

  Expression.ExponentiationOfMinusOne = function (x, y) {
    Expression.Exponentiation.call(this, x, y);
  };
  Expression.ExponentiationOfMinusOne.prototype = Object.create(Expression.Exponentiation.prototype);
  Expression.ExponentiationOfMinusOne.prototype.divideExpression = function (x) {
    return x.multiply(this);
  };

  Expression.ExponentiationOfImaginaryUnit = function (x, y) {
    Expression.Exponentiation.call(this, x, y);
  };
  Expression.ExponentiationOfImaginaryUnit.prototype = Object.create(Expression.Exponentiation.prototype);
  Expression.ExponentiationOfImaginaryUnit.prototype.divideExpression = function (x) {
    return x.multiply(this);
  };

  //!  
  Expression.Division.prototype.negate = function () {
    return new Expression.Division(this.a.negate(), this.b);
  };
  
  Expression.Polynomial = function (polynomial) {
    this.polynomial = polynomial;
  };
  Expression.Polynomial.prototype = Object.create(Expression.prototype);
  Expression.Polynomial.prototype.equals = function (y) {
    return y.equalsPolynomial(this);
  };
  Expression.Polynomial.prototype.equalsPolynomial = function (x) {
    //TODO: test case
    return x.polynomial.equals(this.polynomial);
  };
  Expression.prototype.equalsPolynomial = function (x) {
    return (x.polynomial.equals(Polynomial.ZERO) && this.equals(Expression.ZERO)) || (x.polynomial.getDegree() === 0 && this.equals(x.polynomial.getCoefficient(0)));
  };
  Expression.Polynomial.prototype.multiply = function (p) {
    return p.multiplyPolynomial(this);
  };
  Expression.Polynomial.prototype.multiplyPolynomial = function (x) {
    return new Expression.Polynomial(x.polynomial.multiply(this.polynomial));
  };
  Expression.Division.prototype.multiplyPolynomial = function (p) {
    return this.multiplyExpression(p);
  };
  Expression.Polynomial.prototype.divide = function (l) {
    if (l.equals(Expression.ONE)) {
      return this;
    }
    return l.dividePolynomial(this);
  };
  Expression.Division.prototype.dividePolynomial = function (p) {
    return this.divideExpression(p);
  };
  Expression.Polynomial.prototype.dividePolynomial = function (x) {
    var y = this;
    var a = x.polynomial;
    var b = y.polynomial;
    if (a.getDegree() < 0 && b.getDegree() >= 0) {
      return new Expression.Polynomial(a);
    }
    var t = undefined;
    while (b.getDegree() >= 0) {
      t = a.divideAndRemainder(b).remainder;
      a = b;
      b = t;
    }
    var gcd = a;
    if (y.polynomial.equals(gcd)) {
      return new Expression.Polynomial(x.polynomial.divideAndRemainder(gcd).quotient);
    }
    return new Expression.Division(new Expression.Polynomial(x.polynomial.divideAndRemainder(gcd).quotient), new Expression.Polynomial(y.polynomial.divideAndRemainder(gcd).quotient));
  };
  Expression.Polynomial.prototype.negate = function () {
    return new Expression.Polynomial(this.polynomial.negate());
  };
  Expression.Polynomial.prototype.add = function (y) {
    return y.addPolynomial(this);
  };
  Expression.Polynomial.prototype.addPolynomial = function (x) {
    return new Expression.Polynomial(x.polynomial.add(this.polynomial));
  };
  Expression.prototype.addPolynomial = function () {
    throw new RangeError();
  };
  Expression.Polynomial.prototype.getPrecedence = function () {
    var d = this.polynomial.getDegree();
    var count = 0;
    for (var i = 0; i <= d; i += 1) {
      if (!this.polynomial.getCoefficient(i).equals(Expression.ZERO)) {
        count += 1;
      }
    }
    return (count < 2 ? (this.polynomial.getLeadingCoefficient().equals(Expression.ONE) ? new Expression.Symbol("x") : new Expression.Multiplication(Expression.ONE, Expression.ONE)) : new Expression.Addition(Expression.ONE, Expression.ONE)).getPrecedence();
  };

  Expression.sum = function (array) {
    var count = 0;
    for (var i = 0; i < array.length; i += 1) {
      count += array[i];
    }
    return count;
  };

Expression.Multiplication.prototype.compare4Multiplication = function (y) {
  var x = this;
  if (y instanceof Addition) {//TODO: fix
    return 0 - y.compare4Multiplication(x);
  }
  var i = x.factors();
  var j = y.factors();
  var a = i.next().value;
  var b = j.next().value;
  while (a != null && b != null) {
    var c = a.compare4Multiplication(b);
    if (c !== 0) {
      return c;
    }
    a = i.next().value;
    b = j.next().value;
  }
  return a != null ? +1 : (b != null ? -1 : 0);
};

Expression.Multiplication.compare4Addition = function (x, y) {
  var i = x.factors();
  var j = y.factors();
  var a = i.next().value;
  var b = j.next().value;
  while (a != null && b != null) {
    var c = a.compare4Addition(b);
    if (c !== 0) {
      return c;
    }
    a = i.next().value;
    b = j.next().value;
  }
  return a != null ? +1 : (b != null ? -1 : 0);
};


// cos(2 * x) * cos(x)
Expression.Multiplication.prototype.compare4MultiplicationSymbol = function (x) {
  return 0 - this.compare4Multiplication(x);
};

Expression.Function.prototype.compare4Addition = function (y) {
  if (y instanceof Expression.Function) {
    return this.name < y.name ? -1 : (y.name < this.name ? +1 : this.a.compare4Addition(y.a));
  }
  return +1;
};

Expression.prototype.compare4AdditionSymbol = function (x) {
  var y = this;
  if (y instanceof Expression.Function) {
    return -1;
  }
  return Expression.prototype.compare4Addition.call(x, y);
};

Expression.Symbol.prototype.compare4Addition = function (y) {
  return y.compare4AdditionSymbol(this);
};

Expression.Function.prototype.compare4Multiplication = function (y) {
  if (y instanceof Expression.NthRoot) {
    return +1;
  }
  if (y instanceof Expression.Function) {
    return this.name < y.name ? -1 : (y.name < this.name ? +1 : this.a.compare4Multiplication(y.a));
  }
  return +1;
};

Expression.Function.prototype.compare4MultiplicationInteger = function (x) {
  return 0 - this.compare4Multiplication(x);
};

Expression.Function.prototype.compare4MultiplicationSymbol = function (x) {
  return -1;//?
};

Expression.Function.prototype.compare4MultiplicationNthRoot = function (x) {
  return 0 - this.compare4Multiplication(x);
};

Expression.Function.prototype.pow = function (y) {
  if (this instanceof Expression.NthRoot) {
    return Expression.prototype.pow.call(this, y);
  }
  if (y instanceof Expression.Integer) {
    if (y.compareTo(Expression.ONE) > 0) {
      return new Expression.Exponentiation(this, y);
    }
    return Expression.prototype.pow.call(this, y);
  }
  throw new RangeError("NotSupportedError");
};

  export default Expression;
