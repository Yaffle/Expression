/*jslint plusplus: true, vars: true, indent: 2 */

// Thanks to Eduardo Cavazos
// see also https://github.com/dharmatech/Symbolism/blob/master/Symbolism/Symbolism.cs
// see also "Computer Algebra and Symbolic Computation: Elementary Algorithms" by Joel S. Cohen

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
  import QuadraticInteger from './QuadraticInteger.js';

  import nthRoot from './nthRoot.js';

  var pow = function (x, count, accumulator) {
    if (count < 0) {
      throw new RangeError();
    }
    if (count > 9007199254740991) {
      throw new RangeError("NotSupportedError");
    }
    return (count < 1 ? accumulator : (2 * Math.floor(count / 2) !== count ? pow(x, count - 1, accumulator.multiply(x)) : pow(x.multiply(x), Math.floor(count / 2), accumulator)));
  };

  // https://stackoverflow.com/a/15302448/839199
  var binomialCoefficient = function (n, k) { // binomail coefficient
    return k === 0 ? Expression.ONE : n.multiply(binomialCoefficient(n.subtract(Expression.ONE), k - 1)).divide(Integer.fromNumber(k));
  };

/*
  var powerOfJordanForm = function (J, N) {
    return J.map(function (e, i, j) {
      if (i > j) {
        return Expression.ZERO;
      }
      if (i === j) {
        return J.e(i, i).equals(Expression.ZERO) ? Expression.ZERO : J.e(i, i).pow(N);
      }
      if (J.e(i, i + 1).equals(Expression.ZERO)) {
        return Expression.ZERO;
      }
      var m = j - i;
      for (var k = 0; k < m; k += 1) {
        if (!J.e(j - 1 - k, j - k).equals(Expression.ONE)) { // outside of a block
          return Expression.ZERO;
        }
      }
      return J.e(i, i).equals(Expression.ZERO) ? Expression.ZERO : binomialCoefficient(N, m).multiply(J.e(i, i).pow(N.subtract(Expression.Integer.fromNumber(m))));
    });
  };
*/

  var matrixInN = function (matrix, n) {
    var condition = -1;
    /*
    if (matrix.isDiagonal()) {
      for (var i = 0; i < matrix.cols(); i += 1) {
        if (matrix.e(i, i).equals(Expression.ZERO)) {
          condition = 0;//?
        }
      }
      var result = matrix.map(function (e, i, j) {
        return i === j ? (e.equals(Expression.ZERO) ? Expression.ZERO : e.pow(n)) : Expression.ZERO;
      });
      var an = new Expression.Matrix(result);
      return condition !== -1 ? new ExpressionWithCondition(an, n, '>', condition) : an;
    }
    */
    /*
    if (matrix.isJordanMatrix()) {
      for (var i = 0; i < matrix.cols(); i += 1) {
        if (matrix.e(i, i).equals(Expression.ZERO)) {
          condition = Math.max(condition, 0);
          // should be Jordan block size minus one (?)
          for (var j = 0; i + j + 1 < matrix.cols(); j += 1) {
            if (matrix.e(i + j, i + j + 1).equals(Expression.ONE)) {
              condition = Math.max(condition, j + 1);
            }
          }
        }
      }
      var an = new Expression.Matrix(powerOfJordanForm(matrix, n));
      if (condition > 0) {//TODO: remove(merge)
        var cases = [];
        cases.push(new ExpressionWithCondition(an, n, '>', condition));
        for (var i = 1; i <= condition; i += 1) {
          cases.push(new ExpressionWithCondition(new Expression.Matrix(matrix.pow(i)), n, '=', i));
        }
        return new Expression.Cases(cases);
      }
      return condition !== -1 ? new ExpressionWithCondition(an, n, '>', condition) : an;
    }
    */
    //!
    var D = matrix.map(function (e, i, j) {
      return i === j ? e : Expression.ZERO;
    });
    var N = matrix.subtract(D);
    if (N.isNilpotent()) {//TODO: fix Matrix#isNilpotent
      if (D.multiply(N).eql(N.multiply(D))) {// D and N commute
        //Note: covers diagonal matrices and Jordan matrices
        for (var k = 0; k < D.cols(); k += 1) {
          if (D.e(k, k).equals(Expression.ZERO)) {
            var index = 1;
            while (!N.pow(index).map(function (e, i, j) { return i !== k ? Expression.ZERO : e; }).isZero()) {
              if (index >= N.cols()) {
                throw new TypeError("assertion");
              }
              index += 1;
            }
            condition = Math.max(condition, index - 1);//?
          }
        }
        var result = Matrix.Zero(N.cols(), N.cols());
        for (var k = 0; k < N.cols(); k += 1) {
          var Dnmk = D.map(function (e, i, j) {
            return i === j ? (e.equals(Expression.ZERO) ? Expression.ZERO : e.pow(n.subtract(Expression.Integer.fromNumber(k)))) : Expression.ZERO;
          });
          result = result.add(Dnmk.multiply(N.pow(k)).scale(binomialCoefficient(n, k)));
        }
        var an = new Expression.Matrix(result);
        if (condition > 0) {//TODO: remove(merge)
          var cases = [];
          cases.push(new ExpressionWithCondition(an, n, '>', condition));
          for (var i = 1; i <= condition; i += 1) {
            cases.push(new ExpressionWithCondition(new Expression.Matrix(matrix.pow(i)), n, '=', i));
          }
          return new Expression.Cases(cases);
        }
        return condition !== -1 ? new ExpressionWithCondition(an, n, '>', condition) : an;
      }
    }
    //!
    var canExponentiate = function (k) {
      if (enableAN && k instanceof Expression.Symbol) {
        return true;
      }
      if (k instanceof Exponentiation && getBase(k) instanceof Integer && getBase(k).compareTo(Expression.ONE) > 0 && isIntegerOrN(getExponent(k).inverse())) {//TODO: remove (no need if to change other codes) of fix `isConstant`?
        return true;
      }
      return isConstant(k) || isConstant(k.divide(Expression.E));
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
    var iteration = -1;
    while (!anm1.eql(anm1previous)) {
      iteration += 1;
      anm1previous = anm1;
      // A^(n) = A^(n-1) * A;
      an = anm1.multiply(a);
      anm1 = an.map(function (e, i, j) {
        var isSymbol = anm1.e(i, j) instanceof Expression.Symbol && anm1.e(i, j).symbol.slice(0, symbolName.length) === symbolName;
        if (!isSymbol) {
          return anm1.e(i, j);//?
        }
        // an: {{1,0,0},{0,1,1+aa_23},{0,0,1}}
        // a_n = n + a_(n-1)
        // a_n = k * a_(n-1) + c * k**(n-(m+1)) * choose(n-1, m)
        // =>
        // a_n = c * k**(n-(m+1)) * choose(n, m+1)
        // Note: choose(n-1, m) + choose(n-2, m) + choose(n-3, m) + ... = choose(n, m+1)
        // choose(n-1, m+1) + choose(n-1, m) = choose(n, m+1)
        if (!(e instanceof Integer)) {
          var m = Polynomial.toPolynomial(e.getNumerator(), n).getDegree();
          var previous = anm1.e(i, j);
          var p = Polynomial.toPolynomial(e.getNumerator(), previous);
          var k = p.getLeadingCoefficient().divide(e.getDenominator());
          if (m !== 0 &&
              p.getDegree() === 1 &&
              a.e(i, j).equals(Expression.ZERO) && //TODO: remove
              (k instanceof Integer || k instanceof Expression.Complex || canExponentiate(k))) { //TODO: fix
            var f = k.pow(n).divide(k.pow(Integer.fromNumber(m + 1))).multiply(binomialCoefficient(n.subtract(Expression.ONE), m));
            var c = e.subtract(k.multiply(previous)).divide(f);
            //TODO: remove `k instanceof Integer`
            if (c instanceof Integer) {//?TODO: ?
              console.log("!", e.toString());
              // a.e(i, j).add()
              return c.multiply(k.pow(n).divide(k.pow(Integer.fromNumber(m + 2))).multiply(binomialCoefficient(n.subtract(Expression.ONE), m + 1)));
            }
          }
        }
        // a_n = a_(n-1)
        if (e.equals(anm1.e(i, j))) {
          return a.e(i, j);
        }
        // a_n = k * a_(n-1) + b => a_n = k**(n - 1) * a_1 + b * (1-k**(n-2))/(1-k)
        if (anm1.e(i, j) instanceof Expression.Symbol && anm1.e(i, j).symbol === symbolName + "_(" + i + "," + j + ")" && !e.equals(Expression.ZERO)) {
          var previous = anm1.e(i, j);
          var p = Polynomial.toPolynomial(e.getNumerator(), previous);
          var k = p.getLeadingCoefficient().divide(e.getDenominator());
          var b = p.getCoefficient(0).divide(e.getDenominator());
          if (!Expression.has(b, Expression.Symbol) && //TODO: !!!
              e.equals(k.multiply(previous).add(b))) {
            var s = k.equals(Expression.ONE) ? b.multiply(n.subtract(Expression.TWO)) : b.multiply(Expression.ONE.subtract(k.pow(n.subtract(Expression.TWO))).divide(Expression.ONE.subtract(k)));
            return k.pow(n.subtract(Expression.TWO)).multiply(a.e(i, j)).add(s);
          }
        }
        if (anm1.e(i, j) instanceof Expression.Symbol && anm1.e(i, j).symbol === symbolName + "_(" + i + "," + j + ")" && e.equals(Expression.ZERO)) {
          //!TODO: conditions.push(iteration); //? n > 0 && n <= 3 , n > 3 - ?
          condition = iteration;
          return Expression.ZERO;
        }
        // a_n = a_(n-1) + b => a_n = a_1 + b*(n-1)
        var sub = e.subtract(anm1.e(i, j));
        if (sub instanceof Integer) {
          return a.e(i, j).add(sub.multiply(n.subtract(Expression.TWO)));
        }
        var dpnm1pda = function (k) { // k**(n-1) + k * a_(n-1)
          if (!canExponentiate(k)) {//TODO: remove
            return Expression.ZERO;// cannot do k.pow(n)
          }
          var previous = anm1.e(i, j);
          return k.pow(n.subtract(Expression.ONE)).add(k.multiply(previous));
        };
        // a_n = d**(n-1) + d * a_(n-1)
        // a_n = d**(n-1) + d * a_(n-1) = 2 * d**(n-1) + d**2 * a_(n-2) = ... = n * d**(n-1) + d**n
        if (!e.equals(Expression.ZERO)) {
          var previous = anm1.e(i, j);
          var p = Polynomial.toPolynomial(e.getNumerator(), previous);
          var k = p.getLeadingCoefficient().divide(e.getDenominator());
          var d = k;
          if (e.equals(dpnm1pda(d))) {
            return d.pow(n.subtract(Expression.TWO)).multiply(n.subtract(Expression.TWO).add(a.e(i, j)));
          }
          var d = k.negate();
          if (e.equals(dpnm1pda(d))) {
            return d.pow(n.subtract(Expression.TWO)).multiply(n.subtract(Expression.TWO).add(a.e(i, j)));
          }
        }

        return anm1.e(i, j);
      });
    }
    for (var i = 0; i < anm1.rows(); i += 1) {
      for (var j = 0; j < anm1.cols(); j += 1) {
        var e = anm1.e(i, j);
        if (e instanceof Expression.Symbol && e.symbol.slice(0, symbolName.length) === symbolName) {
          return undefined;
        }
      }
    }
    if (condition > 0) {
      var cases = [];
      cases.push(new ExpressionWithCondition(new Expression.Matrix(an), n, '>', condition));
      for (var i = 1; i <= condition; i += 1) {
        cases.push(new ExpressionWithCondition(new Expression.Matrix(a.pow(i)), n, '=', i));
      }
      return new Expression.Cases(cases);
    }
    var e = new Expression.Matrix(an);
    return condition !== -1 ? new ExpressionWithCondition(e, n, '>', condition) : e;
  };

  var enableEX = true;
  var enable2X = true;
  var enableEC = true;
  var enableAN = true;

  var isPositive = function (x) {
    if (x instanceof Integer) {
      return x.compareTo(Expression.ZERO) > 0;
    }
    if (x instanceof NthRoot) {
      return isPositive(x.a);//?
    }
    if (x instanceof Multiplication) {
      return isPositive(x.a) === isPositive(x.b);
    }
    if (x instanceof Expression.Symbol) {
      return false;
    }
    if (x instanceof Addition && isPositive(x.a) && isPositive(x.b)) { // 1+sqrt(2)
      return true;
    }
    if (x instanceof Addition && !isPositive(x.a) && !isPositive(x.b)) { // -1-sqrt(2)
      return false;
    }
    if (x instanceof Expression.ExponentiationOfMinusOne) {
      return false;
    }
    if (x instanceof Expression.Exponentiation && getExponent(x) instanceof Expression.Integer) {
      return isPositive(x.a); // x**2
    }
    //TODO: tests, fix for algebraic numbers (?)
    throw new TypeError("!" + x);
  };

  var isIntegerOrN = function (e) {
    if (e instanceof Integer) {
      return true;
    }
    if (e instanceof Expression.Symbol && (e.symbol === "n" || e.symbol === "k")) {
      return true;
    }
    if (e instanceof Expression.Addition || e instanceof Expression.Multiplication || e instanceof Expression.Exponentiation) {
      return isIntegerOrN(e.a) && isIntegerOrN(e.b);
    }
    return false;
  };

  Expression.prototype.powExpression = function (x) {
    var y = this;

    if (y instanceof Expression.Symbol && (y.symbol === "t" || y.symbol === "T")) {
      if (Expression.has(x, MatrixSymbol) || Expression.has(x, Expression.Matrix)) {//TODO: fix
        return x.transpose();
      }
    }
    if (y instanceof Expression.Multiplication && y.a instanceof Expression.Integer && y.b instanceof Expression.Symbol && (y.b.symbol === "t" || y.b.symbol === "T")) {
      if (Expression.has(x, MatrixSymbol) || Expression.has(x, Expression.Matrix)) {//TODO: fix
        return x.pow(y.a).transpose();
      }
    }

    //!
    if (y instanceof Division && y.a instanceof Integer && y.b instanceof Integer && x !== Expression.E && !(x instanceof Expression.Symbol) && !Expression.has(x, Expression.Symbol)) {
      if (typeof hit === "function") {
        hit({powExpression: y.toString()});
      }
      var n = y.b.toNumber();
      //if (n >= 2 && n <= 9007199254740991) {//TODO:
        var q = y.a.truncatingDivide(y.b);
        var r = y.a.subtract(q.multiply(y.b));
        if (q.equals(Expression.ZERO)) {// to avoid multiplication event
          return x.pow(r)._nthRoot(n);
        }
        return x.pow(q).multiply(x.pow(r)._nthRoot(n));
      //}
    }
    //!

    if (x instanceof Expression.Integer && y === Expression.CIRCLE) {
      return new Expression.Degrees(x);
    }

    //!new 2017-05-08
    if (enableEX) {
      if (x === Expression.E || (enable2X && x instanceof Integer && x.compareTo(Expression.ONE) > 0 && integerPrimeFactor(x).compareTo(x) === 0)) {
        var isValid = function (y) {
          if (y instanceof Expression.Symbol) {
            return true;
          }
          if (y instanceof Addition) {
            return isValid(y.a) && isValid(y.b);
          }
          if ((y instanceof Integer || y instanceof NthRoot) && (x === Expression.E || (x instanceof Integer && y instanceof NthRoot))) {//TODO: fix
            return true;
          }
          if (y instanceof Multiplication) {
            for (var f of y.factors()) {
              var b = getBase(f);
              if (!(b instanceof Integer || b instanceof NthRoot || b instanceof Expression.Symbol)) {
                return false;
              }
            }
            return true;
          }
          if (y instanceof Exponentiation) {
            return isValid(y.a) && y.b instanceof Integer;
          }
          if ((x === Expression.E || x instanceof Integer && x.compareTo(Expression.ONE) > 0) && y instanceof Division && y.b instanceof Integer) {//!new 2019-08-08
            return isValid(y.a);//?
          }
          return false;
        };
        if (y.getNumerator() instanceof Addition && (y.getNumerator().a.isNegative() || y.getNumerator().b.isNegative())) { // e**(x-y)
          return Expression.ONE.divide(x.pow(y.getNumerator().a.negate().divide(y.getDenominator())).divide(x.pow(y.getNumerator().b.divide(y.getDenominator()))));
        }
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
        if (xf.equals(x)) {
          if (y instanceof Division && y.b instanceof Integer) {
            var n = y.b.toNumber();
            if (n >= 2 && n <= 9007199254740991) {
              return x.pow(y.a)._nthRoot(n);
            }
          }
          //!new 2020-18-01
          if (y instanceof Division && isIntegerOrN(y.a) && isIntegerOrN(y.b)) {
            if (x.compareTo(Expression.ONE) > 0) {
              return new Expression.Exponentiation(x, y);//?
            }
          }
          //!
        } else {
          return xf.pow(y).multiply(x.divide(xf).pow(y));
        }
      }
    }
    //!

    if (enableEX) {
      //TODO: - ?
      if (x instanceof Integer && x.equals(Expression.ONE)) {
        return Expression.ONE;
      }
      if (x instanceof Division || x instanceof Multiplication && (y.getDenominator().equals(Expression.ONE) || isPositive(x.a) || isPositive(x.b))) {
        if (x instanceof Division) {
          return x.a.pow(y).divide(x.b.pow(y));
        }
        if (enable2X) {
          if (x instanceof Multiplication) {
            return x.a.pow(y).multiply(x.b.pow(y));
          }
        }
      }
    }

    var yn = y.getNumerator();
    var yd = y.getDenominator();
    if (x === Expression.E && yn instanceof Multiplication && yn.a instanceof Expression.Complex && yn.a.real.equals(Expression.ZERO) && yn.b instanceof Expression.Symbol) {
      var t = y.multiply(Expression.I.negate());
      return t.cos().add(Expression.I.multiply(t.sin()));
    }
    if (x === Expression.E && getConstant(yn) instanceof Expression.Complex && yd instanceof Expression.Integer) {
      var c = getConstant(yn);
      if (c.real.equals(Expression.ZERO)) {
        var t = y.multiply(Expression.I.negate());
        t = Expression.has(y, Expression.Symbol) ? t : new Expression.Radians(t);
        return t.cos().add(Expression.I.multiply(t.sin()));
      }
      return x.pow(c.real.divide(yd)).multiply(x.pow(c.imaginary.multiply(Expression.I).multiply(yn.divide(c)).divide(yd)));
    }
    if (x === Expression.E && yn instanceof Expression.Addition && yd instanceof Expression.Integer) {
      return x.pow(yn.a.divide(yd)).multiply(x.pow(yn.b.divide(yd)));
    }

    //TODO:
    if (x instanceof Expression.Matrix && y instanceof Expression.Symbol && (y.symbol === "n" || y.symbol === "k")) {
      if (!x.matrix.isSquare()) {
        throw new RangeError("NonSquareMatrixException");
      }
      var an = matrixInN(x.matrix, y);
      if (an != undefined) {
        //?
        var D = x.matrix.map(function (e, i, j) {
          return i === j ? e : Expression.ZERO;
        });
        var N = x.matrix.subtract(D);
        if (x.matrix.isDiagonal()) {
        //  if (Expression.callback != undefined) {
        //    Expression.callback(new Expression.Event("diagonal-matrix-pow", x));
        //  }
        //} else if (x.matrix.isJordanMatrix()) {
        //  if (Expression.callback != undefined) {
        //    Expression.callback(new Expression.Event("Jordan-matrix-pow", x));
        //  }
        } else if (N.isNilpotent() && D.multiply(N).eql(N.multiply(D))) {
          if (Expression.callback != undefined) {
            Expression.callback(new Expression.Event("DpN-matrix-pow", x));
          }
        }

        return an;
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
              return new Expression.Matrix(tmp.T).multiply(SL).multiply(new Expression.Matrix(tmp.T_INVERSED));
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
              return new Expression.Matrix(tmp.P).multiply(JN).multiply(new Expression.Matrix(tmp.P_INVERSED));
            }
          }
        }
      }
      //!
    }

    if (Expression.ExponentiationOfMinusOne != null) {
      if (x instanceof Integer && x.compareTo(Expression.ZERO) < 0 || x.equals(Expression.E.negate())) {
        if (y instanceof Expression.Symbol && (y.symbol === "n" || y.symbol === "k")) {
          return new Expression.ExponentiationOfMinusOne(Expression.ONE.negate(), y).multiply(x.negate().pow(y));
        }
        if (y instanceof Addition && y.a instanceof Expression.Symbol && (y.a.symbol === "n" || y.a.symbol === "k") && y.b instanceof Integer) {
          return new Expression.ExponentiationOfMinusOne(Expression.ONE.negate(), y.a).multiply(Expression.ONE.negate().pow(y.b)).multiply(x.negate().pow(y));
        }
        if (y instanceof Multiplication) {
          return x.pow(y.a).pow(y.b);
        }
        if (y instanceof Addition && y.b instanceof Integer) {
          return x.pow(y.a).multiply(x.pow(y.b));
        }
      }
    }

    if (Expression.ExponentiationOfImaginaryUnit != null) {
      if (x instanceof Expression.Complex && x.equals(Expression.I.negate())) {
        return Expression.ONE.negate().pow(y).multiply(x.negate().pow(y));
      }
      if (x instanceof Expression.Complex && (x.equals(Expression.I) || x.real.compareTo(Expression.ZERO) > 0 && x.primeFactor().equals(x))) {//TODO: -i, other complex numbers - ?
        if (y instanceof Expression.Symbol && (y.symbol === "n" || y.symbol === "k")) {
          return new Expression.ExponentiationOfImaginaryUnit(x, y);
        }
        if (y instanceof Addition && y.a instanceof Expression.Symbol && (y.a.symbol === "n" || y.a.symbol === "k") && y.b instanceof Integer) {
          //var t = x.pow(y.b);
          //return new Expression.ExponentiationOfImaginaryUnit(x, t instanceof Expression.Complex ? y.a.add(Expression.ONE) : y.a).multiply(t instanceof Expression.Complex ? t.divide(x) : t);
          return x.pow(y.a).multiply(x.pow(y.b));
        }
        if (y instanceof Multiplication) {
          return x.pow(y.a).pow(y.b);
        }
        if (y instanceof Addition && y.b instanceof Integer) {
          return x.pow(y.a).multiply(x.pow(y.b));
        }
      }
      if (x instanceof Expression.Complex && x.real.equals(Expression.ZERO) && !x.imaginary.equals(Expression.ONE)) {//TODO: -i, other complex numbers - ?
        if (y instanceof Expression.Symbol && (y.symbol === "n" || y.symbol === "k")) {
          return x.imaginary.pow(y).multiply(Expression.I.pow(y));
        }
        if (y instanceof Multiplication) {
          return x.pow(y.a).pow(y.b);
        }
        if (y instanceof Addition) {
          return x.pow(y.a).multiply(x.pow(y.b));
        }
      }
      if (x instanceof Expression.Complex) {//TODO: ?
        var pf = x.primeFactor();
        return x.divide(pf).pow(y).multiply(pf.pow(y));//TODO: test ?
      }

      //TODO:
      if (x instanceof Expression.Integer && y instanceof Division && y.getDenominator() instanceof Integer && y.getNumerator() instanceof Expression.Symbol &&
          (y.getNumerator().symbol === "n" || y.getNumerator().symbol === "k")) {
        return x.pow(Expression.ONE.divide(y.getDenominator())).pow(y.getNumerator());
      }
      //?
      if (x instanceof Expression.Integer && y instanceof Division && y.getDenominator() instanceof Integer && y.getNumerator() instanceof Addition) {
        return x.pow(Expression.ONE.divide(y.getDenominator())).pow(y.getNumerator());
      }
    }

    if (x === Expression.E && y instanceof Expression.Matrix) {
      if (!y.matrix.isSquare()) {
        throw new RangeError("NonSquareMatrixException");
      }
      // https://en.wikipedia.org/wiki/Matrix_exponential#Using_the_Jordan_canonical_form
      var tmp = Expression.getEigenvalues(y.matrix);
      var eigenvalues = tmp.eigenvalues;
      var multiplicities = tmp.multiplicities;
      if (Expression.sum(multiplicities) === y.matrix.cols()) {
        var tmp = Expression.getFormaDeJordan(y.matrix, eigenvalues, multiplicities);
        // exp(A) = exp(P*J*P^-1) = P*exp(D + N)*P^-1 = P*exp(D)*exp(N)*P^-1
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
            var summand = p.scale(Expression.ONE.divide(Integer.fromNumber(f)));
            s = s.add(summand);
            p = p.multiply(N);
            k += 1;
            f *= k;
          }
          return s;
        };
        if (Expression.callback != undefined) {
          Expression.callback(new Expression.Event("exponential-using-Jordan-canonical-form", y));
        }
        //if (Expression.callback != undefined) {
        //  Expression.callback(new Expression.Event("Jordan-decomposition", y));
        //}
        return new Expression.Matrix(tmp.P.multiply(D.map(function (e, i, j) {
          return i === j ? Expression.E.pow(e) : Expression.ZERO;
        }).multiply(exp(N))).multiply(tmp.P_INVERSED));
      }
    }

    //!2019-04-22
    if (x instanceof NthRoot && x.a instanceof Integer) {
      return x.a.pow(y.divide(Expression.Integer.fromNumber(x.n)));
    }

    if (enableEC) {
      if (x === Expression.E && isConstant(y) && !has(y, Expression.Complex)) {
        return new Expression.Exponentiation(x, y);
      }
      if ((x instanceof Expression.Symbol || Expression.has(x, Expression.Symbol)) && y instanceof Expression.Division && y.getDenominator() instanceof Integer) {
        return x.pow(y.getNumerator())._nthRoot(y.getDenominator().toNumber());
      }
      if (x instanceof Expression.Symbol && y instanceof Expression.Division && isIntegerOrN(y.getNumerator()) && isIntegerOrN(y.getDenominator())) {
        return new Expression.Exponentiation(x, y);
      }
      if (x === Expression.E && y instanceof Expression.Addition) {
        return x.pow(y.a).multiply(x.pow(y.b));
      }
      if (x instanceof Exponentiation && getBase(x) === Expression.E) {//?
        return getBase(x).pow(getExponent(x).multiply(y));
      }
    }

    if (x instanceof Expression.Matrix && y instanceof Expression.Addition) {
      return x.pow(y.a).multiply(x.pow(y.b));
    }
    if (x instanceof Expression.Matrix && y instanceof Expression.Multiplication && y.a instanceof Expression.Integer) {
      return x.pow(y.a).pow(y.b);
    }
    if (x instanceof Expression.Matrix && y instanceof Expression.Division) {
      //?
      if (y.getNumerator().equals(Expression.ONE) && y.getDenominator() instanceof Expression.Symbol && (y.getDenominator().symbol === "n" || y.getDenominator().symbol === "k")) {
        return x._nthRoot(y.getDenominator());//TODO: ?
      }
      if (isIntegerOrN(y.getNumerator()) && isIntegerOrN(y.getDenominator())) {
        return x.pow(y.getNumerator()).pow(Expression.ONE.divide(y.getDenominator()));
      }
    }

    //?
    if (y instanceof Expression.Symbol && (y.symbol === "n" || y.symbol === "k")) {
      var qi = QuadraticInteger.toQuadraticInteger(x);//?
      if (qi != null && qi.equals(qi.primeFactor()) && qi.a > 0 && qi.D > 0 && qi.isValid()) {
        if (qi.b > 0) {
          return new Expression.ExponentiationOfQuadraticInteger(x, y);
        }
        /*
        if (qi.b < 0) {
          var xc = qi.conjugate().toExpression();
          return x.multiply(xc).pow(y).divide(xc.pow(y));
        }
        */
      }
      /*
      if (qi != null) {
        if (x._nthRoot(2) instanceof Expression.Addition) {//TODO: remove
          var t = x._nthRoot(2).pow(y);
          return new Expression.ExponentiationOfQuadraticInteger(x._nthRoot(2), y.multiply(Expression.TWO));
        }
      }
      */
    }
    /*
    if (y instanceof Multiplication && y.a instanceof Integer && y.b instanceof Expression.Symbol) {
      return x.pow(y.a).pow(y.b);
    }
    */
    //?

    if (enableAN) {
      if (x instanceof Expression.Symbol && y instanceof Expression.Symbol && (y.symbol === "n" || y.symbol === "k")) {
        return new Expression.Exponentiation(x, y);
      }
      if (y instanceof Addition && y.a instanceof Expression.Symbol && (y.a.symbol === "n" || y.a.symbol === "k") && y.b instanceof Integer) {
        if (y.b.isNegative()) {
          return x.pow(y.a).multiply(x.pow(y.b));
        }
        return new Expression.Exponentiation(x, y);
      }
      if (y instanceof Multiplication && y.a instanceof Expression.Integer && y.b instanceof Expression.Symbol && (y.b.symbol === "n" || y.b.symbol === "k")) {
        return new Expression.Exponentiation(x, y);
      }
      if (x instanceof Exponentiation && getBase(x) instanceof Integer && getBase(x).compareTo(Expression.ONE) > 0 && isIntegerOrN(getExponent(x).getNumerator()) && isIntegerOrN(getExponent(x).getDenominator())) {//TODO: FIX
        return getBase(x).pow(getExponent(x).multiply(y));
      }
    }

    if (Expression.has(x, Expression.Sin) || Expression.has(x, Expression.Cos)) {
      return Expression._replaceBySinCos(Expression._replaceSinCos(x).pow(y));
    }

    if (x === Expression.E && y instanceof Expression.Logarithm) {
      return y.a;
    }
    if (x === Expression.E && y instanceof Expression.Multiplication && y.a instanceof Integer && y.b instanceof Expression.Logarithm) {
      return x.pow(y.b).pow(y.a);
    }
    if (x === Expression.E && y instanceof Expression.Multiplication && y.a instanceof Expression.Symbol && y.b instanceof Expression.Logarithm) {
      return x.pow(y.b).pow(y.a);
    }
    if (x instanceof Expression.Integer && Expression.has(y, Expression.Logarithm)) {//?
      return Expression.E.pow(x.logarithm().multiply(y));
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
      /*
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
      */
      return 0;
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

    if (x instanceof Expression.Matrix) {
      return -1;
    }
    if (y instanceof Expression.Matrix) {
      return +1;
    }

    //!2019-02-18
    if (x instanceof Integer && y instanceof Expression.Complex) {
      return -1;//?
    }
    if (x instanceof Expression.Complex && y instanceof Integer) {
      return +1;//?
    }
    //!

    if (x.equals(y)) {
      return 0;//!
    }

    if (x instanceof Expression.Exponentiation || y instanceof Expression.Exponentiation) {
      return getBase(x).compare4Addition(getBase(y)) || (0 - getExponent(x).compare4Addition(getExponent(y)));
    }

    //!new 2017-02-10
    if (y instanceof Expression.Symbol) {
      return -1;
    }

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
    //TODO: ?
    //if (x instanceof NthRoot) {
    //  return x.a;
    //}
    return x instanceof Exponentiation ? x.a : x;
  };
  var getExponent = function (x) {
    //TODO: ?
    //if (x instanceof NthRoot) {
    //  return Expression.Integer.fromNumber(x.n);
    //}
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
      for (var y of x.factors()) {
        var t = getConstant(y);
        c = c == undefined ? t : t.multiply(c);
      }
      if (c != undefined) {
        return c;
      }
    } else if (e instanceof Addition) { // -5*x+15
      var c = undefined;
      for (var x of e.summands()) {
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
      for (var y of x.factors()) {
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
    } else if (x instanceof Addition) {
      return x.divide(getConstant(x));
    }
    return x;
  };

  Expression.getConstant = getConstant;

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
    /*
    if (x instanceof IdentityMatrix && y instanceof MatrixSymbol) {
      return y;
    }
    if (y instanceof IdentityMatrix && x instanceof MatrixSymbol) {
      return x;
    }
    */
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
      if (x.a instanceof Integer && y.a instanceof Exponentiation) {//TODO: fix
        return new Multiplication(x, y);
      }
      return x.a.multiply(y.a).squareRoot();
    }
    if (x instanceof CubeRoot && y instanceof CubeRoot) {
      return x.a.multiply(y.a).cubeRoot();
    }
    if (x instanceof NthRoot && y instanceof NthRoot) {
      //if (x.n === y.n) {
      //  return x.a.multiply(y.a)._nthRoot(x.n);
      //}
      var ng = ngcd(x.n, y.n);
      if (!(x.a instanceof Integer) || !x.a.gcd(y.a).equals(Expression.ONE)) {
        return Expression.pow(x.a, y.n / ng).multiply(Expression.pow(y.a, x.n / ng))._nthRoot(x.n / ng * y.n);
      }
      return x.n < y.n ? new Multiplication(x, y) : (x.n > y.n ? new Multiplication(y, x) : x.a.multiply(y.a)._nthRoot(x.n));
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

/*
    // TODO: remove
    if (x instanceof Expression.Complex && y instanceof Expression.ExponentiationOfImaginaryUnit) {
      //!hack
      //TODO: remove
      if (x.real.equals(Expression.ZERO)) {
        if (!x.equals(Expression.I)) {
        //if (x.primeFactor().equals(y.a)) {//TODO: fix
        if (y.a.equals(Expression.I)) {//TODO: fix
        return x.imaginary.multiply(y.multiply(Expression.I));
        }
        //}
        }
      } else {
        if (getBase(y).equals(Expression.I)) {//TODO: remove
          return x.imaginary.multiply(Expression.I).multiply(y).add(x.real.multiply(y));
        }
      }
    }
    if (x instanceof Expression.ExponentiationOfImaginaryUnit && y instanceof Expression.Complex) {
      //!hack
      //TODO: remove
      if (y.real.equals(Expression.ZERO)) {
        if (!y.equals(Expression.I)) {
        return y.imaginary.multiply(x.multiply(Expression.I));
        }
      } else {
        if (getBase(x).equals(Expression.I)) {//TODO: remove
        return y.imaginary.multiply(Expression.I).multiply(x).add(y.real.multiply(x));
        }
      }
    }
*/

    //var cmp = x instanceof Expression.Complex && y instanceof Expression.ExponentiationOfImaginaryUnit && !x.equals(getBase(y)) ? -1 : (x instanceof Expression.ExponentiationOfImaginaryUnit && y instanceof Expression.Complex && !y.equals(getBase(x)) ? +1 : compare4Multiplication(getBase(x), getBase(y)));
    var cmp = x instanceof Expression.Complex && y instanceof Expression.ExponentiationOfImaginaryUnit ? -1 : (x instanceof Expression.ExponentiationOfImaginaryUnit && y instanceof Expression.Complex ? +1 : compare4Multiplication(getBase(x), getBase(y)));
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
  if (typeof globalThis.Symbol === "function") {
    Iterator.prototype[globalThis.Symbol.iterator] = function () {
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
    //!new 2020-02-13
    if (a instanceof Expression.Matrix || b instanceof Expression.Matrix) {
      if (y instanceof Integer && x instanceof Multiplication) {
        return +1;//?
      }
      if (x instanceof Integer && y instanceof Multiplication) {
        return -1;//?
      }
      return 0;
    }
    //!
    return a != null ? +1 : (b != null ? -1 : 0);
  };

  var addSimilar = function (x, y) {
    var c1 = getConstant(x);//TODO: remove (?)
    var c2 = getConstant(y);
    var i = termFactors(getTerm(x));
    var j = termFactors(getTerm(y));
    var a = i.next().value;
    var b = j.next().value;
    var result = Expression.ONE;
    while (a != null || b != null) {
      var f = null;
      if (a instanceof Expression.Matrix || b instanceof Expression.Matrix) {
        f = (a == null ? c1 : a.multiply(c1)).add(b == null ? c2 : b.multiply(c2));
        c1 = Expression.ONE;
        c2 = Expression.ONE;
      } else {
        if (!a.equals(b)) {
          throw new TypeError();
        }
        f = a;
      }
      result = f.multiply(result);//!TODO: depends on the iteration order !!!
      a = i.next().value;
      b = j.next().value;
    }
    return result;
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

    //!new 2019-09-30
    /*
    if (x instanceof Expression.Addition && y instanceof Expression.Matrix) {//TODO:
      return x.a.add(x.b.add(y));
    }
    if (x instanceof Expression.Matrix && y instanceof Expression.Addition) {
      return x.add(y.a).add(y.b);
    }
    if (x instanceof Expression.Addition && y instanceof Expression.Addition) {
      if (x.b instanceof Expression.Matrix && y.b instanceof Expression.Matrix) {
        return x.a.add(x.b.add(y));
      }
    }
    */

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
        if (Expression.has(a, Expression.Matrix) || Expression.has(b, Expression.Matrix)) {
          var last = addSimilar(a, b);
          if (!last.equals(Expression.ZERO)) {
            s.push(last);
          }
        } else {
          var constant = getConstant(a).add(getConstant(b));
          var term = getTerm(a);
          var last = term == undefined ? constant : constant.multiply(term);
          if (!last.equals(Expression.ZERO)) {
            s.push(last);
          }
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

  var divideByInteger = function (e, f) {
    if (f.equals(Expression.ZERO)) {
      throw new RangeError("ArithmeticException");
    }
    var result = Expression.ZERO;
    for (var x of e.summands()) {
      var rest = Expression.ONE;
      var t = undefined;
      // TODO: check, fix?
      for (var y of x.factors()) {
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
    for (var x of e.summands()) {
      var d = Expression.ZERO;
      var c = Expression.ONE;
      for (var y of x.factors()) {
        var t = y;
        for (var ve of getVariableInternal(t)) {
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
      throw new TypeError();
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
      throw new TypeError();
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
    } else if (x instanceof Expression.NthRoot) {//!new 2019-11-30
      value = {v: new Exponentiation(this.v, x), e: Expression.ONE};
    } else if (x instanceof Expression.Exponentiation) {//!new 2019-12-01
      value = {v: new Exponentiation(this.v, x), e: Expression.ONE};
    } else if (x instanceof Multiplication && x.a instanceof Integer) {
      value = {v: new Exponentiation(this.v, x.b), e: x.a};
    } else if (x instanceof Multiplication) {
      value = {v: new Exponentiation(this.v, getTerm(x)), e: getConstant(x)};
    } else if (x instanceof Integer) {
      value = {v: this.v, e: x};
    } else if (x instanceof Expression.Division && x.a instanceof Integer && x.b instanceof Integer) {//!new 2019-06-16
      value = {v: this.v, e: x};
    } else if (x instanceof Expression.Division && x.a instanceof NthRoot && x.b instanceof Integer) {//!new 2019-12-01
      value = {v: new Exponentiation(this.v, x.a), e: x.getDenominator()};
    } else if (x instanceof Expression.Division && x.a instanceof Multiplication && x.a.a instanceof Integer && x.b instanceof Integer) {//!new 2019-06-16
      if (this.v instanceof Integer) {
        value = {v: new Exponentiation(this.v, x.a.b), e: x.divide(x.a.b)};
      } else {
      value = {v: this.v, e: x};//?
      }
    } else {
      // this.v instanceof Integer &&
      // && x.a instanceof Expression.Symbol
      if (x instanceof Division && x.b instanceof Integer) {
        var t = getTerm(x.a);
        value = {v: new Exponentiation(this.v, t), e: x.divide(t)};
      } else {
        if (x instanceof Division && x.a instanceof Integer && x.b instanceof Expression.Symbol && (this.v instanceof Integer || this.v instanceof Expression.Symbol)) {//TODO: fix ? 2**(1/n)
          var t = Expression.ONE.divide(x.b);
          value = {v: new Exponentiation(this.v, t), e: x.divide(t)};
        } else {
      throw new RangeError();
        }
      }
    }
    this.value = value;
    return this;
  };

  function NumeratorSummandsIterator(e) {
    this.internal = e.getNumerator().summands();
    this.denominator = e.getDenominator();
  }
  NumeratorSummandsIterator.prototype = Object.create(Iterator.prototype);
  NumeratorSummandsIterator.prototype.next = function () {
    var next = this.internal.next().value;
    this.value = next == null ? null : next.divide(this.denominator);
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
        var additions = new NumeratorSummandsIterator(e);
        return new VariablesIterator(v, additions);
      }
    }
    //!
    return new VIterator({v: v, e: e});
  };

  var getVariable = function (e) {
    //? square roots at first
    for (var x of e.summands()) {
      for (var y of x.factors()) {
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
    var a = x;
    var b = y;
    while (!b.equals(Expression.ZERO)) {
      var t = a.remainder(b);
      a = b;
      b = t;
    }
    return a;
  };

  //var getIntegerContent = function (x) {
  //  return x instanceof Expression.Complex ? integerGCD(x.real, x.imaginary) : x;
  //};

  var complexGCD = function (a, b) {//?
    // return integerGCD(getIntegerContent(a), getIntegerContent(b));
    var x = integerGCD(a, b);
    if (x instanceof Expression.Complex) {
      //TODO: ?
      if (x.real.compareTo(Expression.ZERO) === 0) {
        return x.imaginary;
      }
    }
    if (x instanceof Expression.Integer) {
      if (x.compareTo(Expression.ZERO) < 0) {
        x = x.negate();
      }
    }
    return x;
  };

  // http://www-troja.fjfi.cvut.cz/~liska/ca/node33.html
  var gcd = function (a, b, v) {
    if (v == undefined) {
      return complexGCD(getConstant(a), getConstant(b));
    }

    var r = getReplacement(a, getReplacement(b, v));
    if (!r.equals(v)) {
      return substitute(substitute(a, v, r, inverseReplacement(r, v)).gcd(substitute(b, v, r, inverseReplacement(r, v))), v, inverseReplacement(r, v), r);
    }

    return Polynomial.polynomialGCD(Polynomial.toPolynomial(a, v), Polynomial.toPolynomial(b, v)).calcAt(v);
  };

  // ! new 21.12.2013 (square roots)

  var getConjugateFactor = function (e) {
    var r = 1 / 0;
    var p = undefined;
    for (var x of e.summands()) {
      for (var y of x.factors()) {
        if (y instanceof NthRoot) {
          var degree = y.getDegree();
          var i = y.a instanceof Integer ? y.a : null;
          if (i == null) {
            i = QuadraticInteger.toQuadraticInteger(y.a);
          }
          if (i == null) {
            throw new TypeError();
          }
          if (r > degree && r === 1 / 0) {
            //TODO: check !!!
            p = i;
            r = degree;
          } else if (r === degree) {
            //TODO: assert(y instanceof Integer)
            if (i != null) {
              var z = integerGCD(p, i);
              if (!z.isUnit()) {
                p = z;//!
              }
              if (z.isUnit() && !(z instanceof Integer) && (p.isUnit() || i.isUnit())) {
                p = z;//?
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
    if (e instanceof Integer) {
      //optimize to not stop the debugger at common code
      return null;
    }
    if (e instanceof SquareRoot) {
      //optimize to not stop the debugger at common code
      return e;
    }
    if (e instanceof CubeRoot) {
      return e._pow(2);
    }
    //!2019-10-20 a workaround
    if (e instanceof Addition &&
        e.a instanceof Multiplication && e.a.a instanceof Integer && e.a.b instanceof CubeRoot &&
        e.b instanceof Multiplication && e.b.a instanceof Integer && e.b.b instanceof CubeRoot) {
      // (aa-ab+bb)
      return e.a._pow(2).subtract(e.a.multiply(e.b)).add(e.b._pow(2));
    }
    //!
    //?
    //if (e instanceof Expression.Exponentiation && getExponent(e).equals(Expression.ONE.divide(Expression.TWO))) {//TODO: FIX
      //return e;//?
    //}
    //?

  //TODO: fix
  //if (true) return undefined;
    var tmp = getConjugateFactor(e);
    var p = tmp.p;
    var r = tmp.r;
    if (p == undefined) {
      return undefined;
    }
    var polynomial = Polynomial.of();
    for (var x of e.summands()) {
      var degree = 0;
      var coefficient = Expression.ONE;
      for (var y of x.factors()) {
        if (y instanceof NthRoot && r === y.getDegree()) {
          var i = y.a instanceof Integer ? y.a : null;
          if (i == null) {
            i = QuadraticInteger.toQuadraticInteger(y.a);
          }
          if (i != null) {
            var j = 0;
            var a = i;
            if (p.isUnit()) {
              if (a.isUnit()) {
                j = 1;
                a = a.truncatingDivide(p);
                // a == 3+sqrt(2), p == 1+sqrt(2)
                while (a.primeFactor().equals(p)) {
                  j += 1;
                  a = a.truncatingDivide(p);
                }
              } else {
                if (!(a instanceof Integer)) {
                //TODO: ?
                //throw new TypeError();
                  var tmp = a.truncatingDivide(p);
                  if (tmp.a * tmp.b >= 0) {
                    j = 1;
                    a = tmp;
                  }
                }
              }
            } else {
              while (a.isDivisibleBy(p) && j < r) {//?
                a = a.truncatingDivide(p);
                j += 1;
              }
            }
            a = a.toExpression();
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
    x = x.toExpression();
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
    for (var x of e.summands()) {
      var v = undefined;
      var c = Expression.ONE;
      var NO_VARIABLE = "";
      for (var y of x.factors()) {
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

  var inverseReplacement = function (e, v) {
    var t = v;
    while (!e.equals(v)) {
      if (e instanceof Expression.Exponentiation && e.b instanceof Multiplication && v instanceof Exponentiation) {
        t = t.pow(e.b.a.inverse());
        e = e.pow(e.b.a.inverse());
      } else if (e instanceof Expression.Exponentiation) {
        t = t.pow(getExponent(e).inverse());
        e = getBase(e);
      } else if (e instanceof Addition) {
        if (!(e.b instanceof Integer)) {
          throw new RangeError();
        }
        t = t.subtract(e.b);
        e = e.a;
      } else if (e instanceof Multiplication) {
        if (!(e.a instanceof Integer)) {
          throw new RangeError();
        }
        t = t.divide(e.a);
        e = e.b;
      } else if (e instanceof Division) {
        if (!(e.b instanceof Integer)) {
          throw new RangeError();
        }
        t = t.multiply(e.b);
        e = e.a;
      } else {
        if ((Expression.E === e) || (e instanceof Integer && e.compareTo(Expression.ONE) > 0) && getBase(v).equals(e)) {//!new 2019-09-23
          t = t.pow(getExponent(v));
          e = v;
        } else {
          throw new TypeError();
        }
      }
    }
    return t;
  };

  var h = function (e, n, q) {
    for (var x of e.summands()) {
      for (var y of x.factors()) {
        if (getBase(q).equals(getBase(y))) {
          n = n.lcm(getExponent(y).getDenominator());
        }
      }
    }
    return n;
  };

  var getReplacement = function (e, v) {
    for (var x of e.summands()) {
      for (var y of x.factors()) {
        if (y instanceof Expression.Exponentiation && (y.a instanceof Expression.Symbol || y.a instanceof Integer && y.a.compareTo(Expression.ONE) > 0) && y.b instanceof Expression.Division) {
          if (getBase(v).equals(y.a)) {
            //v = new Expression.Exponentiation(y.a, y.b.b.lcm(getExponent(v).getNumerator()));
            v = new Expression.Exponentiation(y.a, getConstant(y.b.b).lcm(getExponent(v).getNumerator()).divide(getExponent(v).getDenominator()));
          }
        } else {
          //!TODO: fix
          if (y instanceof Expression.Exponentiation && Expression.has(y.a, Expression.Symbol) && y.b instanceof Expression.Division) {
            if (getBase(v).equals(y.a)) {
              v = new Expression.Exponentiation(y.a, y.b.b.lcm(getExponent(v).getNumerator()));
            } else if (getBase(v) instanceof Expression.Symbol) {
              if ((y.a instanceof Expression.Addition && y.a.a.divide(getBase(v)) instanceof Integer && y.a.b instanceof Integer) && (y.b instanceof Division && y.b.a instanceof Integer && y.b.b instanceof Integer)) {
                var n = y.b.getDenominator();
                n = h(e, n, y);//!
                //TODO: ?
                //debugger;
                // sqrt(y.a.a + y.a.b) = t
                // k * v = y.a.a = t**2 - y.a.b
                // v = (t**2 - y.a.b) / (y.a.a / v)
                var t = getBase(v).pow(n).subtract(y.a.b).divide(y.a.a.divide(getBase(v)));
                return t;
              } else {
                //TODO: test
                //throw new TypeError();
              }
            }
          }
        }
      }
    }
    return v;
  };

  var substitute = function (e, a, b, inv) {
    if (e.equals(a)) {
      return b;
    }
    if (e instanceof Expression.Exponentiation) {
      if (e.equals(inv)) {
        return a;
      }
      if (getBase(e).equals(getBase(inv))) {//!new 2019-09-23
        return a.pow(getExponent(e).divide(getExponent(inv)));
        //TODO: add an assertion below
      }
    }

    if (e instanceof Expression.Addition) {
      return substitute(e.a, a, b, inv).add(substitute(e.b, a, b, inv));
    }
    if (e instanceof Expression.Multiplication) {
      return substitute(e.a, a, b, inv).multiply(substitute(e.b, a, b, inv));
    }
    if (e instanceof Expression.Exponentiation) {
      var x = substitute(e.a, a, b, inv);
      var y = substitute(e.b, a, b, inv);
      //console.log(x + ' ' + y + ' ' + a + ' ' + b + ' ' + inv);
      if (x instanceof Expression.Exponentiation &&
          getBase(x).equals(getBase(inv)) &&
          getExponent(inv).getDenominator().remainder(Expression.TWO).equals(Expression.ZERO)) {
        //TODO: FIX
        return getBase(x).pow(getExponent(x).multiply(y));
      }
      return x.pow(y);
    }
    //return e;

    //! for sin.js:
    if (e instanceof Division) {
      return substitute(e.a, a, b, inv).divide(substitute(e.b, a, b, inv));
    }
    if (e instanceof Expression.Sin) {
      return substitute(e.a, a, b, inv).sin();
    }
    if (e instanceof Expression.Cos) {
      return substitute(e.a, a, b, inv).cos();
    }
    return e;
  };

  Expression._substitute = substitute;

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

if (simplifyIdentityMatrixPower) {
    if (x instanceof Multiplication && x.b instanceof IdentityMatrix) {
      return x.b.equals(y) ? x.a : x.a.divide(y).multiply(x.b);
    } else if (x instanceof IdentityMatrix) {
      return Expression.ONE.divide(y).multiply(x);
    }
    if (y instanceof Multiplication && y.b instanceof IdentityMatrix) {
      return x.divide(y.a).multiply(y.b);
    } else if (y instanceof IdentityMatrix) {
      return x.multiply(y);
    }
}

    if (has(y, MatrixSymbol)) {
      //!?
      var tmp = getBase(y) instanceof MatrixSymbol ? y.pow(Expression.ONE.negate()) : new Expression.Exponentiation(y, Expression.ONE.negate());
      if (x.equals(Expression.ONE)) {
        if (y instanceof Multiplication) {
          //!?
          //TODO: info about properties of the Matrix Inversion
          return x.multiply(Expression.ONE.divide(y.b).multiply(Expression.ONE.divide(y.a)));
        }
        return tmp;
      }
      //return x.multiply(tmp);
      //return new Expression.Multiplication(x, tmp);
      //throw new RangeError("NotSupportedError");
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
        if (e.equals(Expression.ONE)) {
          throw new TypeError(); // "assertion"
        }
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


    //!2019-06-16
    if (v != null) { // e**(1/2)
      v = getVariable(v);
      var r = getReplacement(y, getReplacement(x, v));
      if (!r.equals(v)) {
        var ir = inverseReplacement(r, v);
        var a = substitute(x, v, r, ir);
        var b = substitute(y, v, r, ir);
        //console.log(a + ' ' + b);
        var t = a.divide(b);
        a = substitute(t.getNumerator(), v, ir, r);
        b = substitute(t.getDenominator(), v, ir, r);
        return b.equals(Expression.ONE) ? a : new Expression.Division(a, b);
      }
    }


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
      var g = Polynomial.polynomialGCD(px, py);
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
    if (has(y, MatrixSymbol)) {//?
      return new Expression.Multiplication(x, new Expression.Exponentiation(y, Expression.ONE.negate()));//?
    }//?
    return new Division(x, y);
  };

  function Expression() {
    throw new TypeError("Do not call for better performance");
  }

  Expression.callback = undefined;
  Expression.Event = function (type, data, second) {
    second = second == undefined ? undefined : second;
    this.type = type;
    this.data = data;
    this.second = second;
  };

  Expression.prototype.compare4Multiplication = function (y) {
    throw new TypeError(this.toString());
  };
  Expression.prototype.compare4MultiplicationInteger = function (x) {
    throw new TypeError();
  };
  Expression.prototype.compare4MultiplicationSymbol = function (x) {
    throw new TypeError();
  };
  Expression.prototype.compare4MultiplicationNthRoot = function (x) {
    throw new TypeError();
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
    //!2019-04-22
    if (this.equals(y) && !(y instanceof Expression.Matrix)) { //!TODO: remove - a hack to avoid some exceptions
      //if (this instanceof IdentityMatrix) {
      //  return this;
      //}
      return Expression.ONE;
    }
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
  Expression.prototype.exp = function () {
    return Expression.E.pow(this);
  };


  //TODO: use in Expression#getCoefficients -?
  var variables = function (e) {
    var result = [];
    for (var x of e.summands()) {
      for (var y of x.factors()) {
        for (var ve of getVariableInternal(y)) {
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
    var b = variables(y);
    /*
    for (var i = 0; i < a.length; i += 1) {
      if (a[i] instanceof NthRoot) {
        return a[i];
      }
    }
    for (var i = 0; i < b.length; i += 1) {
      if (b[i] instanceof NthRoot) {
        return b[i];
      }
    }
    */
    for (var i = 0; i < a.length; i += 1) {
      if (a[i] instanceof Addition) {
        //return variables(a[i])[0];//TODO: fix
        a = a.concat(variables(a[i]));
        a[i] = null;
      }
    }
    for (var i = 0; i < b.length; i += 1) {
      if (b[i] instanceof Addition) {
        //return variables(b[i])[0];//TODO: fix
        b = b.concat(variables(b[i]));
        b[i] = null;
      }
    }
    for (var i = 0; i < a.length; i += 1) {
      for (var j = 0; j < b.length; j += 1) {
        if (a[i] != null && b[j] != null) {
        if (a[i].equals(b[j])) {
          return a[i];
        }
        }
      }
    }
    return null;
  };

  // TODO: fix or remove ?
  Expression.prototype.gcd = function (x) {
    //TODO: fix
    //return gcd(this, x, getVariable(this) || getVariable(x));
    //TODO: remove this block (a workaround for buggy gcd)
    var t1 = getTerm(this);
    var t2 = getTerm(x);
    if (t1 != null && t2 != null && t1.equals(t2)) {
      return getConstant(this).gcd(getConstant(x)).multiply(t2);
    }
    //!2020-11-05 more workarounds:
    t1 = getTerm(getFirstAdditionOperand(this));
    t2 = getTerm(getFirstAdditionOperand(x));
    if (t1 != null && t2 != null && t1.equals(t2)) {
      var c1 = getConstant(getFirstAdditionOperand(this));
      var c2 = getConstant(getFirstAdditionOperand(x));
      if (c1 instanceof Integer && c2 instanceof Integer) {
        var alpha = c1.truncatingDivide(c2);
        if (alpha instanceof Expression.Integer && alpha.multiply(c2).equals(c1)) {
          return this.subtract(x.multiply(alpha)).gcd(x);
        }
        var alpha = c2.truncatingDivide(c1);
        if (alpha instanceof Expression.Integer && alpha.multiply(c1).equals(c2)) {
          return this.gcd(x.subtract(this.multiply(alpha)));
        }
      }
    }
    if (this.equals(Expression.ZERO) || x.equals(Expression.ZERO)) {
      return this.add(x);
    }
    //!
    var result = gcd(this, x, getCommonVariable(this, x));
    return result;
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

  var Symbol = null;

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
    if (this.symbol === '\u2147') {
      return 'e';//!
    }
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

  var simplifyIdentityMatrixPower = true; //! TODO:

  function Integer(value) {
    //Expression.call(this);
    this.value = value;
  }

  Integer.prototype = Object.create(Expression.prototype);

  Integer.prototype.powExpression = function (x) {
    var y = this;
    if (x instanceof IdentityMatrix) {
      if (simplifyIdentityMatrixPower) {
        return new IdentityMatrix(x.symbol);
      }
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
    //!new 2019-12-16
    if (x instanceof Exponentiation && getExponent(x) instanceof Integer && y instanceof Integer) {//? (X**2)**(-1)
      return getBase(x).pow(getExponent(x).multiply(y));
    }
    //!
    //!new 2020-03-02
    if (x instanceof Exponentiation && getTerm(getExponent(x)) instanceof Expression.Symbol) {//? (X**2)**(-1)
      return getBase(x).pow(getExponent(x).multiply(y));
    }
    //!
    if (y.compareTo(Expression.ZERO) < 0) {
      return Expression.ONE.divide(x.pow(y.negate()));
    }
    if (x instanceof Expression.Matrix) {
      if (y.compareTo(Expression.ONE) > 0) {
        if (!x.matrix.isDiagonal()) {
          if (Expression.callback != undefined) {
            Expression.callback(new Expression.Event("pow", x, new Expression.Matrix(Matrix.I(1).map(function () { return y; }))));
          }
        }
      }
      var powMatrix = function (matrix, n) {
        if (n.toNumber() > 9007199254740991) {
          return powMatrix(matrix, n.truncatingDivide(Expression.TWO)).pow(2).multiply(matrix.pow(n.remainder(Expression.TWO).toNumber()));
        }
        return matrix.pow(n.toNumber());
      };
      return new Expression.Matrix(powMatrix(x.matrix, y));
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
      var t = x.b.multiply(y);
      if (t.getNumerator() instanceof Integer && t.getDenominator() instanceof Integer) {//TODO: ?
        var i = t.getNumerator().truncatingDivide(t.getDenominator());
        if (i.compareTo(Expression.ZERO) > 0) {
          return x.a.pow(i).multiply(x.a.pow(t.subtract(i)));
        }
      }
      return x.a.pow(x.b.multiply(y));
    }
    if (x instanceof Integer && (x.compareTo(Expression.ZERO) === 0 || x.compareTo(Expression.ONE) === 0 || x.compareTo(Expression.ONE.negate()) === 0)) {
      return y.remainder(Expression.TWO).compareTo(Expression.ZERO) === 0 ? x.multiply(x) : x;
    }
    if (x.equals(Expression.I)) {
      y = y.remainder(Expression.TWO.add(Expression.TWO));
      return Expression.pow(x, y.toNumber());
    }
    // assert(x instanceof Operation || x instanceof Integer);
    return Expression.pow(x, y.toNumber());
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
  Integer.prototype.isUnit = function () {
    return this.equals(Expression.ONE) || this.equals(Expression.ONE.negate());
  };
  Integer.prototype.toExpression = function () {
    return this;
  };
  Integer.prototype.compareTo = function (y) {
    return BigInteger.lessThan(this.value, y.value) ? -1 : (BigInteger.lessThan(y.value, this.value) ? +1 : 0);
  };
  Integer.prototype.abs = function () {
    return BigInteger.lessThan(this.value, BigInteger.BigInt(0)) ? this.negate() : this;
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
    if (gInteger.compareTo(Expression.ZERO) < 0) {
      gInteger = gInteger.negate();
    }
    if (y.compareTo(Expression.ZERO) < 0) {
      gInteger = gInteger.negate();
    }
    x = x.truncatingDivide(gInteger);
    y = y.truncatingDivide(gInteger);
    return y.compareTo(Expression.ONE) === 0 ? x : new Division(x, y);
  };
  Integer.prototype.truncatingDivide = function (y) {
    return y.truncatingDivideInteger(this);
  };
  Integer.prototype.truncatingDivideInteger = function (x) {
    var y = this;
    return new Integer(BigInteger.divide(x.value, y.value));
  };
  Integer.prototype.isDivisibleBy = function (y) {
    return y.isDivisibleByInteger(this);
  };
  Integer.prototype.isDivisibleByInteger = function (x) {
    return x.remainder(this).equals(Expression.ZERO);
  };
  Integer.prototype.remainder = function (y) {
    return y.remainderInteger(this);
  };
  Integer.prototype.remainderInteger = function (x) {
    var y = this;
    return new Integer(BigInteger.remainder(x.value, y.value));
  };
  Integer.prototype.primeFactor = function () {
    return integerPrimeFactor(this);
  };
  Integer.prototype.toNumber = function () {
    return BigInteger.toNumber(this.value);
  };
  Integer.prototype.toBigInt = function () {
    return this.value;
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
  Integer.fromBigInt = function (i) {
    return new Integer(BigInteger.BigInt(i.toString()));
  };

  Expression.ZERO = Integer.fromNumber(0);
  Expression.ONE = Integer.fromNumber(1);
  Expression.TWO = Integer.fromNumber(2);
  Expression.TEN = Integer.fromNumber(10);



  Expression.Matrix = function (matrix) {
    //Expression.call(this);
    this.matrix = matrix;
  };

  Expression.Matrix.fromArray = function (rows) {
    return new Expression.Matrix(Matrix.padRows(rows, null));
  };

  Expression.Matrix.prototype = Object.create(Expression.prototype);

  Expression.Matrix.prototype.equals = function (x) {
    //!new 2019-12-03
    if (x === Expression.ZERO) {
      return this.matrix.isSquare() && this.matrix.eql(this.matrix.map(function (e, i, j) {
        return Expression.ZERO;
      }));
    }
    //!
    if (!(x instanceof Expression.Matrix)) {
      return false;
    }
    return this.matrix.eql(x.matrix);
  };

  Expression.Matrix.prototype.compare4Multiplication = function (y) {
    if (y instanceof Expression.Matrix) {
      return 0;
    }
    if (y instanceof MatrixSymbol) {
      if (this.matrix.isSquare() && this.matrix.isDiagonal()) {//?
        return +1;
      }
      return -1;
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
    if (Expression.callback != undefined) {
      Expression.callback(new Expression.Event("add", x, this));
    }
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
    //if (this instanceof Division && this.isNegative()) {
    //  return '-' + this.negateCarefully().toString(options);
    //}
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

  //TODO: remove - ?
  Exponentiation.prototype.compare4Multiplication = function (y) {
    return y.compare4MultiplicationExponentiation(this);
  };
  Exponentiation.prototype.compare4MultiplicationInteger = function (x) {
    return -1;
  };
  Exponentiation.prototype.compare4MultiplicationExponentiation = function (x) {
    var y = this;
    return getBase(x).compare4Multiplication(getBase(y)) || getExponent(x).compare4Multiplication(getExponent(y));
  };

  function Multiplication(a, b) {
    BinaryOperation.call(this, a, b);
  }

  Multiplication.prototype = Object.create(BinaryOperation.prototype);

  Multiplication.prototype.multiply = function (y) {
    return y.multiplyExpression(this);
  };
  //TODO:
  var compare4Multiplication2 = function (x, y) {//TODO: fix

//  && x.n !== y.n
    if (x instanceof NthRoot && y instanceof NthRoot) {//TODO: fix
      if (x.multiply(y).equals(new Expression.Multiplication(x, y))) {
        return -1;
      }
      if (x.multiply(y).equals(new Expression.Multiplication(y, x))) {
        return +1;
      }
      return 0;
    }
/*
    //!2019-04-22
    if (x instanceof NthRoot && y instanceof NthRoot && x.n === y.n) {//TODO: fix
      if (x.a instanceof Integer && y.a instanceof Integer) {
        return 0;
      }
      if (x.a instanceof Addition && y.a instanceof Integer) {
        return 0;//TODO: fix
      }
      if (x.a instanceof Integer && y.a instanceof Addition) {
        return 0;//TODO: fix
      }
      // -(2^0.5+1)^0.5*(2*2^0.5+2)^0.5
      if (x.a instanceof Addition && y.a instanceof Addition) {
        return 0;//TODO: fix
      }
      // 3 and 3^n
      return compare4Multiplication(x.a, y.a);
    }
  */
    if (x instanceof Integer && y instanceof Exponentiation) {
      return -1;//?
    }
    if (x instanceof Exponentiation && y instanceof Integer) {
      return +1;//?
    }
    if (x instanceof Expression.Complex && y instanceof Exponentiation) {
      return -1;//?
    }
    if (x instanceof Exponentiation && y instanceof Expression.Complex) {
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
  Division.prototype.addDivision = function (x) {
    if (x.b.equals(this.b)) {
      return x.a.add(this.a).divide(this.b);
    }
    return BinaryOperation.prototype.addDivision.call(this, x);
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
  Addition.prototype.compare4MultiplicationMatrixSymbol = function (x) { // (X+{{1}})*X
    return -1;
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
  MatrixSymbol.prototype.transpose = function () {
    // quick solution:
    return new Expression.Exponentiation(this, new Expression.Symbol("T")); // TODO: fix
  };
  //...

  Expression.MatrixSymbol = MatrixSymbol;

  function IdentityMatrix(symbol) {
    MatrixSymbol.call(this, symbol);
  }
  IdentityMatrix.prototype = Object.create(MatrixSymbol.prototype);
  //IdentityMatrix.prototype.multiply = function (y) {
  //  return y.multiplyIdentityMatrix(this);
  //};

  //TODO: move to MatrixSymbol - ?
  IdentityMatrix.prototype.multiplyAddition = function (x) {
    if (isScalar(x)) {
      return new Multiplication(x, this);
    }
    return Expression.prototype.multiplyAddition.call(this, x);
  };

  //Expression.prototype.multiplyIdentityMatrix = function (x) {
  //  return this.multiplyExpression(x);
  //};
  //IdentityMatrix.prototype.multiplyIdentityMatrix = function (x) {
  //  return new IdentityMatrix(this.symbol);
  //};
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

  IdentityMatrix.prototype.compare4MultiplicationMatrixSymbol = function (x) {
    var y = this;
    return x instanceof IdentityMatrix ? (x.symbol < y.symbol ? -1 : (y.symbol < x.symbol ? +1 : 0)) : +1;
  };

  Expression.IdentityMatrix = IdentityMatrix;

  BinaryOperation.prototype.getS = function () {
    throw new TypeError("abstract");
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
    if (x instanceof Expression.Radians) {
      return x.value.isNegative();
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
    var fa = this.a.getPrecedence() <= this.getPrecedence();
    return (fa ? "(" : "") + this.a.toString(setTopLevel(fa || options == undefined || options.isTopLevel, options)) + (fa ? ")" : "") + "^" + (this.n === 2 ? "0.5" : "(1/" + this.n + ")");
  };

  NthRoot.prototype.getDegree = function () {
    return this.n;
  };

  function ngcd(a, b) {
    while (b != 0) {
      var t = a % b;
      a = b;
      b = t;
    }
    return a;
  }
  function isPrime(n) {
    if (typeof n === "bigint") {//TODO: ?
      return n === BigInt(primeFactor(BigInteger.BigInt(n.toString())).toString());
    }
    return n === primeFactor(n);
  }

  function sortKeys(object) {
    var previous = 0;
    var ok = true;
    for (var key in object) {
      if (Object.prototype.hasOwnProperty.call(object, key)) {
        var i = Number(key);
        if (!(i > previous)) {
          ok = false;
        }
        previous = i;
      }
    }
    if (ok) {
      return object;
    }
    var keys = [];
    for (var key in object) {
      if (Object.prototype.hasOwnProperty.call(object, key)) {
        keys.push(Number(key));
      }
    }
    keys.sort(function (a, b) {
      return a - b;
    });
    var result = {};
    for (var i = 0; i < keys.length; i++) {
      result[keys[i]] = object[keys[i]];
    }
    return result;
  }

  function isPerfectCube(n) {
    return BigInteger.equal(BigInteger.exponentiate(nthRoot(n, 3), BigInteger.BigInt(3)), n);
  }

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
          var t = tt instanceof Integer && tt.compareTo(Expression.ZERO) >= 0 ? tt.squareRoot() : undefined;
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
    if (n === 3) {//? new: 2019-08-18
      if (x instanceof Addition) {
        if ((x.a instanceof SquareRoot || x.a instanceof Multiplication && x.a.a instanceof Integer && x.a.b instanceof SquareRoot) && x.b instanceof Integer) {
          // (5^0.5+3)/2 or (-5^0.5+3)/2
          var u = x.a;
          var v = x.b;
          var d = u.multiply(u).subtract(v.multiply(v));
          if (isPerfectCube(d.toBigInt())) {//?
            // (a+b)^3 = aaa+3aab+3bba+bbb = u+v
            // aaa+3bba = v, bb=(v-aaa)/(3a)
            // 3aab+bbb = u, b(3aa+bb)=u, b=u/(3aa+bb), bb=u**2/(3aa+bb)**2
            // (9aaa+(v-aaa))**2*(v-aaa) = 27*aaa*u**2

            // t = aaa
            // (9t+(v-t))**2*(v-t) = 27*t*u**2

            var t = new Expression.Symbol('t');
            //var eq = ExpressionParser.parse('(9t+(v-t))**2*(v-t) - 27*t*u**2', new RPN.Context(function (id) {
            //  return id === 'v' ? v : (id === 'u' ? u : (id === 't' ? t : undefined));
            //})).simplify();
            var NINE = Expression.Integer.fromNumber(9);
            var TWENTY_SEVEN = Expression.Integer.fromNumber(27);
            var eq = Expression.pow(NINE.multiply(t).add(v.subtract(t)), 2).multiply(v.subtract(t)).subtract(TWENTY_SEVEN.multiply(t).multiply(Expression.pow(u, 2)));
            var p = Polynomial.toPolynomial(eq, t);
            p = p.scale(p.getContent().inverse());
            var t = p.doRationalRootTest();
            if (t != null) {
              var a = t._nthRoot(3);
              var b = v.subtract(t).divide(a).divide(Expression.Integer.fromNumber(3))._nthRoot(2);
              return a.add(b);
            }
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
    if (x instanceof Division || x instanceof Multiplication) {
      if (n % 2 !== 0 || x.a instanceof Integer && x.a.compareTo(Expression.ZERO) > 0 || // sqrt(-x) = sqrt(-1 * -1) = i * i = -1
                         x.b instanceof Integer && x.b.compareTo(Expression.ZERO) > 0 || //TODO: fix
                         x.a instanceof Integer && isConstant(x.b) ||
                         isConstant(x) ||
                         isPositive(x.a) ||
                         isPositive(x.b)) { // -sqrt(7)
        if (x instanceof Division) {
          return x.a._nthRoot(n).divide(x.b._nthRoot(n));
        }
        if (x instanceof Multiplication) {
          return x.a._nthRoot(n).multiply(x.b._nthRoot(n));
        }
      }
    }
    var qi = x instanceof Integer ? x : null;
    var qq = QuadraticInteger.toQuadraticInteger(x);
    // sqrt(2sqrt(2)+2) "(2*2^0.5+2)^0.5"
    // sqrt(sqrt(2)+2) "(2^0.5+1)^0.5*2^(1/4)"
    // sqrt(4sqrt(2)+4) 2*(2^0.5+1)^0.5
    // sqrt(2sqrt(2)+4) (2*2^0.5+2)^0.5*2^(1/4)

    //TODO: fix for !isPrime(qq.D)
    //if (qq != null && (isPrime(qq.D) && (qq.a / ngcd(qq.a, qq.b)) % qq.D == 0 || !isPrime(qq.D) && qq.a % qq.D == 0 && ngcd(qq.a, qq.b) % qq.D != 0)) {
      //var D = Expression.Integer.fromNumber(qq.D)._nthRoot(2);
      //return D._nthRoot(n).multiply(x.divide(D)._nthRoot(n));
    //}

    //!
    if (qq == null && x instanceof Expression.Addition) {
      if (x.a instanceof Multiplication && x.b instanceof Multiplication) {
        var g = x.a.pow(Expression.TWO).gcd(x.b.pow(Expression.TWO)).squareRoot();
        if (!g.equals(Expression.ONE)) {
          return g._nthRoot(n).multiply(x.divide(g)._nthRoot(n));
        }
      }
    }
    //!

    if (qq != null) {
      if (n !== 2 && n % 2 === 0) {
        //TODO: check
        return x.squareRoot()._nthRoot(n / 2);
      }
      if ((n === 2 || n === 3) && qq.isValid()) {//Note: isValid should be checked here
      //if (qq.norm() === -1 * Math.pow(ngcd(qq.a, qq.b), 2)) {
        if (qq.isPositive()) {
        qi = qq;
        } else {
          return Expression.ONE.negate()._nthRoot(n).multiply(this.negate()._nthRoot(n));
        }
      //}
      }
    }
    if (qi != null) {
      x = qi;//TODO:
      if (x instanceof Integer && x.compareTo(Expression.ZERO) < 0) {
        if (n % 2 === 0) {
          if (n === 2) {//TODO: ?
            return Expression.I.multiply(this.negate()._nthRoot(n));
          }
          throw new RangeError("NotSupportedError");
        }
        return this.negate()._nthRoot(n).negate();
      }
      if (x.equals(Expression.ZERO)) {
        return this;
      }
      var makeRoot = function (i, n) {
        return n === 1 ? i : (n === 2 ? new SquareRoot(i) : (n === 3 ? new CubeRoot(i) : new NthRoot("n" + "-root", i, n)));
      };
      var roots = {};
      var i = x;
      while (!i.equals(Expression.ONE)) {
        var d = i.primeFactor();
        if (d instanceof QuadraticInteger) {
          if (n !== 2 && n !== 3) {
            throw new TypeError(); // "assertion"
          }
          if (BigInteger.multiply(d.a, d.b) < 0) {
            var a = (n % 2 === 0 ? d.abs() : d);
            return x.toExpression().divide(a.toExpression()._pow(n))._nthRoot(n).multiply(a.toExpression());
          }
          var s = d.norm();
          // https://brownmath.com/alge/nestrad.htm#SurveyDoable
          var isPerfectSquare = function (n) {
            return BigInteger.equal(BigInteger.exponentiate(nthRoot(n, 2), BigInteger.BigInt(2)), n);
          };
          //TODO: s >= 0 - ?
          if (d.b != 0 && d.a != 0 && s >= 0 && isPerfectSquare(s)) {
            if (n === 2) {
            return x.toExpression().divide(d.toExpression())._nthRoot(2).multiply(d.toExpression()._nthRoot(2));
            }
          }
        }
        var e = 0;
        if (i.isUnit()) {
          //TODO:
          // d should be a https://en.wikipedia.org/wiki/Fundamental_unit_(number_theory)
          while (!i.equals(Expression.ONE)) {
            i = i.truncatingDivide(d);
            e += 1;
          }
        } else {
          while (i.isDivisibleBy(d)) {
            i = i.truncatingDivide(d);
            e += 1;
          }
        }
        d = d.toExpression();
        var nn = n;
        if (d instanceof NthRoot) {
          nn *= d.n;
          d = d.a;
        }
        var t = ngcd(nn, e);
        nn /= t;//?
        e /= t;//?

        while (e !== 0) {
          //var g = e;
          //while (nn % g !== 0) {
          //  g -= 1;
          //}
          var g = e >= nn ? nn : 1;

          var e1 = Math.floor(e / g);
          var k = Math.floor(nn / g);
          roots[k] = (roots[k] || Expression.ONE).multiply(Expression.pow(d, e1));

          e = e - g * e1; // e = e % g;
        }
      }
      var y = Expression.ONE;
      roots = sortKeys(roots);
      //for (var j = 1; j <= n; j += 1) {
      //}
      // `for-in` loop is used to have better performance for "a sparse array"
      var f = roots[1];
      for (var jj in roots) {//TODO: fix the iteration order
        if (Object.prototype.hasOwnProperty.call(roots, jj)) {
          var j = Number(jj);
        if (j !== 1) {
          var r = roots[j];
          //y = y.multiply(makeRoot(r, j));
          var x = makeRoot(r, j);
          if (y !== Expression.ONE) {
            y = new Expression.Multiplication(y, x);
          } else {
            y = x;
          }
        }
        }
      }
      if (f != null) {
        y = f.multiply(y);
      }
      return y;
    }
    if (x instanceof Expression.Matrix) {
      if (typeof hit === "function") {
        hit(n === 2 ? {squareRoot: "matrix"} : (n === 3 ? {cubeRoot: "matrix"} : {nthRoot: "Matrix^(1/" + n + ")"}));
      }
      var tmp = Expression.getEigenvalues(x.matrix);
      var eigenvalues = tmp.eigenvalues;
      var multiplicities = tmp.multiplicities;
      var N = typeof n === "number" ? Expression.Integer.fromNumber(n) : n;
      if (Expression.sum(multiplicities) === x.matrix.cols()) {
        var tmp2 = Expression.getEigenvectors(x.matrix, eigenvalues);
        var eigenvectors = tmp2.eigenvectors;
        if (eigenvectors.length === x.matrix.cols()) {
          if (!x.matrix.isDiagonal()) {
            if (Expression.callback != undefined) {
              Expression.callback(new Expression.Event("nth-root-using-diagonalization", x));
            }
            if (Expression.callback != undefined) {//TODO: remove - ?
              Expression.callback(new Expression.Event("diagonalize", x));
            }
          }
          var tmp = Expression.diagonalize(x.matrix, eigenvalues, multiplicities, eigenvectors);
          var L = tmp.L;
          var SL = L.map(function (e, i, j) {
            return i === j ? e.pow(Expression.ONE.divide(N)) : e;
          });
          return new Expression.Matrix(tmp.T.multiply(SL).multiply(tmp.T_INVERSED));
        } else {
          if (!x.matrix.isJordanMatrix()) {
            if (Expression.callback != undefined) {
              Expression.callback(new Expression.Event("nth-root-using-Jordan-normal-form", x));
            }
            if (Expression.callback != undefined) {//TODO: remove - ?
              Expression.callback(new Expression.Event("Jordan-decomposition", x));
            }
          } else {
            if (Expression.callback != undefined) {
              Expression.callback(new Expression.Event("Jordan-matrix-nth-root", x));
            }
          }
          var rootOfJordanForm = function (J, N) {
            var tmp = J.map(function (e, i, j) {
              if (i > j) {
                return Expression.ZERO;
              }
              if (i === j) {
                return J.e(i, j).pow(Expression.ONE.divide(N));
              }
              if (J.e(i, i + 1).equals(Expression.ZERO)) {
                return Expression.ZERO;
              }
              if (!J.e(i, i + 1).equals(Expression.ONE)) {
                throw new TypeError("assertion");
              }
              //if (i + 1 === j) {
                //return J.e(i, i).pow(Expression.ONE.divide(N)).divide(N.multiply(J.e(i, i)));
                //return J.e(i, i + 1).divide(N.multiply(J.e(i, i).pow(Expression.ONE.divide(N)).pow(N.subtract(Expression.ONE))));
              //}
              //return new Expression.Symbol('aa_(' + (j - i) + ',' + j + ')');
              var m = j - i;
              for (var k = 0; k < m; k += 1) {
                if (!J.e(j - 1 - k, j - k).equals(Expression.ONE)) { // outside of a block
                  return Expression.ZERO;
                }
              }
              // 1/n(1/n-1)(1/n-2)(1/n-3)/(4!*lambda**4) * lambda**(1/n)
              var f = Expression.ONE;
              for (var k = 0; k < m; k += 1) {
                f = f.multiply(Expression.ONE.divide(N).subtract(Expression.Integer.fromNumber(k))).divide(Expression.Integer.fromNumber(k + 1));
              }
              return f.divide(J.e(i, i)._pow(m)).multiply(J.e(i, i).pow(Expression.ONE.divide(N)));
            });

            /*
            for (var k = 2; k < J.cols(); k += 1) {
              //var x = tmp.pow(N);
              var x = new Expression.Matrix(tmp).pow(N).matrix;//!?
              tmp = tmp.map(function (e, i, j) {
                if (i + k === j) {
                  if (x.e(i, j).equals(Expression.ZERO)) {
                    return Expression.ZERO;
                  }
                  var s = new Expression.Symbol('aa_(' + (j - i) + ',' + j + ')');
                  var p = Polynomial.toPolynomial(x.e(i, j).getNumerator(), s);
                  if (p.getDegree() === 0) {
                    return x.e(i, j);
                  }
                  if (p.getDegree() !== 1) {
                    throw new TypeError("!");
                  }
                  var y = p.getCoefficient(0).negate().divide(p.getCoefficient(1));
                  return y;
                }
                return e;
              });
            }
            */
            return tmp;
          };
          var tmp = Expression.getFormaDeJordan(x.matrix, eigenvalues, multiplicities);
          var JN = rootOfJordanForm(tmp.J, N);
          //TODO: details - ?
          return new Expression.Matrix(tmp.P.multiply(JN).multiply(tmp.P_INVERSED));
        }
      }
      //TODO: using Jordan normal form -?
    }
    //!2019-04-22
    if (x instanceof Exponentiation && x.a instanceof Integer && x.a.compareTo(Expression.ZERO) > 0) {
      //if (n === 2) {//TODO:
        if (x.b instanceof Expression.Symbol) {
          if (x.a instanceof Expression.Integer && integerPrimeFactor(x.a).equals(x.a)) {
            //return new SquareRoot(x);
            return new Expression.Exponentiation(x.a, x.b.divide(Expression.Integer.fromNumber(n)));
          }
        } else {
          return x.a.pow(x.b.divide(Expression.Integer.fromNumber(n)));
        }
      //}
    }

    //!2019-16-06
    if (x instanceof Exponentiation && x.a === Expression.E) {
      return x.a.pow(x.b.divide(Expression.Integer.fromNumber(n)));
    }
    if (x instanceof IdentityMatrix) {
      if (simplifyIdentityMatrixPower) {
        return x;
      }
    }
    //!2019-17-06
    if (x instanceof Expression.Symbol) {
      return new Expression.Exponentiation(x, Expression.ONE.divide(Expression.Integer.fromNumber(n)));
    }
    if (x instanceof Exponentiation && x.a instanceof Expression.Symbol && (n % 2 !== 0 || (x.b.getNumerator() instanceof Integer && !x.b.getNumerator().remainder(Expression.TWO).equals(Expression.ZERO)))) {
      //TODO: fix condition for n % 2 === 0
      var b = x.b.divide(Expression.Integer.fromNumber(n));
      return b.equals(Expression.ONE) ? x.a : new Expression.Exponentiation(x.a, b);
    }

    //!2019-06-20
    var v = getVariable(x);
    if (v instanceof Expression.Symbol && (n === 2 || n === 3)) {//TODO: other n's
      var p = Polynomial.toPolynomial(x, v);
      if (p.getDegree() === 1 && p.getCoefficient(0) instanceof Integer && !p.getCoefficient(0).equals(Expression.ZERO) && p.getCoefficient(1) instanceof Integer) {
        //TODO:
        var c = p.getContent();
        if (c.isNegative()) {
          c = c.negate();
        }
        if (!c.equals(Expression.ONE)) {
          return x.divide(c)._nthRoot(n).multiply(c._nthRoot(n));
        }
        if (p.getCoefficient(1).compareTo(Expression.ZERO) > 0) {
          return new Expression.Exponentiation(x, Expression.ONE.divide(Expression.Integer.fromNumber(n)));
        } else {
          //!TODO: fix
          if (n % 2 !== 0) {
            return Expression.ONE.negate()._nthRoot(n).multiply(new Expression.Exponentiation(x.negate(), Expression.ONE.divide(Expression.Integer.fromNumber(n))));
          }
        }
      }
      if (p.getDegree() > 1 && !p.getCoefficient(0).equals(Expression.ZERO)) {
        //TODO: check
        var N = p.getDegree();
        var t = v.multiply(p.getCoefficient(N)._nthRoot(N)).add(p.getCoefficient(0)._nthRoot(N));
        if (x.equals(t._pow(N))) {
          //!TODO: remove
          if (N >= n && n % 2 !== 0) {
            return t.pow(Expression.ONE).multiply(t._pow(N - n)._nthRoot(n));
          }
          //!
          if (n % 2 !== 0) {
            return new Expression.Exponentiation(t, Expression.Integer.fromNumber(N).divide(Expression.Integer.fromNumber(n)));
          }
        }
        //TODO: (ax+b)**(n+1)
        if (p.getDegree() > 1 && p.getSquareFreePolynomial().equals(p) && n === 2) {//TODO: fix ?
          return new Expression.Exponentiation(x, Expression.ONE.divide(Expression.Integer.fromNumber(n)));
        }
      }
    }

    if (n % 2 !== 0 && x instanceof Expression.ExponentiationOfMinusOne) {//?
      return getBase(x).pow(getExponent(x).divide(Expression.Integer.fromNumber(n)));
    }

    throw new RangeError("NotSupportedError");
  };

  function SquareRoot(a) {
    NthRoot.call(this, "sqrt", a, 2);
  }

  SquareRoot.prototype = Object.create(NthRoot.prototype);
  //!
  SquareRoot.prototype.divideInteger = function (x) {
    //TODO: check
    return x.multiply(this).divide(this.a);
  };

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
      throw new RangeError("NotSupportedError:matrixArgExpected");//?
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
      throw new RangeError("NotSupportedError:matrixArgExpected");//?
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
      throw new RangeError("NotSupportedError:matrixArgExpected");//?
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
    if (x instanceof Expression.Multiplication) {
      //TODO: info about properties of the Matrix Transpose
      return x.b.transpose().multiply(x.a.transpose());//TODO: ?
    }
    if (x instanceof Expression.Addition) {
      return x.a.transpose().add(x.b.transpose());
    }
    if (isScalar(x)) {
      return x;
    }
    if (!(getBase(x) instanceof MatrixSymbol) && x instanceof Expression.Exponentiation && x.b.equals(Expression.ONE.negate())) {
      //TODO: (X^-2)^T
      return x.a.transpose().pow(x.b);
    }
    if (x instanceof Expression.Exponentiation && x.b.equals(new Expression.Symbol("T"))) {
      //TODO: (X**2)^T
      return x.a;
    }
    if (getBase(x) instanceof MatrixSymbol) {
      var e = getExponent(x).multiply(new Expression.Symbol("T"));
      //TODO: ?
      var p = Polynomial.toPolynomial(e, new Expression.Symbol("T"));
      if (p.getDegree() >= 2) {
        e = e.subtract(p.getCoefficient(2).multiply(new Expression.Symbol("T")._pow(2))).add(p.getCoefficient(2));
      }
      return new Expression.Exponentiation(getBase(x), e);
    }
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
      throw new RangeError("NotSupportedError:matrixArgExpected");//?
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
  Expression.NoAnswerExpression.prototype.add = function () {
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
      //TODO: (!)
      var f = e.primeFactor();
      if (!f.equals(e) && e.divide(f) instanceof Expression.Integer) {
        f = f.multiply(Expression.I);
      }
      return f;
      /*
      var g = integerGCD(e.real, e.imaginary);
      var t = simpleDivisor(g);
      if (t != null) {
        return t;
      }
      if (typeof hit === "function") {
        hit({everySimpleDivisor: e.toString()});
      }
      return e;
      */
    }
    //var v = getVariable(e);
    // To avoid square roots / nth roots:
    var v = getVariableInternal(getLastMultiplicationOperand(getFirstAdditionOperand(e))).next().value.v;
    //if (v instanceof NthRoot || v instanceof Integer || v instanceof Expression.Complex) {
    //  v = undefined;
    //}
    if (v != undefined) {
      var r = getReplacement(e, v);
      if (!r.equals(v)) {
        return substitute(simpleDivisor(substitute(e, v, r, inverseReplacement(r, v))), v, inverseReplacement(r, v), r);
      }

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
        var root = np.doRationalRootTest();
        if (root == null) {
          //TODO: TEST COVERAGE (!)
          var t = np._getFactorByKroneckersMethod();
          if (t != null) {
            //TODO: write a test case
            return simpleDivisor(t.calcAt(v));
          }
        }
        if (root != null) {
          var t = v.multiply(root.getDenominator()).subtract(root.getNumerator());
          return t;
        }
      }

      /*
      if (np.getDegree() >= 2) {
        var roots = np.getroots();
        if (roots.length > 0) {
          var root = roots[0];
          return v.subtract(root);
        }
      }
      */

      var sf = np.getSquareFreePolynomial();
      if (!sf.equals(np)) {
        return sf.calcAt(v);
      }

      e = np.calcAt(v);
      if (e.isNegative()) {//TODO: remove - ?
        e = e.negate();//!?
      }
      return e;
    }
    throw new RangeError();//?
  };
  Expression.simpleDivisor = simpleDivisor;

  Expression.everyDivisor = function (e, callback) {
    if (!callback(Expression.ONE)) {
      return false;
    }
    var divisors = [];
    var rec = function (start, n, s) {
      if (n >= 0) {
        var x = divisors[n];
        for (var i = start; i <= x.e; i += 1) {
          if (!rec(0, n - 1, s.multiply(Expression.pow(x.d, i)))) {
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
    while (!e.equals(Expression.ONE) && !e.equals(Expression.ONE.negate())) {
      var d = simpleDivisor(e);
      if (divisors.length === 0 || !divisors[divisors.length - 1].d.equals(d)) {
        divisors.push({
          d: d,
          e: 0
        });
      }
      divisors[divisors.length - 1].e += 1;
      if (!rec(divisors[divisors.length - 1].e, divisors.length - 1, Expression.ONE)) {
        return false;
      }
      e = e.divide(d);
    }
    return true;
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
  //TODO: remove
  Expression.pow = function (x, count) {
    return x._pow(count);
  };
  Expression.prototype._pow = function (count) {
    return pow(this, count, Expression.ONE);
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
      throw new TypeError();
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
      throw new TypeError();
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
    var xx = undefined;
    for (var x of e.summands()) {
      //TODO: why iteration by additions - (?)
      var state = 0;
      for (var y of x.factors()) {
        var factor = y;
        var factorBase = getBase(y);
        var factorExponent = getExponent(y);
        if ((!(factorBase instanceof Integer) && !(factorBase instanceof Expression.Symbol) && !(factorBase instanceof Expression.Matrix)) ||
            !(factorExponent instanceof Integer)) {//TODO: fix
          throw new RangeError("NotSupportedError");
        }
        var s = factorBase instanceof Expression.Symbol ? factorBase.toString() : "";
        if (s === "X" && state === 0) {
          state = 1;
          xx = factor;
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
    return {s: scalar, l: l, r: r, x: xx};
  };
  Expression.splitX = splitX;
  var groupX = function (a, b) {
    var tmp1 = splitX(a);
    var tmp2 = splitX(b);
    var s1 = tmp1.s;
    var l1 = tmp1.l;
    var r1 = tmp1.r;
    var s2 = tmp2.s;
    var l2 = tmp2.l;
    var r2 = tmp2.r;
    if (r1 == undefined && r2 == undefined && tmp1.x.equals(tmp2.x)) {
      l1 = l1 == undefined ? new IdentityMatrix("I") : l1;
      l2 = l2 == undefined ? new IdentityMatrix("I") : l2;
      return new Multiplication(s1.multiply(l1).add(s2.multiply(l2)), tmp1.x);
    }
    if (l1 == undefined && l2 == undefined && tmp1.x.equals(tmp2.x)) {
      r1 = r1 == undefined ? new IdentityMatrix("I") : r1;
      r2 = r2 == undefined ? new IdentityMatrix("I") : r2;
      return new Multiplication(tmp1.x, s1.multiply(r1).add(s2.multiply(r2)));
    }
    return undefined;
  };

  //TODO: remove (replace with a Condition) - ?
  //?
  var getExpressionWithX = function (e) {
    if (e instanceof Division) {
      if (e.getDenominator() instanceof Expression.Integer) {
        e = e.getNumerator();//!
      } else {
        return {withX: undefined, withoutX: undefined};
      }
    }

    var withX = undefined;
    var withoutX = undefined;
    for (var x of e.summands()) {
      var summand = x;
      var hasX = false;
      for (var y of x.factors()) {
        var factor = y;
        var factorBase = getBase(factor);
        //!new2019-12-03
        if (factorBase instanceof Addition && getExponent(factor).isNegative()) { // (<matrix addition>)**-1 (?)
          //TODO: fix
          var q = null; // e.subtract(x)
          for (var x1 of e.summands()) {
            if (x1 !== x) {
              q = q == null ? x1 : q.add(x1);
            }
          }
          var exponent = getExponent(factor).negate();
          var e1 = q.multiply(factorBase.pow(exponent)).add(x.multiply(new Expression.Exponentiation(factorBase, exponent)));
          var z2 = e1.transformEquality(Expression.ZERO);
          //TODO:
          //var tmp = Polynomial.toPolynomial(factorBase, z2.a).calcAt(z2.b);
          var tmp = Polynomial.toPolynomial(factorBase, getBase(z2.a)).divideAndRemainder(Polynomial.toPolynomial(z2.a.subtract(z2.b), getBase(z2.a))).remainder.calcAt(getBase(z2.a));
          var d = tmp instanceof Expression.Matrix ? tmp.determinant() : null;
          if (!Expression.isConstant(d)) {//TODO: ?
            return {withX: undefined, withoutX: undefined};
          }
          if (d.equals(Expression.ZERO)) {
            return {withX: Expression.ZERO, withoutX: Expression.ONE};//TODO: no solutions
          }
          return getExpressionWithX(e1);//!hack
        }
        //!
        if (!(factorBase instanceof Integer) && !(factorBase instanceof Expression.Symbol)) {
          if (!(factorBase instanceof Expression.Matrix)) {//?
          if (!(factorBase instanceof NthRoot)) {//?TODO: remove - ?
          if (Expression.has(factorBase, Expression.MatrixSymbol)) {//?
            throw new RangeError("NotSupportedError");
          }
          }
          }
        }
        if (factorBase instanceof Expression.Symbol) {
          var s = factorBase.toString();
          if (s === "X") {
            if (hasX) {
              //throw new RangeError("NotSupportedError");
              return {withX: null, withoutX: null};
            }
            hasX = true;
          }
        }
      }
      if (hasX) {
        if (withX != undefined) {
          withX = groupX(withX, summand);
          if (withX == null) {
            //throw new RangeError("NotSupportedError");
            return {withX: null, withoutX: null};
          }
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
    } else if (e instanceof Expression.Division) {
      return isConstant(e.a) && isConstant(e.b);
    } else if (e instanceof Expression.Sin || e instanceof Expression.Cos) {
      return isConstant(e.a);
    } else if (e instanceof Expression.Radians) {
      return isConstant(e.value);
    //TODO:
    //} else if (e instanceof Expression.Logarithm) {
    //  return isConstant(e.a);//TODO: test
    }
    return false;
  };

  Expression.isConstant = isConstant;

  Expression.getMultivariatePolynomial = function (e) {
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
      //throw new TypeError("undefined");
      return undefined;
    }
    //?
    if (v instanceof Expression.Addition) {
      v = getVariableInternal(getLastMultiplicationOperand(getFirstAdditionOperand(v))).next().value.v;
    }
    //?

    //TODO:
    var r = getReplacement(e, v);
    if (!r.equals(v)) {
      e = substitute(e, v, r, inverseReplacement(r, v));
      if (e instanceof Expression.Division && e.b instanceof Expression.Integer) {
        e = e.a;//!
      }
    }

    var p = Polynomial.toPolynomial(e, v);
    //TODO: iteration by sparse coefficients
    for (var i = 0; i <= p.getDegree(); i += 1) {
      var c = p.getCoefficient(i);
      if (!isConstant(c)) {
        var pc = Expression.getMultivariatePolynomial(c);
        if (pc == undefined) {
          return undefined;
        }
      }
    }
    return {p: p.map(function (c) { return substitute(c, v, inverseReplacement(r, v), r); }), v: inverseReplacement(r, v)};
  };
  Expression.isSingleVariablePolynomial = function (e) {
    var tmp = Expression.getMultivariatePolynomial(e);
    if (tmp == null) {
      return false;
    }
    var p = tmp.p;
    //TODO: iteration by sparse coefficients
    for (var i = 0; i <= p.getDegree(); i += 1) {
      var c = p.getCoefficient(i);
      if (!isConstant(c)) {
        return false;
      }
    }
    return true;
  };

  // TODO: NotSupportedError
  Expression.prototype.transformEquality = function (b) {
    var e = this.subtract(b);
    var tmp = getExpressionWithX(e);
    var withX = tmp.withX;
    var withoutX = tmp.withoutX;
    if (withX == undefined) {
      if (e.getDenominator() instanceof Integer &&
          !(e.getNumerator() instanceof Expression.Matrix) &&
          !Expression.has(e, Expression.MatrixSymbol)) {
        //TODO: tests
        var tmp = Expression.getMultivariatePolynomial(e.getNumerator());
        if (tmp != undefined) {
          var p = tmp.p;
          var v = tmp.v;
          var m = Matrix.Zero(1, p.getDegree() + 1).map(function (e, i, j) {
            return p.getCoefficient(p.getDegree() - j);
          });
          return new Expression.NoAnswerExpression(new Expression.Matrix(m), "polynomial-roots", {polynomial: p, variable: v});
        }
      }
      if (e instanceof Expression.Matrix) {
        if (this instanceof Expression.Matrix && (b instanceof Expression.Matrix || b instanceof Expression.IdentityMatrix)) {
          return Expression.SystemOfEquations.from([{left: this, right: b}]);
        }
        //TODO: other things - ?
      }
      if (true) {//!new 2019-11-27
        //TODO: fix
        return Expression.SystemOfEquations.from([{left: this, right: b}]);
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
    for (var y of x.factors()) {
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
    if (left instanceof Expression.Exponentiation && getExponent(left).equals(Expression.ONE.negate())) {//TODO: FIX!!!
      if (right instanceof Expression.Matrix && !right.determinant().equals(Expression.ZERO)) {//TODO: condition - ?
        left = left.inverse();
        right = right.inverse();
        //TODO: add a step (?)
        //console.log(left.toString() + "=" + right.toString());
      }
    }
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
  Expression.E = new Expression.Symbol("\u2147"); // EulerNumber
  Expression.I = new Expression.Symbol("\u2148"); // ImaginaryUnit

  Expression.CIRCLE = new Expression.Symbol("○");

  Expression.prototype.addPosition = function () {
    return this;
  };

  //! 2018-09-30
  Expression.SystemOfEquations = function (equations) {
    throw new TypeError();
    this.equations = equations;
  };
  Expression.SystemOfEquations.from = function (equations) {
    return new Expression.NoAnswerExpression({matrix: null}, "system-of-equations", {equations: equations});
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
    var c = getConjugate(getBase(this)).pow(getExponent(this));
    return x.multiply(c).divide(this.multiply(c));
  };

  Expression.ExponentiationOfQuadraticInteger = function (x, y) {
    Expression.Exponentiation.call(this, x, y);
  };
  Expression.ExponentiationOfQuadraticInteger.prototype = Object.create(Expression.ExponentiationOfImaginaryUnit.prototype);
  Expression.ExponentiationOfQuadraticInteger.prototype.divideExpression = function (x) {
    return Expression.Exponentiation.prototype.divideExpression.call(this, x);
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
    //TODO:
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
  if (y instanceof Multiplication) {
    var x = this;
    return Multiplication.compare4Addition(x, y);
  }
  if (y instanceof Addition) {
    var x = this;
    return Addition.compare4Addition(x, y);
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
  if (y instanceof Expression.MatrixSymbol) {
    return -1;
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

function ExpressionWithCondition(e, variable, operator, value) {
  this.e = e;
  this.variable = variable;
  this.operator = operator;
  this.value = value;
}
ExpressionWithCondition.prototype = Object.create(Expression.prototype);
ExpressionWithCondition.prototype.toString = function () {
  return this.e.toString() + '; ' + this.variable.toString() + ' ' + this.operator + ' ' + this.value;
};
ExpressionWithCondition.prototype._conditionToMathML = function () {
  return this.variable.toMathML() + '<mo>' + this.operator.replace(/</g, '&lt;').replace(/>/g, '&gt;') + '</mo><mn>' + this.value + '</mn>';
};
ExpressionWithCondition.prototype.toMathML = function (options) {
  return this.e.toMathML(options) + '<mtext>; </mtext>' + this._conditionToMathML();
};

ExpressionWithCondition.prototype.multiplyMatrix = function (x) {
  //TODO: apply condition - ?
  return new ExpressionWithCondition(x.multiply(this.e), this.variable, this.operator, this.value);
};
ExpressionWithCondition.prototype.multiply = function (y) {
  return y.multiplyExpressionWithCondition(this);
};
ExpressionWithCondition.prototype.multiplyExpressionWithCondition = function (x) {
  throw new TypeError("TODO:");
};
Expression.prototype.multiplyExpressionWithCondition = function (x) {
  //TODO: apply condition - ?
  return new ExpressionWithCondition(x.e.multiply(this), x.variable, x.operator, x.value);
};


Expression.ExpressionWithCondition = ExpressionWithCondition;

  export default Expression;

  //TODO: ?
  Addition.prototype.compare4MultiplicationNthRoot = function (x) {
    return +1;
  };



// piecewise functions
// https://en.wikipedia.org/wiki/Piecewise
// https://en.wikibooks.org/wiki/LaTeX/Advanced_Mathematics#The_cases_environment

function Cases(cases) {
  this.cases = cases;
}

Cases.prototype = Object.create(Expression.prototype);

Cases.prototype.multiply = function (y) {
  return y.multiplyCases(this);
};
Cases.prototype.multiplyCases = function () {
  throw new TypeError("TODO: ");
};
Expression.prototype.multiplyCases = function (x) {
  var y = this;
  return new Cases(x.cases.map(function (c) {
    return c.multiply(y);
  }));
};

Cases.prototype.multiplyMatrix = function (x) {
  return new Cases(this.cases.map(function (c) {
    return x.multiply(c);
  }));
};

//Cases.prototype.add = function () {
//};

//Cases.prototype.toString = function () {
//};
Cases.prototype.toMathML = function (printOptions) {
  // https://www.w3.org/TR/2006/NOTE-arabic-math-20060131/#Moroccan
  var s = '';
  s += '<mrow>';
  s += '<mo>{</mo>';
  s += '<mtable rowspacing="0ex" columnalign="left">';
  for (var i = 0; i < this.cases.length; i += 1) {
    var x = this.cases[i];
    s += '<mtr>';
    s += '<mtd>';
    s += x.e.toMathML(printOptions);
    s += '</mtd>';
    s += '<mtd>';
    // <mtext> if </mtext> - ?
    s += x._conditionToMathML();
    s += '</mtd>';
    s += '</mtr>';
  }
  s += '</mtable>';
  s += '</mrow>';
  return s;
};

Expression.Cases = Cases;


Expression.Factorial = function (n) {
  this.n = n;
};
Expression.Factorial.prototype = Object.create(Expression.prototype);

Expression.prototype.factorial = function () {
  //a = a.unwrap();
  var n = this;
  if (!(n instanceof Expression.Integer)) {
    throw new TypeError("NotSupportedError");
  }
  if (n.compareTo(Expression.ZERO) < 0) {
    throw new TypeError("NotSupportedError");
  }
  var f = Expression.ONE;
  for (var i = n; i.compareTo(Expression.ONE) >= 0; i = i.subtract(Expression.ONE)) {
    f = f.multiply(i);
  }
  return f;
};

Expression.prototype._abs = function () {
  return this.isNegative() ? this.negate() : this;
};


  //?
  Expression.Comma = function (a, b) {
    BinaryOperation.call(this, a, b);
  };
  Expression.Comma.prototype = Object.create(BinaryOperation.prototype);
  Expression.Comma.prototype.getS = function () {
    // \u200B
    return ",";
  };
  Expression.prototype.transformComma = function (b) {
    var a = this;
    //if (a instanceof Expression.Equality && b instanceof Expression.Equality) {
    //}
    if (a instanceof Expression.NoAnswerExpression && a.name === 'polynomial-roots' && b instanceof Expression.NoAnswerExpression && b.name === 'polynomial-roots') {
      var ae = a.second.polynomial.calcAt(a.second.variable);
      var be = b.second.polynomial.calcAt(b.second.variable);
      //TODO: use original input expression
      return Expression.SystemOfEquations.from([{left: ae, right: Expression.ZERO}, {left: be, right: Expression.ZERO}]);
    }
    if (a instanceof Expression.NoAnswerExpression && a.name === 'system-of-equations' && b instanceof Expression.NoAnswerExpression && b.name === 'polynomial-roots') {
      var be = b.second.polynomial.calcAt(b.second.variable);
      //TODO: do not use NonSimplifiedExpression - ? and systemo-of-equations (?) or change (!)
      return Expression.SystemOfEquations.from(a.second.equations.concat([{left: be, right: Expression.ZERO}]));
    }
    throw new TypeError("NotSupportedError");
  };

Expression.Logarithm = function (argument) {
  Expression.Function.call(this, "log", argument);
};
Expression.Logarithm.prototype = Object.create(Expression.Function.prototype);
Expression.prototype.logarithm = function () {
  var arg = this;
  if (arg instanceof Expression.Integer) {
    if (arg.compareTo(Expression.ZERO) <= 0) {
      throw new RangeError("ArithmeticException");//TODO: better message
    }
    if (arg.compareTo(Expression.ONE) === 0) {
      return Expression.ZERO;
    }
    var p = integerPrimeFactor(arg);
    if (p.equals(arg)) {
      return new Expression.Logarithm(arg);
    }
    return p.logarithm().add(arg.truncatingDivide(p).logarithm());
  }
  if (arg instanceof Expression.Division) {
    var n = arg.getNumerator();
    var d = arg.getDenominator();
    if (d instanceof Integer || isConstant(n)) {
      return n.logarithm().subtract(d.logarithm());
    }
  }
  if (arg instanceof Expression.Multiplication) {
    var c = getConstant(arg);
    if (c instanceof Integer && !c.equals(Integer.ONE)) {
      return c.logarithm().add(arg.divide(c).logarithm());
    }
    if (arg.b instanceof Expression.MatrixSymbol) {//?
      return arg.a.logarithm().add(arg.b.logarithm());
    }
  }
  if (arg instanceof Expression.Symbol) {
    if (arg === Expression.E) {//?
      return Expression.ONE;
    }
    return new Expression.Logarithm(arg);
  }
  if (arg instanceof Expression.NthRoot) {
    var a = arg.a;
    var n = arg.n;
    if (a instanceof Expression.Integer && a.compareTo(Expression.ONE) > 0) {
      return a.logarithm().divide(Expression.Integer.fromNumber(n));
    }
  }
  if (arg instanceof Expression.Exponentiation) {
    var b = getBase(arg);
    var e = getExponent(arg);
    if (b === Expression.E) {//?
      return e;
    }
    if (b instanceof Expression.MatrixSymbol && e instanceof Integer) {//?
      return b.logarithm().multiply(e);
    }
    if (b instanceof Expression.Integer) {
      return b.logarithm().multiply(e);
    }
  }
  if (arg instanceof Expression.Matrix) {
    var matrix = arg.matrix;
    if (matrix.isDiagonal()) {
      return new Expression.Matrix(matrix.map(function (e, i, j) {
        return i === j ? e.logarithm() : Expression.ZERO;
      }));
    }
    var tmp = Expression.getEigenvalues(matrix);
    var eigenvalues = tmp.eigenvalues;
    var multiplicities = tmp.multiplicities;
    if (Expression.sum(multiplicities) === matrix.cols()) {
      var tmp2 = Expression.getEigenvectors(matrix, eigenvalues);
      var eigenvectors = tmp2.eigenvectors;
      if (eigenvectors.length === matrix.cols()) {
        if (Expression.callback != undefined) {
          Expression.callback(new Expression.Event("logarithm-using-diagonalization", new Expression.Matrix(matrix)));//TODO:
        }
        var tmp = Expression.diagonalize(matrix, eigenvalues, multiplicities, eigenvectors);
        // https://en.wikipedia.org/wiki/Logarithm_of_a_matrix#Calculating_the_logarithm_of_a_diagonalizable_matrix
        return new Expression.Matrix(tmp.T).multiply(new Expression.Matrix(tmp.L).logarithm()).multiply(new Expression.Matrix(tmp.T_INVERSED));
      } else {
        var tmp = Expression.getFormaDeJordan(matrix, eigenvalues, multiplicities);
        // https://en.wikipedia.org/wiki/Logarithm_of_a_matrix#The_logarithm_of_a_non-diagonalizable_matrix
        var J = tmp.J;
        var logarithmOfJordanBlockMatrix = function (B) {
          var K = B.map(function (e, i, j) {
            return i !== j ? e.divide(B.e(i, i)) : Expression.ZERO;
          });
          var S = B.map(function (e, i, j) {
            return i === j ? e.logarithm() : Expression.ZERO;
          });
          var n = B.cols();
          for (var i = 1; i < n; i += 1) {
            var x = K.pow(i).scale(Expression.ONE.divide(Expression.Integer.fromNumber(i)));
            S = i % 2 === 1 ? S.add(x) : S.subtract(x);
          }
          return S;
        };
        var LJ = logarithmOfJordanBlockMatrix(J);
        //if (!J.eql(matrix)) {
        //TODO:
        if (Expression.callback != undefined) {
          Expression.callback(new Expression.Event("logarithm-using-Jordan-canonical-form", new Expression.Matrix(matrix)));
        }
        //}
        return new Expression.Matrix(tmp.P).multiply(new Expression.Matrix(LJ)).multiply(new Expression.Matrix(tmp.P_INVERSED));
      }
    }
  }
  throw new TypeError("NotSupportedError");
};
