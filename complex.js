  import Expression from './Expression.js';
  import QuadraticInteger from './QuadraticInteger.js';
  import BigInteger from './BigInteger.js';
  import nthRoot from './nthRoot.compiled.js';

  var Integer = Expression.Integer;

  function Complex(real, imaginary) {
    //Expression.call(this);
    if (!(real instanceof Integer) || !(imaginary instanceof Integer) || imaginary.compareTo(Expression.ZERO) === 0) {
      throw new RangeError();
    }
    this.real = real;
    this.imaginary = imaginary;
  }

  Complex.prototype = Object.create(Expression.prototype);

  Expression.I = new Complex(Expression.ZERO, Expression.ONE);
  Expression.Complex = Complex;

  Complex.prototype.add = function (y) {
    return y.addComplex(this);
  };
  Expression.prototype.addComplex = function (x) {
    return this.addExpression(x);
  };
  Integer.prototype.addComplex = function (x) {
    return new Complex(x.real.add(this), x.imaginary);
  };
  Complex.prototype.addComplex = function (x) {
    var real = x.real.add(this.real);
    var imaginary = x.imaginary.add(this.imaginary);
    return imaginary.compareTo(Expression.ZERO) === 0 ? real : new Complex(real, imaginary);
  };
  Complex.prototype.addInteger = function (x) {
    return new Complex(x.add(this.real), this.imaginary);
  };

  Complex.prototype.equals = function (y) {
    return y instanceof Complex && this.real.equals(y.real) && this.imaginary.equals(y.imaginary) ? true : false;
  };

  Complex.prototype.compare4AdditionSymbol = function (x) {
    return +1;
  };
  Complex.prototype.compare4MultiplicationNthRoot = function (x) {
    return +1;
  };
  Complex.prototype.compare4Addition = function (y) {
    if (y instanceof Complex) {
      if (this.equals(y)) {
        return 0;
      }
      return this.real.compareTo(y.real) || this.imaginary.compareTo(y.imaginary);
    }
    if (y instanceof Integer) {
      return +1;
    }
    if (y instanceof Expression.Division) {
      return Expression.prototype.compare4Addition.call(this, y);
    }
    if (y instanceof Expression.Exponentiation) {
      return Expression.prototype.compare4Addition.call(this, y);
    }
    if (y instanceof Expression.Matrix) {
      return Expression.prototype.compare4Addition.call(this, y);
    }
    return -1;
  };
  // ? zero in case of same "base"
  Complex.prototype.compare4Multiplication = function (y) {
    if (y instanceof Complex) {
      if (y.equals(this)) {
        return 0;
      }
      return this.real.abs().compareTo(y.real.abs()) || this.imaginary.abs().compareTo(y.imaginary.abs());
      //return 0;
      //TODO: fix
      //throw new RangeError("NotSupportedError");//TODO:
    }
    return -1;//?
  };
  Complex.prototype.compare4MultiplicationSymbol = function (y) {
    return +1;
  };
  Complex.prototype.multiply = function (y) {
    return y.multiplyComplex(this);
  };
  Complex.prototype.multiplyComplex = function (x) {
    var real = x.real.multiply(this.real).subtract(x.imaginary.multiply(this.imaginary));
    var imaginary = x.real.multiply(this.imaginary).add(x.imaginary.multiply(this.real));
    return imaginary.compareTo(Expression.ZERO) === 0 ? real : new Complex(real, imaginary);
  };
  Expression.prototype.multiplyComplex = function (x) {
    return this.multiplyExpression(x);
  };
  Integer.prototype.multiplyComplex = function (x) {
    if (this.compareTo(Expression.ZERO) === 0) {
      return this;
    }
    return new Complex(x.real.multiply(this), x.imaginary.multiply(this));
  };
  Complex.prototype.multiplyInteger = function (x) {
    if (x.compareTo(Expression.ZERO) === 0) {
      return x;
    }
    return new Complex(x.multiply(this.real), x.multiply(this.imaginary));
  };

  Complex.prototype.conjugate = function () {
    return new Complex(this.real, this.imaginary.negate());
  };
  //Complex.prototype.divideExpression = function (x) {
  //  var y = this;
  //  return x.multiply(y.conjugate()).divide(y.multiply(y.conjugate()));
  //};
  Complex.prototype.getPrecedence = function () {
    return this.real.equals(Expression.ZERO) ? (this.imaginary.equals(Expression.ONE) ? 1000 : 3) : 2; // precedence.binary['+']
  };

  Complex.prototype.truncatingDivide = function (f) {
    if (f instanceof Integer) {
      return new Complex(this.real.truncatingDivide(f), this.imaginary.truncatingDivide(f));
    }
    return this.multiply(f.conjugate()).truncatingDivide(f.multiply(f.conjugate()));
  };

  Complex.prototype.toStringInternal = function (options, times, i, minus, plus, start, end, toString) {
    if (this.real.equals(Expression.ZERO)) {
      if (this.imaginary.equals(Expression.ONE)) {
        return i;
      }
      if (this.imaginary.equals(Expression.ONE.negate())) {
        return start + minus + i + end;
      }
      return start + toString(this.imaginary, options) + times + i + end;
    }
    var isNegative = this.imaginary.isNegative();
    var imaginary = (isNegative ? this.imaginary.negateCarefully() : this.imaginary);
    var si = (imaginary.equals(Expression.ONE) ? i : start + toString(imaginary, options) + times + i + end);
    var sr = toString(this.real, options);
    return start + sr + (isNegative ? minus : plus) + si + end;
  };

  Complex.prototype.toString = function (options) {
    return this.toStringInternal(options, "", "i", "-", "+", "", "", function (x, options) { return x.toString(options); });
  };

  Expression.getComplexConjugate = function (e) {
    if (!Expression.has(e, Complex)) {
      return undefined;
    }
    var c = Expression.ZERO;
    for (var x of e.summands()) {
      var f = undefined;
      for (var y of x.factors()) {
        if (y instanceof Complex) {
          f = y;
        }
      }
      if (f == undefined) {
        c = c.add(x);
      } else {
        var fc = f.conjugate();
        c = c.add(x.multiply(fc).divide(f.multiply(fc)).multiply(fc));
      }
    }
    //!?
    if (c.equals(e)) {
      return undefined;
    }
    //!?
    return c;
  };

  Complex.prototype.compare4MultiplicationInteger = function (y) {
    return +1;
  };

  Complex.prototype.remainderInteger = function (x) {
    return Complex.prototype.remainder.call(x, this);
  };

  Complex.prototype.remainder = function (y) {
    function norm(x) {
      return x instanceof Expression.Integer ? x.multiply(x) : x.multiply(x.conjugate());
    }
    function roundDivision(a, b) {
      if (b.compareTo(Expression.ZERO) < 0) {
        b = b.negate();
        a = a.negate();
      }
      var sign = a.compareTo(Expression.ONE) < 0 ? Expression.ONE.negate() : Expression.ONE;
      return a.add(b.truncatingDivide(Expression.TWO).multiply(sign)).truncatingDivide(b);
    }
    var x = this;
    var n = y instanceof Expression.Integer ? x : x.multiply(y.conjugate());
    var d = y instanceof Expression.Integer ? y : y.multiply(y.conjugate());
    //TODO: fix
    var q1 = n instanceof Complex ? roundDivision(n.real, d) : roundDivision(n, d);
    var q2 = n instanceof Complex ? roundDivision(n.imaginary, d) : Expression.ZERO;
    var q = q2.compareTo(Expression.ZERO) === 0 ? q1 : new Complex(q1, q2);
    var r =  x.subtract(y.multiply(q));
    if (norm(r).compareTo(norm(y)) >= 0) {
      throw new TypeError();
    }
    return r;
  };

  Complex.prototype.primeFactor = function () {

    function canBeSquare(n) {
      if (typeof n === "object") {
        return true;//TODO:
      }
      // https://www.johndcook.com/blog/2008/11/17/fast-way-to-test-whether-a-number-is-a-square/#comment-15700
      //var bitset = 0;
      //for (var i = 0; i < 32; i += 1) {
      //  bitset |= 1 << ((i * i) % 32);
      //}
      var bitset = 33751571;
      var result = (bitset >> Number(n & n.constructor(31))) & 1;
      return result === 1;
    }
    function norm(a, b) {
      return BigInteger.add(BigInteger.multiply(a, a), BigInteger.multiply(b, b));
    }
    function hasDivisor(r, i, a, b) {
      var d = BigInteger.add(BigInteger.multiply(a, a), BigInteger.multiply(b, b));
      var x = BigInteger.add(BigInteger.multiply(r, a), BigInteger.multiply(i, b));
      var y = BigInteger.subtract(BigInteger.multiply(i, a), BigInteger.multiply(r, b));
      return BigInteger.equal(BigInteger.remainder(x, d), 0) && BigInteger.equal(BigInteger.remainder(y, d), 0);
    }

    var r = this.real.toBigInt();
    var i = this.imaginary.toBigInt();
    var n = norm(r, i);
    //if (n > (9007199254740991 + 1) / 2) {
      //TODO: should not throw (see a call from Polynomial#getroots)
      //throw new RangeError("NotSupportedError");
    //}

    for (var fs = QuadraticInteger._factors(n), p = fs.next().value; p != null; p = fs.next().value) {
      var b = BigInteger.BigInt(0);
      var c = p;
      while (c > 0) {
        if (canBeSquare(c)) {
          var a = nthRoot(c, 2);
          if (BigInteger.equal(BigInteger.exponentiate(a, BigInteger.BigInt(2)), c)) {
            if (norm(a, b) > 1 && hasDivisor(r, i, a, b)) {
              return BigInteger.equal(b, 0) ? new Expression.Complex(Expression.ZERO, new Expression.Integer(a)) : new Complex(new Expression.Integer(a), new Expression.Integer(b));
            }
          }
        }
        b = BigInteger.add(b, BigInteger.BigInt(1));
        c = BigInteger.subtract(p, BigInteger.multiply(b, b));
      }
    }

    if (n > 1) {
      throw new TypeError();
    }
    return this;
  };

  Expression.Complex = Complex;

/*
//!
Expression.Complex.prototype.isValid = function () {
  return true;
};
//!
Expression.Complex.prototype.isPositive = function () {
  return this.imaginary.compareTo(Expression.ZERO) > 0;// || (this.imaginary.compareTo(Expression.ZERO) === 0 && this.real.compareTo(Expression.ZERO) > 0);
};
Expression.Complex.prototype.isUnit = function () {
  return this.multiply(this.conjugate()).equals(Expression.ONE);
};
Expression.Complex.prototype.truncatingDivideInteger = function (x) {
  debugger;
  return x.multiply(this.conjugate()).divide(this.multiply(this.conjugate()));
};

Expression.Complex.prototype.isDivisibleBy = function (y) {
  return !(this.divide(y) instanceof Expression.Division);
};
Expression.Complex.prototype.isDivisibleByInteger = function (x) {
  return !(x.multiply(this.conjugate()).divide(this.multiply(this.conjugate())) instanceof Expression.Division);
};
Expression.Complex.prototype.toExpression = function () {
  return this;
};
*/
