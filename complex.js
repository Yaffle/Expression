/*global Expression*/

(function (global) {
  "use strict";

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
      return 0;
    }
    if (y instanceof Integer) {
      return 0;
    }
    return -1;
  };
  Complex.prototype.compare4Multiplication = function (y) {
    if (y instanceof Complex) {
      if (y.equals(this)) {
        return 0;
      }
      return 0;
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
    // f instancoef Integer
    return new Complex(this.real.truncatingDivide(f), this.imaginary.truncatingDivide(f));
  };

  Complex.prototype.toStringInternal = function (options, times, i, minus, plus, toString) {
    var isNegative = this.imaginary.isNegative();
    var imaginary = (isNegative ? this.imaginary.negateCarefully() : this.imaginary);
    var si = imaginary.equals(Expression.ONE) ? "" : toString(imaginary, options) + times;
    var sr = this.real.equals(Expression.ZERO) ? "" : toString(this.real, options);
    return sr + (isNegative ? minus : (this.real.equals(Expression.ZERO) ? "" : plus)) + si + i;
  };

  Complex.prototype.toString = function (options) {
    return this.toStringInternal(options, "", "i", "-", "+", function (x, options) { return x.toString(options); });
  };

  Expression.getComplexConjugate = function (e) {
    if (!Expression.has(e, Complex)) {
      return undefined;
    }
    var c = Expression.ZERO;
    for (var additions = e.summands(), x = additions.next().value; x != null; x = additions.next().value) {
      var f = undefined;
      for (var multiplications = x.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
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

  Expression.Complex = Complex;
}());
