// https://en.wikipedia.org/wiki/Quadratic_integer
// https://en.wikipedia.org/wiki/Factorization#Unique_factorization_domains


//import primeFactor from './primeFactor.js';
/*

      if (qi instanceof QuadraticInteger && qi.a < 0 || qi.b < 0) {
        var c = pow(qi.conjugate(), n, new QuadraticInteger(1, 0, qi.D));
        return qi.multiply(c).toExpression()._nthRoot(n).divide(new QuadraticInteger(Math.abs(qi.a), Math.abs(qi.b), qi.D).toExpression());
      }

*/

// a + b*sqrt(D)
function QuadraticInteger(a, b, D) {
  //TODO:
  if (Math.abs(a) > 9007199254740991) {throw new Error();}
  if (Math.abs(b) > 9007199254740991) {throw new Error();}
  this.a = a;
  this.b = b;
  this.D = D;
}
QuadraticInteger.prototype.multiply = function (y) {
  var x = this;
  if (x.D !== y.D) {
    throw new Error();
  }
  return new QuadraticInteger(x.a * y.a + x.b * y.b * y.D, x.a * y.b + x.b * y.a, x.D);
};
QuadraticInteger.prototype.conjugate = function (y) {
  return new QuadraticInteger(this.a, -this.b, this.D);
};
QuadraticInteger.prototype.norm = function () {
  return this.a * this.a - this.b * this.b * this.D;
};
QuadraticInteger.prototype.truncatingDivideInteger = function (x) {
  return new QuadraticInteger(x.toNumber(), 0, this.D).truncatingDivide(this);
};
QuadraticInteger.prototype.truncatingDivide = function (y) {
  if (!(y instanceof QuadraticInteger)) {
    y = new QuadraticInteger(y.toNumber(), 0, this.D);
  }
  var x = this;
  if (x.D !== y.D) {
    throw new Error();
  }
  var n = x.multiply(y.conjugate());
  var d = y.norm();
  return n.a % d === 0 && n.b % d === 0 ? new QuadraticInteger(n.a / d, n.b / d, x.D) : null;
};
QuadraticInteger.prototype.negate = function () {
  return new QuadraticInteger(0 - this.a, 0 - this.b, this.D);
};
QuadraticInteger.prototype.primeFactor = function () {
  function primeFactor(n) {
    var i = 2;
    var s = 0;
    var r = Math.floor(Math.sqrt(n + 0.5));
    while (i <= r) {
      if (n % i === 0) {
        return i;
      }
      i += s === 2 ? 2 : s + 1;
      s += 1;
      if (s === 4) {
        s = 2;
      }
    }
    return n;
  }
  function factors(n) {
    var p = n > 1 ? primeFactor(n) : 1;
    var e = 0;
    var t = 1;
    var f = null;
    var fs = null;
    var i = 1;
    return {
      //[Symbol.iterator]: function () {
      //  return this;
      //},
      //get done() {
      //  return this.value == null;
      //},
      value: null,
      next: function () {
        if (p === 1) {
          this.value = null;
          return this;
        }
        if (fs == null) {
          if (n % p === 0) {
            t *= p;
            n /= p;
            e += 1;
            this.value = t;
            return this;
          }
          fs = factors(n);
          i = t;
        }
        if (i === t) {
          i = 1;
          f = fs.next().value;
        } else {
          i *= p;
        }
        this.value = f == null ? null : f * i;
        return this;
      }
    };
  }

  function gcd(a, b) {
    return b === 0 ? a : gcd(b, a % b);
  }
  
  var a = this.a;
  var b = this.b;
  var D = this.D;
  var g = gcd(a, b);
  if (g !== 1 && g !== -1) {
    //TODO:
    return new QuadraticInteger(primeFactor(Math.abs(g)), 0, D);
  }
  var norm = this.norm();
  function quadraticIntegers(norm, D, b) {
    while (true) {
      if (-norm + D * b * b > 9007199254740991 || norm + D * b * b > 9007199254740991) {
        throw new RangeError(norm);
      }
      var a = Math.floor(Math.sqrt(norm + D * b * b));
      if (norm + D * b * b === a * a) { // && Math.abs(gcd(a, b)) === 1
        return new QuadraticInteger(a, b, D);
      }
      var a = Math.floor(Math.sqrt(-norm + D * b * b));
      if (-norm + D * b * b === a * a) { // && Math.abs(gcd(a, b)) === 1
        return new QuadraticInteger(a, b, D);
      }
      b += 1;
    }
  }
  if (a === 0 || b === 0) {
    return this;
  }
  if (norm === 1 || norm === -1) {
    // https://en.wikipedia.org/wiki/Quadratic_field#Orders_of_quadratic_number_fields_of_small_discriminant
    var unit = quadraticIntegers(1, D, 1);
    var uniti = unit.conjugate();
    var x = this;
    if (x.b * x.a < 1) {
      return x.a < 0 ? uniti.negate() : uniti;
    }
    return unit;
  }

  var v = this;
  for (var fs = factors(Math.abs(norm)), p = fs.next().value; p != null; p = fs.next().value) {
    // ? https://www.johndcook.com/blog/2008/11/17/fast-way-to-test-whether-a-number-is-a-square/
    /*if (D === 17) {
      var t = Math.abs(norm);
      while (t % 2 === 0) {
        t /= 2;
      }
      p = t === 1 ? norm : primeFactor(t);
    }*/
    if (p === 2 && (D === 5 || D === 13 || D === 17 || D === 21 || D === 29 || D === 33 || D === 41 || D === 57)) {
      p = 4;
    }
    var i = quadraticIntegers(p, D, 0);
    //console.log(i + '');
    var x = v.truncatingDivide(i);
    if (x != null) {
      return i;
    }
    var x = v.truncatingDivide(i.conjugate());
    if (x != null) {
      return i.conjugate();
    }// 1+9sqrt(2)
  }
  console.log('!');
  return this;
  //throw new Error();
};
QuadraticInteger.prototype.toString = function () {
  return this.a + '+' + this.b + 'sqrt(' + this.D + ')';
};
QuadraticInteger.prototype.isUnit = function () {
  return Math.abs(this.norm()) === 1;
};
QuadraticInteger.prototype.equals = function (y) {
  var x = this;
  if (y === Expression.ZERO) {
    return x.a === 0 && x.b === 0;
  }
  if (y === Expression.ONE) {
    return x.a === 1 && x.b === 0;
  }
  if (!(y instanceof QuadraticInteger)) {
    throw new TypeError();
  }
  return x.a === y.a && x.b === y.b && x.D === y.D;
};
QuadraticInteger.prototype.subtract = function (y) {
  var x = this;
  console.assert(x.D === y.D);
  return new QuadraticInteger(x.a - y.a, x.b - y.b, x.D);
};
QuadraticInteger.prototype.remainder = function (y) {
  if (!(y instanceof QuadraticInteger)) {
    y = new QuadraticInteger(y.toNumber(), 0, this.D);
  }
  var x = this;
  console.assert(x.D === y.D);
  var n = x.multiply(y.conjugate());
  var d = y.norm();
  var q1 = Math.trunc(n.a / d);
  var q2 = Math.trunc(n.b / d);
  var q = new QuadraticInteger(q1, q2, x.D);
  return x.subtract(y.multiply(q));
};
QuadraticInteger.prototype.remainderInteger = function (x) {
  return new QuadraticInteger(x, 0, this.D).remainder(this);
};
QuadraticInteger.prototype.toExpression = function () {
  return Expression.Integer.fromNumber(this.a).add(Expression.Integer.fromNumber(this.b).multiply(Expression.Integer.fromNumber(this.D).squareRoot()));
};

export default QuadraticInteger;


// new QuadraticInteger(1, 1, 2).remainder(new QuadraticInteger(1, 1, 2))
// new QuadraticInteger(2, 2, 2).truncatingDivide(new QuadraticInteger(2, 2, 2))

//!
function toQuadraticInteger(e) {
  if (e instanceof Expression.Addition &&
      e.b instanceof Expression.Integer &&
      e.a instanceof Expression.SquareRoot &&
      e.a.a instanceof Expression.Integer) {
    return new QuadraticInteger(e.b.toNumber(), Expression.ONE.toNumber(), e.a.a.toNumber());
  }
  if (e instanceof Expression.Addition &&
      e.b instanceof Expression.Integer &&
      e.a instanceof Expression.Multiplication &&
      e.a.a instanceof Expression.Integer &&
      e.a.b instanceof Expression.SquareRoot &&
      e.a.b.a instanceof Expression.Integer) {
    return new QuadraticInteger(e.b.toNumber(), e.a.a.toNumber(), e.a.b.a.toNumber());
  }
}
//!

QuadraticInteger.toQuadraticInteger = toQuadraticInteger;

/*
QuadraticInteger.prototype.compareTo = function (e) {
  if (e === Expression.ZERO) {
    var n = this.a * this.a - this.b * this.b * this.D;
    return this.a === 0 && this.b === 0 ? 0 : (this.a < 0 && this.b < 0 || this.a < 0 && n > 0 || this.b < 0 && n < 0 ? -1 : 1);
  }
  if (e === Expression.ONE) {
    return this.a === 1 && this.b === 0 ? 0 : 1;
  }
  throw new TypeError();
};
*/



  //import './QuadraticInteger.js';

/*

    if (n === 2) {
      function gcd(a, b) {
        return b === 0 ? a : gcd(b, a % b);
      }
      var q = isQuadraticInteger(x);
      //TODO: (q.D === 2 || q.D === 3 || q.D === 5 || q.D === 17)
      
      if (q != null && q.D === 2 && Math.abs(q.a * q.a - q.b * q.b * q.D) === Math.pow(gcd(q.a, q.b), 2)) {
        var ff = Expression.ONE;
        if (q.a % q.D === 0) {
          q = {
            a: q.b,
            b: Math.floor(q.a / q.D),
            D: q.D
          };
          ff = new SquareRoot(new Integer(q.D));
          x = x.divide(ff);
        }
        var n = q.a * q.a - q.b * q.b * q.D;
        if (q.a > 0 && q.b > 0 || q.a > 0 && n > 0 || q.b > 0 && n < 0) {
          var t = new QuadraticInteger(q.a, q.b, q.D);
          if (t.primeFactor().equals(t)) {
            return (new SquareRoot(x.multiply(ff)));
          }
          var k1 = new QuadraticInteger(1, 0, q.D);
          var k2 = new QuadraticInteger(1, 0, q.D);
          var i = t;
          var p = null;
          while (!i.equals(Expression.ONE)) {
            var d = i.primeFactor();
            if (p == null) {
              p = d;
            } else {
              if (p.equals(d)) {
                k1 = k1.multiply(d);
                p = null;
              } else {
                k2 = k2.multiply(p);
                p = d;
              }
            }
            i = i.truncatingDivide(d);
          }
          if (p != null) {
            k2 = k2.multiply(p);
          }
          return k1.toExpression().multiply(new Expression.SquareRoot(k2.toExpression().multiply(ff)));
        }
      }
    }
*/


/*

+    var y = evaluateExpression(e.a, context);
+    if (y === "CANNOT_DIVIDE" || y == null) {
+      return y;
+    }
+    //TODO: debug
+    var yy = new Interval(context.nthRoot(y.a, n).a, context.nthRoot(y.b, n).b);
+    var s = context.nthRoot(context.scalingCoefficient, n);
+    yy = context.divide(yy, context.multiply(Interval.degenerate(context.scalingCoefficient), s));
+    return yy;
*/







/*



// +1, -1, +i, -i
// a+bi

// a === 0, i*(a+bi)
// a < 0, -(a+bi)
// b < 0, i*(a+bi)

// a > 0, b > 0


/*
  for (var a = 1; a * a <= n; a += 1) {
    for (var b = 0; b * b <= n - a * a; b += 1) {
      if (norm(a, b) > 1 && hasDivisor(r, i, a, b)) {
        return [a, b];
      }
      if (norm(a, -b) > 1 && hasDivisor(r, i, a, -b)) {
        return [a, -b];
      }
    }
  }
  return [r, i];
*/

/*

function primeFactor(n) {
  var i = 2;
  var s = 0;
  var r = Math.floor(Math.sqrt(n + 0.5));
  while (i <= r) {
    if (n % i === 0) {
      return i;
    }
    i += s === 2 ? 2 : s + 1;
    s += 1;
    if (s === 4) {
      s = 2;
    }
  }
  return n;
}

function isPrime(n) {
  return primeFactor(n) === n;
}

function norm(x) {
  return x instanceof Expression.Integer ? x.multiply(x).value : x.multiply(x.conjugate()).value;
}

function checkFactorization(i) {
  var results = [];
  var x = i;
  while (norm(x) > 1) {
    var p = x.primeFactor();
    results.push(p);

    //A Gaussian integer a + bi is a Gaussian prime if and only if either:
    //  one of a, b is zero and absolute value of the other is a prime number of the form 4n + 3 (with n a nonnegative integer), or
    //  both are nonzero and a**2 + b**2 is a prime number (which will not be of the form 4n + 3).
    var n = norm(p);
    console.assert(isPrime(n) || ((p instanceof Expression.Integer || p.real.equals(Expression.ZERO) || p.imaginary.equals(Expression.ZERO)) && Math.abs(p instanceof Expression.Integer ? p : p.real.add(p.imaginary).value) % 4 === 3), n, p.toString());

    x = x.divide(p);
    if (x instanceof Expression.Integer && norm(x) > 1) {
      x = new Expression.Complex(Expression.ZERO, x);
    }
  }
  console.log(i + '=' + results.map(x => '(' + x + ')').join(''));
}


// 5-5i

checkFactorization(new Complex(new Expression.Integer(3), new Expression.Integer(3)));

checkFactorization(new Complex(new Expression.Integer(0), new Expression.Integer(2)));
checkFactorization(new Complex(new Expression.Integer(0), new Expression.Integer(-2)));
checkFactorization(new Complex(new Expression.Integer(5), new Expression.Integer(1)));

var A = 11;
for (var i = -A; i <= A; i += 1) {
  for (var j = -A; j <= A; j += 1) {
    if (j !== 0) {
      checkFactorization(new Complex(new Expression.Integer(i), new Expression.Integer(j)));
    }
  }
}


*/

/*
          if (i == null && isOnePlusSqrtOf2(y.a)) {
            i = y.a;
          }
          if (i == null) {
            throw new TypeError();
          }
          } else if (isOnePlusSqrtOf2(y.a)) {
            if (!p.equals(y.a)) {
              throw new Error();
            }
            degree += 1;

*/


// http://oeis.org/wiki/Quadratic_integer_rings
// https://oeis.org/A048981
// https://en.wikipedia.org/wiki/Euclidean_domain#Norm-Euclidean_fields
// https://en.wikipedia.org/wiki/Fundamental_unit_(number_theory)
// https://en.wikipedia.org/wiki/Pell%27s_equation
// https://en.wikipedia.org/wiki/Diophantine_equation
// https://ru.wikipedia.org/wiki/Гауссовы_целые_числа#Определение
// https://en.wikipedia.org/wiki/Gaussian_integer
// https://en.wikipedia.org/wiki/Prime_element
