// https://en.wikipedia.org/wiki/Quadratic_integer
// https://en.wikipedia.org/wiki/Factorization#Unique_factorization_domains


// It is possible to use the comparision operators if a is a safe integer or BigIntegerInternal or BigInt and n is a safe integer:
// a < n
// a <= n
// a > n
// a >= n
// a == n
// a != n

import Expression from './Expression.js';
import nthRoot from './nthRoot.js';
import bigIntGCD from './gcd.js';

import primeFactor from './primeFactor.js';

/*

      if (qi instanceof QuadraticInteger && qi.a < 0 || qi.b < 0) {
        var c = pow(qi.conjugate(), n, new QuadraticInteger(1, 0, qi.D));
        return qi.multiply(c).toExpression()._nthRoot(n).divide(new QuadraticInteger(Math.abs(qi.a), Math.abs(qi.b), qi.D).toExpression());
      }

*/

function abs(a) {
  return a < 0n ? -a : a;
}

  function ngcd(a, b) {
    return BigInt(bigIntGCD(a, b));
  }

// a + b*sqrt(D)
function QuadraticInteger(a, b, D) {
  a = typeof a === "number" ? BigInt(a) : a;
  b = typeof b === "number" ? BigInt(b) : b;
  D = Number(D);
  //TODO:
  if (typeof a === "number" && Math.abs(a) > Number.MAX_SAFE_INTEGER) {
    throw new TypeError();
  }
  if (typeof b === "number" && Math.abs(b) > Number.MAX_SAFE_INTEGER) {
    throw new TypeError();
  }
  this.a = a;
  this.b = b;
  this.D = D;
}
QuadraticInteger.prototype.multiply = function (y) {
  var x = this;
  if (!(x.D === y.D)) {
    throw new TypeError();
  }
  var a = x.a * y.a + x.b * y.b * BigInt(y.D);
  var b = x.a * y.b + x.b * y.a;
  return new QuadraticInteger(a, b, x.D);
};
QuadraticInteger.prototype.conjugate = function (y) {
  return new QuadraticInteger(this.a, -this.b, this.D);
};
QuadraticInteger.prototype.norm = function () {
  //var x = this.a * this.a;
  //var y = this.b * this.b;
  //return x % this.D + (((x - x % this.D) / this.D) - y) * this.D;
  var a = this.a;
  var b = this.b;
  var D = this.D;
  var aa = a * a;
  var bb = b * b;
  var norm = aa - bb * BigInt(D);
  if (typeof norm === "number" && Math.abs(norm) > Number.MAX_SAFE_INTEGER) {
    throw new TypeError();
  }
  if (typeof this.a === "number" && norm >= -Number.MAX_SAFE_INTEGER && norm <= +Number.MAX_SAFE_INTEGER) {
    norm = Number(norm);
  }
  return norm;
};
QuadraticInteger.prototype.truncatingDivideInteger = function (x) {
  return new QuadraticInteger(x.toBigInt(), Expression.ZERO.toBigInt(), this.D).truncatingDivide(this);
};
QuadraticInteger.prototype.truncatingDivide = function (y) {
  if (!(y instanceof QuadraticInteger)) {
    if (y instanceof AlmostQuadraticInteger) {
      return null;
    }
    y = new QuadraticInteger(y.toBigInt(), Expression.ZERO.toBigInt(), this.D);
  }
  var x = this;
  if (!(x.D === y.D)) {
    throw new TypeError();
  }
  var n = x.multiply(y.conjugate());
  var d = y.norm();
  return n.a % d === 0n && n.b % d === 0n ? new QuadraticInteger(n.a / d, n.b / d, x.D) : null;
};

QuadraticInteger.prototype.negate = function () {
  return new QuadraticInteger(-this.a, -this.b, this.D);
};

/*
function primeFactor(n) {
  var i = n - n;
  ++i;
  ++i;
  if (n % i == 0) {
    return i;
  }
  ++i;
  while (i * i <= n) {
    if (n % i == 0) {
      return i;
    }
    ++i;
    ++i;
  }
  return n;
}
*/

function factors(n) {
  if (n < 1n) {
    throw new TypeError();
  }
  var p = n > 1n ? BigInt(primeFactor(n)) : 1n;
  var t = 1n;
  var f = 1n;
  var fs = null;
  var i = 1n;
  var result = {
    done: false,
    value: null,
    next: function () {
      if (p === 1n) {
        this.value = null;
        this.done = true;
        return this;
      }
      if (fs == null) {
        if (n % p === 0n) {
          t = t * p;
          n = n / p;
          this.value = t;
          this.done = false;
          return this;
        }

        fs = factors(n);
        i = t;
      }
      if (i === t) {
        i = 1n;
        f = fs.next().value;
      } else {
        i = i * p;
      }
      this.value = f == null ? null : f * i;
      this.done = f == null;
      return this;
    }
  };
  result[Symbol.iterator] = function () {
    return this;
  };
  return result;
}

QuadraticInteger._factors = function (n) {
  return factors(BigInt(n));
};

  
  QuadraticInteger._complexIntegerPrimeFactor = function (r, i) {
    r = BigInt(r);
    i = BigInt(i);

    function canBePerfectSquare(n) {
      // https://www.johndcook.com/blog/2008/11/17/fast-way-to-test-whether-a-number-is-a-square/#comment-15700
      //var bitset = 0;
      //for (var i = 0; i < 32; i += 1) {
      //  bitset |= 1 << ((i * i) % 32);
      //}
      var bitset = 33751571;
      var result = (bitset >> Number(n % 32n)) % 2;
      return result === 1;
    }
    function norm(a, b) {
      return a * a + b * b;
    }
    function hasDivisor(r, i, a, b) {
      var d = a * a + b * b;
      var x = r * a + i * b;
      var y = i * a - r * b;
      return x % d === 0n && y % d === 0n;
    }

    var n = norm(r, i);
    //if (n > (Number.MAX_SAFE_INTEGER + 1) / 2) {
      //TODO: should not throw (see a call from Polynomial#getroots)
      //throw new RangeError("NotSupportedError");
    //}

    for (var p of QuadraticInteger._factors(n)) {
      var b = 0n;
      var c = p;
      while (c > 0n) {
        if (canBePerfectSquare(c)) {
          var a = nthRoot(c, 2);
          if (a * a === c) {
            if (norm(a, b) > 1n && hasDivisor(r, i, a, b)) {
              return b === 0n ? new Expression.Complex(Expression.ZERO, Expression.Integer.fromBigInt(a)) : new Expression.Complex(Expression.Integer.fromBigInt(a), Expression.Integer.fromBigInt(b));
            }
          }
        }
        b = b + 1n;
        c = p - b * b;
      }
    }

    if (n > 1n) {
      throw new TypeError();
    }
    return new Expression.Complex(Expression.Integer.fromBigInt(r), Expression.Integer.fromBigInt(i));
  };

    //  Compute (k n), where k is numerator
    function jacobiSymbol(k, n) {
        if ( k < 0 || n % 2 == 0 ) {
            throw new Error("Invalid value. k = " + k + ", n = " + n);
        }
        k %= n;
        var jacobi = 1;
        while ( k > 0 ) {
            while ( k % 2 == 0 ) {
                k /= 2;
                var r = n % 8;
                if ( r == 3 || r == 5 ) {
                    jacobi = -jacobi;
                }
            }
            var temp = n;
            n = k;
            k = temp;
            if ( k % 4 == 3 && n % 4 == 3 ) {
                jacobi = -jacobi;
            }
            k %= n;
        }
        if ( n == 1 ) {
            return jacobi;
        }
        return 0;
    }

// D**((p - 1) / 2) % p
function LegendreSymbol(a, p) {
  // https://en.wikipedia.org/wiki/Legendre_symbol#:~:text=Special%20formulas
  // https://rosettacode.org/wiki/Jacobi_symbol#Java
  return jacobiSymbol(a, p);
}

QuadraticInteger.prototype.primeFactor = function () {

  var a = this.a;
  var b = this.b;
  var D = BigInt(this.D);
  var g = ngcd(a, D);

  if (g === D) {//TODO: g != 1 - ?
    return new QuadraticInteger(0n, 1n, D);
  }

  var g = ngcd(a, b);
  //!
  //while (g % 2 == 0) {
  //  g = g / 2;
  //}
  //!
  if (g !== 1n) {
    //TODO:
    //return new QuadraticInteger(primeFactor(g), 0, D);
    //return new QuadraticInteger(QuadraticInteger._factors(g).next().value, b - b, D);
  }
  var norm = this.norm();
  function quadraticIntegers(norm, D, b) {
    while (true) {
      var bbD = b * b * D;
      var guess1 = norm + bbD;
      var guess2 = -norm + bbD;
      //if (typeof norm === "number") {//TODO:
      if (Number(guess2) > Number.MAX_SAFE_INTEGER || Number(guess1) > Number.MAX_SAFE_INTEGER) {
        throw new RangeError(norm);
      }
      //}
      var guess = guess1;
      if (guess >= 0n) {
        var a = BigInt(nthRoot(guess, 2));
        if (guess === a**2n) { // && ngcd(a, b) === 1
          return new QuadraticInteger(a, b, D);
        }
      }
      var guess = guess2;
      if (guess >= 0n) {
        var a = BigInt(nthRoot(guess, 2));
        if (guess === a**2n) { // && ngcd(a, b) === 1
          return new QuadraticInteger(a, b, D);
        }
      }
      b = b + 1n;
    }
  }
  //if (a == 0 || b == 0) {
    //return this;
  //}
  if (norm === 1n || norm === -1n) {
    // https://en.wikipedia.org/wiki/Quadratic_field#Orders_of_quadratic_number_fields_of_small_discriminant
    var unit = quadraticIntegers(1n, D, 1n);
    var uniti = unit.conjugate();
    var x = this;
    if (x.b * x.a < 1n) {
      return x.a < 0n ? uniti.negate() : uniti;
    }
    return unit;
  }
  var v = this;
  for (var fs = QuadraticInteger._factors(abs(norm)), p = fs.next().value; p != null; p = fs.next().value) {
  //if (p * p <= norm || abs(norm) === p) {
    // ? https://www.johndcook.com/blog/2008/11/17/fast-way-to-test-whether-a-number-is-a-square/
    /*if (D === 17) {
      var t = Math.abs(norm);
      while (t % 2 === 0) {
        t /= 2;
      }
      p = t === 1 ? norm : primeFactor(t);
    }*/
    if (D % 4n === 1n && p === 2n) {
      continue;
    }
    if (D === 37n && (p === 7n || p === 14n || p === 3n || p === 11n)) {
      continue;
    }
    //?
    // http://oeis.org/wiki/Quadratic_integer_rings
    if (isPrime(p) && D % p !== 0n && p !== 2n && LegendreSymbol(Number(D), Number(p)) !== 1) {
      continue;
    }
    //?
    var i = quadraticIntegers(p, D, 0n);
    //console.log(i + '');
    var x = v.truncatingDivide(i);
    if (x != null) {
      return this.equals(i) ? i : i.primeFactor();
    }
    var x = v.truncatingDivide(i.conjugate());
    if (x != null) {
      return this.equals(i.conjugate()) ? i.conjugate() : i.conjugate().primeFactor();
    }// 1+9sqrt(2)
  //}
  }

  //console.log('!');
  return this;
  //throw new TypeError();
};
QuadraticInteger.prototype.toString = function () {
  return this.a.toString() + '+' + this.b.toString() + 'sqrt(' + this.D + ')';
};
QuadraticInteger.prototype.isUnit = function () {
  var n = this.norm();
  return n === 1n || n === -1n;
};
QuadraticInteger.prototype.equals = function (y) {
  var x = this;
  if (!(y instanceof QuadraticInteger)) {
    if (y.equals(Expression.ZERO)) {
      return x.a === 0n && x.b === 0n;
    }
    if (y.equals(Expression.ONE)) {
      return x.a === 1n && x.b === 0n;
    }
    throw new TypeError();
  }
  return x.a === y.a && x.b === y.b && x.D === y.D;
};
QuadraticInteger.prototype.subtract = function (y) {
  var x = this;
  if (!(x.D === y.D)) {
    throw new TypeError();
  }
  return new QuadraticInteger(x.a - y.a, x.b - y.b, x.D);
};
QuadraticInteger.prototype.isDivisibleBy = function (y) {
  if (y instanceof AlmostQuadraticInteger) {
    var q = this.truncatingDivide(y.qi);
    return q == null ? null : q.truncatingDivide(y.k);
  }
  return this.truncatingDivide(y) != null;
};
QuadraticInteger.prototype.isDivisibleByInteger = function (x) {
  return x.truncatingDivide(this) != null;
};
QuadraticInteger.prototype.remainder = function (y) {
  if (!(y instanceof QuadraticInteger)) {
    if (y instanceof AlmostQuadraticInteger) {
      if (ngcd(this.a, this.b) !== 0n) {
        //?TODO:
      }
      var nk = y.k;
      var remainder = this.remainder(y.qi);
      if (remainder.equals(Expression.ZERO)) {
        return remainder;
      }
      if (remainder.b === 0n) {
        return nk;
      }
      if (ngcd(remainder.a, remainder.b) !== 1n) {
        var i = Expression.Integer.fromBigInt(ngcd(remainder.a, remainder.b));
        nk = nk.multiply(i);
        remainder = remainder.truncatingDivide(i);
      }
      return new AlmostQuadraticInteger(nk, remainder);
    }
    if (y instanceof Expression.Multiplication && y.a instanceof Expression.Integer && y.b instanceof Expression.SquareRoot) {
      return this.remainder(new QuadraticInteger(Expression.ZERO.toBigInt(), y.a.toBigInt(), y.b.a.toBigInt()));
    }
    y = new QuadraticInteger(y.toBigInt(), Expression.ZERO.toBigInt(), this.D);
  }
  var x = this;
  if (!(x.D === y.D)) {
    throw new TypeError();
  }
  var n = x.multiply(y.conjugate());
  var d = y.norm();
  if (d === 1n || d === -1n) { // y.isUnit()
    return x.subtract(x);
  }

  var q1 = (n.a - n.a % d) / d;
  var q2 = (n.b - n.b % d) / d;
  if (q1 === 0n && q2 === 0n) {
    //if (abs(x.norm()) >= abs(y.norm())) {
      //?
      if (x.a > y.a && y.a > 0n) {
        return x.subtract(y.multiply(new QuadraticInteger((x.a - x.a % y.a) / y.a, 0n, x.D)));
      }
      if (x.a > -y.a && y.a < 0n) {
        return x.subtract(y.multiply(new QuadraticInteger((x.a - x.a % -y.a) / -y.a, 0n, x.D)));
      }
      if (y.b === 0n) {
        return new QuadraticInteger(1n, 0n, x.D); //?
      }
    //throw new RangeError("NotSupportedError");//TODO:!!!
    //}
  }
  var q = new QuadraticInteger(q1, q2, x.D);
  var r = x.subtract(y.multiply(q));
  return r;
};
QuadraticInteger.prototype.remainderInteger = function (x) {
  return new QuadraticInteger(x.toBigInt(), Expression.ZERO.toBigInt(), this.D).remainder(this);
};
QuadraticInteger.prototype.toExpression = function () {
  return Expression.Integer.fromBigInt(this.a).add(Expression.Integer.fromBigInt(this.b).multiply(Expression.Integer.fromNumber(this.D).squareRoot()));
};

QuadraticInteger.prototype.abs = function () {
  if (this.a <= 0n && this.b <= 0n ||
      this.a < 0n && this.norm() > 0n ||
      this.b < 0n && this.norm() < 0n) {
    return this.negate();
  }
  return this;
};

//TODO: merge with the QuadraticInteger.toQuadraticInteger
QuadraticInteger.prototype.isValid = function () {
  if (true) {
    //return false;
  }
  if (this.D === 5) {
    return true;//TODO:!!!
  }
  if (this.D === 37) {//TODO: ?
    return true;
  }
  if (this.D === 6) {//TODO: ?
    return false;
  }
  if (this.D % 4 === 1) {//TODO: ?
    return false;
  }
  if ([2, 3, 5, 6, 7, 11, 13, 17, 19, 21, 29, 33, 37, 41, 57, 73].indexOf(this.D) === -1) { // https://oeis.org/A048981
    return false;
  }
  return true;
};
QuadraticInteger.prototype.isPositive = function () {
  var qq = this;
  return qq.a > 0n && qq.b > 0n || qq.a > 0n && qq.norm() > 0n || qq.b > 0n && qq.norm() < 0n;
};

export default QuadraticInteger;

//new QuadraticInteger(-67, 15, 37).primeFactor()

// new QuadraticInteger(1, 1, 2).remainder(new QuadraticInteger(1, 1, 2))
// new QuadraticInteger(2, 2, 2).truncatingDivide(new QuadraticInteger(2, 2, 2))

function toQuadraticInteger(e) {
  //if (e instanceof Expression.Complex) {//!
  //  return e;
  //}
  if (e instanceof Expression.Addition) {
    var g = e.a.gcd(e.b);
    if (!g.equals(Expression.ONE) && g instanceof Expression.Integer) {
      var qi = toQuadraticInteger(e.divide(g));
      return qi == null ? null : new AlmostQuadraticInteger(g, qi);
    }
  }
  // qq.a * qq.a + qq.D * qq.b * qq.b < Number.MAX_SAFE_INTEGER
  if (e instanceof Expression.Addition &&
      e.b instanceof Expression.Integer &&
      e.a instanceof Expression.SquareRoot &&
      e.a.a instanceof Expression.Integer) {
    return new QuadraticInteger(e.b.toBigInt(), Expression.ONE.toBigInt(), e.a.a.toBigInt());
  }
  if (e instanceof Expression.Addition &&
      e.b instanceof Expression.Integer &&
      e.a instanceof Expression.Multiplication &&
      e.a.a instanceof Expression.Integer &&
      e.a.b instanceof Expression.SquareRoot &&
      e.a.b.a instanceof Expression.Integer) {
    return new QuadraticInteger(e.b.toBigInt(), e.a.a.toBigInt(), e.a.b.a.toBigInt());
  }
}
//!

QuadraticInteger.toQuadraticInteger = toQuadraticInteger;

QuadraticInteger.gcd = function (x, y) {
  var a = x;
  var b = y;
  while (!b.equals(Expression.ZERO)) {
    var r = a.remainder(b);
    if (!(abs(r.norm()) <= abs(b.norm()))) {
      throw new TypeError("norm");
    }
    a = b;
    b = r;
  }
  return a;
};

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

/*
// http://oeis.org/wiki/Quadratic_integer_rings#Quadratic_integer_ring_with_discriminant_2
var expected = {
  "2": {
    "1": "",
    "2": "0+1sqrt(2), 0+1sqrt(2)",
    "3": "3+0sqrt(2)",
    "4": "0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2)",
    "5": "5+0sqrt(2)",
    "6": "0+1sqrt(2), 0+1sqrt(2), 3+0sqrt(2)",
    "7": "3+1sqrt(2), 3+-1sqrt(2)",
    "8": "0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2)",
    "9": "3+0sqrt(2), 3+0sqrt(2)",
    "10": "0+1sqrt(2), 0+1sqrt(2), 5+0sqrt(2)",
    "11": "11+0sqrt(2)",
    "12": "0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2), 3+0sqrt(2)",
    "13": "13+0sqrt(2)",
    "14": "0+1sqrt(2), 0+1sqrt(2), 3+1sqrt(2), 3+-1sqrt(2)",
    "15": "3+0sqrt(2), 5+0sqrt(2)",
    "16": "0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2)",
    "17": "5+2sqrt(2), 5+-2sqrt(2)",
    "18": "0+1sqrt(2), 0+1sqrt(2), 3+0sqrt(2), 3+0sqrt(2)",
    "19": "19+0sqrt(2)",
    "20": "0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2), 0+1sqrt(2), 5+0sqrt(2)"
  }
};

var allFactors = function* (qi) {
  while (!qi.isUnit()) {
    var pf = qi.primeFactor();
    yield pf.toString();
    qi = qi.truncatingDivide(pf);
  }
  if (!/^1\+0sqrt\(\d+\)$/.test(qi.toString())) {
    yield qi.toString();
  }
};
for (var D of [2, 3, 5, 6, 7, 11, 13, 17, 19, 21, 29, 33, 37, 41, 57, 73]) {
  expected[D] = expected[D] || {};
  if (D === 3) debugger;
  for (var i = 1; i <= 20; i += 1) {
    expected[D][i] = Array.from(allFactors(new QuadraticInteger(i, 0, D))).join(', ');
    //console.log(i, Array.from(allFactors(new QuadraticInteger(i, 0, D))));
  }
  console.log(D, JSON.stringify(expected[D], null, 2));
}
debugger;
throw new Error();
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
          ff = new SquareRoot(Integer.fromNumber(q.D));
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



// +1, -1, +i, -i
// a+bi

// a === 0, i*(a+bi)
// a < 0, -(a+bi)
// b < 0, i*(a+bi)

// a > 0, b > 0


/*
  for (var a = 1; a * a <= n; a += 1) {
    for (var b = 0; b * b <= n - a * a; b += 1) {
      if (norm(a, b) > 1n && hasDivisor(r, i, a, b)) {
        return [a, b];
      }
      if (norm(a, -b) > 1n && hasDivisor(r, i, a, -b)) {
        return [a, -b];
      }
    }
  }
  return [r, i];
*/

function isPrime(n) {
  return primeFactor(n) === n;
}

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
              throw new TypeError();
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


globalThis.QuadraticInteger = QuadraticInteger;



// ExpressionParser.parse('((17^0.5+7)**3)**(1/3)') + ''


//TODO:
/*

*/

function AlmostQuadraticInteger(k, qi) { // k * qi
  if (!(k instanceof Expression.Integer) || k.compareTo(Expression.TWO) < 0) {
    throw new TypeError();
  }
  if (qi.a === 0n) {
  //  throw new TypeError();
  }
  if (qi.b === 0n) {
    throw new TypeError();
  }
  if (ngcd(qi.a, qi.b) !== 1n) {
    throw new TypeError();
  }
  this.k = k;
  this.qi = qi;
}

//AlmostQuadraticInteger.prototype = Object.create(null);

AlmostQuadraticInteger.prototype.toString = function () {
  return '(' + this.k + ')' + '*' + '(' + this.qi + ')';
};

AlmostQuadraticInteger.prototype.isValid = function () {
  return this.qi.isValid();
};

AlmostQuadraticInteger.prototype.isPositive = function () {
  return this.qi.isPositive();
};

AlmostQuadraticInteger.prototype.norm = function () {
  return BigInt(this.k.abs().toBigInt()) * this.qi.norm();
};

AlmostQuadraticInteger.prototype.equals = function (y) {
  return this.qi.toExpression().multiply(this.k).equals(y.toExpression());
};
AlmostQuadraticInteger.prototype.primeFactor = function () {
  //if (this.k.toBigInt() % this.qi.D == 0) {
    //TODO: ?
    //return new QuadraticInteger(this.k.toBigInt(), 0n, this.qi.D).multiply(this.qi).primeFactor();
  //}
  return !(this.k.toNumber() === this.qi.D) ? this.k.primeFactor() : new QuadraticInteger(0n, 1n, this.qi.D);
};
AlmostQuadraticInteger.prototype.isUnit = function () {
  return false;//!
};
AlmostQuadraticInteger.prototype.isDivisibleBy = function (y) {
  var g = y.toExpression().pow(Expression.TWO).gcd(this.k.pow(Expression.TWO)).squareRoot();
  var nk = this.k.divide(g);

  if (g instanceof Expression.SquareRoot && g.a.toNumber() !== y.D) {
    return this.qi.multiply(new QuadraticInteger(this.k.toBigInt(), 0n, this.qi.D)).isDivisibleBy(y); //?TODO: fix
  }
  //!new 2020-09-27
  if (g instanceof Expression.Multiplication) {
    return this.qi.multiply(new QuadraticInteger(this.k.toBigInt(), 0n, this.qi.D)).isDivisibleBy(y); //?TODO: fix
  }
  //!
  return this.qi.isDivisibleBy(y.truncatingDivide(g instanceof Expression.SquareRoot ? new QuadraticInteger(0n, 1n, g.a.toBigInt()) : g));
};

AlmostQuadraticInteger.prototype.truncatingDivide = function (y) {
  //TODO: Fix
  var g = y.toExpression().pow(Expression.TWO).gcd(this.k.pow(Expression.TWO)).squareRoot();
  var nk = this.k.divide(g);
  if (g instanceof Expression.SquareRoot && g.a.toNumber() !== y.D) {
    return this.qi.multiply(new QuadraticInteger(this.k.toBigInt(), 0n, this.qi.D)).truncatingDivide(y); //?TODO: fix
  }
  var nq = this.qi.truncatingDivide(y.truncatingDivide(g instanceof Expression.SquareRoot ? new QuadraticInteger(0n, 1n, g.a.toBigInt()) : g));
  if (nq.b === 0n) {
    // || nq.a == 0 - 2*sqrt(2)
    return nk.multiply(nq.toExpression());
  }
  if (nq.a === 0n) {
    return nq.multiply(new QuadraticInteger(nk.toBigInt(), 0, nq.D));//?
  }
  if (nk instanceof Expression.Multiplication || nk instanceof Expression.SquareRoot) {
    //TODO: fix
    nq = nq.multiply(new QuadraticInteger(0, 1, nq.D));
    nk = nk.divide(Expression.Integer.fromNumber(nq.D).squareRoot());
    var g = ngcd(nq.a, nq.b);
    nq = nq.truncatingDivide(new QuadraticInteger(g, 0, nq.D));
    nk = nk.multiply(Expression.Integer.fromBigInt(g));
  }
  return nk.equals(Expression.ONE) ? nq : new AlmostQuadraticInteger(nk, nq);
};
AlmostQuadraticInteger.prototype.remainder = function (y) {
  if (y instanceof QuadraticInteger) {
    var x = this.toExpression().divide(y.toExpression());
    var n = x.getNumerator();
    var d = x.getDenominator();
    if (d.equals(Expression.ONE)) {
      return Expression.ZERO;
    }
    if (d instanceof Expression.Integer) {
      d = new QuadraticInteger(d.toBigInt(), 0, this.qi.D);
    } else {
      d = toQuadraticInteger(d);
    }
    if (d == null) debugger;
    if (n instanceof Expression.Multiplication && n.a instanceof Expression.Integer && n.b instanceof Expression.SquareRoot) {
      return new QuadraticInteger(Expression.ZERO.toBigInt(), n.a.toBigInt(), n.b.a.toBigInt()).remainder(this);
    }
    if (!this.equals(toQuadraticInteger(n))) {
      return toQuadraticInteger(n).remainder(d);
    }
  }

  var g = y.toExpression().gcd(this.k);
  //TODO: ?
  var remainder = this.qi.remainder(y.truncatingDivide(g));
  if (remainder.equals(Expression.ZERO)) {
    return remainder;
  }
  if (remainder instanceof AlmostQuadraticInteger) {
    return new AlmostQuadraticInteger(this.k.multiply(remainder.k), remainder.qi);
  }
  if (remainder instanceof Expression.Integer) {
    return this.k.multiply(remainder);//TODO: fix
  }
  var t = Expression.Integer.fromBigInt(ngcd(remainder.a, remainder.b));
  var nk = this.k.multiply(t);
  remainder = remainder.truncatingDivide(t);
  if (remainder.b === 0n) {
    return nk.divide(g);
  }
  return nk.equals(g) ? remainder : new AlmostQuadraticInteger(nk.divide(g), remainder);
};
AlmostQuadraticInteger.prototype.toExpression = function (y) {
  return this.qi.toExpression().multiply(this.k);
};

AlmostQuadraticInteger.prototype.isDivisibleByInteger = function (x) {
  return x.truncatingDivide(this) != null;
};

AlmostQuadraticInteger.prototype.truncatingDivideInteger = function (x) {
  return new QuadraticInteger(x.toBigInt(), Expression.ZERO.toBigInt(), this.qi.D).truncatingDivide(this);
};

AlmostQuadraticInteger.prototype.remainderInteger = function (x) {
  return new QuadraticInteger(x.toBigInt(), Expression.ZERO.toBigInt(), this.qi.D).remainder(this);
};


// new QuadraticInteger(7, 1, 17).remainder(new QuadraticInteger(3, 1, 17))



// http://oeis.org/wiki/Quadratic_integer_rings#Quadratic_integer_ring_with_discriminant_13

/*
//TODO:

function checkFactorization(i) {
  var results = [];
  var x = i;
  while (Math.abs(x.norm()) > 1) {
    debugger;
    var p = x.primeFactor();
    results.push(p);
    x = x.truncatingDivide(p);
  }
  console.log(i + '=' + results.map(x => '(' + x + ')').join(''));
}


for (var i = 3; i <= 20; i += 1) {
  var q = new QuadraticInteger(i, 0, 13);
  checkFactorization(q);
}
*/
