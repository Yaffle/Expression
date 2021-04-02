  const cache = new Array(33);
  function BigIntWrapper() {
  }
  //const BigIntWrapper = typeof BigInt !== 'undefined' ? BigInt : function () {};
  BigIntWrapper.BigInt = function (x) {
    if (typeof x === "number" && Math.abs(x) > Number.MAX_SAFE_INTEGER) {
      throw new RangeError();
    }
    if (typeof x === "number" && Math.abs(x) <= 16) {
      var value = cache[x + 16];
      if (value == null) {
        value = BigInt(x);
        cache[x + 16] = value;
      }
      return value;
    }
    if (typeof x === "bigint") {
      return x;
    }
    return BigInt(x);
  };
  BigIntWrapper.toNumber = function (bigint) {
    return Number(bigint);
  };
  BigIntWrapper.add = function (a, b) {
    return a + b;
  };
  BigIntWrapper.subtract = function (a, b) {
    return a - b;
  };
  BigIntWrapper.multiply = function (a, b) {
    return a * b;
  };
  BigIntWrapper.divide = function (a, b) {
    return a / b;
  };
  BigIntWrapper.remainder = function (a, b) {
    return a % b;
  };
  BigIntWrapper.unaryMinus = function (a) {
    return -a;
  };
  BigIntWrapper.equal = function (a, b) {
    return a === b;
  };
  BigIntWrapper.lessThan = function (a, b) {
    return a < b;
  };
  BigIntWrapper.greaterThan = function (a, b) {
    return a > b;
  };
  BigIntWrapper.notEqual = function (a, b) {
    return a !== b;
  };
  BigIntWrapper.lessThanOrEqual = function (a, b) {
    return a <= b;
  };
  BigIntWrapper.greaterThanOrEqual = function (a, b) {
    return a >= b;
  };
  BigIntWrapper.exponentiate = function (a, b) { // a**b
    if (typeof a !== "bigint" || typeof b !== "bigint") {
      throw new TypeError();
    }
    var n = Number(b);
    if (n < 0) {
      throw new RangeError();
    }
    if (n > 9007199254740991) {
      var y = Number(a);
      if (y === 0 || y === -1 || y === +1) {
        return y === -1 && Number(b % BigInt(2)) === 0 ? -a : a;
      }
      throw new RangeError();
    }
    if (a === BigInt(2)) {
      return BigInt(1) << b;
    }
    if (n === 0) {
      return BigInt(1);
    }
    var x = a;
    while (n % 2 === 0) {
      n = Math.floor(n / 2);
      x *= x;
    }
    var accumulator = x;
    n -= 1;
    if (n >= 2) {
      while (n >= 2) {
        var t = Math.floor(n / 2);
        if (t * 2 !== n) {
          accumulator *= x;
        }
        n = t;
        x *= x;
      }
      accumulator *= x;
    }
    return accumulator;
  };
  BigIntWrapper.signedRightShift = function (a, n) {
    return a >> n;
  };
  BigIntWrapper.leftShift = function (a, n) {
    return a << n;
  };

  var supportsBigInt = typeof BigInt !== "undefined" && BigInt(9007199254740991) + BigInt(2) - BigInt(2) === BigInt(9007199254740991);
  //supportsBigInt = false;//!!!
  var Internal = supportsBigInt ? BigIntWrapper : JSBI;
  if (supportsBigInt) {
    globalThis.JSBI = BigIntWrapper;//!!!
  }

  // noinline
  var n = function (f) {
    return function (x, y) {
      return f(x, y);
    };
  };

  var toNumber = n(function (a) {
    return Internal.toNumber(a);
  });
  var valueOf = function (x) {
    if (typeof x === "number") {
      return Internal.BigInt(x);
    }
    return x;
  };
  var toResult = function (x) {
    var value = Internal.toNumber(x);
    if (value >= -9007199254740991 && value <= +9007199254740991) {
      return value;
    }
    return x;
  };
  var add = n(function (x, y) {
    if (typeof x === "number" && x === 0) {
      return y;
    }
    if (typeof y === "number" && y === 0) {
      return x;
    }
    var a = valueOf(x);
    var b = valueOf(y);
    return toResult(Internal.add(a, b));
  });
  var subtract = n(function (x, y) {
    if (typeof x === "number" && x === 0) {
      return unaryMinus(y);
    }
    // quite good optimization for comparision of big integers
    if (typeof y === "number" && y === 0) {
      return x;
    }
    var a = valueOf(x);
    var b = valueOf(y);
    return toResult(Internal.subtract(a, b));
  });
  var multiply = n(function (x, y) {
    if (x === y) {
      var c = valueOf(x);
      return Internal.multiply(c, c);
    }
    if (typeof x === "number" && x === 0 || typeof y === "number" && y === 0) {
      return 0;
    }
    if (typeof x === "number" && x === 1) {
      return y;
    }
    if (typeof x === "number" && x === -1) {
      return Internal.unaryMinus(y);
    }
    if (typeof y === "number" && y === 1) {
      return x;
    }
    if (typeof y === "number" && y === -1) {
      return Internal.unaryMinus(x);
    }
    var a = valueOf(x);
    var b = valueOf(y);
    return Internal.multiply(a, b);
  });
  var divide = n(function (x, y) {
    if (typeof y === "number" && y === 1) {
      return x;
    }
    var a = valueOf(x);
    var b = valueOf(y);
    return toResult(Internal.divide(a, b));
  });
  var remainder = n(function (x, y) {
    if (typeof x === "number") {
      return x;
    }
    var a = valueOf(x);
    var b = valueOf(y);
    return toResult(Internal.remainder(a, b));
  });
  var exponentiate = n(function (x, y) {
    if (typeof y === "number") {
      if (y === 0) {
        return 1;
      }
      if (y === 1) {
        return x;
      }
      if (y === 2) {
        var c = valueOf(x);
        return Internal.multiply(c, c);
      }
      if (typeof x === "number" && Math.abs(x) > 2 && y >= 0) {
        if (y > 42 && x % 2 === 0) {//TODO: ?
          return multiply(exponentiate(2, y), exponentiate(x / 2, y));
        }
        var k = Math.floor(Math.log(9007199254740991) / Math.log(Math.abs(x) + 0.5));
        if (k >= 2) {
          return multiply(Math.pow(x, y % k), exponentiate(Math.pow(x, k), Math.floor(y / k)));
        }
      }
    }
    var a = valueOf(x);
    var b = valueOf(y);
    return Internal.exponentiate(a, b);
  });
  var unaryMinus = n(function (x) {
    var a = valueOf(x);
    return Internal.unaryMinus(a);
  });
  var equal = n(function (x, y) {
    if (typeof x === "number") {
      return false;
    }
    if (typeof y === "number") {
      return false;
    }
    return Internal.equal(x, y);
  });
  var lessThan = n(function (x, y) {
    if (typeof x === "number") {
      return x < Internal.toNumber(y);
    }
    if (typeof y === "number") {
      return Internal.toNumber(x) < y;
    }
    return Internal.lessThan(x, y);
  });
  var greaterThan = n(function (x, y) {
    if (typeof x === "number") {
      return x > Internal.toNumber(y);
    }
    if (typeof y === "number") {
      return Internal.toNumber(x) > y;
    }
    return Internal.greaterThan(x, y);
  });

  function BigInteger() {
  }

  // Conversion from String:
  // Conversion from Number:
  BigInteger.BigInt = function (x) {
    if (typeof x === "number") {
      return x;
    }
    if (typeof x === "string") {
      var value = 0 + Number(x);
      if (value >= -9007199254740991 && value <= +9007199254740991) {
        return value;
      }
    }
    if (typeof x === "bigint") {
      return toResult(x);
    }
    return toResult(Internal.BigInt(x));
  };
  // Conversion to Number:
  BigInteger.toNumber = function (x) {
    if (typeof x === "number") {
      return x;
    }
    return toNumber(x);
  };

  // Arithmetic:
  BigInteger.add = function (x, y) {
    if (typeof x === "string" || typeof y === "string") {
      return x + y;
    }
    if (typeof x === "number" && typeof y === "number") {
      var value = x + y;
      if (value >= -9007199254740991 && value <= +9007199254740991) {
        return value;
      }
    }
    return add(x, y);
  };
  BigInteger.subtract = function (x, y) {
    if (typeof x === "number" && typeof y === "number") {
      var value = x - y;
      if (value >= -9007199254740991 && value <= +9007199254740991) {
        return value;
      }
    }
    return subtract(x, y);
  };
  BigInteger.multiply = function (x, y) {
    if (typeof x === "number" && typeof y === "number") {
      var value = 0 + x * y;
      if (value >= -9007199254740991 && value <= +9007199254740991) {
        return value;
      }
    }
    return multiply(x, y);
  };
  BigInteger.divide = function (x, y) {
    if (typeof x === "number" && typeof y === "number") {
      if (y !== 0) {
        return x === 0 ? 0 : (x > 0 && y > 0) || (x < 0 && y < 0) ? 0 + Math.floor(x / y) : 0 - Math.floor((0 - x) / y);
      }
    }
    return divide(x, y);
  };
  BigInteger.remainder = function (x, y) {
    if (typeof x === "number" && typeof y === "number") {
      if (y !== 0) {
        return 0 + x % y;
      }
    }
    return remainder(x, y);
  };
  BigInteger.unaryMinus = function (x) {
    if (typeof x === "number") {
      return 0 - x;
    }
    return unaryMinus(x);
  };

  // Comparison:
  BigInteger.equal = function (x, y) {
    if (typeof x === "number" && typeof y === "number") {
      return x === y;
    }
    return equal(x, y);
  };
  BigInteger.lessThan = function (x, y) {
    if (typeof x === "number" && typeof y === "number") {
      return x < y;
    }
    return lessThan(x, y);
  };
  BigInteger.greaterThan = function (x, y) {
    if (typeof x === "number" && typeof y === "number") {
      return x > y;
    }
    return greaterThan(x, y);
  };
  BigInteger.notEqual = function (x, y) {
    return !BigInteger.equal(x, y);
  };
  BigInteger.lessThanOrEqual = function (x, y) {
    return !BigInteger.greaterThan(x, y);
  };
  BigInteger.greaterThanOrEqual = function (x, y) {
    return !BigInteger.lessThan(x, y);
  };

  BigInteger.exponentiate = function (x, y) {
    if (typeof x === "number" && typeof y === "number" && y >= 0 && y < 53) { // 53 === log2(9007199254740991 + 1)
      var value = 0 + Math.pow(x, y);
      if (value >= -9007199254740991 && value <= 9007199254740991) {
        return value;
      }
    }
    return exponentiate(x, y);
  };
  BigInteger.signedRightShift = function (x, n) {
    if (typeof x === "bigint") {
      return toResult(x >> BigInt(n));
    }
    if (n < 0) {
      return BigInteger.leftShift(x, BigInteger.unaryMinus(n));
    }
    if (x < 0) {
      var q = BigInteger.divide(x, BigInteger.exponentiate(2, n));
      if (BigInteger.subtract(x, BigInteger.multiply(q, BigInteger.exponentiate(2, n))) < 0) {
        q = BigInteger.subtract(q, BigInteger.BigInt(1));
      }
      return q;
    }
    return BigInteger.divide(x, BigInteger.exponentiate(2, n));
  };
  BigInteger.leftShift = function (x, n) {
    if (typeof x === "bigint") {
      return x << BigInt(n);
    }
    if (n < 0) {
      return BigInteger.signedRightShift(x, BigInteger.unaryMinus(n));
    }
    return BigInteger.multiply(x, BigInteger.exponentiate(2, n));
  };

  const SmallBigInt = BigInteger;
  export default SmallBigInt;

