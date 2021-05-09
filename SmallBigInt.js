import './node_modules/js-big-integer/BigInteger.js';

  const cache = new Array(16 * 2 + 1).fill(undefined);
  function BigIntWrapper() {
  }
  //const BigIntWrapper = typeof BigInt !== 'undefined' ? BigInt : function () {};
  BigIntWrapper.BigInt = function (x) {
    if (typeof x === "number" && Math.abs(x) <= 16) {
      var value = cache[x + 16];
      if (value == undefined) {
        value = BigInt(x);
        cache[x + 16] = value;
      }
      return value;
    }
    if (typeof x === "number" && Math.abs(x) > Number.MAX_SAFE_INTEGER) { // as Chrome 67 does not support this
      throw new RangeError();
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
  if (supportsBigInt) {
    globalThis.JSBI = BigIntWrapper;//!!!
    globalThis.BigInteger._setInternal(BigIntWrapper);
  } else if (globalThis.JSBI == null) {
    globalThis.JSBI = globalThis.BigInteger._getInternal();
  } else {
    globalThis.BigInteger._setInternal(JSBI);
  }

  const SmallBigInt = globalThis.BigInteger;
  export default SmallBigInt;
  //export default BigIntWrapper;

