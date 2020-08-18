import BigInteger from './BigInteger.js';

// https://en.wikipedia.org/wiki/Find_first_set#CTZ
function countTrailingZeros(x, base) {
  if (BigInteger.equal(x, BigInteger.BigInt(0))) {
    throw new TypeError();
  }

  var k = BigInteger.BigInt(1);

  while (BigInteger.equal(BigInteger.remainder(x, BigInteger.exponentiate(base, k)), BigInteger.BigInt(0))) {
    k = BigInteger.multiply(k, BigInteger.BigInt(2));
  }

  var n = BigInteger.BigInt(0);

  for (var i = BigInteger.divide(k, BigInteger.BigInt(2)); BigInteger.greaterThanOrEqual(i, BigInteger.BigInt(1)); i = BigInteger.divide(i, BigInteger.BigInt(2))) {
    if (BigInteger.equal(BigInteger.remainder(x, BigInteger.exponentiate(base, i)), BigInteger.BigInt(0))) {
      n = BigInteger.add(n, i);
      x = BigInteger.divide(x, BigInteger.exponentiate(base, i));
    }
  }

  return BigInteger.toNumber(n);
}

function abs(a) {
  return BigInteger.lessThan(a, BigInteger.BigInt(0)) ? BigInteger.unaryMinus(a) : a;
}

function gcd(x, y) {
  var a = abs(x);
  var b = abs(y);

  while (BigInteger.notEqual(b, BigInteger.BigInt(0))) {
    var r1 = BigInteger.remainder(a, b);
    var r2 = BigInteger.subtract(b, r1);
    var r = BigInteger.lessThan(r1, r2) ? r1 : r2;
    a = b;
    b = r;
  }

  return a;
}

function naturalLogarithm(n) {
  var number = BigInteger.toNumber(n);

  if (number < 1 / 0) {
    return Math.log(number);
  } // https://github.com/tc39/proposal-bigint/issues/205


  var s = n.toString(16);
  var p = Math.floor(Math.log((Number.MAX_SAFE_INTEGER + 1) / 32 + 0.5) / Math.log(2));
  var l = Math.floor(p / 4);
  return Math.log(Number('0x' + s.slice(0, l)) / Math.pow(2, 4 * l)) + 4 * Math.log(2) * s.length;
}

function nthRootSmall(A, n) {
  var x = Math.exp(Math.log(A) / n); // https://en.wikipedia.org/wiki/Nth_root_algorithm

  x = x + (A / Math.pow(x, n - 1) - x) / n;
  return x;
} // floor(S**(1/n)), S >= 1, n >= 2
// https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Babylonian_method
// https://stackoverflow.com/a/15979957/839199y


function nthRoot(S, n) {
  if (n < 1 || n % 2 === 0 && BigInteger.lessThan(S, BigInteger.BigInt(0))) {
    throw new RangeError();
  }

  if (n === 1) {
    return S;
  }

  if (true) {
    var error = 4 / (9007199254740991 + 1);
    var s = BigInteger.toNumber(S);

    if (s < 1 / (2 * n * error)) {
      // var i = 1; while (Math.floor(Math.sqrt(i**2 - 1 + 0.5)) < i) { i++; } console.log(i**2);
      // var i = 1; while (Math.floor(Math.cbrt(i**3 - 1 + 0.5)) < i) { i++; } console.log(i**3);
      // for (var n = 2; n <= 53; n++) { var i = 1; while (Math.floor(nthRootSmall(i**n - 1 + 0.5, n)) < i) { i++; } console.log(i**n / (2**50 / n)); }
      return BigInteger.BigInt(n === 2 ? Math.floor(Math.sqrt(s + 0.5)) : n === 3 ? Math.floor(Math.cbrt(s + 0.5)) : Math.floor(nthRootSmall(s + 0.5, n)));
    }

    if (s < Math.pow(1 / (2 * error), n)) {
      var y = BigInteger.BigInt(Math.floor((n === 2 ? Math.sqrt(s) : n === 3 ? Math.cbrt(s) : Math.floor(nthRootSmall(s, n))) + 0.5));

      if (BigInteger.lessThan(S, BigInteger.exponentiate(y, BigInteger.BigInt(n)))) {
        y = BigInteger.subtract(y, BigInteger.BigInt(1));
      }

      return y;
    }
  }

  var N = BigInteger.BigInt(n);
  var e = naturalLogarithm(S) / Math.log(2);

  if (e - 1 < n) {
    return BigInteger.BigInt(1);
  }

  var f = BigInteger.BigInt(Math.floor((e + n) / (2 * n)));
  var x = BigInteger.multiply(BigInteger.add(nthRoot(BigInteger.divide(S, BigInteger.exponentiate(BigInteger.BigInt(2), BigInteger.multiply(f, N))), n), BigInteger.BigInt(1)), BigInteger.exponentiate(BigInteger.BigInt(2), f));
  var xprev = BigInteger.add(x, BigInteger.BigInt(1));

  while (BigInteger.lessThan(x, xprev)) {
    xprev = x;
    x = BigInteger.divide(BigInteger.add(BigInteger.multiply(x, BigInteger.subtract(N, BigInteger.BigInt(1))), BigInteger.divide(S, BigInteger.exponentiate(x, BigInteger.subtract(N, BigInteger.BigInt(1))))), N);
  }

  return xprev;
}

nthRoot.countTrailingZeros = countTrailingZeros; //TODO:?

nthRoot.gcd = gcd;
nthRoot.naturalLogarithm = naturalLogarithm;
export default nthRoot;

