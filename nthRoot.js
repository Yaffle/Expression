import BigInteger from './BigInteger.js';

// floor(log(n) / log(b)), n >= 1, base >= 2
// https://programmingpraxis.com/contents/standard-prelude/
function ilogb(n, base) {
  var i = BigInteger.BigInt(0);
  while (!BigInteger.lessThan(n, BigInteger.exponentiate(base, BigInteger.exponentiate(BigInteger.BigInt(2), i)))) {
    i = BigInteger.add(i, BigInteger.BigInt(1));
  }
  var e = BigInteger.BigInt(0);
  var t = BigInteger.BigInt(1);
  while (!BigInteger.lessThan(i, BigInteger.BigInt(0))) {
    var b = BigInteger.exponentiate(BigInteger.BigInt(2), i);
    if (!BigInteger.lessThan(n, BigInteger.multiply(t, BigInteger.exponentiate(base, b)))) {
      t = BigInteger.multiply(t, BigInteger.exponentiate(base, b));
      e = BigInteger.add(e, b);
    }
    i = BigInteger.subtract(i, BigInteger.BigInt(1));
  }
  return e;
}

function ilog2(n) {
  return ilogb(n, BigInteger.BigInt(2));
}

// floor(S**(1/n)), S >= 1, n >= 2
// https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Babylonian_method
// https://stackoverflow.com/a/15979957/839199y
function nthRoot(S, n) {
  if (!BigInteger.lessThan(BigInteger.BigInt(2), n)) {
    var t = BigInteger.toNumber(S);
    if (t < (9007199254740991 + 1) / 2) {
      return BigInteger.BigInt(Math.floor(Math.sqrt(t + 0.5)));
    }
    if (t < (9007199254740991 + 1) / 2 * (9007199254740991 + 1) / 2) {
      var y = BigInteger.BigInt(Math.floor(Math.sqrt(t) + 0.5));
      if (BigInteger.lessThan(S, BigInteger.multiply(y, y))) {
        y = BigInteger.subtract(y, BigInteger.BigInt(1));
      }
      return y;
    }
  }
  var e = ilog2(S);
  if (BigInteger.lessThan(e, n)) {
    return BigInteger.BigInt(1);
  }
  var f = BigInteger.divide(BigInteger.add(e, n), BigInteger.multiply(BigInteger.BigInt(2), n));
  var x = BigInteger.multiply(BigInteger.add(nthRoot(BigInteger.divide(S, BigInteger.exponentiate(BigInteger.BigInt(2), BigInteger.multiply(f, n))), n), BigInteger.BigInt(1)), BigInteger.exponentiate(BigInteger.BigInt(2), f));
  var xprev = BigInteger.add(x, BigInteger.BigInt(1));
  while (BigInteger.lessThan(x, xprev)) {
    xprev = x;
    x = BigInteger.divide(BigInteger.add(BigInteger.multiply(x, BigInteger.subtract(n, BigInteger.BigInt(1))), BigInteger.divide(S, BigInteger.exponentiate(x, BigInteger.subtract(n, BigInteger.BigInt(1))))), n);
  }
  return xprev;
}

export default nthRoot;
