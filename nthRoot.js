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
function nthRoot(S, sn) {
  if (sn < 2 || sn % 2 === 0 && BigInteger.lessThan(S, 0)) {
    throw new RangeError();
  }
  if (sn === 2) {
    // e === 1/2 * 2**-52
    var t = BigInteger.toNumber(S);
    // 1/(2*e)
    if (t < (9007199254740991 + 1) / 2) { // var i = 1; while (Math.floor(Math.sqrt(i**2 - 1 + 0.5)) < i) { i++; } console.log(i**2);
      return BigInteger.BigInt(Math.floor(Math.sqrt(t + 0.5)));
    }
    // (1/(2*e))**2
    if (t < (9007199254740991 + 1) / 2 * (9007199254740991 + 1) / 2) {
      var y = BigInteger.BigInt(Math.floor(Math.sqrt(t) + 0.5));
      if (BigInteger.lessThan(S, BigInteger.exponentiate(y, BigInteger.BigInt(2)))) {
        y = BigInteger.subtract(y, BigInteger.BigInt(1));
      }
      return y;
    }
  }
  //TODO: fix
  if (sn === 3) {
    // e = 2/3 * 2**-52 in some browsers ...
    var t = BigInteger.toNumber(S);
    // 1/(3*e)
    if (t < (9007199254740991 + 1) / 32) { // var i = 1; while (Math.floor(Math.cbrt(i**3 - 1 + 0.5)) < i) { i++; } console.log(i**3);
      return BigInteger.BigInt(Math.floor(Math.cbrt(t + 0.5)));
    }
    // (1/(2*e))**3
    if (t < (9007199254740991 + 1) / 32 * (9007199254740991 + 1) / 32 * (9007199254740991 + 1) / 32) {
      var y = BigInteger.BigInt(Math.floor(Math.cbrt(t) + 0.5));
      if (BigInteger.lessThan(S, BigInteger.exponentiate(y, BigInteger.BigInt(3)))) {
        y = BigInteger.subtract(y, BigInteger.BigInt(1));
      }
      return y;
    }
  }
  if (sn > 3) {
    var t = BigInteger.toNumber(S);
    if (t < (9007199254740991 + 1) / 128) {
      return BigInteger.BigInt(Math.floor(Math.exp(Math.log(t + 0.5) / sn)));
    }
  }
  var n = BigInteger.BigInt(sn);
  var e = ilog2(S);
  if (BigInteger.lessThan(e, n)) {
    return BigInteger.BigInt(1);
  }
  var f = BigInteger.divide(BigInteger.add(e, n), BigInteger.multiply(BigInteger.BigInt(2), n));
  var x = BigInteger.multiply(BigInteger.add(nthRoot(BigInteger.divide(S, BigInteger.exponentiate(BigInteger.BigInt(2), BigInteger.multiply(f, n))), sn), BigInteger.BigInt(1)), BigInteger.exponentiate(BigInteger.BigInt(2), f));
  var xprev = BigInteger.add(x, BigInteger.BigInt(1));
  while (BigInteger.lessThan(x, xprev)) {
    xprev = x;
    x = BigInteger.divide(BigInteger.add(BigInteger.multiply(x, BigInteger.subtract(n, BigInteger.BigInt(1))), BigInteger.divide(S, BigInteger.exponentiate(x, BigInteger.subtract(n, BigInteger.BigInt(1))))), n);
  }
  return xprev;
}

export default nthRoot;
