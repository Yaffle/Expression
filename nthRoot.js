
// https://en.wikipedia.org/wiki/Find_first_set#CTZ
function countTrailingZeros(x, base) {
  if (x === 0n) {
    throw new TypeError();
  }
  var k = 1n;
  while (x % base**k === 0n) {
    k *= 2n;
  }
  var n = 0n;
  for (var i = k / 2n; i >= 1n; i /= 2n) {
    if (x % base**i === 0n) {
      n += i;
      x /= base**i;
    }
  }
  return Number(n);
}

function abs(a) {
  return a < 0n ? -a : a;
}

// https://en.wikipedia.org/wiki/Euclidean_algorithm#Method_of_least_absolute_remainders
function gcd(x, y) {
  var a = abs(x);
  var b = abs(y);
  while (b !== 0n) {
    var r1 = a % b;
    var r2 = b - r1;
    var r = r1 < r2 ? r1 : r2;
    a = b;
    b = r;
  }
  return a;
}

function naturalLogarithm(n) {
  var number = Number(n);
  if (number < 1 / 0) {
    return Math.log(number);
  }
  // https://github.com/tc39/proposal-bigint/issues/205
  var s = n.toString(16);
  var p = Math.floor(Math.log((Number.MAX_SAFE_INTEGER + 1) / 32 + 0.5) / Math.log(2));
  var l = Math.floor(p / 4);
  return Math.log(Number('0x' + s.slice(0, l)) / Math.pow(2, 4 * l)) + 4 * Math.log(2) * s.length;
}

function nthRootSmall(A, n) {
  var x = Math.exp(Math.log(A) / n);
  // https://en.wikipedia.org/wiki/Nth_root_algorithm
  x = x + (A / Math.pow(x, n - 1) - x) / n;
  return x;
}

// floor(S**(1/n)), S >= 1, n >= 2
// https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Babylonian_method
// https://stackoverflow.com/a/15979957/839199y
function nthRoot(S, n) {
  if (n < 1 || n % 2 === 0 && S < 0n) {
    throw new RangeError();
  }
  if (n === 1) {
    return S;
  }
  if (true) {
    var error = 4 / (9007199254740991 + 1);
    var s = Number(S);
    if (s < 1 / (2 * n * error)) {
      // var i = 1; while (Math.floor(Math.sqrt(i**2 - 1 + 0.5)) < i) { i++; } console.log(i**2);
      // var i = 1; while (Math.floor(Math.cbrt(i**3 - 1 + 0.5)) < i) { i++; } console.log(i**3);
      // for (var n = 2; n <= 53; n++) { var i = 1; while (Math.floor(nthRootSmall(i**n - 1 + 0.5, n)) < i) { i++; } console.log(i**n / (2**50 / n)); }
      return BigInt(n === 2 ? Math.floor(Math.sqrt(s + 0.5)) : (n === 3 ? Math.floor(Math.cbrt(s + 0.5)) : Math.floor(nthRootSmall(s + 0.5, n))));
    }
    if (s < Math.pow(1 / (2 * error), n)) {
      var y = BigInt(Math.floor((n === 2 ? Math.sqrt(s) : (n === 3 ? Math.cbrt(s) : Math.floor(nthRootSmall(s, n)))) + 0.5));
      if (S < y**BigInt(n)) {
        y -= 1n;
      }
      return y;
    }
  }
  var N = BigInt(n);
  var e = naturalLogarithm(S) / Math.log(2);
  if (e - 1 < n) {
    return 1n;
  }
  var f = BigInt(Math.floor((e + n) / (2 * n)));
  var x = (nthRoot(S / 2n**(f * N), n) + 1n) * 2n**f;
  var xprev = x + 1n;
  while (x < xprev) {
    xprev = x;
    x = (x * (N - 1n) + S / x**(N - 1n)) / N;
  }
  return xprev;
}

nthRoot.countTrailingZeros = countTrailingZeros;//TODO:?
nthRoot.gcd = gcd;
nthRoot.naturalLogarithm = naturalLogarithm;

export default nthRoot;
