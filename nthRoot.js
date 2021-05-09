import primeFactor from './primeFactor.js';

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
  if (n < 1) {
    throw new RangeError();
  }
  const s = Number(S);
  if (s < 0) {
    if (n % 2 === 0) {
      throw new RangeError();
    }
    return -nthRoot(-S, n);
  }
  if (n === 1) {
    return S;
  }
  if (s < (Number.MAX_SAFE_INTEGER + 1) / (8 * n)) {
    // var i = 1; while (Math.floor(Math.sqrt(i**2 - 1 + 0.5)) < i) { i++; } console.log(i**2);
    // var i = 1; while (Math.floor(Math.cbrt(i**3 - 1 + 0.5)) < i && Math.cbrt(i**3 + 0.5) >= i) { i++; } console.log(i**3);
    // for (var n = 2; n <= 53; n++) { var i = 1; while (Math.floor(nthRootSmall(i**n - 1 + 0.5, n)) < i) { i++; } console.log(i**n); }
    const g = n === 2 ? Math.floor(Math.sqrt(s + 0.5)) : (n === 3 ? Math.floor(Math.cbrt(s + 0.5)) : Math.floor(nthRootSmall(s + 0.5, n)));
    return BigInt(g);
  }
  let g = (n === 2 ? Math.sqrt(s) : (n === 3 ? Math.cbrt(s) : nthRootSmall(s, n)));
  if (g < (Number.MAX_SAFE_INTEGER + 1) / (8)) {
    if (n === 3) {
      g = g + (s / (g * g) - g) / 3;
    }
    if ((n === 2 || n === 3) && Math.floor(g) !== g) {
      return BigInt(Math.floor(g));
    }
    let y = BigInt(Math.floor(g + 0.5));
    if (S < y**BigInt(n)) {
      y -= 1n;
    }
    return y;
  }
  const e = primeFactor._bitLength(S);
  if (e <= n) {
    return 1n;
  }
  const f = Math.floor((e + n) / (2 * n));
  let x = (nthRoot(S >> BigInt(f * n), n) + 1n) << BigInt(f);
  let xprev = x + 1n;
  while (x < xprev) {
    xprev = x;
    x = (BigInt(n - 1) * x + S / x**BigInt(n - 1)) / BigInt(n);
  }
  return xprev;
}

export default nthRoot;
