// floor(S**(1/n)), S >= 1, n >= 2
// https://en.wikipedia.org/wiki/Methods_of_computing_square_roots#Babylonian_method
// https://stackoverflow.com/a/15979957/839199y
function nthRoot(S, n) {
  n = Number(n);
  if (Math.floor(n) !== n || n < 1 || n > Number.MAX_SAFE_INTEGER) {
    throw new RangeError();
  }
  if (n === 1) {
    return S;
  }
  const s = Math.floor(Number(S));
  if (s === 0) {
    return s;
  }
  if (s < 0) {
    if (n % 2 === 0) {
      throw new RangeError();
    }
    return -nthRoot(-S, n);
  }
  const B = Number.MAX_SAFE_INTEGER + 1;
  const E = Math.floor(B / Math.pow(2, n === 2 ? 1 : (n === 3 ? 2 : (1 + Math.ceil(Math.log2(Math.log2(B)))))));
  if (s < E) {
    //var test = function (n, f) { var i = 1; while (f(i**n - 1) === i - 1 && f(i**n) === i) { i += 1; } var a = i**n; while (f(a - 1) === i) { a -= 2**25; } console.log(n, Math.log2(a), a); };
    //test(2, a => Math.floor(Math.sqrt(a + 0.5)));
    //test(3, a => Math.floor(Math.cbrt(a + 0.5)));
    //for (var n = 2; n <= 53; n++) { test(n, a => Math.floor(Math.exp(Math.log(a + 0.5) / n))); }
    const g = n === 2 ? Math.floor(Math.sqrt(s + 0.5)) : (n === 3 ? Math.floor(Math.cbrt(s + 0.5)) : Math.floor(Math.exp(Math.log(s + 0.5) / n)));
    return g;
  }
  const g = (n === 2 ? Math.sqrt(s) : (n === 3 ? Math.cbrt(s) : Math.exp(Math.log(s) / n)));
  if (g < E) {
    if (Math.floor(g - g / E) === Math.floor(g + g / E)) {
      return Math.floor(g);
    }
    let y = Math.floor(g + 0.5);
    if (BigInt(S) < BigInt(y)**BigInt(n)) {
      y -= 1;
    }
    return y;
  }
  const size = BigInt(S).toString(16).length * 4; // TODO: bitLength(BigInt(S))
  if (size <= n) {
    return 1;
  }
  const half = Math.floor((Math.floor(size / n) + 1) / 2);
  let x = (BigInt(nthRoot(BigInt(S) >> BigInt(half * n), n)) + 1n) << BigInt(half);
  let xprev = x + 1n;
  while (x < xprev) {
    xprev = x;
    x = (BigInt(n - 1) * x + BigInt(S) / x**BigInt(n - 1)) / BigInt(n);
  }
  return Number(xprev) <= Number.MAX_SAFE_INTEGER ? Number(xprev) : xprev;
}

export default nthRoot;
