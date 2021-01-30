
// https://en.wikipedia.org/wiki/Euclidean_algorithm#Method_of_least_absolute_remainders
function numbersGCD(a, b) {
  while (b > 0) {
    const r1 = a - Math.floor(a / b) * b;
    const r2 = b - r1;
    const r = r1 < r2 ? r1 : r2;
    a = b;
    b = r;
  }
  return a;
}

function EuclidsGCD(a, b) {
  while (Number(b) > Number.MAX_SAFE_INTEGER) {
    const r1 = a % b;
    const r2 = b - r1;
    const r = r1 < r2 ? r1 : r2;
    a = b;
    b = r;
  }
  if (Number(b) > 0) {
    if (Number(a) > Number.MAX_SAFE_INTEGER) {
      const r = a % BigInt(b);
      a = b;
      b = r;
    }
    return BigInt(numbersGCD(Number(a), Number(b)));
  }
  return a;
}

// https://github.com/tc39/proposal-bigint/issues/205
// https://github.com/tc39/ecma262/issues/1729
// floor(log2(a)) + 1 if a > 0
function bitLength(a) {
  const s = a.toString(16);
  const c = s.charCodeAt(0) - '0'.charCodeAt(0);
  if (c <= 0) {
    throw new RangeError();
  }
  return (s.length - 1) * 4 + (32 - Math.clz32(Math.min(c, 8)));
}

// https://en.wikipedia.org/wiki/Lehmer%27s_GCD_algorithm
// https://www.imsc.res.in/~kapil/crypto/notes/node11.html
// this implementation is good after ~80 bits (?)
function LehmersGCD(a, b) {
  if (a < b) {
    const tmp = a;
    a = b;
    b = tmp;
  }
  while (Number(b) >= Math.sqrt(Math.pow(Number.MAX_SAFE_INTEGER + 1, 3))) {
    console.assert(a >= b);
    const m = BigInt(Math.max(0, bitLength(a) - Math.floor(Math.log2(Number.MAX_SAFE_INTEGER + 1))));
    let x = Number(a >> m);
    let y = Number(b >> m);
    console.assert(x >= (Number.MAX_SAFE_INTEGER + 1) / 2 && x <= Number.MAX_SAFE_INTEGER);
    console.assert(y >= 0 && y <= Number.MAX_SAFE_INTEGER);
    let A = 1;
    let B = 0;
    let C = 0;
    let D = 1;
    let w1 = 0;
    let w2 = 0;
    while (y + C !== 0 && y + D !== 0 && (w1 = Math.floor((x + A) / (y + C))) === (w2 = Math.floor((x + B) / (y + D)))) {
      const w = w1;
      [A, B, x, C, D, y] = [C, D, y, A - w * C, B - w * D, x - w * y];
    }
    if (B === 0) {
      console.assert(A === 1 && B === 0 && C === 0 && D === 1);
      const r = a % b;
      a = b;
      b = r;
    } else {
      [a, b] = [BigInt(A) * a + BigInt(B) * b, BigInt(C) * a + BigInt(D) * b];
    }
  }
  return EuclidsGCD(a, b);
}

function abs(a) {
  return a < 0n ? -a : a;
}

function ctz(a) {
  // https://en.wikipedia.org/wiki/Find_first_set#Properties_and_relations
  //return bitLength(a & (-a)) - 1;
  const s = a.toString(16);
  let n = s.length - 1;
  while (s.charCodeAt(n) === '0'.charCodeAt(0)) {
    n -= 1;
  }
  let x = s.charCodeAt(n);
  if (x < 'a'.charCodeAt(0)) {
    x -= '0'.charCodeAt(0);
  } else {
    x -= 'a'.charCodeAt(0);
    x += 10;
  }
  //var e = (31 - Math.clz32(x & -x));
  var e = 0;
  while (x % 2 === 0) {
    x /= 2;
    e += 1;
  }
  return (s.length - 1 - n) * 4 + e;
}

function bigIntGCD(a, b) {
  if (typeof a === "number" && typeof b === "number") {
    return numbersGCD(Math.abs(a), Math.abs(b));
  }
  a = abs(BigInt(a));
  b = abs(BigInt(b));
  if (Number(a) > Number.MAX_SAFE_INTEGER && Number(b) > Number.MAX_SAFE_INTEGER) {
    const c1 = ctz(a);
    const c2 = ctz(b);
    if (c1 > 4 || c2 > 4) {
      const g = LehmersGCD(c1 === 0 ? a : a >> BigInt(c1), c2 === 0 ? b : b >> BigInt(c2));
      const c = Math.min(c1, c2);
      return c === 0 ? g : (g << BigInt(c));
    }
  }
  return LehmersGCD(a, b);
}

export default bigIntGCD;
