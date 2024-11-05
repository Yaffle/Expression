
// Pollard's rho implementation
// See https://en.algorithmica.org/hpc/algorithms/factorization/#pollard-brent-algorithm
// And https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Variants

//@external("env", "consoleLog")
//declare function consoleLog(a:u64): void;

// ((a * b) >> 64) from https://doc.lagout.org/security/Hackers%20Delight.pdf
function mulhu(a:u64, b:u64):u64 {
  const ahi = a >> 32;
  const alo = a & 0xFFFFFFFF;
  const bhi = b >> 32;
  const blo = b & 0xFFFFFFFF;
  const t = ahi * blo + ((alo * blo) >> 32);
  return ahi * bhi + (t >> 32) + ((alo * bhi + (t & 0xFFFFFFFF)) >> 32);
}

function modInverseMod2pow64(a:u64):u64 {
  let aInv = a;
  for (let e = 2; e < 64; e += e) {
    aInv = u64(aInv * (2 - u64(a * aInv)));
  }
  return aInv;
}

function montgomeryReduce(phi:u64, plo:u64, n:u64, nInv:u64):u64 {
  const x = mulhu(u64(plo * nInv), n);
  return phi - x + (phi < x ? n : 0);
}

function montgomeryMultiply(x:u64, y:u64, n:u64, nInv:u64):u64 {
  return montgomeryReduce(mulhu(x, y), u64(x * y), n, nInv);
}

function mulMod(a:u64, b:u64, m:u64):u64 {
  let s = u64(0);
  while (a !== 0) {
    s = ((a & 1) !== 0 ? s - (m - b) + (s < m - b ? m : 0) : s);
    b = b - (m - b) + (b < m - b ? m : 0);
    a >>= 1;
  }
  return s;
}

function montgomeryTransform(x:u64, n:u64):u64 {
  // (x * 2**64) % n :
  return mulMod(x, u64(-1) % n + 1, n);
}

export function gcdBinary(a:u64, b:u64):u64 {
  if (a === 0) {
    return b;
  }
  if (b === 0) {
    return a;
  }
  const z = i64.ctz(a | b);
  a >>= i64.ctz(a);
  b >>= i64.ctz(b);
  do {
    const absDiff = a < b ? b - a : a - b; // abs(a - b)
    b = a < b ? a : b; // min(a, b)
    a = absDiff;
    a >>= i64.ctz(a);
  } while (a !== 0);
  b <<= z;
  return b;
}

function f(x:u64, c:u64, n:u64, nInv:u64):u64 {
  let y = montgomeryMultiply(x, x, n, nInv);
  y = y - (n - c) + (y < (n - c) ? n : 0);
  return y;
}

function absDiff(a:u64, b:u64):u64 {
  return a < b ? b - a : a - b;
}

function internal(n:u64, x0:u64, c:u64, maxSteps:u32):u64 {
  const nInv = modInverseMod2pow64(n);
  x0 = montgomeryTransform(x0, n);
  c = montgomeryTransform(c, n);
  const one = u64(1); // (2**64)**-1 mod n in Montgomery Form
  // Brent cycle detection:
  let hare = x0;
  let i = u32(0);
  let power = u32(1);
  let product = one;
  let productStartHare = hare;
  let productStartStep = i;
  let cycleFound = false;
  while (i <= maxSteps) {
    const tortoise = hare;
    do {
      hare = f(hare, c, n, nInv);
      i += 1;
      if (i <= ((power * 3) >> 2)) { // see https://maths-people.anu.edu.au/~brent/pd/rpb051i.pdf
        product = one;
        productStartHare = hare;
        productStartStep = i;
      } else {
        const d = absDiff(hare, tortoise);
        if (product === one) {
          product = d;
        } else {
          product = montgomeryMultiply(product, d, n, nInv);
        }
        if (i === power || i % 128 === 0 || cycleFound) {
          const factor = gcdBinary(product, n);
          if (factor !== 1) {
            if (factor !== n) {
              return factor;
            }
            if (cycleFound) {
              return n;
            }
            cycleFound = true;
            hare = productStartHare;
            i = productStartStep;
          }
          product = one;
          productStartHare = hare;
          productStartStep = i;
        }
      }
    } while (i !== power);
    power <<= 1;
  }
  return 0;
}

export function PollardsRho64(n:u64, maxSteps:i32):u64 {
  if (n <= 2) {
    return 0;
  }
  if (n % 2 === 0) {
    return 2;
  }
  const x0 = u64(2);
  let c = u64(0);
  let g = n;
  do {
    // it is much better to switch c, not x0
    // f(x) = x^2 + c, some values are not good with x0=2 (from my tests: 0, -2, -6, -7)
    c += 1;
    if (c > 42) {
      // n is prime?
      return 0;
    }
    g = internal(n, x0, c, maxSteps);
  } while (g === n);
  return g;
}
