

// Pollard's rho implementation
// See https://en.algorithmica.org/hpc/algorithms/factorization/#pollard-brent-algorithm
// And https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm#Variants

function gcd(a, b) {
  if (typeof a !== 'bigint' || typeof b !== 'bigint') {
    throw new RangeError();
  }
  while (b !== 0n) {
    const r = a % b;
    a = b;
    b = r;
  }
  return a;
}

function f(x, c, n) {
  if (typeof x !== 'bigint' || typeof c !== 'bigint' || typeof n !== 'bigint') {
    throw new RangeError();
  }
  let y = (x * x) % n + c;
  y = y >= n ? y - n : y;
  return y;
}

function absDiff(a, b) {
  if (typeof a !== 'bigint' || typeof b !== 'bigint') {
    throw new TypeError();
  }
  return a < b ? b - a : a - b;
}

function internal(n, x0, c, maxSteps) {
  if (typeof n !== 'bigint' || typeof x0 !== 'bigint' || typeof c !== 'bigint') {
    throw new RangeError();
  }
  // Brent cycle detection:
  let hare = x0;
  let i = 0;
  let power = 1;
  let product = 1n;
  let productStartHare = hare;
  let productStartStep = i;
  let cycleFound = false;
  while (i <= maxSteps) {
    const tortoise = hare;
    do {
      hare = f(hare, c, n);
      i += 1;
      if (i <= ((power * 3) >> 2)) { // see https://maths-people.anu.edu.au/~brent/pd/rpb051i.pdf
        product = 1n;
        productStartHare = hare;
        productStartStep = i;
      } else {
        const d = absDiff(hare, tortoise);
        if (product === 1n) {
          product = d;
        } else {
          product = (product * d) % n;
        }
        if (i === power || i % 128 === 0 || cycleFound) {
          const factor = gcd(product, n);
          if (factor !== 1n) {
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
          product = 1n;
          productStartHare = hare;
          productStartStep = i;
        }
      }
    } while (i !== power);
    power <<= 1;
  }
  return 0n;
}

function PollardsRho(n, maxSteps) {
  if (typeof n !== 'bigint') {
    throw new RangeError();
  }
  const x0 = 2n;
  let c = 0n;
  let g = n;
  do {
    // it is much better to switch c, not x0
    // f(x) = x^2 + c, some values are not good with x0=2 (from my tests: 0, -2, -6, -7)
    c += 1n;
    if (c > n) {
      // n is prime?
      return 0n;
    }
    g = internal(n, x0, c, maxSteps);
  } while (g === n);
  return g;
}

export default PollardsRho;
