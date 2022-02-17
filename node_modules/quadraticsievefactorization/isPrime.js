
function log2(x) {
  return BigInt(x.toString(16).length * 4);
}

function modPow(base, exponent, modulus) {
  let accumulator = 1n;
  while (exponent !== 0n) {
    if (exponent % 2n === 0n) {
      exponent /= 2n;
      base = (base * base) % modulus;
    } else {
      exponent -= 1n;
      accumulator = (accumulator * base) % modulus;
    }
  }
  return accumulator;
}

function isPrime(n) {
  n = BigInt(n);
  if (n < 2n) {
    throw new RangeError();
  }

  const s = Number(n % 30n);
  if (s % 2 === 0) {
    return n === 2n;
  }
  if (s % 3 === 0) {
    return n === 3n;
  }
  if (s % 5 === 0) {
    return n === 5n;
  }
  const wheel3 = [0, 4, 6, 10, 12, 16, 22, 24, 24];
  for (let i = 7, max = Math.min(1024, Math.floor(Math.sqrt(Number(n)))); i <= max; i += 30) {
    for (let j = 0; j < wheel3.length; j += 3) {
      const p1 = i + wheel3[j + 0];
      const p2 = i + wheel3[j + 1];
      const p3 = i + wheel3[j + 2];
      const r = Number(n % BigInt((p1 * p2 * p3))) | 0;
      if (r % p1 === 0 || r % p2 === 0 || r % p3 === 0) {
        return false;
      }
    }
  }
  if (n < 1024 * 1024) {
    return true;
  }

  // https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Deterministic_variants
  let r = 0;
  let d = n - 1n;
  while (d % 2n === 0n) {
    d /= 2n;
    r += 1;
  }
  // https://en.wikipedia.org/wiki/Miller–Rabin_primality_test#Testing_against_small_sets_of_bases
  const values = [10, 20, 24, 31, 40, 41, 48, 48, 61, 61, 61, 78, 81];
  const primes = [2, 3, 5, 7, 11, 13, 17, 17, 19, 19, 19, 23, 29];
  let i = 0;
  const x = Math.ceil(Math.log2(Number(n)));
  while (x > values[i] && i < values.length) {
    i += 1;
  }
  let bases = null;
  if (i < values.length) {
    bases = primes.slice(0, i + 1);
  } else {
    const lnN = Number(log2(BigInt(n))) * Math.log(2);
    const max = Math.floor(2 * lnN * Math.log(lnN));
    const range = new Array(max - 2 + 1);
    for (let i = 2; i <= max; i += 1) {
      range[i - 2] = i;
    }
    bases = range;
  }
  for (const a of bases) {
    let x = modPow(BigInt(a), d, n);
    if (x !== 1n) {
      for (let i = r - 1; i > 0 && x !== n - 1n; i -= 1) {
        x = (x * x) % n;
      }
      if (x !== n - 1n) {
        return false;
      }
    }
  }
  return true;
}

export default isPrime;