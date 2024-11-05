
function log2(x) {
  return BigInt(x.toString(16).length * 4);
}

function modPow(base, exponent, modulus) {
  let accumulator = 1n;
  while (exponent !== 0n) {
    if (BigInt.asUintN(1, exponent) === 1n) {
      exponent -= 1n;
      accumulator = (accumulator * base) % modulus;
    }
    exponent >>= 1n;
    base = (base * base) % modulus;
  }
  return accumulator;
}

function primes(MAX) {
  const sieve = new Array(MAX + 1).fill(true);
  const result = [];
  result.push(2);
  for (let i = 3; i <= MAX; i += 2) {
    if (sieve[i]) {
      result.push(i);
      if (i <= Math.floor(MAX / i)) {
        for (let j = i * i; j <= MAX; j += 2 * i) {
          sieve[j] = false;
        }
      }
    }
  }
  return result;
}

function getBases(n) {
  // https://en.wikipedia.org/wiki/Miller–Rabin_primality_test#Testing_against_small_sets_of_bases
  if (typeof n !== 'bigint') {
    throw new TypeError();
  }
  if (n < 2047n) {
    return [2];
  }
  if (n < 1373653n) {
    return [2, 3];
  }
  if (n < 25326001n) {
    return [2, 3, 5];
  }
  if (n < 3215031751n) {
    return [2, 3, 5, 7];
  }
  if (n < 2152302898747n) {
    return [2, 3, 5, 7, 11];
  }
  if (n < 3474749660383n) {
    return [2, 3, 5, 7, 11, 13];
  }
  if (n < 341550071728321n) {
    return [2, 3, 5, 7, 11, 13, 17];
  }
  if (n < 3825123056546413051n) {
    return [2, 3, 5, 7, 11, 13, 17, 19, 23];
  }
  if (n < 318665857834031151167461n) {
    return [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
  }
  if (n < 3317044064679887385961981n) {
    return [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41];
  }
  // https://primes.utm.edu/prove/prove2_3.html
  const lnN = Number(log2(n)) * Math.log(2);
  return primes(Math.floor(1 / Math.log(2) * lnN * Math.log(lnN)));
}

function isPrime(n) {
  if (typeof n !== "bigint") {
    throw new RangeError();
  }
  if (n < 0n) {
    throw new RangeError();
  }
  if (n < 2n) {
    return false;
  }

  const smallPrimes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41];
  const s = Number(n % 304250263527210n);
  for (let i = 0; i < smallPrimes.length; i += 1) {
    const p = smallPrimes[i];
    if (s - Math.floor(s / p) * p === 0) {
      return n === BigInt(p);
    }
  }
  if (n < 43 * 43) {
    return true;
  }

  // https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test#Deterministic_variants
  let r = 0;
  let d = n - 1n;
  while (d % 2n === 0n) {
    d /= 2n;
    r += 1;
  }
  const bases = getBases(n);
  for (const a of bases) {
    let x = modPow(BigInt(a), d, n);
    if (x !== 1n) {
      for (let i = r - 1; i > 0 && x !== n - 1n; i -= 1) {
        x = (x * x) % n;
      }
      if (x !== n - 1n) {
        //if (bases.indexOf(base) >= 1) console.log(n, bases.indexOf(base));
        return false;
      }
    }
  }
  return true;
}

export default isPrime;
