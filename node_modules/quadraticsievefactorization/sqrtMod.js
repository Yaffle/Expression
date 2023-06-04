
// https://rosettacode.org/wiki/Tonelli-Shanks_algorithm#Python

function modPowSmall(base, exponent, modulus) {
  let accumulator = 1;
  while (exponent !== 0) {
    if (exponent % 2 === 0) {
      exponent /= 2;
      base = (base * base) % modulus;
    } else {
      exponent -= 1;
      accumulator = (accumulator * base) % modulus;
    }
  }
  return accumulator;
}

function legendre(a, p) {
  console.assert((p - 1) % 2 === 0);
  return modPowSmall(a, (p - 1) / 2, p);
}

function sqrtMod(n, p) {
  // Tonelli–Shanks algorithm:
  console.assert(p * p <= Number.MAX_SAFE_INTEGER);
  if (legendre(n, p) !== 1) {
    // not a square (mod p)
    return -1;
  }
  let q = p - 1;
  let s = 0;
  while (q % 2 === 0) {
    q /= 2;
    s += 1;
  }
  if (s === 1) {
    return modPowSmall(n, (p + 1) / 4, p);
  }
  let z = 2;
  while (z <= p && p - 1 !== legendre(z, p)) {
    z += 1;
  }
  let c = modPowSmall(z, q, p);
  let r = modPowSmall(n, (q + 1) / 2, p);
  let t = modPowSmall(n, q, p);
  let m = s;
  let t2 = 0;
  while ((t - 1) % p !== 0) {
    t2 = (t * t) % p;
    let i = 1;
    while (i <= m && (t2 - 1) % p !== 0) {
      t2 = (t2 * t2) % p;
      i += 1;
    }
    const b = modPowSmall(c, 1 << (m - i - 1), p);
    r = (r * b) % p;
    c = (b * b) % p;
    t = (t * c) % p;
    m = i;
  }
  return r;
}

export default sqrtMod;
