"use strict";

// https://rosettacode.org/wiki/Tonelli-Shanks_algorithm#Python

function modPowSmall(base, exponent, modulus) {
  let accumulator = 1;
  while (exponent !== 0) {
    if (exponent % 2 === 0) {
      exponent /= 2;
      base = (base * base);
      base = base - Math.floor(base / modulus) * modulus;
    } else {
      exponent -= 1;
      accumulator = (accumulator * base);
      accumulator = accumulator - Math.floor(accumulator / modulus) * modulus;
    }
  }
  return accumulator;
}


//function legendre(a, p) {
//  console.assert((p - 1) % 2 === 0);
//  return (modPowSmall(a, (p - 1) / 2, p) + 1) % p - 1;
//}

// a/n is represented as (a,n)
function legendre(a, n) {
    if (typeof a !== 'number' || typeof n !== 'number') {
      throw new TypeError();
    }
    // from https://en.wikipedia.org/wiki/Jacobi_symbol#Implementation_in_C++ :
    a = a | 0;
    n = n | 0;
    //console.assert(n > 0 && n%2 == 1);
    //step 1
    a = (a | 0) % (n | 0);
    var t = 1;
    var r = 0;
    //step 3
    while (a !== 0) {
        //step 2
        while ((a & 1) === 0) {
            a >>= 1;
            r = n & 7;
            if (r === 3 || r === 5) {
                t = -t;
            }
        }
        //step 4
        r = n;
        n = a;
        a = r;
        if ((a & 3) === 3 && (n & 3) === 3) {
            t = -t;
        }
        a = (a | 0) % (n | 0);
    }
    return n === 1 ? t : 0;
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
  while (z <= p && legendre(z, p) !== -1) {
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
