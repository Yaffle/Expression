/*jshint esversion:11, bitwise:false*/

import solve from './solve.js';
import sqrtMod from './sqrtMod.js';
import PollardsRho64 from './libs/PollardsRho64.js';
import isProbablyPrime64 from './libs/isProbablyPrime64.js';

function modInverse(a, m) {
  if (typeof a !== 'bigint' || typeof m !== 'bigint') {
    throw new TypeError();
  }
  if (Number(m) <= (2**30 - 1)) {
    return BigInt(modInverseSmall(Number(a), Number(m)));
  }
  if (a < 0n || a >= m || m <= 0n) {
    throw new RangeError();
  }
  // We use the extended Euclidean algorithm:
  let oldR = a;
  let r = m;
  let oldS = 1n;
  let s = 0n;
  while (r !== 0n) {
    const q = (oldR - oldR % r) / r; // floor(oldR / r)
    const newR = oldR - q * r;
    oldR = r;
    r = newR;
    const newS = oldS - q * s;
    oldS = s;
    s = newS;
  }
  if (oldR !== 1n) {
    return 0n;
  }
  return oldS < 0n ? oldS + m : oldS;
}

function modInverseSmall(a, m) {
  if (typeof a !== 'number' || typeof m !== 'number') {
    throw new TypeError();
  }
  const maxSMI = (~(-1 << 30));
  if (a < 0 || a >= m || m <= 0 || m > maxSMI) {
    throw new RangeError();
  }
  // We use the extended Euclidean algorithm:
  let oldR = a & maxSMI;
  let r = m & maxSMI;
  let oldS = 1;
  let s = 0;
  while (r !== 0) {
    const q = Math.floor(oldR / r);
    const newR = oldR % r;
    oldR = r;
    r = newR;
    const newS = oldS - q * s;
    oldS = s;
    s = newS;
  }
  if (oldR !== 1) {
    return 0;
  }
  return oldS < 0 ? oldS + m : oldS;
}

function ChineseRemainderTheorem(r1, r2, m1, m2) {
  if (typeof r1 !== 'bigint' || typeof r2 !== 'bigint' || typeof m1 !== 'bigint' || typeof m2 !== 'bigint') {
    throw new TypeError();
  }
  // https://en.wikipedia.org/wiki/Chinese_remainder_theorem#Case_of_two_moduli
  // x = r1 (mod m1)
  // x = r2 (mod m2)
  const c = modInverse(m1 % m2, m2);
  return r1 + (((r2 - r1) * c) % m2) * m1;
}

function squareRootsModuloOddPrimesProduct(n, primes, e = 1) {
  // Chinese Remainder Theorem idea from https://en.wikipedia.org/wiki/Quadratic_residue#Complexity_of_finding_square_roots
  let result = [];
  result.push(0n);
  let P = 1n;
  for (let i = 0; i < primes.length; i += 1) {
    const p = BigInt(Math.pow(primes[i], e));
    if (Number(p) > Number.MAX_SAFE_INTEGER) {
      throw new RangeError();
    }
    const x2 = BigInt(squareRootModuloOddPrime(Number(n % p), primes[i], e));
    const result2 = [];
    for (let j = 0; j < result.length; j += 1) {
      const x1 = result[j];
      result2.push(ChineseRemainderTheorem(x1, x2, P, p));
      result2.push(ChineseRemainderTheorem(x1, -x2, P, p));
    }
    P *= p;
    result = result2;
  }
  return result;
}

function squareRootsModuloTwo(n, e = 1) {
  if (e >= 3) {
    if (n % 8 === 1) { // from Cohen H.
      const m = Math.pow(2, e);
      const candidate = +squareRootsModuloTwo(n, e - 1)[0];
      const candidate2 = m / 4 - candidate;
      const r = (candidate * candidate) % m !== n % m ? candidate2 : candidate;
      return [r, m / 2 - r, m / 2 + r, m - r];
    }
    return [];
  }
  if (e >= 2) {
    if (n % 4 === 1) {
      return [1, 3];
    }
    return [];
  }
  if (e >= 1) {
    return [1];
  }
  return [];
}

function squareRootModuloOddPrime(n, p, e = 1) { // slow for non-small p
  if (typeof n !== 'number' || typeof p !== 'number' || typeof e !== 'number') {
    throw new TypeError();
  }
  const m = Math.pow(p, e);
  if (!(n > 0 && n < m && p > 0 && p % 2 !== 0 && e >= 1 && +n % +p !== 0 && m < Math.floor(Math.sqrt(Number.MAX_SAFE_INTEGER * 4)))) { // + p is a prime number
    throw new RangeError();
  }
  // r**2 == n (mod p)
  if (e > 1) {
    const x = squareRootModuloOddPrime(n % Math.pow(p, e - 1), p, e - 1);
    // x**2 = n mod p**(e - 1)
    // x1 = x + a * p**(e-1)
    // x1**2 = x**2 + (a * p**(e-1))**2 + 2*x*a*p**(e-1) = n mod p**e
    // a*p**(e-1) = (n - x**2) * (2*x)**-1 mod p**e
    let inv = modInverseSmall(2 * x, m) % m;
    let v = (n - x * x) % m;
    inv = inv > m / 2 ? inv - m : inv; // use sign bit
    v = v > m / 2 ? v - m : v; // use sign bit
    let x1 = x + (v * inv) % m;
    if (x1 >= m) {
      x1 -= m;
    }
    if (x1 < 0) {
      x1 += m;
    }
    return Math.min(x1, m - x1);
  }
  const r = sqrtMod(n, p) | 0;
  return Math.min(r, p - r);
}

function bitLength(x) {
  return x.toString(16).length * 4;
}

function sqrt(x) {
  if (x < BigInt((Number.MAX_SAFE_INTEGER + 1) / 2)) {
    return BigInt(Math.floor(Math.sqrt(Number(BigInt(x)) + 0.5)));
  }
  const q = BigInt(bitLength(x) >> 2);
  const initialGuess = ((sqrt(x >> (q * 2n)) + 1n) << q);
  let a = initialGuess;
  let b = a + 1n;
  while (a < b) {
    b = a;
    a = (b + x / b) >> 1n;
  }
  return b;
}

function getSmoothFactorization(a, base) {
  let value = BigInt(a);
  if (value === 0n) {
    return [0];
  }
  const result = [];
  if (value < 0n) {
    result.push(-1);
    value = -value;
  }
  let i = 0;
  while (i < base.length) {
    const p = base[i];
    while (value % BigInt(p) === 0n) {
      value /= BigInt(p);
      result.push(p);
    }
    i += 1;
  }
  return value === 1n ? result : null;
}

// (X**2 - Y) % N === 0, where Y is a smooth number
function CongruenceOfsquareOfXminusYmoduloN(X, Y, N) {
  this.X = X;
  this.Y = Y;
  this.N = N;
}
CongruenceOfsquareOfXminusYmoduloN.prototype.toString = function () {
  return 'X**2 ≡ Y (mod N)'.replaceAll('X', this.X)
                           .replaceAll('N', this.N)
                           .replaceAll('Y', this.Y.join(' * '));
};

function isQuadraticResidueModuloPrime(a, p) {
  if (typeof a !== 'bigint' || typeof p !== 'number') {
    throw new TypeError();
  }
  if (p === 2) {
    // "Modulo 2, every integer is a quadratic residue." - https://en.wikipedia.org/wiki/Quadratic_residue#Prime_modulus
    return true;
  }
  // https://en.wikipedia.org/wiki/Euler%27s_criterion
  const amodp = Number(BigInt(a) % BigInt(p));
  if (amodp === 0) {
    return true;
  }
  return legendre(amodp, p) === 1;
}

function significand(a) {
  const s1 = Number(a);
  if (Math.abs(s1) < 1/0) {
    return s1 / 2**Math.floor(Math.log2(Math.abs(s1)));
  }
  const e = Math.max(0, bitLength(a) - 1023);
  const s = Number(a >> BigInt(e));
  return s / 2**Math.floor(Math.log2(Math.abs(s)));
}

function exponent(a) {
  const s1 = Number(a);
  if (Math.abs(s1) < 1/0) {
    return Math.floor(Math.log2(Math.abs(s1)));
  }
  const e = Math.max(0, bitLength(a) - 1023);
  const s = Number(a >> BigInt(e));
  return e + Math.floor(Math.log2(Math.abs(s)));
}

function log2(a) {
  return Math.log2(significand(a)) + exponent(a);
}

function L(N) {  // exp(sqrt(log(n)*log(log(n))))
  const lnn = log2(N) * Math.log(2);
  return Math.exp(Math.sqrt(lnn * Math.log(lnn)));
}

function modPow(base, exponent, modulus) {
  if (typeof base !== 'bigint' || typeof exponent !== 'bigint' || typeof modulus !== 'bigint') {
    throw new TypeError();
  }
  let e = exponent;
  let b = base;
  let accumulator = 1n;
  while (e !== 0n) {
    if (BigInt.asUintN(1, e) === 1n) {
      e -= 1n;
      accumulator = (accumulator * b) % modulus;
    }
    e >>= 1n;
    b = (b * b) % modulus;
  }
  return accumulator;
}

function primes(MAX) {
  // Note: it is slow in Chrome to create array this way when MAX >= 2**25
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

//!copy-paste

function FastModBigInt(a) {
  const array = [];
  while (a !== 0n) {
    const x = Number(BigInt.asUintN(51, a));
    array.push(x);
    a >>= 51n;
  }
  return new Float64Array(array);
}
function FastMod(array, m, mInv) { // mInv === (1 + 2**-52) / m
  m = -0 + m;
  mInv = -0 + mInv;
  const n = array.length - 1;
  let result = array[n];
  result = result - Math.floor(result * mInv) * m;
  if (n > 0) {
    const x = 2**51 - Math.floor(2**51 * mInv) * m;
    let i = n;
    do {
      i -= 1;
      result = result * x + array[i];
      result = result - Math.floor(result * mInv) * m;
    } while (i !== 0);
  }
  return -0 + result;
}

//squareRootModuloOddPrime(4865648, 9749, 2)  // huge size of p**e

function exp2(x) {
  let y = Math.floor(Math.pow(2, Math.floor(x)) * Math.exp(Math.LN2 * (x - Math.floor(x))));
  if (y % 2 === 0) {
    y += 1;
  }
  return y;
}

const useMultiplePolynomials = true;


// Q=((ax+b)**2-n)/a
// a~=sqrt(2n)/M
// max value = n/a = sqrt(n)*M * 1/sqrt(2)

// Q2=((2ax+b)**2-n)/(4a)
// a=sqrt(n/2)/M
// max value = n/(4a) = sqrt(n)*M * 1/(2*sqrt(2)) - max is two times smaller

// (A * x + B)**2 - N = A * (A * x**2 + 2 * B * x + C), A * C = B**2 - N
function QuadraticPolynomial(A, B, N, AFactors, useQ2Form) {
  if (typeof A !== 'bigint' || typeof B !== 'bigint' || typeof N !== 'bigint') {
    throw new TypeError();
  }
  const BBmN = (B * B - N);
  if (useQ2Form && N % 4n !== 1n) {
    throw new RangeError();
  }
  if (BBmN % ((useQ2Form ? 4n : 1n) * A) !== 0n) {
    throw new TypeError();
  }
  const C = BBmN / (useQ2Form ? 4n * A : A);
  this.A = A;
  this.B = B;
  this.C = C;
  this.AFactors = AFactors;
  const u = (significand(B) / significand(A)) * Math.pow(2, exponent(B) - exponent(A)); // B / A
  const v = Math.sqrt(significand(N) / Math.pow(significand(A), 2)) * Math.pow(2, exponent(N) / 2 - exponent(A)); // N**(1/2) / A
  this.x1 = (useQ2Form ? 1/2 : 1) * (-u - v);
  this.x2 = (useQ2Form ? 1/2 : 1) * (-u + v);
  this.log2a = log2(A);
  this.useQ2Form = useQ2Form;
}
QuadraticPolynomial.generator = function (M, primes, N) {
  // see https://www.cs.virginia.edu/crab/QFS_Simple.pdf for multiple primes optimization
  const getCombinations = function (elements, k) {
    if (elements.length === 0) {
      return [];
    }
    if (k === 0) {
      return [[]];
    }
    if (k === 1) {
      return elements.map(e => [e]);
    }
    return getCombinations(elements.slice(1), k - 1).map(c => [elements[0]].concat(c)).concat(getCombinations(elements.slice(1), k));
  };
  const nthRootApprox = function (A, n) {
    return Math.round(Math.pow(significand(A), 1 / n) * Math.pow(2, exponent(A) / n));
  };
  const max = function (a, b) {
    return a < b ? b : a;
  };
  const useQ2Form = N % 8n === 1n;
  const S = max(1n, useQ2Form ? BigInt(sqrt(N / 2n)) / BigInt(M) : BigInt(sqrt(2n * N)) / BigInt(M));
  const e = log2(S);
  const max1 = Math.log2(primes[primes.length - 1]);
  
  // see https://www.cecm.sfu.ca/~mmonagan/papers/NT4.pdf
  // "... Contini [5] recommends minimizing the number of duplicate relations found by requiring that the sets {qi} differ by at least two primes ..."
  const elementPrimes = 2;
  const k = Math.max(elementPrimes, Math.ceil(e / Math.min(e < 180 ? 14.5 : max1, max1) / elementPrimes) * elementPrimes); // number of small primes
  console.debug('k: ', k, 'useQ2Form: ', useQ2Form);
  const p = nthRootApprox(S, k);
  let s = 0;
  const nextPrime = function () {
    let p3 = 0;
    do {
      p3 = p - p % 2 + 1 + (s % 2 === 0 ? s : (-1 - s));
      s += 1;
      if (p3 > primes[primes.length - 1]) {
        throw new RangeError();
      }
    } while (indexOf(primes, p3) === -1);
    return p3;
  };
  let combinations = [];
  const polynomials = [];
  const elements = [];
  QuadraticSieveFactorization.polynomialsCounter = 0;
  return {
    next: function generator() {
      while (polynomials.length === 0) {
        // There must be at least two different primes from previous selections. - from https://www.rieselprime.de/ziki/Self-initializing_quadratic_sieve
        while (combinations.length === 0) {
          const element = [];
          for (let i = 0; i < elementPrimes; i += 1) {
            element.push(nextPrime());
          }
          console.assert(k % elementPrimes === 0);
          combinations = getCombinations(elements, k / elementPrimes - 1).map(c => [element].concat(c));
          elements.push(element);
          //console.log(elements.length, combinations.length, p**k / Number(S));
        }
        const qPrimes = combinations.pop().reduce((array, pair) => array.concat(pair), []);
        //console.debug('qPrimes', qPrimes.join(' '));
        const q = BigInt(qPrimes.reduce((p, a) => p * BigInt(a), 1n));
        const qInv = modInverse(q % N, N);
        if (qInv === 0n) {
          //TODO: what to do here - ?
          return this.next();
        }
        const A = q;
        let Bs = squareRootsModuloOddPrimesProduct(N, qPrimes, 1);
        for (let i = 0; i < Bs.length; i += 1) {
          Bs[i] = Bs[i] < 0n ? Bs[i] + A : Bs[i];
          if (Bs[i] < 0n || Bs[i] >= A) throw new Error();
        }
        if (useQ2Form) {
          Bs = Bs.map(B => B % 2n === 0n ? B - A : B);
        }
        Bs = Bs.slice(0, Bs.length / 2);
        //DO NOT SORT!!!
        //Bs.sort((a, b) => Number(BigInt(a) - BigInt(b)));
        const aFactors = useQ2Form ? [2, 2].concat(qPrimes) : qPrimes.slice(0);
        aFactors.sort((a, b) => a - b);
        for (let i = 0; i < Bs.length; i += 1) {
          const B = Bs[i];
          polynomials.push(new QuadraticPolynomial(A, B, N, aFactors, useQ2Form));
        }
      }
      QuadraticSieveFactorization.polynomialsCounter += 1;
      return polynomials.shift();
    }
  };
};
QuadraticPolynomial.prototype.X = function (x) {
  return (this.A * BigInt((this.useQ2Form ? 2 : 1) * x) + this.B);
};
QuadraticPolynomial.prototype.Y = function (x, s, primes) {
  if (typeof x !== 'number') {
    throw new TypeError();
  }
  const Y = this.A * (x * x >= 2**53 ? BigInt(x) * BigInt(x) : BigInt(x * x)) + this.B * BigInt((this.useQ2Form ? 1 : 2) * x) + this.C;
  if (s === undefined && primes === undefined) {
    return Y;// for debugging
  }
  if (Y % s !== 0n) {
    return null;
  }
  const YFactors = getSmoothFactorization(Y / s, primes);
  if (YFactors == null) {
    return null;
  }
  if (YFactors.length === 1 && YFactors[0] === 0) {
    return YFactors;
  }
  return this.AFactors.concat(YFactors);
};
QuadraticPolynomial.prototype.log2AbsY = function (x) {
  if (typeof x !== 'number') {
    throw new TypeError();
  }
  //const v1 = Math.log2(Math.abs(Number(this.Y(x))));
  const v2 =  Math.log2(Math.abs((x - this.x1) * (x - this.x2))) + this.log2a;
  return v2;
};

function thresholdApproximationInterval(polynomial, x, threshold, sieveSize) {
  let w = sieveSize > 2048 ? (sieveSize > 2**18 ? 1024 : 256) : 1;
  while (w >= 2 && Math.abs(polynomial.log2AbsY(x + w) - threshold) > 0.5) {
    w = Math.floor(w / 2);
  }
  return x + w;
}

// https://ru.wikipedia.org/wiki/Алгоритм_Диксона
// https://www.youtube.com/watch?v=TvbQVj2tvgc
// https://www.rieselprime.de/ziki/Self-initializing_quadratic_sieve

const wasmCode = new Uint8Array([0, 97, 115, 109, 1, 0, 0, 0, 1, 28, 3, 96, 8, 127, 127, 127, 127, 127, 127, 127, 127, 1, 127, 96, 2, 127, 127, 1, 127, 96, 6, 127, 127, 127, 127, 127, 127, 0, 2, 15, 1, 3, 101, 110, 118, 6, 109, 101, 109, 111, 114, 121, 2, 0, 0, 3, 5, 4, 0, 1, 2, 0, 7, 85, 4, 16, 115, 105, 110, 103, 108, 101, 66, 108, 111, 99, 107, 83, 105, 101, 118, 101, 0, 0, 15, 102, 105, 110, 100, 83, 109, 111, 111, 116, 104, 69, 110, 116, 114, 121, 0, 1, 24, 117, 112, 100, 97, 116, 101, 87, 104, 101, 101, 108, 115, 73, 110, 116, 101, 114, 110, 97, 108, 78, 101, 120, 116, 0, 2, 17, 104, 97, 110, 100, 108, 101, 83, 109, 97, 108, 108, 87, 104, 101, 101, 108, 115, 0, 3, 10, 154, 8, 4, 226, 1, 1, 4, 127, 32, 2, 32, 5, 106, 33, 8, 32, 0, 32, 4, 106, 33, 0, 32, 1, 32, 4, 106, 33, 1, 32, 2, 32, 4, 106, 33, 2, 32, 3, 32, 4, 65, 2, 117, 106, 33, 3, 3, 64, 32, 2, 32, 8, 71, 4, 64, 32, 0, 40, 2, 0, 33, 5, 32, 1, 40, 2, 0, 33, 4, 32, 2, 40, 2, 0, 33, 10, 32, 3, 45, 0, 0, 33, 11, 3, 64, 32, 4, 32, 6, 72, 4, 64, 32, 4, 45, 0, 0, 32, 11, 106, 33, 9, 32, 5, 32, 5, 45, 0, 0, 32, 11, 106, 58, 0, 0, 32, 4, 32, 9, 58, 0, 0, 32, 5, 32, 10, 106, 33, 5, 32, 4, 32, 10, 106, 33, 4, 12, 1, 11, 11, 32, 5, 32, 6, 72, 4, 64, 32, 5, 32, 5, 45, 0, 0, 32, 11, 106, 58, 0, 0, 32, 5, 32, 10, 106, 33, 9, 32, 4, 33, 5, 32, 9, 33, 4, 11, 32, 0, 32, 5, 32, 7, 107, 54, 2, 0, 32, 1, 32, 4, 32, 7, 107, 54, 2, 0, 32, 0, 65, 4, 106, 33, 0, 32, 1, 65, 4, 106, 33, 1, 32, 2, 65, 4, 106, 33, 2, 32, 3, 65, 1, 106, 33, 3, 12, 1, 11, 11, 65, 0, 11, 42, 1, 1, 123, 32, 0, 253, 15, 33, 2, 32, 1, 65, 16, 107, 33, 1, 3, 64, 32, 1, 65, 16, 106, 34, 1, 253, 0, 4, 0, 32, 2, 253, 44, 253, 100, 69, 13, 0, 11, 32, 1, 11, 173, 1, 2, 5, 123, 1, 127, 32, 5, 253, 17, 33, 6, 3, 64, 32, 11, 32, 0, 65, 2, 116, 72, 4, 64, 32, 3, 32, 3, 253, 0, 4, 0, 34, 8, 32, 1, 253, 0, 4, 0, 34, 7, 32, 6, 253, 177, 1, 32, 7, 32, 6, 253, 55, 32, 2, 253, 0, 4, 0, 34, 9, 253, 78, 253, 174, 1, 34, 7, 253, 177, 1, 32, 8, 32, 7, 253, 57, 32, 9, 253, 78, 253, 174, 1, 34, 8, 32, 4, 253, 0, 4, 0, 34, 10, 32, 7, 253, 177, 1, 32, 10, 32, 7, 253, 57, 32, 9, 253, 78, 253, 174, 1, 34, 7, 253, 182, 1, 253, 11, 4, 0, 32, 4, 32, 8, 32, 7, 253, 184, 1, 253, 11, 4, 0, 32, 1, 65, 16, 106, 33, 1, 32, 2, 65, 16, 106, 33, 2, 32, 3, 65, 16, 106, 33, 3, 32, 4, 65, 16, 106, 33, 4, 32, 11, 65, 16, 106, 33, 11, 12, 1, 11, 11, 11, 217, 4, 2, 10, 127, 5, 123, 32, 6, 32, 7, 70, 4, 64, 0, 11, 32, 5, 33, 8, 3, 64, 32, 13, 32, 0, 65, 2, 116, 72, 4, 64, 32, 1, 253, 0, 4, 0, 65, 16, 253, 171, 1, 65, 16, 253, 173, 1, 32, 1, 253, 0, 4, 16, 65, 16, 253, 171, 1, 253, 80, 33, 20, 32, 2, 253, 0, 4, 0, 65, 16, 253, 171, 1, 65, 16, 253, 173, 1, 32, 2, 253, 0, 4, 16, 65, 16, 253, 171, 1, 253, 80, 33, 21, 32, 3, 253, 0, 4, 0, 65, 16, 253, 171, 1, 65, 16, 253, 173, 1, 32, 3, 253, 0, 4, 16, 65, 16, 253, 171, 1, 253, 80, 33, 19, 32, 4, 253, 0, 4, 0, 65, 16, 253, 171, 1, 65, 16, 253, 173, 1, 32, 4, 253, 0, 4, 16, 65, 16, 253, 171, 1, 253, 80, 33, 22, 32, 6, 65, 4, 107, 33, 10, 253, 12, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 33, 18, 3, 64, 32, 18, 32, 20, 32, 10, 65, 4, 106, 34, 10, 253, 8, 1, 0, 34, 18, 253, 145, 1, 32, 19, 253, 149, 1, 32, 21, 32, 18, 253, 145, 1, 32, 19, 253, 149, 1, 253, 151, 1, 253, 151, 1, 33, 18, 32, 7, 32, 10, 74, 13, 0, 11, 32, 18, 32, 22, 253, 52, 253, 100, 4, 64, 32, 6, 33, 10, 3, 64, 32, 7, 32, 10, 74, 4, 64, 32, 20, 32, 10, 253, 8, 1, 0, 34, 18, 253, 145, 1, 32, 19, 253, 149, 1, 32, 21, 32, 18, 253, 145, 1, 32, 19, 253, 149, 1, 253, 151, 1, 32, 22, 253, 52, 253, 100, 34, 11, 4, 64, 65, 0, 33, 9, 3, 64, 32, 11, 65, 3, 32, 11, 104, 34, 12, 116, 115, 33, 11, 32, 9, 65, 1, 32, 12, 65, 1, 118, 34, 9, 65, 1, 113, 65, 2, 116, 32, 9, 65, 1, 117, 114, 116, 115, 33, 9, 32, 11, 13, 0, 11, 32, 9, 33, 11, 32, 10, 40, 2, 0, 33, 14, 3, 64, 32, 11, 104, 34, 9, 65, 2, 116, 33, 15, 32, 11, 65, 1, 32, 9, 116, 115, 33, 11, 32, 4, 32, 15, 106, 40, 2, 0, 34, 12, 32, 3, 32, 15, 106, 40, 2, 0, 34, 9, 32, 2, 32, 15, 106, 40, 2, 0, 34, 16, 32, 14, 107, 108, 79, 32, 12, 32, 9, 32, 1, 32, 15, 106, 40, 2, 0, 34, 17, 32, 14, 107, 108, 79, 114, 65, 0, 65, 1, 32, 16, 32, 17, 27, 27, 4, 64, 32, 8, 65, 2, 116, 32, 13, 32, 15, 106, 65, 2, 117, 34, 12, 54, 2, 0, 32, 8, 65, 1, 106, 65, 2, 116, 32, 10, 32, 6, 107, 65, 2, 117, 34, 9, 54, 2, 0, 32, 8, 65, 2, 106, 33, 8, 32, 16, 32, 17, 70, 4, 64, 32, 8, 65, 2, 116, 32, 12, 54, 2, 0, 32, 8, 65, 1, 106, 65, 2, 116, 32, 9, 54, 2, 0, 32, 8, 65, 2, 106, 33, 8, 11, 11, 32, 11, 13, 0, 11, 11, 32, 10, 65, 4, 106, 33, 10, 12, 1, 11, 11, 11, 32, 1, 65, 32, 106, 33, 1, 32, 2, 65, 32, 106, 33, 2, 32, 3, 65, 32, 106, 33, 3, 32, 4, 65, 32, 106, 33, 4, 32, 13, 65, 32, 106, 33, 13, 12, 1, 11, 11, 32, 8, 32, 5, 107, 11]);

let wasmModule = null;
function instantiateWasm(memorySize) {
  if (wasmModule == null) {
    wasmModule = new WebAssembly.Module(wasmCode);
  }
  const pages = Math.ceil(memorySize / 2**16);
  const memory = new WebAssembly.Memory({
    initial: pages,
    maximum: pages
  });
  const buffer = memory.buffer;
  const exports = new WebAssembly.Instance(wasmModule, {env: { memory: memory }}).exports;
  return Object.assign({}, exports, {memory: {buffer: buffer}});
}

// TOWO: WebAssembly (~17% faster)


function congruencesUsingQuadraticSieve(primes, N, sieveSize0) {
  if (typeof N !== 'bigint') {
    throw new TypeError();
  }
  let sieveSize1 = Number(sieveSize0 || 0);
  if (sieveSize1 === 0) {
    sieveSize1 = 3 * 2**14;
    sieveSize1 = Math.min(sieveSize1, Math.ceil(Math.pow(+primes[primes.length - 1], 1.15)));
    sieveSize1 = Math.max(sieveSize1, primes[primes.length - 1] + 1);
    if (Number(N) > 2**285) {
      sieveSize1 = Math.floor(sieveSize1 / 3.2 / 1.5);
    } else if (Number(N) > 2**240) {
      sieveSize1 = Math.floor(sieveSize1 / 1.6 / 1.5);
    } else {
      sieveSize1 = Math.floor(sieveSize1 / 1.5);
    }
  }

  const q = Math.ceil(sieveSize1 / (typeof navigator !== 'undefined' && navigator.hardwareConcurrency === 12 ? 2.75 * 2**20 : 6 * 2**20));
  console.debug('q', q);
  const segmentSize = Math.ceil(Math.ceil(sieveSize1 / q) / 2**15) * 2**15;
  const sieveSize = segmentSize * q;
  console.debug('sieveSize', sieveSize);

  const SHIFT = 0;
  const MAX = 255;
  const SCALE = 2**0;//TODO:

  const log2B = Math.log2(primes.length === 0 ? Math.sqrt(2) : +primes[primes.length - 1]);
  const largePrimesThreshold = log2B + Math.min(Math.log2(N < 2**240 ? 50 : (N < 2**285 ? 400 : 1000)), log2B);
  console.debug('largePrimesThreshold', largePrimesThreshold);
  const largePrimes = new Map(); // faster (?)

  const doubleLargePrimes = Number(N) > 2**285;
  const doubleLargePrimesThreshold = Math.min(2 * log2B + Math.min(Math.log2(200), log2B), 64);
  if (doubleLargePrimes) {
    console.log('doubleLargePrimesThreshold', doubleLargePrimesThreshold);
  }
  

  // see https://www.youtube.com/watch?v=TvbQVj2tvgc
  const wheels0 = [];
  for (let i = 0; i < primes.length; i += 1) {
    const p = +primes[i];
    for (let beta = 1, pInBeta = p; pInBeta <= sieveSize || beta === 1; beta += 1, pInBeta *= p) {
      const nmodpInBeta = Number(N % BigInt(pInBeta));
      if (nmodpInBeta % p === 0) {
        //console.warn('N has a factor in prime base', N, p);
        if (beta === 1) {
          wheels0.push({step: pInBeta, p: p, root: 0});
        }
      } else {
        if (p === 2) {
          const roots = squareRootsModuloTwo(nmodpInBeta, beta);
          for (let j = 0; j < Math.ceil(roots.length / 2); j += 1) {
            wheels0.push({step: pInBeta, p: p, root: roots[j] | 0});
          }
        } else {
          const root = squareRootModuloOddPrime(nmodpInBeta, p, beta);
          wheels0.push({step: pInBeta, p: p, root: root | 0});
        }
      }
    }
  }
  wheels0.sort((a, b) => +a.step - +b.step);
  const realWheels = wheels0.length;
  while (wheels0.length % 16 !== 0) {
    wheels0.push(wheels0.at(-1));
  }
  console.debug('wheels', wheels0.length);
  console.debug('max wheel', wheels0.at(-1).step);

  const wheelLogs0 = new Float64Array(wheels0.length);
  let invCacheKey = 0n;
  const zeroInvs = [];
  
  const wheelsCount = wheels0.length;
  console.assert(segmentSize % 4 === 0);

  let memorySize = 0;
  memorySize += segmentSize * 1;

  const wheelRoots1 = memorySize >> 2;
  memorySize += wheelsCount * 4;
  const wheelRoots2 = memorySize >> 2;
  memorySize += wheelsCount * 4;
  const wheelSteps = memorySize >> 2;
  memorySize += wheelsCount * 4;
  const wheelLogs = memorySize >> 0;
  memorySize += wheelsCount;

  const storage = memorySize >> 2;
  memorySize += (16 * 1024) * 4;//TODO: what size to use?

  const wheelRoots = memorySize >> 2;
  memorySize += wheelsCount * 4;
  const invCache = memorySize >> 2;
  memorySize += wheelsCount * 4;
  const smoothEntriesX = memorySize >> 2;
  memorySize += 512 * 4;//TODO: what size?

  const alpha = memorySize >> 2;
  memorySize += 30 * wheelsCount * 4; //TODO: what size to use?

  const divTestA = memorySize >> 2;
  memorySize += wheelsCount * 4;
  const divTestB = memorySize >> 2;
  memorySize += wheelsCount * 4;

  const exports = instantiateWasm(memorySize);
  const singleBlockSieve = exports.singleBlockSieve;
  const findSmoothEntry = exports.findSmoothEntry;
  const updateWheelsInternalNext = exports.updateWheelsInternalNext;
  const handleSmallWheels = exports.handleSmallWheels;

  const arrayBuffer = exports.memory.buffer;
  const SIEVE_SEGMENT = new Uint8Array(arrayBuffer);
  const heap32 = new Int32Array(arrayBuffer);
  const heap8 = new Uint8Array(arrayBuffer);

  for (let i = 0; i < wheelsCount; i += 1) {
    const w = wheels0[i];
    const wheelLog = i >= realWheels ? 0 : Math.log2(w.p) * (w.step === 2 || w.root === 0 ? 0.5 : 1);
    const log = Math.round(wheelLog * SCALE) | 0;
    if (log >= 2**8 || w.step >= 2**32) {
      throw new RangeError();
    }

    heap32[wheelSteps + i] = w.step;
    heap8[wheelLogs + i] = log;
    heap32[wheelRoots + i] = w.root;
    wheelLogs0[i] = wheelLog;
  }

  const lpStrategy = function (p, polynomial, x, pb) {
    // https://ru.wikipedia.org/wiki/Алгоритм_Диксона#Стратегия_LP
    const lp = largePrimes.get(p);
    if (lp == undefined) {
      // storing polynomial + x has smaller memory usage
      largePrimes.set(p, {polynomial: polynomial, x: x, pb: pb.slice(0)});
      QuadraticSieveFactorization.lpRels += 1;
    } else {
      const s = BigInt(p);
      const sInverse = modInverse(s, N);
      if (sInverse === 0n) {
        return new CongruenceOfsquareOfXminusYmoduloN(s, [0], N);//?
      } else {
        const X = polynomial.X(x);
        const Y = polynomial.Y(x, s, pb);
        const lpX = lp.polynomial.X(lp.x);
        const lpY = lp.polynomial.Y(lp.x, s, lp.pb);
        const X1 = (sInverse * lpX * X) % N;
        if (Y != null && lpY != null) {
          const Y1 = Y.concat(lpY);
          return new CongruenceOfsquareOfXminusYmoduloN(X1, Y1, N);
        } else {
          console.count('bad lp');
        }
      }
    }
    return null;
  };

  const polynomialGenerator = useMultiplePolynomials ? QuadraticPolynomial.generator(sieveSize / 2, primes, N) : null;
  let polynomial = null;
  let baseOffsets = null;
  if (!useMultiplePolynomials) {
    const baseOffset = BigInt(sqrt(N)) + 1n;
    polynomial = new QuadraticPolynomial(1n, baseOffset, N, []);
    baseOffsets = new Float64Array(wheelsCount);
    // - Number(baseOffset % BigInt(pInBeta))
    for (let i = 0; i < wheelsCount; i += 1) {
      baseOffsets[i] = Number(baseOffset % BigInt(heap32[wheelSteps + i])) | 0;
    }
  }

  function checkWheels(offset) {
    for (let k = 0; k < wheelsCount; k += 1) {
      const p = heap32[wheelSteps + k];
      for (let v = 0; v <= 1; v += 1) {
        const root = (v === 0 ? heap32[wheelRoots1 + k] : heap32[wheelRoots2 + k]);
        if (root !== sieveSize) {
          const x = BigInt(+root + offset);
          const X = ((polynomial.useQ2Form ? 2n : 1n) * polynomial.A * x + polynomial.B); // polynomial.X(x)
          const Y = X * X - N;
          if (root < 0 || root >= p || Y % ((polynomial.useQ2Form ? 4n : 1n) * polynomial.A) !== 0n || (Y / ((polynomial.useQ2Form ? 4n : 1n) * polynomial.A)) % BigInt(p) !== 0n) {
            throw new Error();
          }
        }
      }
    }
  }

  /*

  export function singleBlockSieve(wheelRoots1:i32, wheelRoots2:i32, wheelSteps:i32, wheelLogs:i32, startWheel:i32, endWheel:i32, subsegmentEnd:i32, s:i32):i32 {
    const end = wheelSteps + endWheel;
    wheelRoots1 += startWheel;
    wheelRoots2 += startWheel;
    wheelSteps += startWheel;
    wheelLogs += (startWheel >> 2);
    while (wheelSteps !== end) {
      let kpplusr = i32.load(wheelRoots1);
      let kpplusr2 = i32.load(wheelRoots2);
      const step = i32.load(wheelSteps);
      const log2p = i32.load8_u(wheelLogs);
      while (kpplusr2 < subsegmentEnd) {
        const log = i32.load8_u(kpplusr) + log2p;
        const log2 = i32.load8_u(kpplusr2) + log2p;
        i32.store8(kpplusr, log);
        i32.store8(kpplusr2, log2);
        kpplusr += step;
        kpplusr2 += step;
      }
      if (kpplusr < subsegmentEnd) {
        const log = i32.load8_u(kpplusr) + log2p;
        i32.store8(kpplusr, log);
        kpplusr += step;
        const tmp = kpplusr;
        kpplusr = kpplusr2;
        kpplusr2 = tmp;
      }
      i32.store(wheelRoots1, kpplusr - s);
      i32.store(wheelRoots2, kpplusr2 - s);
      wheelRoots1 += 4;
      wheelRoots2 += 4;
      wheelSteps += 4;
      wheelLogs += 1;
    }
    return 0;
  }
  export function findSmoothEntry(thresholdApproximation:i32, i:i32):i32 {
    const t = i8x16.splat(i8(thresholdApproximation));
    i -= 16;
    do {
      i += 16;
    } while (i8x16.bitmask(i8x16.ge_u(v128.load(i), t)) == 0);
    return i;
  }
  export function updateWheelsInternalNext(wheelsCount:i32, alphav:i32, wheelSteps:i32, wheelRoots1:i32, wheelRoots2:i32, e:i32):void {
    const ev = i32x4.splat(e);
    for (let i = 0; i < (wheelsCount << 2); i += 16) {
      let a = v128.load(alphav);
      const p = v128.load(wheelSteps);
      let r1 = v128.load(wheelRoots1);
      let r2 = v128.load(wheelRoots2);
      a = i32x4.add(i32x4.sub(a, ev), v128.and(i32x4.eq(a, ev), p));
      r1 = i32x4.add(i32x4.sub(r1, a), v128.and(i32x4.lt_s(r1, a), p)); // r1 mod p
      r2 = i32x4.add(i32x4.sub(r2, a), v128.and(i32x4.lt_s(r2, a), p)); // r2 mod p
      v128.store(wheelRoots1, i32x4.min_s(r1, r2));
      v128.store(wheelRoots2, i32x4.max_s(r1, r2));
      alphav += 16;
      wheelSteps += 16;
      wheelRoots1 += 16;
      wheelRoots2 += 16;
    }
  }
  @inline
  function reorderBitmask(bm:i32):i32 {
    let bitmask = 0;
    do {
      const i = i32.ctz(bm);
      bm ^= (3 << i);
      const sortedIndex = ((i >> 1) >> 1) | (((i >> 1) & 1) << 2);
      bitmask ^= (1 << sortedIndex);
    } while (bm != 0);
    return bitmask;
  }
  export function handleSmallWheels(wheelsCount:i32, wheelRoots1:i32, wheelRoots2:i32, divTestA:i32, divTestB:i32, storage:i32, smoothEntriesXStart:i32, smoothEntriesXEnd:i32):i32 {
    if (smoothEntriesXStart === smoothEntriesXEnd) {
      unreachable();
    }
    const max = i32x4.splat(-1);
    let k = storage;
    for (let j = 0; j < (wheelsCount << 2); j += 32) {
      const proot1 = v128.or(i32x4.shr_u(i32x4.shl(v128.load(wheelRoots1, 0), 16), 16), i32x4.shl(v128.load(wheelRoots1, 16), 16));
      const proot2 = v128.or(i32x4.shr_u(i32x4.shl(v128.load(wheelRoots2, 0), 16), 16), i32x4.shl(v128.load(wheelRoots2, 16), 16));
      const ta = v128.or(i32x4.shr_u(i32x4.shl(v128.load(divTestA, 0), 16), 16), i32x4.shl(v128.load(divTestA, 16), 16));
      const tb = v128.or(i32x4.shr_u(i32x4.shl(v128.load(divTestB, 0), 16), 16), i32x4.shl(v128.load(divTestB, 16), 16));
      // divisibility test is based on https://math.stackexchange.com/a/1251328 and https://lomont.org/posts/2017/divisibility-testing/ :
      // +working with 16 bit integers to make the test giving false positives, but faster
      let i = smoothEntriesXStart - 4;
      let m = max;
      do {
        i += 4;
        const e = v128.load16_splat(i);
        m = i16x8.min_u(m, i16x8.min_u(i16x8.mul(i16x8.sub(proot1, e), ta), i16x8.mul(i16x8.sub(proot2, e), ta)));
      } while (i < smoothEntriesXEnd);
      if (i8x16.bitmask(i16x8.le_u(m, tb)) != 0) {
        for (let i = smoothEntriesXStart; i < smoothEntriesXEnd; i += 4) {
          const e = v128.load16_splat(i);
          // i32x4.bitmask is faster than v128.any_true
          // i8x16.bitmask is faster than i16x8.bitmask
          let bitmask = i8x16.bitmask(i16x8.le_u(i16x8.min_u(i16x8.mul(i16x8.sub(proot1, e), ta), i16x8.mul(i16x8.sub(proot2, e), ta)), tb));
          if (bitmask != 0) {
            bitmask = reorderBitmask(bitmask);
            const e = i32.load(i);
            do {
              const index = i32.ctz(bitmask);
              const j1 = (index << 2);
              bitmask ^= (1 << index);
              const proot1 = i32.load(wheelRoots1 + j1);
              const proot2 = i32.load(wheelRoots2 + j1);
              const ta = i32.load(divTestA + j1);
              const tb = i32.load(divTestB + j1);
              if (u32(((proot1 - e) * ta) & 0xFFFFFFFF) <= u32(tb) ||
                  u32(((proot2 - e) * ta) & 0xFFFFFFFF) <= u32(tb)) {
                if (proot1 !== 0 || proot2 !== 0) {
                  i32.store(k << 2, (j + j1) >> 2);
                  i32.store((k + 1) << 2, (i - smoothEntriesXStart) >> 2);
                  k += 2;
                  if (proot1 === proot2) {
                    i32.store(k << 2, (j + j1) >> 2);
                    i32.store((k + 1) << 2, (i - smoothEntriesXStart) >> 2);
                    k += 2;
                  }
                }
              }
            } while (bitmask != 0);
          }
        }
      }
      wheelRoots1 += 32;
      wheelRoots2 += 32;
      divTestA += 32;
      divTestB += 32;
    }
    return k - storage;
  }

  */

  function modInverseMod2pow32(a) {
    let aInv = a;
    for (let e = 2; e < 32; e += e) {
      aInv = Math.imul(aInv, 2 - Math.imul(a, aInv));
    }
    return aInv;
  }

  function computeDivTestAB(d, Nmax) {
    // see https://math.stackexchange.com/a/1251328 and https://lomont.org/posts/2017/divisibility-testing/
    if (d % 2 === 0) {
      if (d !== 2**Math.round(Math.log2(d))) {
        throw new RangeError();
      }
      // power of two
      return {a: (2**32 / d) | 0, b: 0};
    }
    const a = modInverseMod2pow32(d);
    let b = Math.floor(Nmax / d);
    b = b < 2**16 ? b : (b | 0xFFFF);
    console.assert(b <= Math.floor((2**32 - 1) / d));
    return {a: a, b: b};
  }

  //console.log(computeDivTestAB(5)); - {a: -858993459, b: 858993459}

  for (let i = 0; i < wheelsCount; i += 1) {
    const p = heap32[wheelSteps + i];
    const tmp = computeDivTestAB(p, wheels0.at(-1).step + sieveSize);
    heap32[divTestA + i] = tmp.a;
    heap32[divTestB + i] = tmp.b;
  }

  let BPrev = 0n;
  let Bdistances = [];
  let counternext = 0;

  const updateWheels = function (polynomial, offset) {
    offset = -0 + offset;
    //recalculate roots based on the formula:
    //proot = ((-B + root) * modInv(A, p)) % p;
    //+some optimizations to minimize bigint usage and modInverseSmall calls
    const useCache = BigInt(polynomial.A) === BigInt(invCacheKey);
    if (!useCache) {
      const AA = FastModBigInt((polynomial.useQ2Form ? 2n : 1n) * polynomial.A);
      const bsign = polynomial.B < 0n ? -1 : +1;
      const BB = FastModBigInt(polynomial.B < 0n ? -polynomial.B : polynomial.B);
      //  first one: ((-B + root) * modInv(A, p) - offset) % p
      //  next proots: (-(B - Bprev) * modInv(A, p) + prootPrev) % p, where (-(B - Bprev) * modInv(A, p)) mod p is cached
      zeroInvs.length = 0;
      for (let i = 0; i < wheelsCount; i += 1) {
        const p = -0 + heap32[wheelSteps + i];
        const pInv = (1 + 2**-52) / p;
        //const a = Number(polynomial.A % BigInt(p));
        const a = -0 + FastMod(AA, p, pInv);
        const invA = modInverseSmall(a, p) | 0;
        heap32[invCache + i] = invA;
        if (invA === 0) {
          zeroInvs.push(i);
        }
        // find roots for the first polynomial:
        let b = -0 + FastMod(BB, p, pInv);
        b = bsign < 0 ? -b : b;
        const root = -0 + heap32[wheelRoots + i];
        let x1 = ((p - b) + (p - root)) * invA - offset;
        let x2 = ((p - b) + root) * invA - offset;
        x1 = x1 - Math.floor(x1 * pInv) * p; // x1 mod p
        x2 = x2 - Math.floor(x2 * pInv) * p; // x2 mod p
        const r1 = Math.min(x1, x2);
        const r2 = Math.max(x1, x2);
        heap32[wheelRoots1 + i] = r1;
        heap32[wheelRoots2 + i] = r2;
      }
      BPrev = 0n;
      //console.log('Bdistances.length', Bdistances.length, counternext);
      Bdistances.length = 0;
      counternext = 0;
    }
    if (BPrev === 0n) {
      BPrev = polynomial.B;
    } else {
      let d = polynomial.B - BPrev;
      BPrev = polynomial.B;
      let e = 0;
      if (d < 0n) {
        d += (polynomial.useQ2Form ? 2n : 1n) * polynomial.A;
        e = 1;
      }
      if (d < 0n || d >= polynomial.A * 2n) {
        throw new RangeError();
      }
      let v = Bdistances.indexOf(d);
      if (v === -1) {
        Bdistances.push(d);
        const dd = FastModBigInt(d < 0n ? -d : d);
        v = Bdistances.length - 1;
        const alphav = alpha + v * wheelsCount;
        for (let i = 0; i < wheelsCount; i += 1) {
          const p = -0 + heap32[wheelSteps + i];
          const invA = -0 + heap32[invCache + i];
          const pInv = (1 + 2**-52) / p;
          const d = FastMod(dd, p, pInv);
          let a = (p - d) * invA;
          a = a + sieveSize; // correction
          a = a - Math.floor(a * pInv) * p;
          a = p - a;
          heap32[alphav + i] = a;
        }
      }
      counternext += 1;
      updateWheelsInternalNext(wheelsCount, (alpha + v * wheelsCount) << 2, wheelSteps << 2, wheelRoots1 << 2, wheelRoots2 << 2, e);
    }
    const a = Number(polynomial.A & BigInt(2**29 - 1));
    const b = Number(polynomial.B & BigInt(2**29 - 1));
    const c = Number(polynomial.C & BigInt(2**29 - 1));
    const aInv = modInverseSmall(a, 2**29);
    for (let j = 0; j < zeroInvs.length; j += 1) {
      const i = zeroInvs[j];
      // single root:
      // x = (2B)^-1*(-C) (mod p)
      // skip as the performance is not better
      heap32[wheelRoots1 + i] = sieveSize;
      heap32[wheelRoots2 + i] = sieveSize;
      
      if (polynomial.useQ2Form) {
        const p = heap32[wheelSteps + i];
        if (p % 2 === 0) {
          const r0 = heap32[wheelRoots + i];
          for (let k = 0; k < (p === 2 ? 1 : 2); k += 1) {
            const r = k === 0 ? r0 : p - r0;
            console.assert((-b - r) % 2 === 0);
            console.assert((-b + r) % 2 === 0);
            let r1 = (Math.imul((-b - r) / 2, aInv) - offset) & (p - 1);
            let r2 = (Math.imul((-b + r) / 2, aInv) - offset) & (p - 1);
            let y1 = (Math.imul(Math.imul(a, r1 + offset) + b, r1 + offset) + c) & (p - 1);
            let y2 = (Math.imul(Math.imul(a, r2 + offset) + b, r2 + offset) + c) & (p - 1);
            if (y1 === 0 && y2 === 0) {
              heap32[wheelRoots1 + i] = Math.min(r1, r2);
              heap32[wheelRoots2 + i] = Math.max(r1, r2);
            }
            if (p === 2) {
              // for (let i = 0; i < p; i++) { if (polynomial.Y(i) % BigInt(p) === 0n) { console.log('root', i+offset%p); } } console.log('done');
              if (r1 !== r2) {
                heap8[wheelLogs + i] = Math.round(1 * SCALE) | 0;
                wheelLogs0[i] = 1;
              } else {
                throw new Error();
              }
            }
          }
        }
      }
    }
    //...
    invCacheKey = polynomial.A;
    //checkWheels(offset);
  };

  const gcd = function (a, b) {
    while (b !== 0) {
      const r = +a % +b;
      a = b;
      b = r;
    }
    return a;
  };
  const lcm = function (a, b) {
    return Math.floor(a / gcd(a, b)) * b;
  };

  const getSmallWheels = function () {
    let p = 1;
    let i = 0;
    while (i < wheelsCount && lcm(p, heap32[wheelSteps + i]) <= segmentSize / 4) {
      p = lcm(p, heap32[wheelSteps + i]);
      i += 1;
    }
    return i;
  };
  const smallWheels = getSmallWheels();

  const copyCycle = function (array, cycleLength, limit) {
    if (typeof limit !== 'number' || typeof cycleLength !== 'number') {
      throw new TypeError();
    }
    if (limit > array.length) {
      throw new RangeError();
    }
    for (let i = cycleLength; i < limit; i += cycleLength) {
      array.copyWithin(i, 0, Math.min(limit - i, cycleLength));
    }
  };

  QuadraticSieveFactorization.smallSegmentTime = 0;
  QuadraticSieveFactorization.largeSegmentTime = 0;

  const updateSieveSegment = function (segmentStart) {
    let cycleLength = 1;
    SIEVE_SEGMENT[0] = SHIFT;
    for (let j = 0; j < smallWheels; j += 1) {
      const newCycleLength = +lcm(cycleLength, heap32[wheelSteps + j]);
      copyCycle(SIEVE_SEGMENT, cycleLength, newCycleLength);
      cycleLength = newCycleLength;
      const p = heap32[wheelSteps + j];
      const log2p = heap8[wheelLogs + j];
      const r1 = (heap32[wheelRoots1 + j] | 0);
      if (r1 !== sieveSize) {
        for (let k = (r1 + newCycleLength - segmentStart % newCycleLength) % p; k < newCycleLength; k += p) {
          SIEVE_SEGMENT[k] = (SIEVE_SEGMENT[k] + log2p) | 0;
        }
      }
      const r2 = (heap32[wheelRoots2 + j] | 0);
      if (r2 !== sieveSize) {
        for (let k = (r2 + newCycleLength - segmentStart % newCycleLength) % p; k < newCycleLength; k += p) {
          SIEVE_SEGMENT[k] = (SIEVE_SEGMENT[k] + log2p) | 0;
        }
      }
    }
    copyCycle(SIEVE_SEGMENT, cycleLength, segmentSize);
    //for (let j = 0; j < segmentSize; j += 1) {
    //  SIEVE_SEGMENT[j] = SHIFT;
    //}
    // "Block Sieving Algorithms" by Georg Wambach and Hannes Wettig May 1995
    const m = (typeof navigator !== 'undefined' && navigator.hardwareConcurrency === 12 ? 1 : 1.5);
    const V = Math.min(0 + wheelsCount - smallWheels, Math.floor(64 * 3 * m * (N > 2**285 ? 4 : (N > 2**240 ? 2 : 1))));
    const S = Math.floor((1 << 15) * m);
    const t1 = performance.now();
    let subsegmentEnd = 0;
    while (subsegmentEnd + S <= segmentSize) {
      subsegmentEnd += S;
      singleBlockSieve(wheelRoots1 * 4, wheelRoots2 * 4, wheelSteps * 4, wheelLogs * 1, smallWheels * 4, smallWheels * 4 + V * 4, subsegmentEnd, 0);
    }
    QuadraticSieveFactorization.smallSegmentTime += performance.now() - t1;
    const t2 = performance.now();
    singleBlockSieve(wheelRoots1 * 4, wheelRoots2 * 4, wheelSteps * 4, wheelLogs * 1, smallWheels * 4, wheelsCount * 4, segmentSize, segmentSize);
    QuadraticSieveFactorization.largeSegmentTime += performance.now() - t2;
  };

  const smoothEntryPositions = [];
  const smoothEntryValues = [];
  const smoothEntryFactorBases = [];

  const findSmoothEntries = function (offset, polynomial) {
    if (typeof offset !== "number") {
      throw new TypeError();
    }
    
    const commonThreshold = Math.max(SHIFT, (Math.round((polynomial.log2AbsY(0) - (doubleLargePrimes ? doubleLargePrimesThreshold : largePrimesThreshold)) * SCALE + SHIFT) | 0) - 9);

    let i = 0;
    let thresholdApproximation = 0;
    while (i < segmentSize) {
      // it is slow to compute the threshold on every iteration, so trying to optimize:

      //TODO: the threshold calculation is much more simple in the Youtube videos (?)
      thresholdApproximation = useMultiplePolynomials ? commonThreshold : Math.floor((polynomial.log2AbsY(i + offset) - (doubleLargePrimes ? doubleLargePrimesThreshold : largePrimesThreshold)) * SCALE + SHIFT + 0.5) | 0;
      const j = useMultiplePolynomials ? segmentSize : Math.min(segmentSize, thresholdApproximationInterval(polynomial, i + offset, (thresholdApproximation - SHIFT) * (1 / SCALE) + (doubleLargePrimes ? doubleLargePrimesThreshold : largePrimesThreshold), sieveSize) - offset);

      while (i < j) {
        if (i < j - 1) {
          const tmp = SIEVE_SEGMENT[j - 1];
          SIEVE_SEGMENT[j - 1] = MAX;
          while ((i & 0xF) !== 0 && thresholdApproximation > SIEVE_SEGMENT[i]) {
            i += 1;
          }
          if ((i & 0xF) === 0) {
            i = findSmoothEntry(thresholdApproximation, i);
            var c = 0;
            while (thresholdApproximation > SIEVE_SEGMENT[i]) {
              i += 1;
              c += 1;
            }
            if (c >= 16) {
              console.error('too big c', c);
            }
          }
          SIEVE_SEGMENT[j - 1] = tmp;
        }
        const r = polynomial.log2AbsY(i + offset) - (SIEVE_SEGMENT[i] - SHIFT) * (1 / SCALE);
        if (r < largePrimesThreshold || doubleLargePrimes && r >= log2B * 2 && r < doubleLargePrimesThreshold) {
          smoothEntryPositions.push(i + offset);
          smoothEntryValues.push(-0 + (SIEVE_SEGMENT[i] - SHIFT) * (1 / SCALE));
        }
        i += 1;
      }
    }
  };

  function checkFactorization(x) {
    let p = 0;
    for (let n = 0; n < wheelsCount; n += 1) {
      const log2p = heap8[wheelLogs + n];
      const step = heap32[wheelSteps + n];
      for (let v = 0; v <= 1; v += 1) {
        if ((x - (v === 0 ? (heap32[wheelRoots1 + n] | 0) : (heap32[wheelRoots2 + n] | 0)) - (n < smallWheels ? 0 : segmentSize)) % step === 0) {
          if (polynomial.AFactors.indexOf(step) === -1) {
            console.log(step);
            p += log2p;
          }
        }
      }
    }
    return p;
  }

  function applyOffset(offset) {
    for (let j = 0; j < wheelsCount; j += 1) {
      const step = heap32[wheelSteps + j];
      let r1 = (0 + (heap32[wheelRoots + j] | 0) - baseOffsets[j] - offset) % step;
      r1 += (r1 < 0 ? step : 0);
      let r2 = (0 - (heap32[wheelRoots + j] | 0) - baseOffsets[j] - offset) % step;
      r2 += (r2 < 0 ? step : 0);
      heap32[wheelRoots1 + j] = Math.min(r1, r2);
      heap32[wheelRoots2 + j] = Math.max(r1, r2);
    }
  }

  function findPreciseSmoothEntries(offset) {
    if (smoothEntryPositions.length > 512) {
      console.warn('too many smooth entries: ' + smoothEntryPositions.length + ' N = ' + N);
      smoothEntryPositions.length = 512;
    }
    for (let i = 0; i < smoothEntryPositions.length; i += 1) {
      heap32[smoothEntriesX + i] = (smoothEntryPositions[i] - offset) - sieveSize;
    }
    for (let j = 0; j < smallWheels; j += 1) {
      const step = heap32[wheelSteps + j];
      if (heap32[wheelRoots1 + j] !== sieveSize) {
        let x = (heap32[wheelRoots1 + j] + (0 - sieveSize) % step);
        heap32[wheelRoots1 + j] = x <= 0 ? x + step : x;
      } else {
        heap32[wheelRoots1 + j] = 0;
      }
      if (heap32[wheelRoots2 + j] !== sieveSize) {
        let x = (heap32[wheelRoots2 + j] + (0 - sieveSize) % step);
        heap32[wheelRoots2 + j] = x <= 0 ? x + step : x;
      } else {
        heap32[wheelRoots2 + j] = 0;
      }
    }
    //console.log('smoothEntryPositions.length', smoothEntryPositions.length);
    const k = smoothEntryPositions.length === 0 ? 0 : handleSmallWheels(wheelsCount, wheelRoots1 << 2, wheelRoots2 << 2, divTestA << 2, divTestB << 2, storage, smoothEntriesX << 2, (smoothEntriesX + smoothEntryPositions.length) << 2);

    const preciseValues = new Float64Array(smoothEntryPositions.length);
    for (let v = 0; v < k; v += 2) {
      const j = heap32[storage + v];
      const i = heap32[storage + v + 1];
      const step = heap32[wheelSteps + j];
      preciseValues[i] += +wheelLogs0[j];
      smoothEntryFactorBases[i].push(step);
    }
    for (let i = 0; i < smoothEntryPositions.length; i += 1) {
      const e = Math.abs(smoothEntryValues[i] - preciseValues[i]);
      if (e >= 9 && e < 100) {
        console.error(e);
      }
      smoothEntryValues[i] = preciseValues[i];
    }
  }

  // Double Large Primes:
  
  function Queue() {
    // Queue using two arrays - https://stackoverflow.com/a/73957258
    this.a = [];
    this.b = [];
    this.length = 0;
  }
  Queue.prototype.push = function (e) {
    this.a.push(e);
    this.length = this.a.length + this.b.length;
  };
  Queue.prototype.shift = function (e) {
    if (this.b.length === 0) {
      while (this.a.length > 0) {
        this.b.push(this.a.pop());
      }
    }
    var e = this.b.pop();
    this.length = this.a.length + this.b.length;
    return e;
  };

  function Graph() {
    this.edges = 0;
    this.components = 0;
    this.vertices = 0;
    this._edges = new Map();
    this._edgesArray = [];
    this._g = new Map();
  }
  Graph.prototype._insertVertex = function (p) {
    if (!this._g.has(p)) {
      this._g.set(p, p);
      this.vertices += 1;
      this.components += 1;
      return p;
    }
    let root = this._g.get(p);
    while (this._g.get(root) !== root) {
      root = this._g.get(root);
    }
    let x = p;
    while (x !== root) {
      const next = this._g.get(x);
      this._g.set(x, root);
      x = next;
    }
    return root;
  };
  const key = function (p1, p2) {
    return BigInt.asUintN(64, (BigInt(Math.min(p1, p2)) << 32n) + BigInt(Math.max(p1, p2)));
  };
  Graph.prototype.insertEdge = function (p1, p2, data) {
    const p1p2 = key(p1, p2);
    if (this._edges.get(p1p2) !== undefined) {
      console.count('same p1p2');
      return;
    }
    this._edges.set(p1p2, data);
    this._edgesArray.push({p1: p1, p2: p2});
    if (p1 !== p2) {
      this._edgesArray.push({p1: p2, p2: p1});
    }
    this.edges += 1;
    const p1root = this._insertVertex(p1);
    const p2root = this._insertVertex(p2);
    if (p1root < p2root) {
      this.components -= 1;
      this._g.set(p2root, p1root);
    }
    if (p1root > p2root) {
      this.components -= 1;
      this._g.set(p1root, p2root);
    }
  };
  Graph.prototype.iterateCycles = function (onCycle) {
    const nodes = new Map();
    this._edgesArray.sort((a, b) => a.p1 - b.p1 || a.p2 - b.p2);

    const roots = [];
    for (let v of this._g.keys()) {
      if (this._g.get(v) === v) {
        roots.push(v);
      }
      nodes.set(v, {explored: false, finished: false, parent: 0});
    }
    roots.sort((a, b) => a - b);
    const p1Array = this._edgesArray.map(edge => edge.p1);

    const path = function (w) {
      const p = [];
      while (w !== 0) {
        p.push(w);
        w = nodes.get(w).parent;
      }
      return p;
    };

    for (const r of roots) {
      // https://en.wikipedia.org/wiki/Breadth-first_search#Pseudocode
      const Q = new Queue();
      const root = nodes.get(r);
      console.assert(root.explored === false && root.parent === 0);
      root.explored = true;
      root.parent = 0;
      Q.push(r);
      while (Q.length > 0) {
        const v = Q.shift();
        const start = indexOf(p1Array, v);
        //console.assert(p1Array.lastIndexOf(v) === start);
        for (let i = start; i >= 0 && this._edgesArray[i].p1 === v; i -= 1) {
          const edge = this._edgesArray[i];
          const p2 = nodes.get(edge.p2);
          if (p2.explored) {
            if (p2.finished) {
              continue;
            }
            // cycle
            const path1 = path(edge.p2);
            const path2 = path(edge.p1);
            let j = 0;
            while (path1.length - 1 - j >= 0 && path2.length - 1 - j >= 0 && path1[path1.length - 1 - j] === path2[path2.length - 1 - j]) {
              j += 1;
            }
            console.assert(j > 0);
            const vertices = path1.slice(0, path1.length - (j - 1)).reverse().concat(path2.slice(0, path2.length - (j - 1)));
            const edges = [];
            for (let i = 1; i < vertices.length; i += 1) {
              edges.push({p1: vertices[i - 1], p2: vertices[i]});
            }
            const edgesMap = this._edges;
            const data = edges.map(e => edgesMap.get(key(e.p1, e.p2)));
            onCycle(edges, data);
          } else {
            p2.explored = true;
            p2.parent = edge.p1; // v
            Q.push(edge.p2);
          }
        }
        nodes.get(v).finished = true;
      }
    }
  };

  const onCycle = function (edges, data) {
    let s = 1n;
    let X = 1n;
    let Y = [];
    for (let j = 0; j < edges.length; j += 1) {
      const lp = edges[j];
      const lpData = data[j];
      const lpX = lpData.polynomial.X(lpData.x);
      const lpY = lpData.polynomial.Y(lpData.x, BigInt(lp.p1) * BigInt(lp.p2), lpData.pb);
      X = (lpX * X) % N;
      if (Y == null || lpY == null) {
        Y = null;
      } else {
        Y = Y.concat(lpY);
      }
      s = (s * BigInt(lp.p2)) % N;
    }
    const sInverse = modInverse(s, N);
    if (sInverse === 0n) {
      foundGraphRelations.push(new CongruenceOfsquareOfXminusYmoduloN(s, [0], N));//?
    } else if (Y != null) {
      X = (sInverse * X) % N;
      const c = new CongruenceOfsquareOfXminusYmoduloN(X, Y, N);
      if ((X**2n - Y.reduce((p, a) => p * BigInt(a), 1n)) % N !== 0n) {
        throw new Error();
      }
      foundGraphRelations.push(c);
    } else {
      console.count('BAD CYCLE');
    }
  };

  let graph = new Graph();
  const foundGraphRelations = [];
  let total = 0;

  function handleDoubleLargePrimeNext(p1, p2, x, polynomial, pb) {
    if (graph == null) {
      return;
    }
    graph.insertEdge(p1, p2, {polynomial: polynomial, x: x, pb: pb.slice(0)});
    const cyclesCount = graph.edges + graph.components - graph.vertices;
    if (graph.edges % 10000 === 0) {
      console.debug('graph:', graph.edges, graph.components, graph.vertices, cyclesCount);
    }
    if (cyclesCount > primes.length + 64 - total) {//TODO: !?
      console.debug('graph:', graph.edges, graph.components, graph.vertices, cyclesCount);
      console.time('collectGraphRelations');
      graph.iterateCycles(onCycle);
      console.timeEnd('collectGraphRelations');
      console.log(foundGraphRelations.length);
      graph = null;
    }
  }

  function handleDoubleLargePrime(x, polynomial, pb) {
    let Y = abs(polynomial.Y(x));
    for (let j = pb.length - 1; j >= 0; j -= 1) { // backwards is slighly faster
      const p = pb[j];
      while (Y % BigInt(p) === 0n) {
        Y /= BigInt(p);
      }
    }
    const r = BigInt(Y);
    if (Number(r) >= 2**64) {
      //console.count('2**64');
      return;
    }
    if (isProbablyPrime64(r)) {
      //console.count('prime');
      return;
    }
    const f = Number(PollardsRho64(r, 2**18));
    if (f === 0 || f === 1) {
      console.count('cannot factor');
      return;      
    }
    if (BigInt(r) % BigInt(f) !== BigInt(0) || !(f > 1 && BigInt(f) < BigInt(r))) {
      throw new Error();
    }
    const f2 = Math.floor(Number(BigInt(r) / BigInt(f)));
    const p1 = Math.min(f2, f);
    const p2 = Math.max(f2, f);
    if (!(Math.log2(p1) > log2B && p1 <= 1073741823 && p2 <= 1073741823)) {
      //console.count('invalid values');
      return;
    }
    handleDoubleLargePrimeNext(p1 | 0, p2 | 0, x, polynomial, pb);
  }

  QuadraticSieveFactorization.lpCounter = 0;
  QuadraticSieveFactorization.lpRels = 0;
  let i1 = -1;
  let k = 0;
  const iterator = {
    next: function congruencesUsingQuadraticSieve() {
      total += 1;
      while (k < 32 || (useMultiplePolynomials ? 2 : 1/16) * k * sieveSize <= Math.pow(primes[primes.length - 1], 2)) {
        if (foundGraphRelations.length !== 0) {
          return {value: foundGraphRelations.pop(), done: false};
        }

        if (i1 === sieveSize) {
          k += 1;
          i1 = -1;
        }
        const offset = useMultiplePolynomials ? -sieveSize / 2 : (k % 2 === 0 ? 1 : -1) * Math.floor((k + 1) / 2) * sieveSize;
        if (i1 === -1) {

          if (useMultiplePolynomials) {
            polynomial = polynomialGenerator.next();
            updateWheels(polynomial, offset);
          } else {
            applyOffset(offset);
          }

          smoothEntryPositions.length = 0;
          smoothEntryValues.length = 0;

          for (let segmentStart = 0; segmentStart < sieveSize; segmentStart += segmentSize) {
            updateSieveSegment(segmentStart);
            findSmoothEntries(offset + segmentStart, polynomial);
          }

          for (let i = 0; i < smoothEntryPositions.length; i += 1) {
            if (i < smoothEntryFactorBases.length) {
              smoothEntryFactorBases[i].length = 0;
            } else {
              smoothEntryFactorBases.push([]);
            }
          }

          findPreciseSmoothEntries(offset);
        }

          //Note: separate loop over "smooth entries" is better for performance, seems
          for (let i = i1 + 1; i < smoothEntryPositions.length; i += 1) {
            const x = smoothEntryPositions[i];
            const value = +smoothEntryValues[i];
            const pb = smoothEntryFactorBases[i];
            const threshold = +polynomial.log2AbsY(x);
            let c = null;
            if (threshold - value <= log2B) {
              const X = polynomial.X(x);
              let Y = polynomial.Y(x, 1n, pb);
              if (Y == null) {
                // this may happen because of prime powers
                // or primes used to make "polynomial.A"
                Y = polynomial.Y(x, 1n, pb.concat(zeroInvs.map(i => heap32[wheelSteps + i])));
              }
              if (Y != null) {
                c = new CongruenceOfsquareOfXminusYmoduloN(X, Y, N);
              } else {
                // may happen when exp2(threshold - value) is a multiplier
                console.count('wrong entry', exp2(threshold - value));
                //console.log(threshold, value, checkFactorization(x - offset));
              }
            } else if (threshold - value < 2 * log2B) {
              const p = exp2(threshold - value);
              //if (!isProbablyPrime(p)) {
                //console.debug('wrong large prime?', p);
              //}
              c = lpStrategy(p, polynomial, x, pb);
              if (c != null) {
                QuadraticSieveFactorization.lpCounter += 1;
              } else {
                if (doubleLargePrimes && p <= 1073741823) {
                  if (polynomial.Y(x, BigInt(p), pb) != null) {
                    handleDoubleLargePrimeNext(1, p | 0, x, polynomial, pb);
                  } else {
                    console.count('wrong?');
                    console.log(polynomial.Y(x), p);
                  }
                }
              }
            } else if (doubleLargePrimes && threshold - value < 3 * log2B) {
              handleDoubleLargePrime(x, polynomial, pb);
            } else {
              console.count('too big', (threshold - value) / log2B);
            }
            if (c != null) {
              i1 = i;
              return {value: c, done: false};
            }
          }
        i1 = sieveSize;
      }
      return {value: undefined, done: true};
    }
  };
  iterator[globalThis.Symbol.iterator] = function () {
    return this;
  };
  return iterator;
}

function gcd(a, b) {
  while (b !== 0n) {
    const r = BigInt(a) % BigInt(b);
    a = b;
    b = r;
  }
  return a;
}

function abs(x) {
  if (typeof x !== 'bigint') {
    throw new TypeError();
  }
  return x < 0n ? -x : x;
}

function indexOf(sortedArray, x) {
  if (typeof x !== 'number' || (x | 0) !== x) {
    throw new TypeError();
  }
  let min = 0;
  let max = sortedArray.length - 1;
  while (min < max) {
    const mid = Math.ceil((min + max) / 2);
    if ((sortedArray[mid] | 0) > (x | 0)) {
      max = mid - 1;
    } else {
      min = mid;
    }
  }
  if ((sortedArray[min] | 0) === (x | 0)) {
    return min;
  }
  return -1;
}

function computeY(primeBase, solution, N) {
  const Y = new Array(primeBase.length + 1).fill(0);
  for (let i = 0; i < solution.length; i += 1) {
    const v = solution[i].v;
    for (let j = 0; j < v.length; j += 1) {
      Y[v[j]] += 1;
    }
  }
  let y = 1n;
  for (let i = 0; i < Y.length; i += 1) {
    if (Y[i] % 2 !== 0) {
      throw new RangeError();
    }
    if (i !== 0) {
      const p = primeBase[i - 1];
      const e = Y[i] / 2;
      if (e > 0) {
        if (e <= 2) {
          y = (y * BigInt(Math.pow(p, e))) % N;
        } else {
          y = (y * modPow(BigInt(p), BigInt(e), N)) % N;
        }
      }
    }
  }
  return y;
}


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

function getBestMultiplier(n, primesList) {
  // https://ir.cwi.nl/pub/27839/27839.pdf

  // f (m, n) = --2 - + L g(p, mn) logp,
  const scores = new Array(150).fill(0);
  for (let m = 1; m < scores.length; m += 1) {
    scores[m] = -Math.log(m) / 2;
  }

  for (let i = 0; i < primesList.length && i < 300; i += 1) {
    const p = primesList[i];
    if (p === 2) {
      for (let m = 1; m < scores.length; m += 1) {
        const q2formExtraScore = 1;
        scores[m] += [0, 2 + q2formExtraScore, 0, 0.5, 0, 1, 0, 0.5][(m * Number(n % 8n)) % 8] * Math.log(2);
      }
    } else {
      const lnp = legendre(Number(n % BigInt(p)), p);
      const cp = 2 / (p - 1) * Math.log(p);
      for (let m = 1; m < scores.length; m += 1) {
        scores[m] += (lnp * legendre(m, p) === 1 ? cp : 0);
      }
    }
  }

  let max = 0;
  let best = 1;
  for (let m = 1; m <= scores.length; m += 1) {
    var y = +scores[m];
    if (+y > +max) {
      max = y;
      best = m;
    }
  }

  //console.log('best: ', best, 'scores: ', scores);
  return best;
}

function QuadraticSieveFactorization(N) { // N - is not a prime
  if (typeof N !== 'bigint') {
    throw new TypeError();
  }
  // https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html#:~:text=optimal%20value :
  // to limit memory usage during "solve" to 2GB:
  const memoryBasedLimit = Math.floor(((performance.memory || {jsHeapSizeLimit: 0}).jsHeapSizeLimit || 0) / 2**32 < 0.99 ? 2**23.5 : 2**23.75);
  const limit = Math.min(memoryBasedLimit, (1 << 25) - 1);
  const dlpMultiplier = 2.5;//TODO: when using double large primes
  const B = Math.max(Math.min(Math.floor(Math.sqrt(L(N) / (Number(N) > 2**285 ? 24 * dlpMultiplier : (Number(N) > 2**240 ? 12 : 6)))), limit), 512);
  const primesList = primes(B);
  let k = 1n;
  k = Number(N) > 2**64 ? BigInt(getBestMultiplier(N, primesList)) : 1n;
  for (;; k += 1n) {
    if (k !== 1n) {
      console.debug('multiplier', k);
    }
    if (k !== 1n && N % k === 0n) {
      console.error('k');
      return BigInt(k);
    }
    const kN = k * N;
    const primeBase = primesList.filter(p => isQuadraticResidueModuloPrime(kN, p));
    for (let i = 0; i < primeBase.length; i += 1) {
      if (Number(N % BigInt(primeBase[i])) === 0) {
        return BigInt(primeBase[i]);
      }
    }
    const congruences = congruencesUsingQuadraticSieve(primeBase, kN); // congruences X_k^2 = Y_k mod N, where Y_k is smooth over the prime base
    const solutions = solve.sparseSolve(1 + primeBase.length); // find products of Y_k = Y, so that Y is a perfect square
    solutions.next();
    let c = null;
    const start = performance.now();
    let congruencesFound = 0;
    let last = start;
    while ((c = congruences.next().value) != undefined) {
      if (c.Y.length === 1 && c.Y[0] === 0) {
        const g = BigInt(gcd(abs(c.X), N));
        if (g !== 1n && g !== N) {
          return g;
        }
      } else {
        const t = function () {
          throw new TypeError(N);
        };
        const v = c.Y.map(p => (p === -1 ? 0 : 1 + indexOf(primeBase, p) || t()));
        const solution = solutions.next([v, {c: c, v: v}]).value;
        if (true) {
          const now = performance.now();
          congruencesFound += 1;
          if (false && congruencesFound % 50 === 0) {
            console.debug('smallSegmentTime: ' + QuadraticSieveFactorization.smallSegmentTime,
                          'largeSegmentTime: ' + QuadraticSieveFactorization.largeSegmentTime);
            return 1n;
          }
          if (now - last > 5000 || solution != null) {
            console.debug('congruences found: ', congruencesFound, '/', primeBase.length,
                          'expected time: ', Math.round((now - start) / congruencesFound * primeBase.length),
                          'large prime congruences: ', QuadraticSieveFactorization.lpCounter + '(' + QuadraticSieveFactorization.lpRels + ')',
                          'polynomials used: ', QuadraticSieveFactorization.polynomialsCounter);
            last = now;
          }
        }
        if (solution != null) {
          let x = 1n;
          for (let i = 0; i < solution.length; i += 1) {
            x = (x * solution[i].c.X) % N;
          }
          // we cannot just compute product as it is larger 2**(2**20) (max BigInt in Firefox)
          let y = computeY(primeBase, solution, N); // Y mod N === X^2 mod N
          const g = BigInt(gcd(abs(x + y), N));
          if (g !== 1n && g !== N) {
            return g;
          }
        }
      }
    }
  }
}

QuadraticSieveFactorization.testables = {
  congruencesUsingQuadraticSieve: congruencesUsingQuadraticSieve,
  squareRootModuloOddPrime: squareRootModuloOddPrime,
  isQuadraticResidueModuloPrime: isQuadraticResidueModuloPrime,
  solve: solve,
  QuadraticPolynomial: QuadraticPolynomial,
  thresholdApproximationInterval: thresholdApproximationInterval
};

export default QuadraticSieveFactorization;