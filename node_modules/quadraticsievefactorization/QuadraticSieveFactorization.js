/*jshint esversion:6*/
import isPrime from './isPrime.js';

function modInverse(a, m) {
  a = BigInt(a);
  m = BigInt(m);
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
  a = Number(a);
  m = Number(m);
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

function getSquareRootsModuloPrimeBig(n, p, e = 1) {
  n = BigInt(n);
  p = BigInt(p);
  e = Number(e);
  const m = e === 1 ? p : (e === 2 ? p * p : 0n); // p**BigInt(e)
  n %= m;
  if ((p + 1n) % 4n === 0n) {
    if (e !== 1) {
      const x = getSquareRootsModuloPrimeBig(n, p, e - 1)[0];
      let x1 = (x + (n - (x * x) % m) * modInverse(x + x, m)) % m;
      if (x1 < 0n) {
        x1 += m;
      }
      if (x1 > m - x1) {
        x1 = m - x1;
      }
      return [x1, m - x1];
    }
    // from https://en.wikipedia.org/wiki/Quadratic_residue#Prime_or_prime_power_modulus :
    let r = modPow(n, (p + 1n) / 4n, p);
    if ((r * r) % p === n) {
      if (r > p - r) {
        r = p - r;
      }
      return [r, p - r];
    }
    return [];
  }
  throw new RangeError('implemented only for primes of the form 4k+3');
}

function getSquareRootsModuloPrime(n, p, e = 1) { // slow for non-small p
  n = Number(n);
  p = Number(p);
  e = Number(e);
  const m = Math.pow(p, e);
  n = n % m;
  if (!(n > 0 && p > 0 && e >= 1 && n % p !== 0 && m < Math.floor(Math.sqrt(Number.MAX_SAFE_INTEGER)))) { // + p is a prime number
    throw new RangeError();
  }
  // r**2 == n (mod p)
  if (e > 1) {
    if (p === 2) {
      if (e >= 3) {
        if (n % 8 === 1) { // from Cohen H.
          const candidates = getSquareRootsModuloPrime(n, p, e - 1);
          let i = 0;
          while ((candidates[i] * candidates[i]) % m !== n) {
            i += 1;
          }
          const r = candidates[i];
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
    } else {
      const roots = getSquareRootsModuloPrime(n, p, e - 1);
      if (roots.length === 0) {
        return [];
      }
      if (roots.length !== 2) {
        throw new Error();
      }
      const x = roots[0];
      // x**2 = n mod p**(e - 1)
      // x1 = x + a * p**(e-1)
      // x1**2 = x**2 + (a * p**(e-1))**2 + 2*x*a*p**(e-1) = n mod p**e
      // a*p**(e-1) = (n - x**2) * (2*x)**-1 mod p**e
      let x1 = x + ((n - x * x) % m * modInverseSmall(2 * x, m)) % m;
      if (x1 >= m) {
        x1 -= m;
      }
      if (x1 < 0) {
        x1 += m;
      }
      if (x1 > m - x1) {
        x1 = m - x1;
      }
      return [x1, m - x1];
    }
  }
  if (p % 2 === 0) {
    return [1];
  }
  let rrmnmodp = 1 - n; // r**2 % p - n
  for (let tworm1 = -1; tworm1 < p; tworm1 += 2) {
    rrmnmodp += tworm1;
    if (rrmnmodp >= p) {
      rrmnmodp -= p;
    }
    if (rrmnmodp === 0) {
      const r = Math.floor((tworm1 + 1) / 2);
      return [r, p - r];
    }
  }
  return [];
}

function log2(x) {
  return BigInt(x.toString(16).length * 4);
}

function sqrt(x) {
  if (x < BigInt((Number.MAX_SAFE_INTEGER + 1) / 2)) {
    return BigInt(Math.floor(Math.sqrt(Number(x) + 0.5)));
  }
  const q = (log2(x) / 4n);
  const initialGuess = ((sqrt(x >> (q * 2n)) + 1n) << q);
  let a = initialGuess;
  let b = a + 1n;
  while (a < b) {
    b = a;
    a = (b + x / b) / 2n;
  }
  return b;
}

function getSmoothFactorization(a, base) {
  let value = BigInt(a);
  if (value === 0n) {
    return [];
  }
  const result = [];
  if (value < 0n) {
    result.push(-1);
    value = -value;
  }
  let i = 0;

  while (i < base.length && Number(value) > Number.MAX_SAFE_INTEGER) {
    let primesProduct = 1;
    const productStart = i;
    while (i < base.length && primesProduct * base[i] <= Number.MAX_SAFE_INTEGER) {
      primesProduct *= base[i];
      i += 1;
    }
    const valueModPrimesProduct = Number(value % BigInt(primesProduct));
    for (let j = productStart; j < i; j += 1) {
      const p = base[j];
      if (valueModPrimesProduct - Math.floor(valueModPrimesProduct / p) * p === 0) {
        while (value % BigInt(p) === 0n) {
          value /= BigInt(p);
          result.push(p);
        }
      }
    }
  }

  let n = Number(value);
  while (i < base.length) {
    const p = base[i];
    while (n - Math.floor(n / p) * p === 0) {
      n /= p;
      result.push(p);
    }
    if (n !== 1 && n < p * p) {
      // n should be prime (?)
      const index = indexOf(base, n);
      if (index === -1) {
        return null;
      }
      result.push(n);
      return result;
    }
    i += 1;
  }
  return n === 1 ? result : null;
}

// (X**2 - Y) % N === 0, where Y is a smooth number
function CongruenceOfsquareOfXminusYmoduloN(X, Y, N, factorization) {
  this.X = X;
  this.Y = Y;
  this.N = N;
  this.factorization = factorization;
}
CongruenceOfsquareOfXminusYmoduloN.prototype.toString = function () {
  const X = this.X;
  const Y = this.Y;
  const N = this.N;
  return 'X**2 ≡ Y (mod N), Y = F'.replaceAll('X', X).replaceAll('Y', Y).replaceAll('N', N).replaceAll('F', this.factorization.join(' * '));
};

function isQuadraticResidueModuloPrime(a, p) {
  a = BigInt(a);
  p = Number(p);
  if (p === 2) {
    // "Modulo 2, every integer is a quadratic residue." - https://en.wikipedia.org/wiki/Quadratic_residue#Prime_modulus
    return true;
  }
  // https://en.wikipedia.org/wiki/Euler%27s_criterion
  const amodp = Number(BigInt(a) % BigInt(p));
  if (amodp === 0) {
    return true;
  }
  console.assert(p % 2 === 1);
  const value = modPowSmall(amodp, (p - 1) / 2, p);
  console.assert(value === 1 || value === p - 1);
  return value === 1;
}

function isQuadraticResidueModuloPrimeBig(a, p) {
  a = BigInt(a);
  p = BigInt(p);
  if (p === 2n) {
    // "Modulo 2, every integer is a quadratic residue." - https://en.wikipedia.org/wiki/Quadratic_residue#Prime_modulus
    return true;
  }
  // https://en.wikipedia.org/wiki/Euler%27s_criterion
  const amodp = a % p;
  if (amodp === 0n) {
    return true;
  }
  console.assert(p % 2n === 1n);
  const value = modPow(amodp, (p - 1n) / 2n, p);
  console.assert(value === 1n || value === p - 1n);
  return value === 1n;
}

function log(N) {
  const e = Math.max(N.toString(16).length * 4 - 4 * 12, 0);
  const lnn = Math.log(Number(N >> BigInt(e))) + Math.log(2) * e;
  return lnn;
}

function L(N) {  // exp(sqrt(log(n)*log(log(n))))
  const lnn = log(N);
  return Math.exp(Math.sqrt(lnn * Math.log(lnn)));
}

function product(array) {
  const n = array.length;
  const m = Math.floor(n / 2);
  return n === 0 ? 1n : (n === 1 ? BigInt(array[0]) : product(array.slice(0, m)) * product(array.slice(m)));
}

function modPowSmall(base, exponent, modulus) {
  base = Number(base);
  exponent = Number(exponent);
  modulus = Number(modulus);
  if (Math.max(Math.pow(modulus, 2), Math.pow(base, 2)) > Number.MAX_SAFE_INTEGER) {
    throw new RangeError();
  }
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

function modPow(base, exponent, modulus) {
  base = BigInt(base);
  exponent = BigInt(exponent);
  modulus = BigInt(modulus);
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

function primes(MAX) {
  const sieve = new Array(MAX + 1).fill(true);
  const result = [];
  result.push(2);
  for (let i = 3; i <= MAX; i += 2) {
    if (sieve[i]) {
      result.push(i);
      for (let j = i * i; j <= MAX; j += 2 * i) {
        sieve[j] = false;
      }
    }
  }
  return result;
}

const BitSetWordSize = 31; // see https://v8.dev/blog/pointer-compression

function packedArray(n) {
  // see https://v8.dev/blog/elements-kinds
  const array = [];
  for (let i = 0; i < n; i += 1) {
    array.push(0);
  }
  return array;
}
function BitSet(size) {
  const n = Math.ceil(size / BitSetWordSize);
  this.data = packedArray(n);
  this.size = size;
}
BitSet.prototype.nextSetBit = function (index) {
  if (index >= this.size) {
    return -1;
  }
  const data = this.data;
  let q = Math.floor(index / BitSetWordSize);
  let r = index % BitSetWordSize;
  let x = data[q] >> r;
  while (x === 0) {
    q += 1;
    if (q === data.length) {
      return -1;
    }
    x = data[q];
    r = 0;
  }
  // https://stackoverflow.com/questions/61442006/whats-the-most-efficient-way-of-getting-position-of-least-significant-bit-of-a
  r += 31 - Math.clz32(x & -x);
  return q * BitSetWordSize + r;
};
BitSet.prototype.toggle = function (index) {
  if (index >= this.size) {
    throw new RangeError();
  }
  const q = Math.floor(index / BitSetWordSize);
  const r = index % BitSetWordSize;
  this.data[q] ^= (r === BitSetWordSize - 1 ? ((-1) << r) : (1 << r));
};
BitSet.prototype.xor = function (other) {
  const a = this.data;
  const b = other.data;
  const n = a.length;
  if (n !== b.length) {
    throw new RangeError();
  }
  for (let i = 0; i < n; i += 1) {
    a[i] ^= b[i];
  }
};
BitSet.prototype.toString = function () {
  return this.data.map(x => (x >>> 0).toString(2).padStart(BitSetWordSize, '0').split('').reverse().join('')).join('').slice(0, this.size);
};

// pass factorizations with associated values to the next call
// returns linear combinations of vectors which result in zero vector by modulo 2
// (basis of the kernel of the matrix)
function solve(matrixSize) {
  // We build the augmented matrix in row-echelon form with permuted rows, which can grow up to matrixSize rows:
  // The augmented matrix is stored in the lower triangle!
  const M = new Array(matrixSize).fill(null); // We will fill the matrix so pivot elements will be placed on the diagonal
  const associatedValues = new Array(matrixSize).fill(undefined);
  let nextSolution = null;
  let state = 1;
  const iterator = {
    next: function solve(tmp) {
      while (true) {
        if (state === 1) {
          state = 0;
          return {value: nextSolution, done: false};
        }
        state = 1;
        const [rawRow, associatedValue] = tmp;
        let row = new BitSet(matrixSize);
        const reverseColumns = true; // makes it much faster when the data is more dense from the beginning (?)
        for (let j = 0; j < rawRow.length; j += 1) {
          const unitaryColumn = rawRow[j];
          const c = reverseColumns ? matrixSize - 1 - unitaryColumn : unitaryColumn;
          row.toggle(c);
        }
        // add row to the matrix maintaining it to be in row-echelon form:
        for (let pivotColumn = row.nextSetBit(0); pivotColumn !== -1 && row != null; pivotColumn = row == null ? -1 : row.nextSetBit(pivotColumn + 1)) {
          const pivotRow = M[pivotColumn];
          if (pivotRow != null) {
            // row-reduction:
            row.xor(pivotRow);
            console.assert(row.nextSetBit(pivotColumn) > pivotColumn || row.nextSetBit(pivotColumn) === -1);
            row.toggle(pivotColumn);
          } else {
            //row.toggle(matrixSize + pivotColumn);
            associatedValues[pivotColumn] = associatedValue;
            M[pivotColumn] = row;
            row = null;
          }
        }
        if (row != null) {
          // row has a solution
          // extract solution from the augmented part of the matrix:
          const solution = [];
          for (let i = row.nextSetBit(0); i !== -1; i = row.nextSetBit(i + 1)) {
            solution.push(associatedValues[i]);
          }
          solution.push(associatedValue);
          nextSolution = solution;
        } else {
          nextSolution = null;
        }
      }
      //console.log(M.filter(x => x != null).map(x => x.toString()).join('\n').replaceAll('0', ' '))
    }
  };
  iterator[globalThis.Symbol.iterator] = function () {
    return this;
  };
  return iterator;
}

//!copy-paste




function exp2(x) {
  return Math.pow(2, Math.floor(x)) * Math.exp(Math.LN2 * (x - Math.floor(x)));
}

const useMultiplePolynomials = true;

// A*x**2 + 2*B*x + C, A = q**2, qInv = q**-1 mod N
function QuadraticPolynomial(A, B, C, q, qInv, N) {
  this.A = A;
  this.B = B;
  this.C = C;
  this.q = q;
  this.qInv = qInv;
  this.N = N;
}

QuadraticPolynomial.prototype.calculateNewPolynomial = function (M, primes, N) {
  let q = this.q;
  if (q === 1n) {
    q = BigInt(sqrt(BigInt(sqrt(2n * N)) / BigInt(M)));//TODO: !?
    const B = BigInt(primes[primes.length - 1]);
    if (q <= B) {
      q = B + 2n;
    }
    q += 3n - q % 4n;
  } else {
    q += 4n;
  }
  while (!isPrime(q) || !isQuadraticResidueModuloPrimeBig(N, q)) {
    q += 4n;
  }
  const qInv = modInverse(q % N, N);
  if (qInv === 0n) {
    //TODO: what to do here - ?
    return new QuadraticPolynomial(0n, 0n, 0n, q, 0n, 0n).calculateNewPolynomial(M, primes, N);
  }
  const A = q * q;
  const B = getSquareRootsModuloPrimeBig(N, q, 2)[0];
  const AC = (B * B - N);
  if (AC % A !== 0n) {
    throw new Error();
  }
  const C = AC / A;
  return new QuadraticPolynomial(A, B, C, q, qInv, N);
};
QuadraticPolynomial.prototype.X = function (x) {
  return (this.A * BigInt(x) + this.B) * this.qInv;
};
QuadraticPolynomial.prototype.Y = function (x) {
  return this.A * (x * x <= Number.MAX_SAFE_INTEGER ? BigInt(x * x) : (a => a * a)(BigInt(x))) + this.B * BigInt(2 * x) + this.C;
};


function thresholdApproximationInterval(polynomial, x, threshold) {
  let w = 256;
  while (w >= 2 && Math.abs(Math.log2(Math.abs(Number(polynomial.Y(x + w)))) - threshold) > 0.5) {
    w /= 2;
  }
  return x + w;
}

// https://ru.wikipedia.org/wiki/Алгоритм_Диксона
// https://www.youtube.com/watch?v=TvbQVj2tvgc

function congruencesUsingQuadraticSieve(primes, N, sieveSize) {
  if (sieveSize == undefined) {
    sieveSize = Math.pow(2, 18);
    sieveSize = Math.min(sieveSize, Math.ceil(Math.pow(primes[primes.length - 1], 1.15)));
    sieveSize = Math.max(sieveSize, primes[primes.length - 1] + 1);
  }
  sieveSize += sieveSize % 2;

  N = BigInt(N);
  const SIEVE = new Array(sieveSize).fill(-0);

  const twoB = 2 * Math.log2(primes.length === 0 ? Math.sqrt(2) : Number(primes[primes.length - 1]));
  const largePrimes = Object.create(null); //TODO: new Map(); // faster (?)
  
  const xxx = [];

  // see https://www.youtube.com/watch?v=TvbQVj2tvgc
  const wheels = [];
  for (let i = 0; i < primes.length; i += 1) {
    const p = primes[i];
    for (let beta = 1, pInBeta = p; pInBeta <= sieveSize; beta += 1, pInBeta *= p) {
      const nmodpInBeta = Number(N % BigInt(pInBeta));
      if (nmodpInBeta % p === 0) {
        xxx.push(new CongruenceOfsquareOfXminusYmoduloN(BigInt(p), 0n, N, null));//?
      } else {
        const roots = getSquareRootsModuloPrime(nmodpInBeta, p, beta);
        for (let j = 0; j < roots.length; j += 1) {
          wheels.push({root: roots[j], log2p: Math.log2(p), step: pInBeta});
        }
      }
    }
  }

  const lpStrategy = function (p, X, Y) {
    // https://ru.wikipedia.org/wiki/Алгоритм_Диксона#Стратегия_LP
    const lp = largePrimes[p];
    if (lp == undefined) {
      largePrimes[p] = {X: X, Y: Y};
    } else {
      const s = BigInt(p);
      const sInverse = modInverse(s, N);
      if (sInverse === 0n) {
        return new CongruenceOfsquareOfXminusYmoduloN(s, 0n, N, null);//?
      } else {
        const X1 = (sInverse * lp.X * X) % N;
        if (Y % s === 0n && lp.Y % s === 0n) {
          const a = lp.Y / s;
          const b = Y / s;
          const Y1 = a * b;
          const fa = getSmoothFactorization(a, primes);
          const fb = getSmoothFactorization(b, primes);
          if (fa != null && fb != null) {
            const factorization = fa.concat(fb).sort((a, b) => a - b);
            return new CongruenceOfsquareOfXminusYmoduloN(X1, Y1, N, factorization);
          }
        }
      }
    }
    return null;
  };

  const wheelsRoots = wheels.map(wheel => wheel.root);

  let polynomial = null;
  if (!useMultiplePolynomials) {
    // - Number(baseOffset % BigInt(pInBeta))
    const baseOffset = BigInt(sqrt(N)) + 1n;
    polynomial = new QuadraticPolynomial(1n, baseOffset, baseOffset * baseOffset - N, 1n, 1n, N);
    for (let i = 0; i < wheels.length; i += 1) {
      const wheel = wheels[i];
      const pInBeta = wheel.step;
      wheelsRoots[i] = wheel.root - Number(baseOffset % BigInt(pInBeta));
    }
  } else {
    polynomial = new QuadraticPolynomial(1n, 0n, -N, 1n, 1n, N);
  }

  const updateWheels = function (polynomial) {
    const fmod = function (a, b) {
      return (a - Math.floor(a / b) * b) | 0;
    };
    //recalculate roots based on the formulat:
    //proot = ((-B + root) * modInv(A, pInBeta)) % pInBeta;
    //+some optimizations to minimize bigint usage and modInverseSmall calls
    let smallProduct = 1;
    let k = wheels.length - 1;
    for (let i = wheels.length - 1; i >= -1; i -= 1) {
      const pInBeta = i === -1 ? Number.MAX_SAFE_INTEGER : wheels[i].step;
      if (fmod(smallProduct, pInBeta) !== 0) {
        if (smallProduct * pInBeta >= Number.MAX_SAFE_INTEGER) {
          const a1 = Number(polynomial.A % BigInt(smallProduct));
          const b1 = Number(polynomial.B % BigInt(smallProduct));
          let invA = 0;
          let m = 1;
          while (k > i) {
            const wheel = wheels[k];
            const prevm = m;
            m = wheel.step;
            invA = prevm % m === 0 ? invA % m : modInverseSmall(fmod(a1, m), m);
            if (invA === 0) {
              throw new Error('unsupported A');
            }
            const b = fmod(b1, m);
            const proot = fmod((wheel.root - b) * invA, m);
            wheelsRoots[k] = proot;
            k -= 1;
          }
          smallProduct = 1;
        }
        smallProduct *= pInBeta;
      }
    }
    if (false) {
      for (let k = 0; k < wheels.length; k += 1) {
        const x = BigInt(wheelsRoots[k]);
        const X = (polynomial.A * x + polynomial.B);
        const Y = X * X - N;
        if (Y % BigInt(wheels[k].step) !== 0n) {
          throw new Error();
        }
      }
    }
  };
  
  const updateSieve = function (offset) {
    for (let j = 0; j < sieveSize; j += 1) {
      SIEVE[j] = -0;
    }
    for (let j = 0; j < wheels.length; j += 1) {
      const w = wheels[j];
      const root = wheelsRoots[j];
      const log2p = w.log2p;
      const step = w.step;
      const start = (root - offset) - Math.floor((root - offset) / step) * step;
      for (let kpplusr = start; kpplusr < sieveSize; kpplusr += step) {
        SIEVE[kpplusr] += log2p;
      }
    }
  };

  let i = -1;
  let k = 0;
  const iterator = {
    next: function congruencesUsingQuadraticSieve() {
      if (xxx.length > 0) {
        return {value: xxx.pop(), done: false};
      }
      while (2 * k * sieveSize <= Math.pow(primes[primes.length - 1], 2)) {
        if (i === sieveSize) {
          k += 1;
          i = -1;
        }
        const offset = useMultiplePolynomials ? -sieveSize / 2 : (k % 2 === 0 ? 1 : -1) * Math.floor((k + 1) / 2) * sieveSize;
        if (i === -1) {

          if (useMultiplePolynomials) {
            polynomial = polynomial.calculateNewPolynomial(sieveSize / 2, primes, N);
            updateWheels(polynomial);
          }

          updateSieve(offset);
        }

        let j = -1;
        let thresholdApproximation = 0.5;
        while (++i < sieveSize) {
          // it is slow to compute the threshold on every iteration, so trying to optimize:

          //TODO: the threshold calculation is much more simple in the Youtube videos (?)
          if (i >= j && sieveSize >= 2048) {
            const Y = polynomial.Y(i + offset);
            thresholdApproximation = Math.log2(Math.abs(Number(Y))) - twoB;
            j = thresholdApproximationInterval(polynomial, i + offset, thresholdApproximation + twoB) - offset;
          }

          const value = SIEVE[i];
          if (thresholdApproximation < value) {
            const X = polynomial.X(i + offset);
            const Y = polynomial.Y(i + offset);
            const threshold = Math.log2(Math.abs(Number(Y))) - twoB;
            if (threshold + twoB - 1 < value) {
              const factorization = getSmoothFactorization(Y, primes);
              if (factorization != null) {
                return {value: new CongruenceOfsquareOfXminusYmoduloN(X, Y, N, factorization), done: false};
              } else {
                console.count('?');
                /*let p = 1n;
                for (let n = 0; n < wheels.length; n += 1) {
                  const w = wheels[n];
                  if ((i + offset - wheelsRoots[n]) % w.step === 0) {
                    console.log(w);
                    p *= BigInt(w.step);
                  }
                }*/
              }
            } else {
              if (threshold < value) {
                const p = Math.round(exp2(threshold + twoB - value));
                const c = lpStrategy(p, X, Y);
                if (c != null) {
                  return {value: c, done: false};
                }
              }
            }
          }
        }
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
    const r = a % b;
    a = b;
    b = r;
  }
  return a;
}

function abs(x) {
  return x < 0n ? -x : x;
}

function indexOf(sortedArray, x) {
  let min = 0;
  let max = sortedArray.length - 1;
  while (min < max) {
    const mid = Math.ceil((min + max) / 2);
    if (sortedArray[mid] > x) {
      max = mid - 1;
    } else {
      min = mid;
    }
  }
  if (sortedArray[min] === x) {
    return min;
  }
  return -1;
}

function QuadraticSieveFactorization(N) { // N - is not a prime
  N = BigInt(N);
  for (let k = 1n;; k += 1n) {
    const kN = k * N;
    // https://trizenx.blogspot.com/2018/10/continued-fraction-factorization-method.html#:~:text=optimal%20value :
    const B = Math.min(Math.floor(Math.sqrt(L(kN) / 1.5)), (1 << 25) - 1);
    const primeBase = primes(B).filter(p => isQuadraticResidueModuloPrime(kN, p));
    const congruences = congruencesUsingQuadraticSieve(primeBase, kN); // congruences X_k^2 = Y_k mod N, where Y_k is smooth over the prime base
    const solutions = solve(1 + primeBase.length); // find products of Y_k = Y, so that Y is a perfect square
    solutions.next();
    let c = null;
    while ((c = congruences.next().value) != undefined) {
      const solution = c.Y === 0n ? [c] : solutions.next([c.factorization.map(p => (p === -1 ? 0 : 1 + indexOf(primeBase, p))), c]).value;
      if (solution != null) {
        const X = product(solution.map(c => c.X));
        const Y = product(solution.map(c => c.Y)); // = sqrt(X**2 % N)
        const x = X;
        const y = BigInt(sqrt(Y));
        console.assert(y * y === BigInt(Y));
        const g = gcd(abs(x + y), N);
        if (g !== 1n && g !== N) {
          return g;
        }
      }
    }
  }
}

QuadraticSieveFactorization.testables = {
  isPrime: isPrime,
  congruencesUsingQuadraticSieve: congruencesUsingQuadraticSieve,
  getSquareRootsModuloPrime: getSquareRootsModuloPrime,
  getSquareRootsModuloPrimeBig: getSquareRootsModuloPrimeBig,
  isQuadraticResidueModuloPrime: isQuadraticResidueModuloPrime,
  isQuadraticResidueModuloPrimeBig: isQuadraticResidueModuloPrimeBig,
  solve: solve,
  QuadraticPolynomial: QuadraticPolynomial,
  thresholdApproximationInterval: thresholdApproximationInterval
};

export default QuadraticSieveFactorization;
