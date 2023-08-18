import primeFactor from './primeFactor.js';
import Polynomial from './Polynomial.js';
import Expression from './Expression.js';
import './node_modules/seedrandom/seedrandom.js';
import combinations from './combinations.js';
import IntPolynomial from './IntPolynomial.js';//?

function toInt(c, p) {
  if (p instanceof Expression.Integer) {
    return Expression.Integer.fromBigInt(c);
  }
  return c;
}

function fromInt(c) {
  if (c instanceof Expression.Integer) {
    return c;
  }
  return Expression.Integer.fromBigInt(c);
}

function toIntPolynomial(f, p) {
  if (p instanceof Expression.Integer) {
    return f.mod(p);
  }
  var ep = Expression.Integer.fromBigInt(p);
  var coefficients = new Array(f.getDegree() + 1);
  for (var i = 0; i < coefficients.length; i += 1) {
    coefficients[i] = f.getCoefficient(i).modulo(ep).toBigInt();
  }
  return IntPolynomial.from(coefficients);
}

function fromIntPolynomial(f) {
  var coefficients = new Array(f.getDegree() + 1);
  for (var i = 0; i < coefficients.length; i += 1) {
    coefficients[i] = Expression.Integer.fromBigInt(f.getCoefficient(i));
  }
  return Polynomial.from(coefficients);
}

// Books:
// Henri Cohen "A Course in Computational Algebraic Number Theory"
// "Computer algebra and symbolic computation Mathematical Methods" Joel S. Cohen
// "The art of computer programming. Vol.2: Seminumerical algorithms" Donald E. Knuth
// https://en.wikipedia.org/wiki/Finite_field
// https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
// https://en.wikipedia.org/wiki/Factorization_of_polynomials#Factoring_univariate_polynomials_over_the_integers


const isPrime = primeFactor._isPrime;

function ExtendedEuclideanAlgorithm(A, B, p) {
  // U * A + V * B = gcd(A, B) (mod p)
  A = A.mod(p);
  B = B.mod(p);
  // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Pseudocode
  let [old_r, r] = [A, B];
  const ZERO = Polynomial.of(Expression.ZERO);
  const ONE = Polynomial.of(p.divide(p));//TODO: ?
  let [old_s, s] = [ONE, ZERO];
  let [old_t, t] = [ZERO, ONE];
  while (r.getDegree() >= 0) {
    const multiplier = r.getLeadingCoefficient().modInverse(p);
    const tmp = old_r.scale(multiplier).mod(p).divideAndRemainderModP(r.scale(multiplier).mod(p), p);
    const quotient = tmp.quotient;
    //[old_r, r] = [r, old_r.subtract(quotient.multiply(r)).mod(p)];
    [old_r, r] = [r, tmp.remainder.scale(r.getLeadingCoefficient()).mod(p)];
    [old_s, s] = [s, old_s.subtract(quotient.multiply(s)).mod(p)];
    [old_t, t] = [t, old_t.subtract(quotient.multiply(t)).mod(p)];
  }
  const k = old_r.getLeadingCoefficient().modInverse(p);
  const gcd = old_r.scale(k).mod(p);
  const U = old_s.scale(k).mod(p);
  const V = old_t.scale(k).mod(p);
  return {
    U: U,
    V: V,
    gcd: gcd
  };
}

// https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Equal-degree_factorization

function modPow(polynomial, n, m, q) {
  let p = polynomial;
  let accumulator = null;
  while (!n.equals(Expression.ZERO)) {
    if (!n.remainder(Expression.TWO).equals(Expression.ZERO)) {
      n = n.subtract(Expression.ONE);
      if (accumulator == null) {
        accumulator = p;
      } else {
        accumulator = accumulator.multiply(p).mod(q).divideAndRemainderModP(m, q).remainder;
      }
    } else {
      n = n.truncatingDivide(Expression.TWO);
      p = p.multiply(p).mod(q).divideAndRemainderModP(m, q).remainder;
    }
  }
  return accumulator;
}

function distinctDegreeFactorization(f, p) {
  f = f.mod(p);
  // copy-paste of pseudo code from Wikipedia - https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Distinct-degree_factorization
  let i = 1;
  const S = [];
  let fStar = f;
  fStar = toMonic(fStar, p);
  const q = p;
  const x = f.constructor.from([toInt(0, p), toInt(1, p)]);//TODO: ?
  let xInQInI = modPow(x, fromInt(q), fStar, q); // x**(q**i)
  while (fStar.getDegree() >= 2 * i) {
    //TODO: see the Wikipedia page for some optimizations - ?
    const h = xInQInI.subtract(x).divideAndRemainderModP(fStar, q).remainder;
    const g = gcdOfPolynomialsOverFiniteField(fStar, h, q);
    if (g.getDegree() !== 0) {
      S.push({factor: g, degree: i});
      fStar = fStar.divideAndRemainderModP(g, q).quotient;
    }
    i = i + 1;
    xInQInI = modPow(xInQInI, fromInt(q), fStar, q);
  }
  if (fStar.getDegree() > 0) {
    S.push({factor: fStar, degree: fStar.getDegree()});
  }
  if (S.length === 0) {
    S.push({factor: f, degree: 1});
  }
  return S;
}

function randomBigInt0(size, random = Math.random) {
  console.assert(Math.floor(size) === size);
  if (size <= 52) {
    return Math.floor(random() * 2**size);
  }
  const q = Math.ceil(size / (2 * 52)) * 52;
  return (BigInt(randomBigInt0(size - q, random)) << BigInt(q)) + BigInt(randomBigInt0(q, random));
}

// [0; max - 1]
function randomBigInt(max, random = Math.random) {
  if (Number(max) <= 2**52) {
    return Math.floor(random() * Number(max));
  }
  var size = Expression.Integer.fromBigInt(max).bitLength();
  return (BigInt(randomBigInt0(size, random)) * BigInt(max)) >> BigInt(size);
}

function CantorZassenhausAlgorithm(f, p, factorsDegree) {
  f = f.mod(p);
  f = toMonic(f, p);//TODO: is it needed here, test - ?
  // copy-paste of pseudo code from Wikipedia - https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Cantor–Zassenhaus_algorithm
  const n = f.getDegree();
  const d = factorsDegree;
  const r = n / d;
  const q = p;
  let Factors = [];
  Factors.push(f);
  const random = Math.seedrandom != null ? new Math.seedrandom('hello.') : Math.random;
  while (Factors.length < r) {
    const h = f.constructor.from(new Array(n).fill(q).map(q => toInt(randomBigInt(fromInt(q).toBigInt(), random), q)));
    const ONE = f.constructor.from([toInt(1, p)]);
    const qInDminusOneOverTwo = (fromInt(q)._pow(d).subtract(Expression.ONE)).truncatingDivide(Expression.TWO);
    const g = modPow(h, qInDminusOneOverTwo, f, q).subtract(ONE).divideAndRemainderModP(f, q).remainder;
    const updatedFactors = [];
    for (const u of Factors) {
      const gcd = gcdOfPolynomialsOverFiniteField(g, u, q);
      if (gcd.getDegree() !== 0 && gcd.getDegree() !== u.getDegree()) {
        updatedFactors.push(gcd);
        updatedFactors.push(u.divideAndRemainderModP(gcd, q).quotient);
      } else {
        updatedFactors.push(u);
      }
    }
    Factors = updatedFactors;
  }
  return Factors;
}

function factorizeOverTheFiniteField(f, p) {
  if (!isFactorizationOverZpSquareFree(f, p)) {
    throw new RangeError("implemented only for square-free polynomials");
  }
  let factorization = [];
  const distinctDegreeFactors = distinctDegreeFactorization(f, p);
  for (const ddf of distinctDegreeFactors) {
    const equalDegreeFactors = CantorZassenhausAlgorithm(ddf.factor, p, ddf.degree);
    factorization = factorization.concat(equalDegreeFactors);
  }
  return factorization;
}

function toMonic(f, p) {
  if (f.getDegree() < 0) {
    return f;
  }
  const scale = toInt(fromInt(f.getLeadingCoefficient()).modInverse(fromInt(p)).toBigInt(), p);
  return f.scale(scale).mod(p);
}



function gcdOfPolynomialsOverFiniteField(a, b, p) {
  a = a.mod(p);//?
  b = b.mod(p);//?
  b = toMonic(b, p);
  while (b.getDegree() >= 0) {
    let r = a.divideAndRemainderModP(b, p).remainder;
    r = toMonic(r, p);
    a = b;
    b = r;
  }
  return a;
}

function _gcdOfPolynomialsOverFiniteField0(a, b, p) {
  return fromIntPolynomial(gcdOfPolynomialsOverFiniteField(toIntPolynomial(a, p), toIntPolynomial(b, p), p));
}

function isFactorizationOverZpSquareFree(u, prime) {
  const f = u;
  return gcdOfPolynomialsOverFiniteField(f, f.derive().mod(prime), prime).getDegree() === 0;
}

// The art of computer programming. Vol.2: Seminumerical algorithms, page 452
// returns some factors, not necessary irreducible
function factorizeOverTheIntegers(u, useHenselLifting = true) {
  const polynomial = u;
  let factorsIterator1 = null;
  let factorsIterator2 = null;
  let results1 = null;
  let checkBothVariants = true;
  let prime = 0;
  let factors = null;
  let q = 0;
  let c = 0;
  let countOfFactors = 1;
  let combinationsIterator = null;
  const iterator = {
    next: function () {
      if (u.getDegree() === 0) {
        return {value: null, done: true};
      }
      while (u.getDegree() > 1 && u.getCoefficient(0).equals(Expression.ZERO)) {
        const factor = Polynomial.of(Expression.ZERO, Expression.ONE);
        u = u.divideAndRemainder(factor, "undefined").quotient;
        return {value: factor, done: false};
      }
      if (u.getCoefficient(0).abs().bitLength() - u.getLeadingCoefficient().abs().bitLength() < -42) {//?
        if (factorsIterator1 == null) {
          factorsIterator1 = factorizeOverTheIntegers(u._exponentiateRoots(-1), useHenselLifting);
        }
        const factor = factorsIterator1.next().value;
        if (factor != null) {
          const f = factor._exponentiateRoots(-1);
          return {value: f, done: false};
        }
        return {value: null, done: true};
      }
      if (u.isEven() && polynomial.getDegree() >= 4) { // https://math.stackexchange.com/a/2894104
        if (factorsIterator2 == null) {
          factorsIterator2 = factorizeOverTheIntegers(u._exponentiateRoots(2), useHenselLifting);
        }
        const factor = factorsIterator2.next().value;
        if (factor != null) {
          const f = factor._exponentiateRoots(1 / 2);
          u = u.divideAndRemainder(f, "undefined").quotient;
          return {value: f, done: false};
        }
        // see below
      }
      
      // from https://mathoverflow.net/a/449002 ( https://arxiv.org/pdf/0904.3057.pdf )
      function centralBinomialCoefficientBound(n) {
        return (n - Math.log2(Math.sqrt(Math.PI * Math.ceil(n / 2))));
      }
      function MignotteFactorBound(p) {
        return (centralBinomialCoefficientBound(p.getDegree()) + p._log2hypot());
      }
      function MignotteSingleFactorBound(p) {
        return MignotteFactorBound(p) * 0.5;
      }
      const getBound = function (p) {
        if (checkBothVariants) {
          return MignotteSingleFactorBound(p);
        }
        return MignotteFactorBound(p);
      };

      const nextGoodPrime = function (integer) {
        let p = integer;
        let pp = p;
        do {
          do {
            p = Expression.Integer.fromBigInt(p).add(Expression.TWO).toBigInt();
          } while (!isPrime(p));
          if (!useHenselLifting) {
            pp = Expression.Integer.fromBigInt(p);
          } else {
            pp = p;
          }
        } while (u.getLeadingCoefficient().remainder(Expression.Integer.fromBigInt(p)).equals(Expression.ZERO) || !isFactorizationOverZpSquareFree(toIntPolynomial(u, pp), pp));
        return pp;
      };

      if (factors == null) {
        const B = 1 + getBound(u.scale(u.getLeadingCoefficient().abs()));
        //const useHenselLifting = true;//TODO: ?
        if (!useHenselLifting) {
          prime = nextGoodPrime(Expression.TWO._pow(Math.ceil(B)).add(Expression.ONE).toBigInt());
        } else {
          prime = nextGoodPrime(1);
        }
        const tryMultiplePrimes = !useHenselLifting ? 0 : 2;
        if (tryMultiplePrimes !== 0) {
          let best = prime;
          let bestFactorsNumber = 1 / 0;
          for (let tries = 0; tries < tryMultiplePrimes && bestFactorsNumber > 1; tries += 1) {
            let factorsNumber = 0;
            const ddfs = distinctDegreeFactorization(toIntPolynomial(u, prime), prime);
            for (const entry of ddfs) {
              console.assert((entry.factor.getDegree() % entry.degree) === 0);
              factorsNumber += (entry.factor.getDegree() / entry.degree);
            }
            if (bestFactorsNumber > factorsNumber) {
              best = prime;
              bestFactorsNumber = factorsNumber;
            }
            prime = nextGoodPrime(prime);
          }
          prime = best;
        }
        factors = factorizeOverTheFiniteField(toIntPolynomial(u, prime), prime).map(factor => useHenselLifting ? fromIntPolynomial(factor) : factor);
        if (useHenselLifting && factors.length > 15) {
          checkBothVariants = false;
        }

        if (useHenselLifting) {
          //This gives smaller value:
          const e = Math.ceil(B / Math.log2(prime));
          //if (useQuadraticHenselLift) {
          //  e = Math.pow(2, Math.ceil(Math.log2(e)));
          //}
          const p = Expression.Integer.fromNumber(prime);
          factors[factors.length - 1] = factors[factors.length - 1].scale(u.getLeadingCoefficient().modulo(p)).mod(p);
          factors = HenselLifting(u, factors, p, e);
          q = p._pow(e);
          factors = factors.map(factor => factor.scale(factor.getLeadingCoefficient().modulo(q).modInverse(q)).mod(q));//TODO: ?
        } else {
          q = prime;
        }
      }

      //!new 2022-07-27
      if (u.isEven() && polynomial.getDegree() >= 4) { // https://math.stackexchange.com/a/2894104
        // see above
        if (results1 == null) {
          for (let i = 0; i < factors.length; i += 1) {
            if (factors[i] != null) {
              const f1 = factors[i].mod(q);
              const f2 = factors[i]._scaleRoots(Expression.ONE.negate()).mod(q);
              let found = false;
              for (let j = i + 1; j < factors.length && !found; j += 1) {
                if (factors[j] != null && factors[j].mod(q).equals(f2)) {
                  factors[j] = null;
                  found = true;
                }
              }
              if (!found) {
                return {value: null, done: true};
              }
            }
          }
          factors = factors.filter(f => f != null);
          for (let i = 0; i < Math.pow(2, factors.length - 1) && results1 == null; i += 1) {
            const combination = factors.map((f, index) => Math.floor(i / 2**index) % 2 === 0 ? f : f._scaleRoots(Expression.ONE.negate()));
            let ok = true;
            if (true) { //TODO: !?
              const lc = u.getLeadingCoefficient();
              // 2 = 2 + 0 = 1 + 1 = 0 + 2
              // 1 = 1 + 0 = 0 + 1
              // 0 = 0 + 0
              let t2 = Expression.ZERO;
              let t1 = Expression.ZERO;
              let t0 = lc;
              for (const f of combination) {
                const f2 = f.getCoefficient(2);
                const f1 = f.getCoefficient(1);
                const f0 = f.getCoefficient(0);
                t2 = t2.multiply(f0).add(t1.multiply(f1)).add(t0.multiply(f2)).modulo(q);
                t1 = t1.multiply(f0).add(t0.multiply(f1)).modulo(q);
                t0 = t0.multiply(f0).modulo(q);
              }
              let tail = Polynomial.of(t0, t1, t2).mod2(q);
              tail = tail.multiply(tail._scaleRoots(Expression.ONE.negate()));
              ok = tail.getCoefficient(2).abs().equals(u.getCoefficient(2).multiply(u.getLeadingCoefficient()).abs());
            }
            if (ok) {
              const candidate = productModQ(combination, q).scale(u.getLeadingCoefficient()).mod2(q).primitivePart();
              const candidate2 = candidate._scaleRoots(Expression.ONE.negate());//TODO: return candidate2
              if (candidate.multiply(candidate2).equals(u)) {
                results1 = [];
                results1.push(candidate);
                results1.push(candidate2);
              }
            }
          }
        }
        if (results1 != null) {
          while (results1.length > 0) {
            return {value: results1.pop(), done: false};
          }
        }
        return {value: null, done: true};
      }

      //!!! (number of factors depends on the choise of prime numbers)
      //TODO: how to reduce number of iterations (?) (see Donald Knuth's book)
      for (; countOfFactors <= (checkBothVariants ? factors.length - 1 : Math.floor(factors.length / 2)); countOfFactors += 1) {
        combinationsIterator = combinations(factors, countOfFactors);
        let combination = null;
        while ((combination = combinationsIterator.next().value) != null) {
          c += 1;
          const lc = u.getLeadingCoefficient();
          // an optimization from the Donald Knuth's book, page 452
          let productTrailingCoefficient = lc;
          for (const f of combination) {
            productTrailingCoefficient = productTrailingCoefficient.multiply(f.getCoefficient(0)).modulo(q);
          }
          productTrailingCoefficient = Polynomial.of(productTrailingCoefficient).mod2(q).getCoefficient(0);
          if (u.getCoefficient(0).multiply(lc).remainder(productTrailingCoefficient).equals(Expression.ZERO)) {
            let v = productModQ(combination, q);
            console.assert(v.getLeadingCoefficient().equals(Expression.ONE));
            v = v.scale(lc);
            v = v.mod2(q);
            //TODO: test, we need to try w(x) = productModQ(factors.filter(factor => combination.indexOf(factor) === -1), q) as well (see Donald Knuth's book) - ?
            console.assert(v.getDegree() < u.getDegree());
            //if (v.getDegree() <= u.getDegree() / 2 || v.getDegree() < u.getDegree()) {
              //v = v.primitivePart();
              const tmp = u.scale(lc).divideAndRemainder(v, "undefined");
              if (tmp != undefined && tmp.remainder.getDegree() < 0) {
                v = v.primitivePart();
                factors = factors.filter(factor => combination.indexOf(factor) === -1);
                combinationsIterator = combinations(factors, countOfFactors);//!?
                u = tmp.quotient;
                if (c > 16) {
                  console.debug(c);
                }
                return {value: v, done: false};
              }
            //}
          }
        }
      }
      if (c > 16) {
        console.debug(c);
      }
      if (polynomial.getDegree() > u.getDegree()) {
        const f = u.primitivePart();//?
        u = Polynomial.of(u.getContent());
        return {value: f, done: false};
      }
      return {value: null, done: true};
    }
  };
  iterator[globalThis.Symbol.iterator] = function () {
    return this;
  };
  return iterator;
}

// if q == p, then C = A * B (mod p) -> A1 * B1 (mod p**2), A1 = A (mod p) and B1 = B (mod p)
function HenselLift(C, A, B, U, V, q, p) { // q -> q * p
  // C = A * B mod p
  /*
  // https://www.csd.uwo.ca/~mmorenom/CS874/Lectures/Newton2Hensel.html/node17.html#eq:FactorizationLiftingProblem
  //C = C.map(c => new IntegerModuloPrimeNumber(c, Expression.Integer.fromBigInt(p**2)))
  //TODO: ?
  const e = C.subtract(A.multiply(B));
  const A1 = A.add(U.multiply(e));
  const B1 = B.add(V.multiply(e));
  return {A1: A1, B1: B1};
  */
  // http://tomlr.free.fr/Math%E9matiques/Math%20Complete/Number%20theory/A%20course%20in%20computational%20algebraic%20number%20theory%20-%20Cohen%20H..pdf
  // Algorithm 3.5.5 (Hensel Lift).
  // A,B,C are polynomials over Integers:
  //console.assert(A.hasIntegerCoefficients());
  //console.assert(B.hasIntegerCoefficients());
  //console.assert(C.hasIntegerCoefficients());
  console.assert(q.isDivisibleBy(p));
  //C = C.mod(q.multiply(p));//TODO: ???
  const f = C.subtract(A.multiply(B)).scale(q.inverse()).mod(p);
  const tmp = V.multiply(f).mod(p).divideAndRemainderModP(A.mod(p), p);
  const t = tmp.quotient;
  //const A0 = V.multiply(f).subtract(A.multiply(t)).mod(p);
  const A0 = tmp.remainder;
  const B0 = U.multiply(f).add(B.mod(p).multiply(t)).mod(p);
  const A1 = A.add(A0.scale(q));
  const B1 = B.add(B0.scale(q));
  //console.assert(A1.mod(q.multiply(p)).equals(A1));
  //console.assert(B1.mod(q.multiply(p)).equals(B1));
  return [A1, B1];
}
function QuadraticHenselLift(A1, B1, U, V, p) {
  // http://tomlr.free.fr/Math%E9matiques/Math%20Complete/Number%20theory/A%20course%20in%20computational%20algebraic%20number%20theory%20-%20Cohen%20H..pdf
  // Algorithm 3.5.6 
  const one = p.divide(p);
  const g = Polynomial.of(one).subtract(U.multiply(A1)).subtract(V.multiply(B1)).scale(p.inverse()).mod(p);
  const tmp = V.multiply(g).mod(p).divideAndRemainderModP(A1.mod(p), p);
  const t = tmp.quotient;
  const U0 = U.multiply(g).add(B1.mod(p).multiply(t)).mod(p);
  //const V0 = V.multiply(g).subtract(A1.multiply(t)).mod(p);
  const V0 = tmp.remainder;
  const U1 = U.add(U0.scale(p));
  const V1 = V.add(V0.scale(p));
  //console.assert(U1.mod(p.multiply(p)).equals(U1));
  //console.assert(V1.mod(p.multiply(p)).equals(V1));
  return [U1, V1];
}
function HenselLiftingOfTwoFactors(C, A, B, p, k) {
  const useQuadraticHenselLift = true;
  const tmp1 = ExtendedEuclideanAlgorithm(A, B, p);
  console.assert(tmp1.gcd.getDegree() === 0);
  let U = tmp1.U;
  let V = tmp1.V;
  const ok = !(p instanceof Expression.Polynomial); // somehow the quadratic hensel lifting is slower in other case
  if (useQuadraticHenselLift && ok) { // TODO: any degree
    const originalP = p;
    let e = 1;
    while (e < k / 4) {
      [A, B] = HenselLift(C, A, B, U, V, p, p);
      [U, V] = QuadraticHenselLift(A, B, U, V, p);
      p = p.multiply(p);
      e *= 2;
      if (true) {
        let c = 1;
        while (e * c < k) {
          c *= 2;
        }
        if ((e - 1) * c >= k) {
          e -= 1;
          p = p.divide(originalP);
        }
      }
    }
    let e0 = e;
    let q = p;
    while (e < k) {
      [A, B] = HenselLift(C, A, B, U, V, q, p);
      q = q.multiply(p);
      e += e0;
    }
    if (e !== k) {
      const pInK = originalP._pow(k);
      A = A.mod(pInK);
      B = B.mod(pInK);
    }
    return [A, B];
  }
  //TODO: ?
  let q = p;
  for (let i = 1; i < k; i += 1) {
    [A, B] = HenselLift(C, A, B, U, V, q, p);
    //console.assert(U.multiply(A).add(V.multiply(B)).mod(p).toString() === '1');
    q = q.multiply(p);
  }
  return [A, B];
}
function productModQ(factors, q) {
  console.assert(factors.length > 0);
  return factors.length > 1 ? productModQ(factors.slice(0, Math.ceil(factors.length / 2)), q).multiply(productModQ(factors.slice(Math.ceil(factors.length / 2)), q)).mod(q) : factors[0];
}
function HenselLifting(C, factors, p, e) {
  // https://scholar.rose-hulman.edu/cgi/viewcontent.cgi?article=1163&context=math_mstr
  // "2.3 Factoring mod p e: Hensel Lifting"
  if (factors.length === 1) {
    return [C];
  }
  // divide and conquer
  const s = Math.ceil(factors.length / 2);
  const A = factors.slice(0, s);
  const B = factors.slice(s);
  const [A1, B1] = HenselLiftingOfTwoFactors(C, productModQ(A, p), productModQ(B, p), p, e);
  return HenselLifting(A1, A, p, e).concat(HenselLifting(B1, B, p, e));
}

factorizeOverTheIntegers._gcdOfPolynomialsOverFiniteField0 = _gcdOfPolynomialsOverFiniteField0; //TODO: ?

function factorizeMultivariateIntegerPolynomial(p) {
  // see "Art of Computer Programming, Volume 2: Seminumerical Algorithms"
  function factorizeInternal(p) {//TODO: REMOVE
    var factors = [];
    p = p.primitivePart();
    var f = p.getDegree() > 1 ? (!p.hasIntegerCoefficients() ? p.factorize() : factorizeOverTheIntegers(p).next().value) : null;//TODO: ?
    if (f != null) {
      factors = factors.concat(factorizeInternal(f));
      factors = factors.concat(factorizeInternal(p.divideAndRemainder(f, "throw").quotient));
    } else {
      factors.push(p);
    }
    return factors;
  }
  var degreeByY = 0;
  for (var i = 0; i <= p.getDegree(); i += 1) {
    if (!(p.getCoefficient(i).equals(Expression.ZERO))) {
      degreeByY = Math.max(degreeByY, p.getCoefficient(i).polynomial.getDegree());
    }
  }
  //if (degreeByY > p.getDegree()) {
  //  return toPolynomialByAnotherVar(factorizeMultivariateIntegerPolynomial(toPolynomialByAnotherVar(p)));
  //}
  //degreeByY = Math.pow(2, Math.ceil(Math.log2(degreeByY + 1))) - 1;//!?
  for (var y = 0;; y += 1) {
    var p_r = p.map(c => c.polynomial.calcAt(Expression.Integer.fromNumber(y)));
    if (p_r.getDegree() === p.getDegree() && p_r.isSquareFreePolynomial()) {
      var factors = Array.from(factorizeInternal(p_r)).map(f => f.map(c => new Expression.Polynomial(Polynomial.of(c))));
      if (factors.length < 2) {
        return null; // primitive (?)
      }
      var s = new Expression.Polynomial(Polynomial.of(p_r.getContent()));
      factors[factors.length - 1] = factors[factors.length - 1].scale(s);
      const r = new Expression.Polynomial(Polynomial.of(Expression.Integer.fromNumber(0 - y), Expression.ONE));
      var q = r._pow(degreeByY + 1);
      factors = HenselLifting(p, factors, r, degreeByY + 1);
      console.assert(p.subtract(productModQ(factors, q)).mod(q).toString() === '0');
      for (var number = 1; number <= factors.length - 1; number += 1) {
        for (var c of combinations(factors, number)) {
          var candidate = productModQ(c, q);
          candidate = candidate.scale(p.getLeadingCoefficient()).mod(q);
          if (candidate._hasIntegerLikeCoefficients()) {
            candidate = candidate.primitivePart();
            if (p.scale(p.getLeadingCoefficient()._pow(p.getDegree() - candidate.getDegree() + 1)).isDivisibleBy(candidate)) {
              return candidate;
            }
          }
        }
      }
    }
  }
}

factorizeOverTheIntegers._factorizeMultivariateIntegerPolynomial = factorizeMultivariateIntegerPolynomial;

export default factorizeOverTheIntegers;

factorizeOverTheIntegers.testables = {
  gcdOfPolynomialsOverFiniteField: gcdOfPolynomialsOverFiniteField,
  distinctDegreeFactorization: distinctDegreeFactorization,
  CantorZassenhausAlgorithm: CantorZassenhausAlgorithm,
  isFactorizationOverZpSquareFree: isFactorizationOverZpSquareFree,
  factorizeOverTheFiniteField: factorizeOverTheFiniteField,
  modPow: modPow,

  randomBigInt: randomBigInt,

  ExtendedEuclideanAlgorithm: ExtendedEuclideanAlgorithm,
  HenselLift: HenselLift,
  QuadraticHenselLift: QuadraticHenselLift,
  HenselLiftingOfTwoFactors: HenselLiftingOfTwoFactors,
  HenselLifting: HenselLifting,
  productModQ: productModQ
};
