import primeFactor from './primeFactor.js';
import Polynomial from './Polynomial.js';
import Expression from './Expression.js';
import './node_modules/seedrandom/seedrandom.js';

// Books:
// Henri Cohen "A Course in Computational Algebraic Number Theory"
// "Computer algebra and symbolic computation Mathematical Methods" Joel S. Cohen
// "The art of computer programming. Vol.2: Seminumerical algorithms" Donald E. Knuth
// https://en.wikipedia.org/wiki/Finite_field
// https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields
// https://en.wikipedia.org/wiki/Factorization_of_polynomials#Factoring_univariate_polynomials_over_the_integers


const modInverse = primeFactor._modInverse;
const nextPrime = primeFactor._nextPrime;

function ExtendedEuclideanAlgorithm(A, B, p) {
  // U * A + V * B = gcd(A, B) (mod p)
  p = Expression.Integer.fromBigInt(p);
  A = A.mod(p);
  B = B.mod(p);
  // https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Pseudocode
  let [old_r, r] = [A, B];
  const ONE = Polynomial.of(Expression.ONE);//TODO: ?
  let [old_s, s] = [ONE, Polynomial.ZERO];
  let [old_t, t] = [Polynomial.ZERO, ONE];
  while (!r.equals(Polynomial.ZERO)) {
    const multiplier = r.getLeadingCoefficient().modInverse(p);
    const quotient = old_r.scale(multiplier).mod(p).divideAndRemainderModP(r.scale(multiplier).mod(p), "throw", p).quotient;
    [old_r, r] = [r, old_r.subtract(quotient.multiply(r)).mod(p)];
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

Expression.Integer.prototype.modulo = function modulo(p) {
  const r = this.remainder(p);
  return r.compareTo(Expression.ZERO) < 0 ? r.add(p) : r;
};

Polynomial.prototype.modPow = function (n, m, q) {
  let p = this;
  let accumulator = null;
  while (n.compareTo(Expression.ZERO) > 0) {
    if (n.remainder(Expression.TWO).compareTo(Expression.ZERO) !== 0) {
      n = n.subtract(Expression.ONE);
      accumulator = (accumulator == null ? p : accumulator.multiply(p).mod(q)).divideAndRemainderModP(m, "throw", q).remainder;
    } else {
      n = n.truncatingDivide(Expression.TWO);
      p = p.multiply(p).mod(q).divideAndRemainderModP(m, "throw", q).remainder;
    }
  }
  return accumulator;
};

function distinctDegreeFactorization(f, p) {
  p = Expression.Integer.fromBigInt(p);
  f = f.mod(p);
  // copy-paste of pseudo code from Wikipedia - https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Distinct-degree_factorization
  let i = 1;
  const S = [];
  let fStar = f;
  fStar = toMonic(fStar, p);
  const q = p;
  const x = Polynomial.of(Expression.ONE).shift(1);//TODO: ?
  let xInQInI = x.modPow(q, fStar, q); // x**(q**i)
  while (fStar.getDegree() >= 2 * i) {
    //TODO: see the Wikipedia page for some optimizations - ?
    const h = xInQInI.subtract(x).divideAndRemainderModP(fStar, "throw", q).remainder;
    const g = gcdOfPolynomialsOverFiniteField(fStar, h, q);
    if (g.getDegree() !== 0) {
      S.push({factor: g, degree: i});
      fStar = fStar.divideAndRemainderModP(g, "throw", q).quotient;
    }
    i = i + 1;
    xInQInI = xInQInI.modPow(q, fStar, q);
  }
  if (fStar.getDegree() > 0) {
    S.push({factor: fStar, degree: fStar.getDegree()});
  }
  if (S.length === 0) {
    S.push({factor: f, degree: 1});
  }
  return S;
}

function randomBigInt(max, random) {
  if (Number(max) > 2**53) {
    let result = BigInt(0);
    let i = BigInt(1);
    while (i <= BigInt(max)) {
      result *= BigInt(2**53);
      result += BigInt(Math.floor(random() * 2**53));
      i *= BigInt(2**53);
    }
    return result * BigInt(max) / i;
  }
  return Math.floor(random() * Number(max));
}
function randomPolynomial(maxCoefficient, maxDegree, random) {
  const coefficients = new Array(maxDegree);
  for (let i = 0; i < maxDegree; i += 1) {
    coefficients[i] = Expression.Integer.fromBigInt(randomBigInt(maxCoefficient, random));
  }
  return Polynomial.from(coefficients);
}
Polynomial.random = randomPolynomial;//TODO: remove - ?

function CantorZassenhausAlgorithm(f, p, factorsDegree) {
  p = Expression.Integer.fromBigInt(p);
  f = f.mod(p);
  f = toMonic(f, p);//TODO: is it needed here, test - ?
  // copy-paste of pseudo code from Wikipedia - https://en.wikipedia.org/wiki/Factorization_of_polynomials_over_finite_fields#Cantor–Zassenhaus_algorithm
  const n = f.getDegree();
  const d = factorsDegree;
  const r = n / d;
  const q = p;
  let Factors = [];
  Factors.push(f);
  const random = new Math.seedrandom('hello.');
  while (Factors.length < r) {
    const h = randomPolynomial(q.toBigInt(), n, random);
    const ONE = Polynomial.of(Expression.ONE);
    const g = h.modPow((q._pow(d).subtract(Expression.ONE)).truncatingDivide(Expression.TWO), f, q).subtract(ONE).divideAndRemainderModP(f, "throw", q).remainder;
    const updatedFactors = [];
    for (const u of Factors) {
      const gcd = gcdOfPolynomialsOverFiniteField(g, u, q);
      if (gcd.getDegree() !== 0 && gcd.getDegree() !== u.getDegree()) {
        updatedFactors.push(gcd);
        updatedFactors.push(u.divideAndRemainderModP(gcd, "throw", q).quotient);
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


Polynomial.prototype.mod = function (m) {
  return this.map(c => c.modulo(m));
};
Polynomial.prototype.mod2 = function (m) {
  return this.mod(m).map(c => c.subtract(m).add(c).compareTo(Expression.ZERO) < 0 ? c : c.subtract(m));
};

Expression.Integer.prototype.modInverse = function (p) {
  return Expression.Integer.fromBigInt(modInverse(this.toBigInt(), p.toBigInt()));
};

function toMonic(f, p) {
  if (f.getLeadingCoefficient().compareTo(Expression.ZERO) === 0) {
    return f;
  }
  const scale = f.getLeadingCoefficient().modInverse(p);
  return f.map(c => c.multiply(scale).modulo(p));
}

//@private
Polynomial.prototype.divideAndRemainderModP = function (divisor, w, p) {
  if (w !== "throw") {
    throw new RangeError();
  }
  //const tmp = this.divideAndRemainder(divisor, "throw");
  //return {quotient: tmp.quotient.mod(p), remainder: tmp.remainder.mod(p)};
  const dividend = this;
  const divisorLeadingCoefficient = divisor.getLeadingCoefficient();
  if (divisorLeadingCoefficient.compareTo(Expression.ONE) !== 0) {
    throw new RangeError();
  }
  const minusDivisor = divisor.negate().mod(p);
  const divisorDegree = divisor.getDegree();
  //let remainder = dividend;
  let remainder = new Array(dividend.getDegree() + 1).fill(Expression.ZERO);
  for (let i = 0; i < dividend.a.size; i += 1) {
    remainder[dividend.a.degree(i)] = dividend.a.coefficient(i);
  }
  //let remainderDegree = remainder.getDegree();
  let remainderDegree = remainder.length - 1;
  let quotient = new Array(Math.max(remainderDegree - divisorDegree, 0)).fill(Expression.ZERO);
  while (remainderDegree >= divisorDegree) {
    const n = remainderDegree - divisorDegree;
    const q = remainder[remainderDegree];
    //const q = remainder.getLeadingCoefficient();
    quotient[n] = q;
    //TODO: optimize
    //remainder = remainder.add(minusDivisor.shift(n).scale(q)).mod(p);
    for (let j = 0; j < minusDivisor.a.size; j += 1) {
      const degree = minusDivisor.a.degree(j);
      const coefficient = minusDivisor.a.coefficient(j);
      remainder[degree + n] = remainder[degree + n].add(q.multiply(coefficient)).modulo(p);
    }
    while (remainderDegree >= 0 && remainder[remainderDegree].compareTo(Expression.ZERO) === 0) {
      remainderDegree -= 1;
    }
  }
  return {quotient: Polynomial.from(quotient), remainder: Polynomial.from(remainder)};
};

function gcdOfPolynomialsOverFiniteField(a, b, p) {
  a = a.mod(p);
  b = b.mod(p);
  b = toMonic(b, p);
  //TODO: make monic - ?
  while (!b.equals(Polynomial.ZERO)) {
    let r = a.divideAndRemainderModP(b, "throw", p).remainder;
    r = toMonic(r, p);
    a = b;
    b = r;
  }
  return a;
}

function isFactorizationOverZpSquareFree(u, prime) {
  prime = Expression.Integer.fromBigInt(prime);
  const f = u;
  return gcdOfPolynomialsOverFiniteField(f, f.derive(), prime).getDegree() === 0;
}

// The art of computer programming. Vol.2: Seminumerical algorithms, page 452
function factorizeOverTheIntegers(u) {  
  const polynomial = u;
  const B = u._log2OfBoundForCoefficientsOfFactor();
  let prime = undefined;
  const useHenselLifting = true;//TODO: ?
  if (!useHenselLifting) {
    prime = nextPrime(BigInt(2) * BigInt(Math.pow(2, Math.ceil(B))));
  } else {
    prime = 3;
    while (u.getLeadingCoefficient().remainder(Expression.Integer.fromBigInt(prime)).compareTo(Expression.ZERO) === 0 || !isFactorizationOverZpSquareFree(u, prime)) {
      prime = Number(nextPrime(prime));
    }
  }
  const tryMultiplePrimes = 2;
  if (tryMultiplePrimes !== 0) {
    let best = prime;
    let bestFactorsNumber = 1 / 0;
    for (let tries = 0; tries < tryMultiplePrimes; tries += 1) {
      let factorsNumber = 0;
      for (const entry of distinctDegreeFactorization(u, prime)) {
        factorsNumber += (entry.factor.getDegree() / entry.degree);
      }
      if (bestFactorsNumber > factorsNumber) {
        best = prime;
        bestFactorsNumber = factorsNumber;
      }
      prime = nextPrime(prime);
      while (u.getLeadingCoefficient().remainder(Expression.Integer.fromBigInt(prime)).compareTo(Expression.ZERO) === 0 || !isFactorizationOverZpSquareFree(u, prime)) {
        prime = Number(nextPrime(prime));
      }
    }
    prime = best;
  }
  let factors = factorizeOverTheFiniteField(u, prime);
  let q = prime;
  if (useHenselLifting) {
    //TODO: bitLength(u_n)
    let e = Math.ceil((1 + u.getLeadingCoefficient().abs().toBigInt().toString(16).length * 4 + B) / Math.log2(Number(prime)));
    if (useQuadraticHenselLift) {
      e = Math.pow(2, Math.ceil(Math.log2(e)));
    }
    factors = HenselLifting(u, factors, prime, e);
    q = Expression.Integer.fromBigInt(prime)._pow(e).toBigInt();
  }
  //!!! (number of factors depends on the choise of prime numbers)
  //TODO: how to reduce number of iterations (?) (see Donald Knuth's book)
  let c = 0;
for (let countOfFactors = 1; countOfFactors <= Math.floor(factors.length / 2); countOfFactors += 1) {
  let combinations = Polynomial._combinations(factors, countOfFactors);
  let combination = null;
  while ((combination = combinations.next().value) != null) {
    c += 1;
    let v = product(combination);
    console.assert(v.getLeadingCoefficient().compareTo(Expression.ONE) === 0);
    v = v.scale(u.getLeadingCoefficient());
    v = v.mod2(Expression.Integer.fromBigInt(q));
    if (v.getDegree() <= u.getDegree() / 2 || v.getDegree() < u.getDegree()) {
      //v = v.primitivePart();
      const tmp = u.scale(u.getLeadingCoefficient()).divideAndRemainder(v, "undefined");
      if (tmp != undefined && tmp.remainder.equals(Polynomial.ZERO)) {
        v = v.primitivePart();
        factors = factors.filter(factor => combination.indexOf(factor) === -1);
        combinations = Polynomial._combinations(factors, countOfFactors);//!?
        u = tmp.quotient;
        return v;
      }
    }
  }
}
if (c > 16) {
  console.debug(c);
}
  if (!polynomial.subtract(u).equals(Polynomial.ZERO)) {
    u = u.primitivePart();//?
    return u;
  }
  return null;
}


function HenselLift(C, A, B, U, V, q, r) { // q -> q * r
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
  q = Expression.Integer.fromBigInt(q);
  r = Expression.Integer.fromBigInt(r);
  console.assert(A.hasIntegerCoefficients());
  console.assert(B.hasIntegerCoefficients());
  console.assert(C.hasIntegerCoefficients());
  const f = C.subtract(A.multiply(B)).scale(q.inverse()).mod(r);
  const t = V.multiply(f).mod(r).divideAndRemainderModP(A, "throw", r).quotient;
  const A0 = V.multiply(f).subtract(A.multiply(t)).mod(r);
  const B0 = U.multiply(f).add(B.multiply(t)).mod(r);
  const A1 = A.add(A0.scale(q));
  const B1 = B.add(B0.scale(q));
  return {A1: A1.mod(q.multiply(r)), B1: B1.mod(q.multiply(r))};
}
function QuadraticHenselLift(C, A, B, U, V, q, p) {
  // http://tomlr.free.fr/Math%E9matiques/Math%20Complete/Number%20theory/A%20course%20in%20computational%20algebraic%20number%20theory%20-%20Cohen%20H..pdf
  // Algorithm 3.5.6 
  const tmp = HenselLift(C, A, B, U, V, q, p);
  const A1 = tmp.A1;
  const B1 = tmp.B1;
  p = Expression.Integer.fromBigInt(p);
  q = Expression.Integer.fromBigInt(q);
  const g = Polynomial.of(Expression.ONE).subtract(U.multiply(A1)).subtract(V.multiply(B1)).scale(p.inverse()).mod(p);
  const t = V.multiply(g).mod(q).divideAndRemainderModP(A1, "throw", q).quotient;
  const U0 = U.multiply(g).add(B1.multiply(t)).mod(q);
  const V0 = V.multiply(g).subtract(A1.multiply(t)).mod(q);
  const U1 = U.add(U0.scale(p));
  const V1 = V.add(V0.scale(p));
  //return HenselLift(C, A1, B1, U1, V1, BigInt(q) * BigInt(p.toBigInt()), BigInt(p.toBigInt()) * BigInt(p.toBigInt()));
  return {A1: A1, B1: B1, U1: U1.mod(q.multiply(p)), V1: V1.mod(q.multiply(p))};
}
const useQuadraticHenselLift = true;
function HenselLiftingOfTwoFactors(C, A, B, p, k) {
  if (useQuadraticHenselLift && k === Math.pow(2, Math.ceil(Math.log2(k)))) { // TODO: any degree
    //TODO: ???
    let q = p;
    const tmp1 = ExtendedEuclideanAlgorithm(A, B, p);
    console.assert(tmp1.gcd.getDegree() === 0);
    let U = tmp1.U;
    let V = tmp1.V;
    for (let i = 1; i < k / 2; i *= 2) {
      const tmp = QuadraticHenselLift(C, A, B, U, V, q, p);
      A = tmp.A1;
      B = tmp.B1;
      U = tmp.U1;
      V = tmp.V1;
      q = Expression.Integer.fromBigInt(q).multiply(Expression.Integer.fromBigInt(q)).toBigInt();
      p = Expression.Integer.fromBigInt(p).multiply(Expression.Integer.fromBigInt(p)).toBigInt();
    }
    const tmp = HenselLift(C, A, B, U, V, q, p);
    A = tmp.A1;
    B = tmp.B1;
    return {A1: A, B1: B};
  }
  //TODO: ?
  let q = p;
  for (let i = 1; i < k; i += 1) {
    const tmp1 = ExtendedEuclideanAlgorithm(A, B, p);
    console.assert(tmp1.gcd.getDegree() === 0);
    const tmp = HenselLift(C, A, B, tmp1.U, tmp1.V, q, p);
    A = tmp.A1;
    B = tmp.B1;
    q = Expression.Integer.fromBigInt(q).multiply(Expression.Integer.fromBigInt(p)).toBigInt();
  }
  return {A1: A, B1: B};
}
function product(factors) {
  console.assert(factors.length > 0);
  return factors.length > 1 ? product(factors.slice(0, Math.ceil(factors.length / 2))).multiply(product(factors.slice(Math.ceil(factors.length / 2)))) : factors[0];
}
function HenselLifting(f, factors, p, e) {
  // https://scholar.rose-hulman.edu/cgi/viewcontent.cgi?article=1163&context=math_mstr
  // "2.3 Factoring mod p e: Hensel Lifting"
  let C = f;
  let newFactors = [];
  const c = C.getLeadingCoefficient().modulo(Expression.Integer.fromBigInt(p));
  if (c.compareTo(Expression.ONE) !== 0) {
    factors = factors.concat([Polynomial.of(c)]);
  }
  if (true && factors.length > 1) {
    //TODO: check this code !!!
    // divide and conquer
    let A = factors.slice(0, Math.ceil(factors.length / 2));
    let B = factors.slice(Math.ceil(factors.length / 2));
    const tmp = HenselLiftingOfTwoFactors(C, product(A), product(B), p, e);
    if (c.compareTo(Expression.ONE) !== 0) {
      B = B.slice(0, -1);
    }
    return HenselLifting(tmp.A1, A, p, e).concat(HenselLifting(tmp.B1, B, p, e));
  }
  for (let i = 0; i < factors.length - 1; i += 1) {
    const A = factors[i];
    const tmp = HenselLiftingOfTwoFactors(C, A, product(factors.slice(i + 1)), p, e);
    newFactors.push(tmp.A1);
    C = tmp.B1;
  }
  if (C.getDegree() > 0) {
    newFactors.push(C);
  } else {
    console.assert(c.compareTo(Expression.ONE) !== 0);
  }
  const pInE = Expression.Integer.fromBigInt(p)._pow(e);
  newFactors = newFactors.map(factor => toMonic(factor.mod(pInE), pInE));//TODO: ?
  return newFactors;
}

//  import factorizeOverTheIntegers from './polynomialFactorization.js';
Polynomial.prototype._factorizeOverTheIntegers = function () {
  //return factorizeOverTheIntegers(this).next().value;
  return factorizeOverTheIntegers(this);
};
Polynomial._gcdOfPolynomialsOverFiniteField = gcdOfPolynomialsOverFiniteField;//TODO: ?

export default factorizeOverTheIntegers;

factorizeOverTheIntegers.testables = {
  gcdOfPolynomialsOverFiniteField: gcdOfPolynomialsOverFiniteField,
  distinctDegreeFactorization: distinctDegreeFactorization,
  CantorZassenhausAlgorithm: CantorZassenhausAlgorithm,
  factorizeOverTheFiniteField: factorizeOverTheFiniteField,
  ExtendedEuclideanAlgorithm: ExtendedEuclideanAlgorithm,
  HenselLift: HenselLift,
  HenselLifting: HenselLifting,
  randomBigInt: randomBigInt
};
