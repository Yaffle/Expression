//@external("env", "consoleLog")
//declare function consoleLog(a:u64): void;

// mulhu from https://doc.lagout.org/security/Hackers%20Delight.pdf
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

function montgomeryPow(b:u64, e:u64, m:u64, mInv:u64):u64 {
  let a = montgomeryTransform(1, m);
  while (e !== 0) {
    if (e % 2 !== 0) {
      a = montgomeryMultiply(a, b, m, mInv);
    }
    e = e >> 1;
    b = montgomeryMultiply(b, b, m, mInv);
  }
  return a;
}

// https://en.wikipedia.org/wiki/Miller–Rabin_primality_test#Miller–Rabin_test
export function isProbablyPrime64(n:u64):boolean {
  // only Miller-Rabin test for base = 2
  if (n < 2) {
    return false;
  }
  if (n % 2 === 0) {
    return n === 2;
  }
  const s = i64.ctz(n - 1);
  const d = (n - 1) >> s;
  const nInv = modInverseMod2pow64(n);
  const one = montgomeryTransform(1, n);
  if (true) {
    const base = 2;
    let x = montgomeryPow(montgomeryTransform(base, n), d, n, nInv);//!!!
    //consoleLog(montgomeryReduce(0, x, n, nInv));
    let y = u64(0);
    for (let i = 0; i < s; i += 1) {
      y = montgomeryMultiply(x, x, n, nInv);
      if (y === one && x !== one && x !== n - one) {
        return false;
      }
      x = y;
      //consoleLog(montgomeryReduce(0, x, n, nInv));
    }
    if (y !== one) {
      return false;
    }
  }
  return true;
}
