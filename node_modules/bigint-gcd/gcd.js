/*jshint esversion:11*/

const SUBQUADRATIC_GCD_THRESHOLD = 32 * 1024;
const SUBQUADRATIC_HALFGCD_THRESHOLD = 4096;
const DOUBLE_DIGIT_METHOD = true;
const USE_HALF_EXTENDED = true;

let SMALL_GCD_MAX = 0n;
let DIGITSIZE = 0;

let smallgcd = null;
let smallgcdext = null;

let wasmHelper = null;
let jsHelper = null;

if (true && typeof WebAssembly !== 'undefined') {
  const wasmCode = new Uint8Array([0,97,115,109,1,0,0,0,1,20,3,96,0,1,126,96,2,126,126,1,126,96,5,126,126,126,126,127,1,127,3,8,7,1,1,2,0,0,0,0,6,21,4,126,1,66,0,11,126,1,66,0,11,126,1,66,0,11,126,1,66,0,11,7,47,7,6,117,54,52,103,99,100,0,0,9,117,54,52,103,99,100,101,120,116,0,1,6,104,101,108,112,101,114,0,2,1,65,0,3,1,66,0,4,1,67,0,5,1,68,0,6,10,141,4,7,37,1,1,126,32,1,66,0,82,4,64,3,64,32,0,32,1,130,33,2,32,1,33,0,32,2,34,1,66,0,82,13,0,11,11,32,0,11,97,1,6,126,66,1,33,6,66,1,33,3,32,1,66,0,82,4,64,3,64,32,0,32,0,32,1,128,34,7,32,1,126,125,33,4,32,1,33,0,32,6,32,2,32,7,126,125,33,1,32,5,32,3,32,7,126,125,33,7,32,2,33,6,32,3,33,5,32,1,33,2,32,7,33,3,32,4,34,1,66,0,82,13,0,11,11,32,6,36,0,32,5,36,1,32,0,11,238,2,2,7,126,2,127,66,1,33,9,66,1,33,7,32,0,66,127,82,32,2,66,127,82,113,4,64,3,64,32,0,32,6,124,33,8,32,0,32,9,124,33,6,32,2,32,5,124,33,5,32,2,32,7,124,33,7,32,12,65,1,113,4,64,32,8,33,9,32,6,33,8,32,9,33,6,32,5,33,9,32,7,33,5,32,9,33,7,11,3,64,32,5,32,6,32,8,32,7,128,34,11,32,5,126,125,34,9,86,4,64,32,12,65,1,106,33,12,32,8,32,7,32,11,126,125,33,10,32,0,32,2,32,11,126,125,33,11,32,5,33,8,32,7,33,6,32,2,33,0,32,10,33,5,32,9,33,7,32,11,33,2,12,1,11,11,32,6,32,0,125,33,9,32,8,32,0,125,33,6,32,5,32,2,125,33,5,32,7,32,2,125,33,7,32,12,65,1,113,4,64,32,9,33,8,32,6,33,9,32,8,33,6,32,5,33,8,32,7,33,5,32,8,33,7,11,32,4,32,0,32,9,124,34,8,32,0,32,6,124,34,10,32,8,32,10,86,27,121,167,34,13,32,4,32,13,72,27,34,13,4,64,32,1,32,1,32,4,32,13,107,34,4,172,34,8,136,34,10,32,8,134,125,33,1,32,3,32,3,32,8,136,34,11,32,8,134,125,33,3,32,9,32,10,126,32,6,32,11,126,124,32,0,32,13,172,34,8,134,124,33,0,32,5,32,10,126,32,7,32,11,126,124,32,2,32,8,134,124,33,2,11,32,13,13,0,11,11,32,9,36,0,32,6,36,1,32,5,36,2,32,7,36,3,65,0,11,4,0,35,0,11,4,0,35,1,11,4,0,35,2,11,4,0,35,3,11]);
  try {
    const exports = new WebAssembly.Instance(new WebAssembly.Module(wasmCode)).exports;
    if (exports.helper(1n, 0n, 1n, 0n) != null) {
      wasmHelper = function (x, xlo, y, ylo, lobits) {
        exports.helper(x, xlo, y, ylo, lobits);
        return [exports.A(), exports.B(), exports.C(), exports.D()];
      };
      DIGITSIZE = 64;
      SMALL_GCD_MAX = BigInt.asUintN(64, -1n);
    }
    if (exports.u64gcd(0n, 0n) === 0n) {
      smallgcd = exports.u64gcd;
    }
    if (exports.u64gcdext(0n, 0n) === 0n) {
      smallgcdext = function (a, b) {
        const g = exports.u64gcdext(a, b);
        return [exports.A(), exports.B(), g];
      };
    }
  } catch (error) {
    console.log(error);
  }
}


function AsmModule(stdlib) {
  "use asm";
  
  var floor = stdlib.Math.floor;
  var max = stdlib.Math.max;
  var clz32 = stdlib.Math.clz32;

  var gA = -0.0;
  var gB = -0.0;
  var gC = -0.0;
  var gD = -0.0;

function f64gcd(a, b) {
  a = +a;
  b = +b;
  var b1 = -0.0;
  var q = -0.0;
  while (b > -0.0) {
    q = +floor(a / b);
    b1 = a - q * b;
    a = b;
    b = b1;
  }
  return +a;
}

function f64gcdext(a, b) {
  a = +a;
  b = +b;
  var b1 = -0.0;
  var q = -0.0;
  var A = 1.0;
  var B = -0.0;
  var C = -0.0;
  var D = 1.0;
  var C1 = -0.0;
  var D1 = -0.0;
  while (b > -0.0) {
    q = +floor(a / b);
    b1 = a - q * b;
    a = b;
    b = b1;
    C1 = A - q * C;
    D1 = B - q * D;
    A = C;
    B = D;
    C = C1;
    D = D1;
  }
  gA = A;
  gB = B;
  return +a;
}

// 1 + floor(log2(x))
function log2(x) {
  x = +x;
  var e = 0;
  while (x >= 4294967296.0) {
    x = x * 2.3283064365386963e-10;
    e = (e + 32) | 0;
  }
  e = (e + (32 - (clz32(~~x) | 0))) | 0;
  return e | 0;
}

// 2**n
function exp2(n) {
  n = n | 0;
  var result = 1.0;
  while ((n | 0) < 0) {
    n = (n + 32) | 0;
    result = result * 2.3283064365386963e-10; // * 2**-32
  }
  while ((n | 0) >= 32) {
    n = (n - 32) | 0;
    result = result * 4294967296.0; // * 2**32
  }
  result = result * +((1 << n) >>> 0);
  return +result;
}

// @Deprecated, see helper64.js
function jsHelper(x, xlo, y, ylo, lobits) {
  x = +x;
  xlo = +xlo;
  y = +y;
  ylo = +ylo;
  lobits = lobits | 0;

  var A = 1.0;
  var B = -0.0;
  var C = -0.0;
  var D = 1.0;

  var bits = 0;
  var sameQuotient = 0;
  var q = -0.0;
  var y1 = -0.0;
  var A1 = -0.0;
  var B1 = -0.0;
  var C1 = -0.0;
  var D1 = -0.0;
  var b = 0;
  var d = -0.0;
  var dInv = -0.0;
  var xlo1 = -0.0;
  var ylo1 = -0.0;
  var p = -0.0;
  if (y != -0.0) {
    do {
      do {
        q = floor(x / y);
        y1 = x - q * y;
        A1 = C;
        B1 = D;
        C1 = A - q * C;
        D1 = B - q * D;

        // The quotient for a point (x_initial + alpha, y_initial + beta), where 0 <= alpha < 1 and 0 <= beta < 1:
        // floor((x + A * alpha + B * beta) / (y + C * alpha + D * beta))
        // As the sign(A) === -sign(B) === -sign(C) === sign(D) (ignoring zero entries) the maximum and minimum values are floor((x + A) / (y + C)) and floor((x + B) / (y + D))

        // floor((x + A) / (y + C)) === q  <=>  0 <= (x + A) - q * (y + C) < (y + C)  <=>  0 <= y1 + C1 < y + C
        // floor((x + B) / (y + D)) === q  <=>  0 <= (x + B) - q * (y + D) < (y + D)  <=>  0 <= y1 + D1 < y + D
        sameQuotient = (-0.0 <= y1 + C1) & (y1 + C1 < y + C) &
                       (-0.0 <= y1 + D1) & (y1 + D1 < y + D);
        if (sameQuotient) {
          x = y;
          y = y1;
          A = A1;
          B = B1;
          C = C1;
          D = D1;
        }
      } while (sameQuotient);

      b = (53 - (log2(x + max(A, B)) | 0)) | 0;
      bits = (b | 0) < 0 ? b : ((b | 0) > (lobits | 0) ? lobits : b);
      if ((b | 0) != 0) {
        d = +exp2((lobits - bits) | 0);
        dInv = +exp2((bits - lobits) | 0);
        xlo1 = +floor(xlo * dInv);
        ylo1 = +floor(ylo * dInv);
        xlo = xlo - xlo1 * d;
        ylo = ylo - ylo1 * d;
        lobits = (lobits - bits) | 0;
        p = +exp2(bits);
        x = A * xlo1 + B * ylo1 + x * p;
        y = C * xlo1 + D * ylo1 + y * p;
      }

    } while ((bits | 0) != 0);
  }
  gA = A;
  gB = B;
  gC = C;
  gD = D;
  return 0;
}

  function A() {
    return gA;
  }
  function B() {
    return gB;
  }
  function C() {
    return gC;
  }
  function D() {
    return gD;
  }

  return {f64gcdext: f64gcdext, f64gcd: f64gcd, helper: jsHelper, A: A, B: B, C: C, D: D};
}

if (DIGITSIZE === 0) {
  const asmExports = AsmModule(globalThis);
  jsHelper = function (x, xlo, y, ylo, lobits) {
    asmExports.helper(x, xlo, y, ylo, lobits);
    return [BigInt(asmExports.A()), BigInt(asmExports.B()), BigInt(asmExports.C()), BigInt(asmExports.D())];
  };
  DIGITSIZE = 53;
  SMALL_GCD_MAX = BigInt.asUintN(53, -1n);
  smallgcd = function (a, b) {
    return BigInt(asmExports.f64gcd(Number(BigInt(a)), Number(BigInt(b))));
  };
  smallgcdext =  function (a, b) {
    const g = asmExports.f64gcdext(Number(BigInt(a)), Number(BigInt(b)));
    return [BigInt(asmExports.A()), BigInt(asmExports.B()), BigInt(g)];
  };
}

// https://github.com/tc39/proposal-bigint/issues/205
// https://github.com/tc39/ecma262/issues/1729
// floor(log2(a)) + 1 if a > 0
function bitLength(a) {
  const s = a.toString(16);
  const c = s.charCodeAt(0) - 0 - '0'.charCodeAt(0);
  if (c <= 0) {
    throw new RangeError();
  }
  return (s.length - 1) * 4 + (32 - Math.clz32(Math.min(c, 8)));
}

const frexpf64 = typeof Float64Array !== 'undefined' ? new Float64Array(1) : null;
const frexpi32 = typeof Float64Array !== 'undefined' ? new Int32Array(frexpf64.buffer) : null;

let previousValue = 0;
// some terrible optimization as bitLength is slow
function bitLength2(a) {
  if (previousValue <= 1024) {
    const n = -0.0 + Number(BigInt(a));
    if (frexpf64 != null) {
      frexpf64[0] = n;
      const e = (frexpi32[1] >> 20) - 1023;
      if (e < 1024 && frexpi32[0] !== 0 || (frexpi32[1] & 0xFFFFF) !== 0) {
        previousValue = e + 1;
        return previousValue;
      }
    }
    const x = Math.log2(n) + 1024 * 4 - 1024 * 4;
    const y = Math.ceil(x);
    if (x !== y) {
      previousValue = y;
      return previousValue;
    }
  }
  if (previousValue < DIGITSIZE) {
    previousValue = DIGITSIZE;
  }
  const n = -0.0 + Number(a >> BigInt(previousValue - DIGITSIZE));
  if (n >= 1 && n <= 9007199254740992) { // 2**53
    let x = +n;
    let e = 0;
    while (x > +(1 << 30)) {
      x = Math.floor(x / +(1 << 30));
      e += 30;
    }
    e += (32 - Math.clz32(x));
    previousValue = previousValue - DIGITSIZE + e;
    return previousValue;
  }
  previousValue = bitLength(a);
  return previousValue;
}


function helper(X, Y) {
  if (typeof X !== 'bigint' || typeof Y !== 'bigint') {
    throw new TypeError();
  }
  if (wasmHelper != null) {
    if (DIGITSIZE !== 64) {
      throw new RangeError();
    }
    if (!DOUBLE_DIGIT_METHOD) {
      return wasmHelper(X, 0n, Y, 0n, 0);
    }
    const x = BigInt.asUintN(64, X >> 64n);
    const xlo = BigInt.asUintN(64, X);
    const y = BigInt.asUintN(64, Y >> 64n);
    const ylo = BigInt.asUintN(64, Y);
    return wasmHelper(x, xlo, y, ylo, 64);
  } else {
    if (DIGITSIZE !== 53) {
      throw new RangeError();
    }
    if (!DOUBLE_DIGIT_METHOD) {
      return jsHelper(-0.0 + Number(X), 0, -0.0 + Number(Y), 0, 0);
    }
    const x = -0.0 + Number(X >> 53n);
    const xlo = -0.0 + Number(BigInt.asUintN(53, X));
    const y = -0.0 + Number(Y >> 53n);
    const ylo = -0.0 + Number(BigInt.asUintN(53, Y));
    return jsHelper(x, xlo, y, ylo, 53);
  }
}

//TODO: remove those comments:
  // floor((a + 1) / b) < q = floor(a / b) < floor(a / (b + 1))
  // ([A, B], [C, D]) * (a + x, b + y) = (A*(a+x)+B*(b+y), C*(a+x)+D*(b+y)) = (A*a+B*b, C*a+D*b) + (A*x+B*y, C*x+D*y)

  //console.assert(A * D >= 0 && B * C >= 0 && A * B <= 0 && D * C <= 0);//TODO: why - ?
  //const [a1, b1] = [a + A, b + C]; // T * (a_initial + 1n, b_initial);
  //const [a2, b2] = [a + B, b + D]; // T * (a_initial, b_initial + 1n);
  //if (!isSmall && n <= size * (2 / 3)) { // TODO: ?, the constant is based on some testing with some example
  //  return [A, B, C, D, a, b];
  //}

function abs(a) {
  if (typeof a !== 'bigint') {
    throw new TypeError();
  }
  return a < 0n ? -a : a;
}
function max(a, b) {
  if (typeof a !== 'bigint' || typeof b !== 'bigint') {
    throw new TypeError();
  }
  return a < b ? b : a;
}

function halfgcd(a, b, extended = true, reallyhalfgcd = true, wrapper = false) {
  if (typeof a !== 'bigint' || typeof b !== 'bigint') {
    throw new TypeError();
  }
  extended = extended || reallyhalfgcd;

  // 2x2 matrix:
  let A = 1n;
  let B = 0n;
  let C = 0n;
  let D = 1n;
  let step = 0;

  //TODO: TEST
  if (a < 0n) {
    a = -a;
    if (extended) {
      A = -A;
      B = -B;
      step += 1;
    }
  }
  //TODO: TEST
  if (b < 0n) {
    b = -b;
    if (extended) {
      C = -C;
      D = -D;
      step += 1;
    }
  }
  if (a < b) {
    const tmp = a;
    a = b;
    b = tmp;
    if (extended) {
      const C1 = A;
      A = C;
      C = C1;
      const D1 = B;
      B = D;
      D = D1;
      step += 1;
    }
  }

  // The function calculates the transformation matrix for numbers (x, y), where a <= x < a + 1 and b <= y < b + 1
  // Seems, this definition is not the same as in https://mathworld.wolfram.com/Half-GCD.html
  // Note: for debugging it is useful to compare the quotients in the simple Euclidean algorithm vs the quotients there
  // see helper64.js for the properties used

  let isSmall = false;
  if (reallyhalfgcd) {
    if (BigInt.asUintN(SUBQUADRATIC_HALFGCD_THRESHOLD, a) === a) {
      isSmall = true;
    }
  }

  let isVerySmall = false;
  let sizea0 = 0;

  while ((reallyhalfgcd || a > SMALL_GCD_MAX) && b !== 0n) {
    //console.assert(a >= b);
    step += 1;

    if (!isSmall) {
      // Subquadratic Lehmer's algorithm:
      if (!reallyhalfgcd) {
        if (BigInt.asUintN(SUBQUADRATIC_GCD_THRESHOLD * (extended ? 1 / 16 : 1), b) === b) {
          isSmall = true;
          continue;
        }
      }
      if (reallyhalfgcd && step === 1) {
        sizea0 = bitLength(a);
      }
      const s1 = reallyhalfgcd ? bitLength(max(abs(C), abs(D))) : 0;
      let m = reallyhalfgcd ? Math.max(0, Math.ceil((sizea0 - s1 - s1) * (1 / 2))) : (extended ? 0 : Math.floor(bitLength(a) * 2 / 3)); // 2/3 is somehow faster
      let M = BigInt(m);
      if (step !== 1 && reallyhalfgcd) {
        if (m < 256) {
          //if (s1 / sizea0 < 0.46) {
          //  console.debug(s1 / sizea0);
          //}
          isSmall = true;
          continue;
        }
        //if (s1 / sizea0 > 0.251) {
        //  console.debug(s1 / sizea0);
        //}
        let k = 8;
        while (((a + A) >> M) !== ((a + B) >> M) ||
               ((b + C) >> M) !== ((b + D) >> M)) {
          // min and max values have different high bits
          //isSmall = true;
          //continue;
          m += k;
          M = BigInt(m);
          k += k;
        }
      }
      const [$A1, $B1, $C1, $D1, $ahi1, $bhi1] = halfgcd(a >> M, b >> M);
      const [A1, B1, C1, D1, ahi1, bhi1] = [BigInt($A1), BigInt($B1), BigInt($C1), BigInt($D1), BigInt($ahi1), BigInt($bhi1)];
      //if (typeof A1 !== 'bigint' || typeof B1 !== 'bigint' || typeof C1 !== 'bigint' || typeof D1 !== 'bigint' || typeof ahi1 !== 'bigint' || typeof bhi1 !== 'bigint') {
        //throw new TypeError();
      //}
      if (B1 !== 0n) {
        if (extended) {
          // T := T1 * T:
          if (step === 1) {
            A = A1;
            B = B1;
            C = C1;
            D = D1;
          } else {
            const B2 = A1 * B + B1 * D;
            const D2 = C1 * B + D1 * D;
            B = B2;
            D = D2;
            if (!USE_HALF_EXTENDED || reallyhalfgcd) {
              const A2 = A1 * A + B1 * C;
              const C2 = C1 * A + D1 * C;
              A = A2;
              C = C2;
            }
          }
        }
        const alo = BigInt.asUintN(m, a);
        const blo = BigInt.asUintN(m, b);
        // (a, b) := T1 * (alo, blo) + T1 * (ahi, bhi) * 2**m:
        const a1 = (A1 * alo + B1 * blo) + (ahi1 << M);
        const b1 = (C1 * alo + D1 * blo) + (bhi1 << M);
        a = a1;
        b = b1;
        if (a < 0n || b < 0n) {
          throw new TypeError("assertion");
        }
        continue;
      } else {
        //console.assert(A1 === 1n && B1 === 0n && C1 === 0n && D1 === 1n);
      }
    }
    if (isSmall && !isVerySmall) {
      // Lehmer's algorithm:
      let m = Math.max(0, bitLength2(a) - DIGITSIZE * (DOUBLE_DIGIT_METHOD ? 2 : 1));
      let M = m === 0 ? 0n : BigInt(m);
      if (step !== 1 && reallyhalfgcd) {
        if (((a + A) >> M) !== ((a + B) >> M) ||
            ((b + C) >> M) !== ((b + D) >> M)) {
          // min and max values have different high bits
          if (!wrapper) {
            break;
          }
          do {
            m += 8;
            M = BigInt(m);
          } while (((a + A) >> M) !== ((a + B) >> M) ||
                   ((b + C) >> M) !== ((b + D) >> M));
          if ((b >> M) === 0n) {
            isVerySmall = true;
            continue;
          }
        }
      }
      const [$A1, $B1, $C1, $D1] = helper((m === 0 ? a : a >> M), (m === 0 ? b : b >> M));
      const [A1, B1, C1, D1] = [BigInt($A1), BigInt($B1), BigInt($C1), BigInt($D1)];
      //if (typeof A1 !== 'bigint' || typeof B1 !== 'bigint' || typeof C1 !== 'bigint' || typeof D1 !== 'bigint') {
        //throw new TypeError();
      //}
      if (B1 !== 0n) {
        if (extended) {
          // T := T1 * T:
          if (step === 1) {
            A = A1;
            B = B1;
            C = C1;
            D = D1;
          } else {
            const B2 = A1 * B + B1 * D;
            const D2 = C1 * B + D1 * D;
            B = B2;
            D = D2;
            if (!USE_HALF_EXTENDED || reallyhalfgcd) {
              const A2 = A1 * A + B1 * C;
              const C2 = C1 * A + D1 * C;
              A = A2;
              C = C2;
            }
          }
        }
        // (a, b) := T1 * (a, b):
        const a1 = A1 * a + B1 * b;
        const b1 = C1 * a + D1 * b;
        a = a1;
        b = b1;
        if (a < 0n || b < 0n) {
          throw new TypeError("assertion");
        }
        continue;
      } else {
        //console.assert(A1 === 1n && B1 === 0n && C1 === 0n && D1 === 1n);
      }
    }

    //console.debug('%');
    const q = a / b;
    const b1 = a - q * b;
    if (extended) {
      const C1 = A - q * C;
      const D1 = B - q * D;
      if (reallyhalfgcd) {
        const sameQuotient = 0n <= b1 + C1 && b1 + C1 < b + C &&
                             0n <= b1 + D1 && b1 + D1 < b + D;
        if (!sameQuotient) {
          break;
        }
      }
      // T := {{0, 1}, {1, -q}} * T:
      A = C;
      B = D;
      C = C1;
      D = D1;
    }
    // (a, b) := {{0, 1}, {1, -q}} * (a, b)
    a = b;
    b = b1;
    //gcd.debug(q);
  }
  // see "2. General structure of subquadratic gcd algorithms" in “On Schönhage’s algorithm and subquadratic integer GCD computation” by Möller
  return [A, B, C, D, a, b]; // for performance transformedA and transformedB are returned
}

// https://en.wikipedia.org/wiki/Lehmer%27s_GCD_algorithm
// https://www.imsc.res.in/~kapil/crypto/notes/node11.html
function LehmersGCD(a, b) {
  const [A1, B1, C1, D1, a1, b1] = halfgcd(a, b, false, false);
  a = BigInt(a1);
  b = BigInt(b1);
  if (b !== 0n) {
    a = BigInt.asUintN(64, smallgcd(a, b));
  }
  return a;
}

function LehmersGCDExt(a, b) {
  const [A1, B1, C1, D1, a1, b1] = halfgcd(a, b, true, false);
  const a0 = a;
  const b0 = b;
  let A = BigInt(A1);
  let B = BigInt(B1);
  let C = BigInt(C1);
  let D = BigInt(D1);
  a = BigInt(a1);
  b = BigInt(b1);
  if (b1 !== 0n) {
    const [A1, B1, g] = smallgcdext(a, b);
    a = BigInt.asUintN(64, g);
    b = 0n;
    B = A1 * B + B1 * D;
    if (!USE_HALF_EXTENDED) {
      A = A1 * A + B1 * C;
    }
  }
  if (USE_HALF_EXTENDED) {
    // A*a + B*b = g
    A = a0 === 0n ? 0n : (a - B * b0) / a0; // exact division
  }
  return [A, B, a];
}

function gcd(a, b) {
  return LehmersGCD(BigInt(a), BigInt(b));
}

function gcdext(a, b) {
  return LehmersGCDExt(BigInt(a), BigInt(b));
}

function halfgcdWrapper(a, b) {
  // reduce numbers as much as possible:
  return halfgcd(a, b, true, true, true);
}

export default gcd;

gcd.halfgcd = halfgcdWrapper;//TODO:?
gcd.gcdext = gcdext;//TODO:?
