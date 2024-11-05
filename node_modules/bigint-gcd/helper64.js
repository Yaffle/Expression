
export function u64gcd(a:u64, b:u64):u64 {
  if (b !== 0) {
    do {
      const r = a % b;
      a = b;
      b = r;
    } while (b !== 0);
  }
  return a;
}

export function u64gcdext(a:u64, b:u64):u64 {
  let A = i64(1);
  let B = i64(0);
  let C = i64(0);
  let D = i64(1);
  if (b !== 0) {
    do {
      const q = a / b;
      const b1 = a - q * b;
      a = b;
      b = b1;
      const C1 = i64(A - i64(q * C));
      const D1 = i64(B - i64(q * D));
      A = C;
      B = D;
      C = C1;
      D = D1;
    } while (b !== 0);
  }
  gA = A;
  gB = B;
  return a;
}

export function helper(x:u64, xlo:u64, y:u64, ylo:u64, lobits:i32):i32 {
  // computes the transformation matrix, which is the product of all {{0, 1}, {1, -q_i}} matrices,
  // where q_i are the quotients produced by Euclidean algorithm for any pair of integers (a, b),
  // where a within [x, x + 1] and b within [y, y + 1]

  // 2x2-matrix transformation matrix of (x_initial, y_initial) into (x, y):
  let A = i64(1);
  let B = i64(0);
  let C = i64(0);
  let D = i64(1);

  let i = 0;

  let bits = 0;
  if (u64(y) != u64(-1) && u64(x) != u64(-1)) { // overflow
    do {
      // switch from a matrix to pairs of (xmin,xmax) and (ymin,ymax):
      // any pair of initial integers looks like: (x_initial + alpha, y_initial + beta), where 0 <= alpha < 1 and 0 <= beta < 1
      // if we multiply transformation by it we get (x + alpha * A + beta * B, y + alpha * C + beta * D)
      // as A,B and C,D have different signs, and the signs are changed after every Euclidean step, then:
      //console.assert(sign(A) !== sign(B) && sign(C) !== sign(D));
      //console.assert(sign(A) !== sign(C) && sign(B) !== sign(D));
      let xmin = u64(x + B);
      let xmax = u64(x + A);
      let ymin = u64(y + C);
      let ymax = u64(y + D);
      if ((i & 1) != 0) {
        const xmin0 = xmin;
        xmin = xmax;
        xmax = xmin0;
        const ymin0 = ymin;
        ymin = ymax;
        ymax = ymin0;
      }
      do {
        const q = u64(xmin / ymax);
        // The quotient for any pair (x,y) is within floor(xmin / ymax) and floor(xmax / ymin) as x > 0 and y > 0
        // So we need to check that u64(xmax / ymin) == q if q === u64(xmax / ymin):
        // 0 <= xmax - q * ymin < ymin
        if (u64(xmax - u64(q * ymin)) >= u64(ymin)) {
          // not same quotient
          break;
        }
        // continue Euclidean step:
        i = i32(i + 1);
        const ymin1 = u64(xmin - u64(q * ymax));
        const ymax1 = u64(xmax - u64(q * ymin));
        const y1 = u64(x - u64(q * y));
        xmin = ymin;
        xmax = ymax;
        x = y;
        ymin = ymin1;
        ymax = ymax1;
        y = y1;
        //console.debug(q);
      } while (true);
      // switch back to the matrix:
      A = i64(xmax - x);
      B = i64(xmin - x);
      C = i64(ymin - y);
      D = i64(ymax - y);
      if ((i & 1) != 0) {
        const A0 = A;
        A = B;
        B = A0;
        const C0 = C;
        C = D;
        D = C0;
      }

      // add more bits from xlo and ylo:
      bits = i32(i64.clz(i64(u64((u64(x + A) > u64(x + B) ? u64(x + A) : u64(x + B)))))); // xmax ?
      bits = i32(lobits) < i32(bits) ? lobits : bits;
      if (i32(bits) != i32(0)) {
        const s = i32(lobits - bits);
        const xlo1 = u64(u64(xlo) >>> u64(s));
        const ylo1 = u64(u64(ylo) >>> u64(s));
        xlo = u64(xlo - u64(xlo1 << u64(s)));
        ylo = u64(ylo - u64(ylo1 << u64(s)));
        x = u64(u64(u64(A * xlo1) + u64(B * ylo1)) + u64(x << u64(bits)));
        y = u64(u64(u64(C * xlo1) + u64(D * ylo1)) + u64(y << u64(bits)));
        lobits = s;
      }
    } while (i32(bits) != i32(0));
  }
  // AssemblyScript does not support multi-value return.
  // Performance of it is similar to the used method. Perhaps it would be faster to use the memory.
  gA = A;
  gB = B;
  gC = C;
  gD = D;
  return 0;
}

let gA = i64(0);
let gB = i64(0);
let gC = i64(0);
let gD = i64(0);
export function A():i64 {
  return gA;
}
export function B():i64 {
  return gB;
}
export function C():i64 {
  return gC;
}
export function D():i64 {
  return gD;
}
