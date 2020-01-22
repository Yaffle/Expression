  import Expression from './Expression.js';

  var Integer = Expression.Integer;
  var Addition = Expression.Addition;
  var Multiplication = Expression.Multiplication;
  var Division = Expression.Division;
  var Exponentiation = Expression.Exponentiation;
  var BinaryOperation = Expression.BinaryOperation;

var separateSinCos = function (e) {
  if (!(e instanceof Multiplication)) {
    throw new TypeError();
  }
  var sinCos = undefined;
  var other = undefined;
  var x = e;
  for (var multiplications = x.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
    var v = y;
    if (v instanceof Sin || v instanceof Cos || 
        (v instanceof Exponentiation && (v.a instanceof Sin || v.a instanceof Cos))) {
      sinCos = sinCos == undefined ? v : sinCos.multiply(v);
    } else {
      other = other == undefined ? v : other.multiply(v);
    }
  }
  return {
    sinCos: sinCos == undefined ? Expression.ONE : sinCos,
    other: other == undefined ? Expression.ONE : other
  };
};

var expandMainOp = function (u) {
  return u;
};

var contractTrigonometryInternal = function (a, b) {
  // sin(a) * sin(b) = (cos(a - b) - cos(a + b)) / 2
  // sin(a) * cos(b) = (sin(a + b) + sin(a - b)) / 2
  // cos(a) * sin(b) = (sin(a + b) - sin(a - b)) / 2
  // cos(a) * cos(b) = (cos(a - b) + cos(a + b)) / 2
  var ax = a.a;
  var bx = b.a;
  if (a instanceof Sin && b instanceof Sin) {
    return ax.subtract(bx).cos().divide(Expression.TWO).subtract(ax.add(bx).cos().divide(Expression.TWO));
  }
  if (a instanceof Sin && b instanceof Cos) {
    return ax.add(bx).sin().divide(Expression.TWO).add(ax.subtract(bx).sin().divide(Expression.TWO));
  }
  if (a instanceof Cos && b instanceof Sin) {
    return ax.add(bx).sin().divide(Expression.TWO).subtract(ax.subtract(bx).sin().divide(Expression.TWO));
  }
  if (a instanceof Cos && b instanceof Cos) {
    return ax.subtract(bx).cos().divide(Expression.TWO).add(ax.add(bx).cos().divide(Expression.TWO));
  }
  throw new TypeError();
};

// page 306
var contractTrigonometryPower = function (u) {
  var b = u.a;
  if (!(b instanceof Sin) && !(b instanceof Cos)) {
    return u;
  }
  var e = contractTrigonometryInternal(b, b).multiply(u.divide(b.multiply(b)));
  return contractTrigonometryRules(e.getNumerator()).divide(e.getDenominator());
};

// page 318
var contractTrigonometryProduct = function (u) {
  var i = u.factors();
  var a = i.next().value;
  var b = i.next().value;
  var rest = Expression.ONE;
  var y = i.next().value;
  while (y != null) {
    rest = y.multiply(rest);//TODO: fix
    y = i.next().value;
  }

  if (a instanceof Exponentiation) {
    a = contractTrigonometryPower(a);
    return contractTrigonometryRules(a.multiply(b).multiply(rest));
  }
  if (b instanceof Exponentiation) {
    b = contractTrigonometryPower(b);
    return contractTrigonometryRules(a.multiply(b).multiply(rest));
  }
  // (a instanceof Sin || a instanceof Cos) && (b instanceof Sin || b instanceof Cos)
  var c = contractTrigonometryInternal(a, b);

  return contractTrigonometryRules(c.multiply(rest));
};

// page 317
var contractTrigonometryRules = function (u) {
  var v = expandMainOp(u);
  if (v instanceof Exponentiation) {
    return contractTrigonometryPower(v);
  }
  if (v instanceof Multiplication) {
    var tmp = separateSinCos(v);
    var c = tmp.other;
    var d = tmp.sinCos;
    if (d.equals(Expression.ONE)) {
      return v;
    }
    if (d instanceof Sin || d instanceof Cos) {
      return v;
    }
    if (d instanceof Exponentiation) {
      return expandMainOp(c.multiply(contractTrigonometryPower(d)));
    }
    if (d instanceof Multiplication) {
      return expandMainOp(c.multiply(contractTrigonometryProduct(d)));
    }
    throw new TypeError();
  }
  if (v instanceof Addition) {
    var s = Expression.ZERO;
    var e = v;
    for (var additions = e.summands(), x = additions.next().value; x != null; x = additions.next().value) {
      if (x instanceof Multiplication || x instanceof Exponentiation) {
        s = s.add(contractTrigonometryRules(x));
      } else {
        s = s.add(x);
      }
    }
    return s;
  }
  return v;
};

var map = function (f, u) {
  if (u instanceof Integer) {
    return f(u);
  }
  if (u instanceof Addition) {
    return f(map(f, u.a).add(map(f, u.b)));
  }
  if (u instanceof Multiplication) {
    return f(map(f, u.a).multiply(map(f, u.b)));
  }
  if (u instanceof Division) {
    return f(map(f, u.a).divide(map(f, u.b)));
  }
  if (u instanceof Exponentiation) {
    return f(map(f, u.a).pow(map(f, u.b)));
  }
  if (u instanceof Sin) {
    return f(map(f, u.a).sin());
  }
  if (u instanceof Cos) {
    return f(map(f, u.a).cos());
  }
  if (u instanceof Expression.Matrix) {
    return new Expression.Matrix(u.matrix.map(function (e, i, j) {
      return map(f, e);
    }));
  }
  if (u instanceof Expression.Polynomial) {//TODO: test case
    return new Expression.Polynomial(u.polynomial.map(function (c, d) {
      return map(f, c);
    }));
  }
  if (u instanceof Expression.GF2Value) {
    return u;
  }
  if (u instanceof Expression.NthRoot) {
    return u;
  }
  if (u instanceof Expression.Negation) {
    return u;//?
  }
  if (u instanceof Expression.Complex) {
    return u;//?
  }
  if (u instanceof Expression.NonSimplifiedExpression) {
    //TODO: fix
    return u;//?
  }
  if (u instanceof Expression.Degrees) {
    return u;//?
  }
  if (u instanceof Expression.Radians) {
    return u;//?
  }
  if (u instanceof Expression.Symbol) {
    return u;//?
  }
  throw new TypeError();
};

Expression._map = map;

// page 303

var expandTrigonometryRulesInternal = function (a, b, type) {
  if (type === "cos") {
    // cos(a + b) = cos(a) * cos(b) - sin(a) * sin(b)
    return expandTrigonometryRules(a, "cos").multiply(expandTrigonometryRules(b, "cos")).subtract(expandTrigonometryRules(a, "sin").multiply(expandTrigonometryRules(b, "sin")));
  }
  if (type === "sin") {
    // sin(a + b) = sin(a) * cos(b) + cos(a) * sin(b)
    return expandTrigonometryRules(a, "sin").multiply(expandTrigonometryRules(b, "cos")).add(expandTrigonometryRules(a, "cos").multiply(expandTrigonometryRules(b, "sin")));
  }
  throw new TypeError(type);
};

var expandTrigonometryRules = function (A, type) {
  if (A instanceof Addition) {
    return expandTrigonometryRulesInternal(A.a, A.b, type);
  } else if (A instanceof Multiplication) {
    var i = Expression.ONE;
    for (var multiplications = A.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
      if (y instanceof Expression.Integer) {
        i = i.multiply(y);
      }
    }
    var a = i;
    var b = A.divide(i);
    //var a = A.a;
    //var b = A.b;
    if (!(a instanceof Integer)) {
      throw new TypeError();
    }
    if (a.equals(Expression.ONE)) {
      if (type === "cos") {
        return A.cos();
      }
      if (type === "sin") {
        return A.sin();
      }
    }
    if (a.compareTo(Expression.ONE.negate()) === 0) {
      if (type === "cos") {
        return expandTrigonometryRules(b, type);
      }
      if (type === "sin") {
        return expandTrigonometryRules(b, type).negate();
      }
    }
    var c = a.compareTo(Expression.ZERO) > 0 ? Expression.ONE : Expression.ONE.negate();
    return expandTrigonometryRulesInternal(c.multiply(b), a.subtract(c).multiply(b), type);
  } else if (A instanceof Division) {
    var t = simplifyConstantValue(A, type);
    if (t != null) {
      return t;
    }
    var a = A.a;
    var b = A.b;
    if (a instanceof Addition) {
      return expandTrigonometryRulesInternal(a.a.divide(b), a.b.divide(b), type);
    }
  }
  if (A instanceof Expression.Symbol || A instanceof Expression.Degrees || A instanceof Expression.Radians) {
    if (type === "cos") {
      return A.cos();
    }
    if (type === "sin") {
      return A.sin();
    }
  }
  throw new TypeError();
};

// CA and SC, EA, p. 303

var expandTrigonometry = function (u) {
  return map(function (v) {
    if (v instanceof Sin) {
      return expandTrigonometryRules(v.a, "sin");
    }
    if (v instanceof Cos) {
      return expandTrigonometryRules(v.a, "cos");
    }
    return v;
  }, u);
};

Expression._expandTrigonometry = expandTrigonometry;//!

var contractTrigonometry = function (u) {
  return map(function (v) {
    if (v instanceof Multiplication || v instanceof Exponentiation || v instanceof Addition) {//! Addition - ?
      return contractTrigonometryRules(v);
    }
    if (v instanceof Division) {
      return contractTrigonometryRules(v.getNumerator()).divide(v.getDenominator());
    }
    return v;
  }, u);
};

// page 323

var hasTrigonometry = function (e) {//TODO: remove
  if (e instanceof BinaryOperation) {
    return hasTrigonometry(e.a) || hasTrigonometry(e.b);
  }
  return e instanceof Cos || e instanceof Sin;
};

var simplifyTrigonometry = function (u) {
  if (!hasTrigonometry(u)) {
    return u;
  }
  //!new
  var v = null;
  var r = null;
  Expression._map(function (e) {
    // sin(x/2) -> sin(t), t = 2x
    if (e instanceof Expression.Sin || e instanceof Expression.Cos) {
      var a = e.a;
      if (a instanceof Division) {
        var n = a.getNumerator();
        var d = a.getDenominator();
        if (!d.equals(Expression.ONE)) {
          for (var additions = n.summands(), x = additions.next().value; x != null; x = additions.next().value) {
            var g = x.gcd(d);
            if (!g.equals(d)) {
              for (var multiplications = x.factors(), y = multiplications.next().value; y != null; y = multiplications.next().value) {
                if (y instanceof Expression.Symbol && y !== Expression.PI) {
                  r = d.divide(g);
                  v = y;
                }
              }
            }
          }
        }
      }
    }
    return e;
  }, u);
  if (v != null && r != null) {
    // sin(x/2) -> sin(x)
    u = Expression._substitute(u, v, v.multiply(r), v.divide(r));
    u = simplifyTrigonometry(u);
    // sin(x) -> sin(x/2)
    u = Expression._substitute(u, v, v.divide(r), v.multiply(r));
    return u;
  }

  //!
  // TODO: https://en.wikipedia.org/wiki/Euler%27s_formula#Relationship_to_trigonometry - is it possible to do this with anoher method?
  var n = u.getNumerator();
  n = expandTrigonometry(n);
  n = contractTrigonometry(n);
  var d = u.getDenominator();
  d = expandTrigonometry(d);
  d = contractTrigonometry(d);
  return n.divide(d);
};

Expression.Division.prototype.compare4Multiplication = function (y) {
  if (y instanceof Division) {
    return this.a.compare4Multiplication(y.a) || this.b.compare4Multiplication(y.b);
  }
  return -1;//TODO:
};
Expression.Division.prototype.compare4MultiplicationSymbol = function () {
  return +1;//TODO:
};


Expression.simplifyTrigonometry = simplifyTrigonometry;//?


function Sin(x) {
  Expression.Function.call(this, "sin", x);
}
Sin.prototype = Object.create(Expression.Function.prototype);

//TODO: new 2017-04-26
var simplifyConstantValueInternal = function (d) {
  if (d >= +360 || d <= -360) {
    throw new RangeError();
  }
  if (d < 0) {
    d += 360;
  }
  if (d >= 180) {
    d = d - 180;
    var tmp = simplifyConstantValueInternal(d);
    return tmp == null ? null : tmp.negate();
  }
  if (d > 90) {
    d = 180 - d;
    var tmp = simplifyConstantValueInternal(d);
    return tmp == null ? null : tmp.negate();
  }

  function f(x) {
    // https://en.wikipedia.org/wiki/Trigonometric_constants_expressed_in_real_radicals#Calculated_trigonometric_values_for_sine_and_cosine
    // cos(x) = sqrt(2+2cos(2x))/2 - sign - ?
    var y = simplifyConstantValueInternal(x * 2);
    return y == null ? null : Expression.TWO.add(Expression.TWO.multiply(y)).squareRoot().divide(Expression.TWO);
  }
  function cosapb(a, b) { // cos(d) = cos(a + b)
    var cosa = simplifyConstantValueInternal(a);
    var cosb = simplifyConstantValueInternal(b);
    var sina = simplifyConstantValueInternal(90 - a);
    var sinb = simplifyConstantValueInternal(90 - b);
    return cosa.multiply(cosb).subtract(sina.multiply(sinb));
  }

  //function cos2x(d) { // cos(d) = cos(a + b)
    // cos(2x) = 2cos^2(x)-1
    //var y = simplifyConstantValueInternal(d / 2);
    //return Expression.TWO.multiply(y.multiply(y)).subtract(Expression.ONE);
  //}


  //if (d === 24) {
    // https://en.wikipedia.org/wiki/Trigonometric_constants_expressed_in_real_radicals#24°:_sum_12°_+_12°
    //return RPN('(sqrt(6*(5-sqrt(5)))+sqrt(5)+1)/8');
    //return cosapb(60, -36);
  //}
  //if (d === 42) {
    // cos(42) = sin(48) = 2*sin(24)*cos(24)
    //return RPN('2*sin(24)*cos(24)');
  //}

  // 0, 15, 30, 36, 45, 60, 72, 75, 90 - more simple

  if (d === 0) {
    return Expression.ONE;
  }
  if (d === 30) {
    return Expression.ONE.add(Expression.TWO).squareRoot().divide(Expression.TWO);
  }
  if (d === 45) {
    return Expression.ONE.divide(Expression.TWO.squareRoot());
  }
  if (d === 60) {
    return Expression.ONE.divide(Expression.TWO);
  }
  if (d === 90) {
    return Expression.ZERO;
  }

  if (d === 15) {
    return f(d);
  }
  if (d === 75) {
    return f(d);
  }

  function phi() {
    return Expression.ONE.add(Expression.TWO.add(Expression.TWO).add(Expression.ONE).squareRoot()).divide(Expression.TWO);
  }

  if (d === 36) {
    return phi().divide(Expression.TWO);
  }
  if (d === 72) {
    return phi().subtract(Expression.ONE).divide(Expression.TWO);
  }

  if (d === 18) {
    return Expression.TWO.add(phi()).squareRoot().divide(Expression.TWO);
  }
  // http://www.maths.surrey.ac.uk/hosted-sites/R.Knott/Fibonacci/simpleTrig.html#section4.2
  if (d === 54) {
    //return Expression.TWO.subtract(Expression.TWO.subtract(phi()).squareRoot()).squareRoot().divide(Expression.TWO);
    // https://www.cut-the-knot.org/pythagoras/cos36.shtml
    return Expression.TWO.add(Expression.ONE).subtract(phi()).squareRoot().divide(Expression.TWO);
  }

  if (d === 7.5) {
    return cosapb(30, -22.5);
  }
  if (d === 22.5) {
    return f(d);
  }
  if (d === 37.5) {
    return cosapb(60, -22.5);
  }
  if (d === 52.5) {
    return cosapb(30, +22.5);
  }
  if (d === 67.5) {
    return f(d);
  }
  if (d === 82.5) {
    return cosapb(60, +22.5);
  }

  // https://en.wikipedia.org/wiki/Sine

/*
  for (var d = 3; d <= 90; d += 3) {
    if (d % 15 !== 0 && d % 18 !== 0) {
      for (var a = 75; a > 0; a -= 15) {
        var b = -(a - d);
        if (b % 18 === 0) {
          console.log(`cos(${d})=cos(${a}+${b})`);
        }
      }
    }
  }
*/

  //if (d === 3) {
    // https://en.wikipedia.org/wiki/Trigonometric_constants_expressed_in_real_radicals#3°:_regular_hexacontagon_(60-sided_polygon)
    //return cosapb(18, -15);
  //}

  if (d % 3 === 0) {
    var a = 90 - (Math.floor(d / 3) % 6) * 15;
    var b = -(a - d);
    return cosapb(a, b);
  }

  //TODO: sin(1.5)
  return undefined;
};

var simplifyConstantValue = function (x, type) {
  var a = undefined;
  var b = undefined;
  if (x instanceof Integer && x.compareTo(Expression.ZERO) === 0) {
    a = Expression.ZERO;
    b = Expression.ONE;
  } else if (x === Expression.PI) {
    a = Expression.ONE;
    b = Expression.ONE;
  } else if (x instanceof Multiplication && x.a instanceof Integer && x.b === Expression.PI) {
    a = x.a;
    b = Expression.ONE;
  } else if (x instanceof Division && x.b instanceof Integer && x.a === Expression.PI) {
    a = Expression.ONE;
    b = x.b;
  } else if (x instanceof Division && x.b instanceof Integer && x.a instanceof Multiplication && x.a.a instanceof Integer && x.a.b === Expression.PI) {
    a = x.a.a;
    b = x.b;
  } else if (x instanceof Expression.Degrees) {
    var t = x.value.simplify();
    return simplifyConstantValue(t.multiply(Expression.PI).divide(Integer.fromNumber(180)), type);
  }
  if (a != undefined && b != undefined) {
    b = b.toNumber();
    var k = 2;
    if (b >= 1 && b <= 180 && (180 * k) % b === 0) {
      var d = a.multiply(Integer.fromNumber(Math.floor((180 * k) / b))).remainder(Integer.fromNumber(360 * k)).toNumber();
      d /= k;
      if (type === "sin") {
        d = 90 - d;
        if (d >= 360 - 90) {
          d -= 360;
        }
      } else if (type !== "cos") {
        throw new TypeError();
      }
      return simplifyConstantValueInternal(d);
    }
  }
  if (x instanceof Expression.Radians && x.value.equals(Expression.ZERO)) {
    return simplifyConstantValue(x.value, type);
  }
  var c = Expression.getConstant(x);
  if (c instanceof Expression.Complex && c.real.equals(Expression.ZERO)) {
    if (type === "sin") {
      return Expression.I.multiply(x.divide(Expression.I).sinh());
    }
    if (type === "cos") {
      return x.divide(Expression.I).cosh();
    }
  }
  return undefined;
};

Expression.prototype.sinh = function () {
  var a = this;
  return a.exp().subtract(a.negate().exp()).divide(Expression.TWO);
};

Expression.prototype.cosh = function () {
  var a = this;
  return a.exp().add(a.negate().exp()).divide(Expression.TWO);
};

var isArgumentValid = function (x, type) {
  if (x instanceof Expression.Radians) {
    var isAlgebraicInteger = function (x) {
      return x instanceof Expression.Integer ||
             x instanceof Expression.NthRoot && typeof x.n === "number" && x.n % 1 === 0 && isAlgebraicInteger(x.a) ||
             x instanceof Expression.Division && isAlgebraicInteger(x.getNumerator()) && isAlgebraicInteger(x.getDenominator()) ||
             x instanceof Expression.Addition && isAlgebraicInteger(x.a) && isAlgebraicInteger(x.b) ||
             x instanceof Expression.Multiplication && isAlgebraicInteger(x.a) && isAlgebraicInteger(x.b);
    };
    // https://ru.wikipedia.org/wiki/Трансцендентное_число#Примеры_трансцендентных_чисел
    return isAlgebraicInteger(x.value);
  }
  if (x instanceof Expression.Degrees) {
    return simplifyConstantValue(x, type) != undefined;
  }
  if (x instanceof Expression.Symbol) {
    return Expression.isScalar(x);
  }
  if (x instanceof Addition) {
    return isArgumentValid(x.a, type) && isArgumentValid(x.b, type);
  }
  if (x instanceof Multiplication) {
    if (x.a instanceof Integer && Expression.isScalar(x.b) && x.b instanceof Expression.Symbol) {
      return true;
    }
    if (Expression.isScalar(x.b) && x.b instanceof Expression.Symbol) {
      if (x.a instanceof Expression.NthRoot && x.a.a instanceof Integer) {
        return true;
      }
      if (x.a instanceof Multiplication && x.a.a instanceof Integer && x.a.b instanceof Expression.NthRoot && x.a.b.a instanceof Integer) {
        return true;
      }
    }
  }
  if (x instanceof Division) {
    if (x.b instanceof Integer) {
      return isArgumentValid(x.a);
    }
  }
  return false;
};

Expression.prototype.sin = function () {
  var x = this;
  var t = simplifyConstantValue(x, "sin");
  if (t != undefined) {
    return t;
  }
  if (x.isNegative()) {
    return x.negate().sin().negate();
  }
  if (!isArgumentValid(x, "sin")) {
    throw new RangeError("NotSupportedError");
  }
  return new Sin(x);
};

function Cos(x) {
  Expression.Function.call(this, "cos", x);
}
Cos.prototype = Object.create(Expression.Function.prototype);

Expression.prototype.cos = function () {
  var x = this;
  var t = simplifyConstantValue(x, "cos");
  if (t != undefined) {
    return t;
  }
  if (x.isNegative()) {
    return x.negate().cos();
  }
  if (!isArgumentValid(x, "cos")) {
    throw new RangeError("NotSupportedError");
  }
  return new Cos(x);
};

Expression.simplifications.push(Expression.simplifyTrigonometry);

Expression.Sin = Sin;
Expression.Cos = Cos;

//Expression.Negation.prototype.compare4Multiplication = function (y) {
//TODO: fix, more tests
//  return new Expression.Multiplication(Expression.ONE.negate(), this.a).compare4Multiplication(y);
//};

Expression.Addition.prototype.compare4Addition = function (y) {
  // cos(a + b) + cos(a + b)
  var x = this;
  return Expression.Addition.compare4Addition(x, y);
};

/*
Expression.Multiplication.prototype.compare4MultiplicationInteger = function (x) {
  return -1;
};
Expression.MatrixSymbol.prototype.compare4MultiplicationExponentiation = function (x) {
  return -1;//?
};
*/

//!!!
Expression.Addition.prototype.compare4Multiplication = function (y) {
  if (y instanceof Integer) {
    return -1;
  }
  if (y instanceof Expression.MatrixSymbol) {
    return +1;
  }
  //TODO: fix

  var x = this;
  var i = x.summands();
  var j = y.summands();
  var a = i.next().value;
  var b = j.next().value;
  while (a != null && b != null) {
    var c = a.compare4Multiplication(b);
    if (c !== 0) {
      return c;
    }
    a = i.next().value;
    b = j.next().value;
  }
  return a != null ? +1 : (b != null ? -1 : 0);
};

Expression.Addition.prototype.compare4MultiplicationSymbol = function (x) {
  return 0 - this.compare4Multiplication(x);
};

Expression.Addition.prototype.compare4MultiplicationInteger = function (x) {
  return +1;
};

Expression.Addition.compare4Addition = function (x, y) {
  var i = x.summands();
  var j = y.summands();
  var a = i.next().value;
  var b = j.next().value;
  while (a != null && b != null) {
    var c = a.compare4Addition(b);
    if (c !== 0) {
      return c;
    }
    a = i.next().value;
    b = j.next().value;
  }
  return a != null ? +1 : (b != null ? -1 : 0);
};

//!!!


//!new 2017-04-26
Expression.Degrees = function (value) {
  this.value = value;
};
Expression.Degrees.prototype = Object.create(Expression.prototype);
Expression.Degrees.prototype.toString = function (options) {
  return this.value.toString(options) + "\u00B0";
};
Expression.Degrees.prototype.equals = function (y) {
  return y instanceof Expression.Degrees && this.value.equals(y.value);
};
Expression.Degrees.prototype.compare4AdditionSymbol = function (x) {
  return -1;
};


//!new 2019-12-27
Expression.Radians = function (value) {
  this.value = value;
};
Expression.Radians.prototype = Object.create(Expression.prototype);
Expression.Radians.prototype.toString = function (options) {
  var b = this.value;
  var fb = b instanceof Expression.Integer ? false : true;
  return (fb ? "(" : "") + b.toString(options) + (fb ? ")" : "") + " rad";
};
Expression.Radians.prototype.equals = function (y) {
  return y instanceof Expression.Radians && this.value.equals(y.value);
};
Expression.Radians.prototype.compare4AdditionSymbol = function (x) {
  return x.compare4Addition(this.value);
};
Expression.Radians.prototype.compare4Addition = function (y) {
  return this.value.compare4Addition(y instanceof Expression.Radians ? y.value : y);
};
Expression.Radians.prototype.compare4Multiplication = function (y) {
  return this.value.compare4Multiplication(y instanceof Expression.Radians ? y.value : y);
};
Expression.Radians.prototype.compare4MultiplicationSymbol = function (x) {
  return x.compare4Multiplication(this.value);
};
Expression.Radians.prototype.compare4MultiplicationInteger = function (x) {
  return x.compare4Multiplication(this.value);
};

Expression.Radians.prototype.negate = function () {
  return new Expression.Radians(this.value.negate());
};

Expression.Radians.prototype.add = function (y) {
  if (y instanceof Expression.Radians) {
    var x = this.value.add(y.value);
    if (x.equals(Expression.ZERO)) {
      return x;
    }
    return new Expression.Radians(x);//!?
  }
  return Expression.prototype.add.call(this, y);
};


//!new 2019-11-23
// https://en.wikipedia.org/wiki/Trigonometric_functions_of_matrices#cite_note-3
Expression.Matrix.prototype.sin = function () {
  var X = this;
  var i = Expression.I;
  var TWO = Expression.TWO;
  return i.multiply(X).exp().subtract(i.negate().multiply(X).exp()).divide(TWO.multiply(i));
};
Expression.Matrix.prototype.cos = function () {
  var X = this;
  var i = Expression.I;
  var TWO = Expression.TWO;
  return i.multiply(X).exp().add(i.negate().multiply(X).exp()).divide(TWO);
};

