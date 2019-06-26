import Expression from './Expression.js';
import Polynomial from './Polynomial.js';

function Condition(array) {
  this.array = array;
}

Condition.NEZ = " != 0";
Condition.EQZ = " == 0";

Condition.prototype._and = function (operator, e) {
  if (operator !== Condition.NEZ && operator !== Condition.EQZ) {
    throw new Error();
  }
  if (e == undefined) {
    throw new RangeError();
  }
  if (this === Condition.FALSE) {
    return this;
  }
  if (e instanceof Expression.GF2Value) {
    return this._and(operator, e.value);//?
  }
  if (e instanceof Expression.NthRoot) {
    return this._and(operator, e.a);
  }
  if (e instanceof Expression.Integer || e instanceof Expression.Complex) {
    if (e.equals(Expression.ZERO)) {
      if (operator === Condition.NEZ) {
        return Condition.FALSE;
      }
      return this;
    }
    if (operator === Condition.EQZ) {
      return Condition.FALSE;
    }
    return this;
  }
  if (operator === Condition.NEZ) {
    if (e instanceof Expression.Multiplication) {
      return this._and(Condition.NEZ, e.a)._and(Condition.NEZ, e.b);
    }
  }
  if (e instanceof Expression.Division) {
    return this._and(operator, e.a)._and(Condition.NEZ, e.b);
  }
  if (Expression.isConstant(e)) {
    if (operator === Condition.NEZ) {
      return this;
    }
    if (operator === Condition.EQZ) {
      return Condition.FALSE;
    }
  }
  var contains = function (array, operator, e) {
    for (var i = 0; i < array.length; i += 1) {
      if (array[i].operator === operator && array[i].expression.equals(e)) {
        return true;
      }
    }
    return false;
  };
  if (contains(this.array, operator, e)) {
    return this;
  }
  if (contains(this.array, operator === Condition.EQZ ? Condition.NEZ : Condition.EQZ, e)) {
    return Condition.FALSE;
  }
  if (e instanceof Expression.Exponentiation && e.b instanceof Expression.Integer && e.b.compareTo(Expression.ZERO) > 0 && e.a instanceof Expression.Symbol) {
    return this._and(operator, e.a);
  }

  //!new
  if (e.isNegative()) {
    e = e.negate();
  }
  if (e instanceof Expression.Addition && e.a instanceof Expression.Exponentiation && e.a.b.inverse().equals(Expression.TWO) && e.b instanceof Expression.Integer) {
    if (e.b.compareTo(Expression.ZERO) > 0) {
      return this._and(operator, Expression.ONE);
    }
    //?
    //TODO: fix
    return this._and(operator, e.a.a.subtract(e.b.pow(Expression.TWO)));
  }

  var add = function (oldArray, y) {
    if (y.expression instanceof Expression.Multiplication && y.expression.b instanceof Expression.IdentityMatrix) {
      return add(oldArray, {expression: y.expression.a, operator: y.operator});
    }
    if (y.expression instanceof Expression.Division) {
      var tmp = oldArray;
      tmp = add(tmp, {expression: y.expression.a, operator: y.operator});
      if (tmp == null) {
        return null;
      }
      tmp = add(tmp, {expression: y.expression.b, operator: Condition.NEZ});
      return tmp;
    }
    /*if (y.expression instanceof Expression.Division && Expression.isConstant(y.expression.b)) {
      y = {
        expression: y.expression.a,
        operator: y.operator
      };
    }*/
    if (y.expression instanceof Expression.Integer || y.expression instanceof Expression.Complex) {
      if (y.operator === Condition.NEZ && y.expression.equals(Expression.ZERO) ||
          y.operator === Condition.EQZ && !y.expression.equals(Expression.ZERO)) {
        return null;
      }
      return oldArray;
    }
    //TODO: check code coverage, remove extra branches
    if (Expression.isConstant(y.expression) && !y.expression.equals(Expression.ZERO)) {
      if (y.operator === Condition.NEZ) {
        return oldArray;
      }
      if (y.operator === Condition.EQZ) {
        return null;
      }
    }
    if (y.expression instanceof Expression.NthRoot) {
      return add(oldArray, {expression: y.expression.a, operator: y.operator});
    }
    if (y.expression instanceof Expression.Multiplication) {
      if (y.operator === Condition.EQZ) {
        if (y.expression.a instanceof Expression.Integer && !y.expression.a.equals(Expression.ZERO)) {
          //TODO: fix - ?
          y = {
            expression: y.expression.b,
            operator: y.operator
          };
        }
      }
      if (y.operator === Condition.NEZ) {
        var tmp = oldArray;
        tmp = add(tmp, {expression: y.expression.a, operator: Condition.NEZ});
        if (tmp == null) {
          return null;
        }
        tmp = add(tmp, {expression: y.expression.b, operator: Condition.NEZ});
        return tmp;
      }
    }

    var p = Expression.getMultivariatePolynomial(y.expression);
    if (p != null) {
      var content = p.p.getContent();
      if (!content.equals(Expression.ONE) && !content.equals(Expression.ONE.negate())) {
        y = {
          expression: y.expression.divide(content),
          operator: y.operator
        };
      }
    }

    var newArray = [];
    for (var i = 0; i < oldArray.length; i += 1) {
      var x = oldArray[i];
      if (x.operator === Condition.NEZ && y.operator === Condition.EQZ ||
          x.operator === Condition.EQZ && y.operator === Condition.NEZ) {
        var g = x.expression.gcd(y.expression);
        while (!g.equals(Expression.ONE) && !g.equals(Expression.ONE.negate())) {
          if (x.operator === Condition.EQZ) {
            x = {
              operator: x.operator,
              expression: x.expression.divide(g)
            };
          } else {
            y = {
              operator: y.operator,
              expression: y.expression.divide(g)
            };
          }
          g = x.expression.gcd(y.expression);
        }
        if (x.operator === Condition.EQZ) {
          var tmp = y;
          y = x;
          x = tmp;
        }
        if (Expression.isConstant(y.expression)) {
          return null;
        }
        //if (!Expression.isSingleVariablePolynomial(x.expression.multiply(y.expression))) {
        //  newArray.push(x);
        //}
      }
      if (x.operator === Condition.NEZ && y.operator === Condition.EQZ && Expression.isSingleVariablePolynomial(x.expression.multiply(y.expression))) {
        y = y;
      } else if (x.operator === Condition.EQZ && y.operator === Condition.EQZ && Expression.isSingleVariablePolynomial(x.expression.multiply(y.expression))) {
        var g = x.expression.gcd(y.expression);
        if (g instanceof Expression.Integer) {
          return null;
        }
        y = {
          operator: y.operator,
          expression: g
        };
      } else if (Expression.isSingleVariablePolynomial(x.expression.multiply(y.expression))) {
        // TODO: both Condition.NEZ - ?
        var g = x.expression.gcd(y.expression);
        x = {
          operator: x.operator,
          expression: x.expression.divide(g)
        };
        if (!Expression.isConstant(x.expression)) {
          newArray.push(x);
        }
      } else { // !isSingleVariablePolynomial
        var p = null;
        var pOperator = null;
        var pp = null;
        var other = null;
        var px = Expression.getMultivariatePolynomial(x.expression);
        var py = Expression.getMultivariatePolynomial(y.expression);
        //var xy = Expression.getMultivariatePolynomial(x.expression.multiply(y.expression));

        //console.assert(px != null && py != null);

        if (//xy != null &&
            //x.operator === Condition.EQZ &&
            //y.operator === Condition.EQZ &&
            px != null && py != null) {

          //if (px != null && px.p.getDegree() !== 1 && py == null) {
            //py = {p: Polynomial.toPolynomial(y.expression, px.v), v: px.v};
            //if (xy.v === py.v) {
            //  py = null;
            //}
          //}
          //if (px == null && py != null && py.p.getDegree() !== 1) {
            //px = {p: Polynomial.toPolynomial(x.expression, py.v), v: py.v};
            //if (xy.v === px.v) {
            //  px = null;
            //}
          //}
          //if (px == null && py == null) {//?TODO:
          //  px = {p: Polynomial.toPolynomial(x.expression, xy.v), v: xy.v};
          //  py = {p: Polynomial.toPolynomial(y.expression, xy.v), v: xy.v};
          //}

          if (x.operator === Condition.EQZ && px != null && px.p.getDegree() === 1 && y.operator === Condition.EQZ && py != null && py.p.getDegree() === 1) {
            if (px.v.symbol < py.v.symbol) {//!
              px = null;
            }
          }

          if (y.operator === Condition.EQZ && py != null && py.p.getDegree() === 1) {
            pp = py;
            p = x.expression;
            pOperator = x.operator;
            other = y;
          }
          if (x.operator === Condition.EQZ && px != null && px.p.getDegree() === 1) {
            pp = px;
            p = y.expression;
            pOperator = y.operator;
            other = x;
          }
        }
        if (pp != null) {
          var polynomial = Polynomial.toPolynomial(p, pp.v);
          var a = polynomial.calcAt(pp.p.getCoefficient(0).negate().divide(pp.p.getCoefficient(1)));
          if (!a.equals(p) && (a.getDenominator() instanceof Expression.Integer)) {
            var tmp = {
              operator: pOperator,
              expression: a
            };
            newArray = add(newArray, tmp);
            if (newArray == null) {
              return null;
            }
            if (true) {
              newArray = add(newArray, other);
              if (newArray == null) {
                return null;
              }
              for (var j = i + 1; j < oldArray.length; j += 1) {
                newArray = add(newArray, oldArray[j]);
                if (newArray == null) {
                  return null;
                }
              }
              return newArray;
            } else {
              y = other;
            }
          } else {
            newArray.push(x);
          }
        } else {
          newArray.push(x);
        }
      }
    }
    newArray.push(y);
    return newArray;
  };
  var newArray = add(this.array, {
    operator: operator,
    expression: e
  });
  if (newArray == null) {
    return Condition.FALSE;
  }
  
  return new Condition(newArray);
};

Condition.prototype.andNotZero = function (e) {
  return this._and(Condition.NEZ, e);
};
Condition.prototype.andZero = function (e) {
  return this._and(Condition.EQZ, e);
};
Condition.prototype.isFalse = function () {
  return this === Condition.FALSE;
};
Condition.prototype.isTrue = function () {
  return this === Condition.TRUE;
};
Condition.prototype.toString = function (options) {
  if (this === Condition.FALSE || this === Condition.TRUE || this.array.length === 0) {
    //throw new RangeError();
    return "";
  }
  var s = '';
  for (var i = 0; i < this.array.length; i += 1) {
    s += (i !== 0 ? ', ' : '') + this.array[i].expression.toString(options) + (this.array[i].operator === Condition.NEZ ? ' != 0' : '') + (this.array[i].operator === Condition.EQZ ? ' == 0' : '');
  }
  return s;
};

Condition.TRUE = new Condition(new Array(0));
Condition.FALSE = new Condition(undefined);

export default Condition;
