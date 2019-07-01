
import toDecimalStringInternal from './toDecimalString.js';
import Polynomial from './Polynomial.js';//TODO:
import Expression from './Expression.js';//TODO:
import NonSimplifiedExpression from './NonSimplifiedExpression.js';//TODO:

//TODO: remove - ?
import Condition from './Condition.js';


//!TODO: remove
Polynomial.prototype.toString = function (options) {
  options = options || {};
  return this.toExpression(options.polynomialVariable || (new Expression.Symbol("x"))).toString(options);
};
Polynomial.prototype.toMathML = function (options) {
  options = options || {};
  return this.toExpression(options.polynomialVariable || (new Expression.Symbol("x"))).toMathML(options);
};
Expression.Polynomial.prototype.toString = function (options) {
  return this.polynomial.toString(options);
};
Expression.Polynomial.prototype.toMathML = function (options) {
  return this.polynomial.toMathML(options);
};

// coefficient - Expression
// variable - Expression
var printPartOfAddition = function (isLast, isFirst, coefficient, variable, options) {
  if (coefficient.equals(Expression.ZERO)) {
    return (isLast && isFirst ? "<mn>0</mn>" : "");
  }
  var isNegative = false;
  if (coefficient.isNegative() && !isFirst) {
    isNegative = true;
    coefficient = coefficient.negateCarefully();//?
  }
  var precedenceOfMultiptication = new Expression.Multiplication(Expression.ZERO, Expression.ZERO).getPrecedence();
  var areBracketsRequired = coefficient.getPrecedence() < precedenceOfMultiptication; //?
  var c = coefficient.equals(Expression.ONE);
  //TODO: fix
  return (isFirst ? '' : '') +
         (!isFirst && isNegative ? '<mo form="infix">&minus;</mo>' : '') +
         (!isFirst && !isNegative ? '<mo form="infix">+</mo>' : '') +
         (c ? '' : '<mrow>') +
         (c || !areBracketsRequired ? '' : '<mrow><mo>(</mo>') +
         (c ? '' : coefficient.toMathML(options)) +
         (c || !areBracketsRequired ? '' : '<mo>)</mo></mrow>') +
         (c ? '' : '<mo>&times;</mo>') +
         variable.toMathML(options) +
         (c ? '' : '</mrow>');
};


var decimalToMathML = function (sign, number) {
  return (sign < 0 ? "<mrow>" : "") + (sign < 0 ? "<mo>&minus;</mo>" : "") + "<mn>" + number + "</mn>" + (sign < 0 ? "</mrow>" : "");
};
var complexToMathML = function (real, imaginarySign, imaginaryAbs) {
  return '<mrow>' + real + (imaginarySign > 0 ? '<mo>+</mo>' : '<mo>&minus;</mo>') + (imaginaryAbs !== '' ? imaginaryAbs + '<mo>&it;</mo>' : '') + '<mi>&ii;</mi>' + '</mrow>';
};

//TODO: move
Expression.toDecimalString = function (x, options) {
  
  function isConstant(e) {
    if (e instanceof Expression.Symbol) {
      return false;
    }
    if (e instanceof Expression.NonSimplifiedExpression) {
      return false;
    }
    if (e instanceof Expression.Matrix) {
      return false;
    }
    if (e instanceof Expression.Polynomial) {
      return false;
    }
    if (e instanceof Expression.BinaryOperation) {
      return (e.a === Expression.E || isConstant(e.a)) && isConstant(e.b);
    }
    if (e instanceof Expression.Negation) {
      return isConstant(e.b);
    }
    if (e instanceof Expression.Function) {
      return isConstant(e.a);
    }
    return true;
  }

  var fractionDigits = options != null ? options.fractionDigits : -1;
  if (fractionDigits >= 0 && isConstant(x)) {
    return toDecimalStringInternal(x, fractionDigits, decimalToMathML, complexToMathML);
  }
  return undefined;
};

Expression.idCounter = 0;
Expression.id = function () {
  return (Expression.idCounter += 1).toString();
};

//TODO: ?
Expression.escapeHTML = function (s) {
  return s.replace(/&/g, "&amp;")
          .replace(/"/g, "&quot;")
          .replace(/</g, "&lt;")
          .replace(/>/g, "&gt;");
};

Expression.Matrix.prototype.toMathML = function (options) {
  var x = this.matrix;
  options = Expression.setTopLevel(true, options);

  var useMatrixContainer = options.useMatrixContainer == undefined ? true : options.useMatrixContainer;
  //TODO: fix!
  var braces = options.useBraces == undefined ? undefined : options.useBraces;
  var columnlines = options.columnlines == undefined ? 0 : options.columnlines;
  var variableNames = options.variableNames == undefined ? undefined : options.variableNames;

  var verticalStrike = options.verticalStrike == undefined ? -1 : options.verticalStrike;
  var horizontalStrike = options.horizontalStrike == undefined ? -1 : options.horizontalStrike;

  var cellIdGenerator = options.cellIdGenerator == undefined ? undefined : options.cellIdGenerator;
  var pivotCell = options.pivotCell == undefined ? undefined : options.pivotCell;

  var isLUDecomposition2 = options.isLUDecomposition2 == undefined ? undefined : options.isLUDecomposition2;
  var highlightRow = options.highlightRow == undefined ? -1 : options.highlightRow;
  var highlightCol = options.highlightCol == undefined ? -1 : options.highlightCol;

  options = Object.assign({}, options, {
    useBraces: undefined,
    columnlines: undefined,
    variableNames: undefined,
    verticalStrike: undefined,
    horizontalStrike: undefined,
    cellIdGenerator: undefined,
    pivotCell: undefined,
    isLUDecomposition2: undefined,
    highlightRow: undefined,
    highlightCol: undefined
  });

  var result = "";
  var rows = x.rows();
  var cols = x.cols();
  var i = -1;

  //TODO: remove `Expression.id()`
  var containerId = options.idPrefix + "-" + Expression.id();
  if (useMatrixContainer) {
    result += "<munder>";
    result += "<menclose notation=\"none\" href=\"#\" id=\"" + containerId + "\" data-matrix=\"" + Expression.escapeHTML(x.toString()) + "\" draggable=\"true\" tabindex=\"0\" contextmenu=\"matrix-menu\">";
  }

  result += braces == undefined ? '<mrow><mo>(</mo>' : '<mrow>' + (braces[0] === ' ' ? '' : '<mo>' + braces[0] + '</mo>');
  var columnlinesAttribute = "";
  if (columnlines !== 0 && cols - 1 > 0) {
    var k = -1;
    while (++k < cols - 1) {
      columnlinesAttribute += (cols - 1 + columnlines === k ? "solid " : "none ");
    }
    // whitespace
    columnlinesAttribute = columnlinesAttribute.slice(0, -1);
  }
  //! 2017-07-06 rowspacing="0ex" was added to make it look better with Native MathML (when it is supported) and to have the same style as in mathml.css
  var useColumnspacing = verticalStrike !== -1 || horizontalStrike !== -1 || pivotCell != undefined || cellIdGenerator != undefined;
  result += "<mtable" + (useColumnspacing ? " rowspacing=\"0ex\"" + " columnspacing=\"0em\"" : " rowspacing=\"0ex\"") + (variableNames != undefined ? " columnalign=\"right\"" : "") + (columnlinesAttribute !== "" ? " columnlines=\"" + columnlinesAttribute + "\"" : "") + ">";
  while (++i < rows) {
    var j = -1;
    if (variableNames != undefined) {// TODO: fix?
      //TODO: use code from polynomialToExpression3 (shared)
      var row = "";
      var wasNotZero = false;
      while (++j < cols - 1) {
        // TODO: fix `new Expression.Symbol()`
        row += "<mtd>";
        row += printPartOfAddition(j === cols - 2, !wasNotZero, x.e(i, j), new Expression.Symbol(variableNames[j]), options);
        row += "</mtd>";
        wasNotZero = wasNotZero || !x.e(i, j).equals(Expression.ZERO);
      }
      row += "<mtd><mo>=</mo></mtd><mtd>" + x.e(i, cols - 1).toMathML(options) + "</mtd>";
      if (wasNotZero || !x.e(i, cols - 1).equals(Expression.ZERO)) {
        result += "<mtr>";
        result += row;
        result += "</mtr>";
      }
    } else {
      result += "<mtr>";
      while (++j < cols) {
        result += "<mtd" + (cellIdGenerator != undefined ? " id=\"" + cellIdGenerator(i, j) + "\"" : "") + ">";
        if (pivotCell != undefined && i === pivotCell.i && j === pivotCell.j) {
          result += "<mstyle mathvariant=\"bold\">";
          result += "<menclose notation=\"circle\">";
        }
        if (verticalStrike === j || horizontalStrike === i) {
          var notation = ((verticalStrike === j ? " " + "verticalstrike" : "") +
                          (horizontalStrike === i ? " " + "horizontalstrike" : "")).slice(1);
          result += "<menclose notation=\"" + notation + "\">";
        }
        if (useColumnspacing) {
          result += "<mpadded width=\"+0.8em\" lspace=\"+0.4em\">";
        }
        var highlight = j < i && isLUDecomposition2 ||
                        highlightRow === i && (columnlines === 0 || j <= cols - 1 + columnlines) || highlightCol === j;
        if (highlight) {
          result += "<mrow mathbackground=\"#80FF80\">";
        }
        result += x.e(i, j).toMathML(options);
        if (highlight) {
          result += "</mrow>";
        }
        if (useColumnspacing) {
          result += "</mpadded>";
        }
        if (verticalStrike === j || horizontalStrike === i) {
          result += "</menclose>";
        }
        if (pivotCell != undefined && i === pivotCell.i && j === pivotCell.j) {
          result += "</menclose>";
          result += "</mstyle>";
        }
        result += "</mtd>";
      }
      result += "</mtr>";
    }
  }
  result += "</mtable>";
  result += braces == undefined ? '<mo>)</mo></mrow>' : (braces[1] === ' ' ? '' : '<mo>' + braces[1] + '</mo>') + '</mrow>';

  if (useMatrixContainer) {
    result += "</menclose>";

    result += "<mtext>";
    result += "<button type=\"button\" class=\"matrix-menu-show matrix-menu-show-new\" data-for-matrix=\"" + containerId + "\" aria-haspopup=\"true\">&#x2630;</button>";
    result += "</mtext>";

    result += "</munder>";
  }

  return result;
};

Expression.Determinant.prototype.toMathML = function (options) {
  var x = this;
  if (x.a instanceof Expression.Matrix || (x.a instanceof NonSimplifiedExpression && x.a.e instanceof Expression.Matrix)) {
    options = Object.assign({}, options, {
      useBraces: ["|", "|"]
    });
    //TODO: fix
    return x.a.toMathML(options);
  }
  return "<mrow><mo>|</mo>" + x.a.toMathML(options) + "<mo>|</mo></mrow>";
};
Expression.Transpose.prototype.toMathML = function (options) {
  var x = this;
  //TODO: ^T ?
  // https://www.w3.org/TR/MathML3/chapter4.html#contm.transpose
  return "<msup>" +
         x.a.toMathML(options) +
         "<mi>T</mi>" +
         "</msup>";
};
Expression.SquareRoot.prototype.toMathML = function (options) {
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
  return "<msqrt>" +
         this.a.toMathML(Expression.setTopLevel(true, options)) +
         "</msqrt>";
};
Expression.CubeRoot.prototype.toMathML = function (options) {
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
  return "<mroot>" +
         this.a.toMathML(Expression.setTopLevel(true, options)) +
         "<mi>" + 3 + "</mi>" +
         "</mroot>";
};
Expression.NthRoot.prototype.toMathML = function (options) {
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
  return "<mroot>" +
         this.a.toMathML(Expression.setTopLevel(true, options)) +
         "<mi>" + this.n + "</mi>" +
         "</mroot>";
};
Expression.Function.prototype.toMathML = function (options) {
  var x = this;
  var fa = !(x.a instanceof Expression.Matrix) && !(x.a instanceof NonSimplifiedExpression && x.a.e instanceof Expression.Matrix);//?
  //TODO: fix
  return "<mrow>" +
         (typeof i18n !== "undefined" ? "<mi>" + (x.name === "rank" ? i18n.rankDenotation : (x.name === "sin" ? i18n.sinDenotation : (x.name === "tan" ? i18n.tanDenotation : x.name))) + "</mi>" : "<mi>" + x.name + "</mi>") +
         "<mo>&af;</mo>" +
         (fa ? "<mrow><mo>(</mo>" : "") +
         x.a.toMathML(Expression.setTopLevel(true, options)) +
         (fa ? "<mo>)</mo></mrow>" : "") +
         "</mrow>";
};
Expression.Division.prototype.toMathML = function (options) {
  if (options != null && options.nofractions) {
    return Expression.BinaryOperation.prototype.toMathML.call(this, options);
  }
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
  var x = this;
  var denominator = x.getDenominator();
  var numerator = x.getNumerator();
  //???
  //if (numerator.isNegative()) {
  //  return "<mrow><mo>&minus;</mo>" + x.negateCarefully().toMathML(options) + "</mrow>";
  //}
  return "<mfrac>" +
         numerator.toMathML(Expression.setTopLevel(true, options)) + 
         denominator.toMathML(Expression.setTopLevel(true, options)) +
         "</mfrac>";
};

Expression.Integer.prototype.toMathML = function (options) {
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
  var x = this;
  var tmp = x.toString();
  return tmp.slice(0, 1) === "-" ? "<mrow>" + "<mo>&minus;</mo>" + "<mn>" + tmp.slice(1) + "</mn>" + "</mrow>" : "<mn>" + tmp + "</mn>";
};
Expression.BinaryOperation.prototype.toMathML = function (options) {
  if (options != null &&
      options.fractionDigits >= 0 &&
      this.unwrap() instanceof Expression.Exponentiation &&
      this.unwrap().a.unwrap() instanceof Expression.Symbol &&
      this.unwrap().a.unwrap() !== Expression.E &&
      (this.unwrap().b.unwrap() instanceof Expression.Integer || this.unwrap().b.unwrap() instanceof Expression.Negation && this.unwrap().b.unwrap().b.unwrap() instanceof Expression.Integer)) {
    options = Object.assign({}, options, {fractionDigits: -1});
  }
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }

  //!2019-05-16
  if (this instanceof Expression.Addition) {
    var s = [];
    var b = null;
    for (var additions = this.summands(), x = additions.next().value; x != null; x = additions.next().value) {
      if (b != null) {
        var n = false;
        if (b.isNegative()) {
          n = true;
          b = b.negateCarefully();
        }
        var fence = this.getPrecedence() >= b.getPrecedence();
        fence = fence || b.isUnaryPlusMinus();
        s.push((fence ? '<mrow><mo>(</mo>' : '') + b.toMathML(Expression.setTopLevel(fence, options)) + (fence ? '<mo>)</mo></mrow>' : ''));
        s.push(n ? '<mo>&minus;</mo>' : '<mo>+</mo>');
      }
      b = x;
    }
    s = s.reverse().join('');
    var a = b;
    var fence = a.getPrecedence() + (a.isRightToLeftAssociative() ? -1 : 0) < this.getPrecedence();
    if (options != undefined && options.isTopLevel != undefined && options.isTopLevel === false) {
      fence = fence || a.isUnaryPlusMinus();
    }
    s = (fence ? "<mrow><mo>(</mo>" : "") + a.toMathML(Expression.setTopLevel(fence || options == undefined || options.isTopLevel, options)) + (fence ? "<mo>)</mo></mrow>" : "") + s;
    return '<mrow>' + s + '</mrow>';
  }

  //!
  var a = this.a;
  var b = this.b;
  var isSubtraction = false;
  // TODO: check
  if (this instanceof Expression.Addition && b.isNegative()) {
    isSubtraction = true;
    b = b.negateCarefully();//?
  }

  var fa = a.getPrecedence() + (a.isRightToLeftAssociative() ? -1 : 0) < this.getPrecedence();
  var fb = this.getPrecedence() + (this.isRightToLeftAssociative() ? -1 : 0) >= b.getPrecedence();
  if (options != undefined && options.isTopLevel != undefined && options.isTopLevel === false) {
    fa = fa || a.isUnaryPlusMinus();
  }
  fb = fb || b.isUnaryPlusMinus();
  fb = fb || (this.unwrap() instanceof Expression.Exponentiation && b.unwrap() instanceof Expression.Exponentiation);// 2^3^4
  fa = fa || (this.unwrap() instanceof Expression.Exponentiation && a.unwrap() instanceof Expression.Function); // cos(x)^(2+3)
  var s = isSubtraction ? "-" : this.getS();

  if (this instanceof Expression.Exponentiation) {
    if (a.unwrap() === Expression.E && b.unwrap() instanceof Expression.Matrix) {
      return '<mrow><mi>exp</mi><mo>&af;</mo>' + b.toMathML(options) + '</mrow>';
    }
    var boptions = options;
    if (!(a.unwrap() instanceof Expression.Matrix)) {
      boptions = Object.assign({}, options || {}, {nofractions: true});
    }

      return "<msup>" + 
             (fa ? "<mrow><mo>(</mo>" : "") + a.toMathML(Expression.setTopLevel(fa || options == undefined || options.isTopLevel, options)) + (fa ? "<mo>)</mo></mrow>" : "") +
             (fb ? "<mrow><mo>(</mo>" : "") + b.toMathML(Expression.setTopLevel(fb, boptions)) + (fb ? "<mo>)</mo></mrow>" : "") + 
             "</msup>";
  }
  if (this.isNegation()) {
    // assert(fa === false);
      return "<mrow><mo>&minus;</mo>" + (fb ? "<mrow><mo>(</mo>" : "") + b.toMathML(Expression.setTopLevel(fb, options)) + (fb ? "<mo>)</mo></mrow>" : "") + "</mrow>";
  }  
  //TODO: fix spaces (matrix parsing)
  return "<mrow>" + 
         (fa ? "<mrow><mo>(</mo>" : "") + a.toMathML(Expression.setTopLevel(fa || options == undefined || options.isTopLevel, options)) + (fa ? "<mo>)</mo></mrow>" : "") +
         (s === '*' ? '<mo>&times;</mo>' : (s === '-' ? '<mo>&minus;</mo>' : (s === '/' ? '<mo>&#x2215;</mo>' : '<mo>' + s + '</mo>'))) +
         (fb ? "<mrow><mo>(</mo>" : "") + b.toMathML(Expression.setTopLevel(fb, options)) + (fb ? "<mo>)</mo></mrow>" : "") + 
         "</mrow>";
};

Expression.Symbol.prototype.toMathML = function (options) {
  var x = this;
  var s = x.symbol;
  var i = s.indexOf("_");
  if (i !== -1) {
    var indexes = s.slice(i + 1).replace(/^\(|\)$/g, "").split(",");
    var indexesMathML = "";
    for (var j = 0; j < indexes.length; j += 1) {
      indexesMathML += (j !== 0 ? "<mo>,</mo>" : "") +
                       (/^\d+$/.exec(indexes[j]) != undefined ? "<mn>" : "<mi>") +
                       indexes[j] +
                       (/^\d+$/.exec(indexes[j]) != undefined ? "</mn>" : "</mi>");
    }
    if (indexes.length > 1) {
      indexesMathML = "<mrow>" + indexesMathML + "</mrow>";
    }
    return "<msub>" + 
           "<mi>" + s.slice(0, i) + "</mi>" +
           indexesMathML +
           "</msub>";
  }
  return "<mi>" + s + "</mi>";
};
Expression.Negation.prototype.toMathML = function (options) {
  var b = this.b;
  var fb = this.getPrecedence() + (this.isRightToLeftAssociative() ? -1 : 0) >= b.getPrecedence();
  fb = fb || b.isUnaryPlusMinus();
  // assert(fa === false);
  return "<mrow><mo>&minus;</mo>" + (fb ? "<mrow><mo>(</mo>" : "") + b.toMathML(Expression.setTopLevel(fb, options)) + (fb ? "<mo>)</mo></mrow>" : "") + "</mrow>";
};

Condition.prototype.toMathML = function (options) {
  if (this === Condition.FALSE || this === Condition.TRUE || this.array.length === 0) {
    //throw new RangeError();
    return "";
  }
  var s = '';
  for (var i = 0; i < this.array.length; i += 1) {
    s += (i !== 0 ? '<mo>,</mo>' : '');
    s += '<mrow>';
    s += this.array[i].expression.toMathML(options) + (this.array[i].operator === Condition.NEZ ? '<mo>&ne;</mo><mn>0</mn>' : '') + (this.array[i].operator === Condition.EQZ ? '<mo>=</mo><mn>0</mn>' : '');
    s += '</mrow>';
  }
  return this.array.length === 1 ? s : '<mrow>' + s + '</mrow>';
};

Expression.Complex.prototype.toMathML = function (options) {
  return "<mrow>" + this.toStringInternal(options, "<mo>&it;</mo>", "<mi>&ii;</mi>", "<mo>&minus;</mo>", "<mo>+</mo>", function (x, options) { return x.toMathML(options); }) + "</mrow>";
};

Expression.GF2.prototype.toMathML = function (options) {
  //TODO: fix
  return this.a.toMathML(options);
};
Expression.GF2Value.prototype.toMathML = function (options) {
  return "<mrow>" + "<mn>" + this.value.toString() + "</mn>" + "</mrow>";
};
Expression.Degrees.prototype.toMathML = function (options) {
  return "<mrow>" + this.value.toMathML(options) + "<mo>&it;</mo><mi>&deg;</mi></mrow>";
};

NonSimplifiedExpression.prototype.toMathML = function (options) {
  //?
  //options = options.fractionDigits >= 0 ? Object.assign({}, options, {fractionDigits: -1}) : options;
  if (options != null && options.printId != undefined) {
    return "<mrow id=\"" + this.getId() + "\">" + this.e.toMathML(options) + "</mrow>";
  }
  return this.e.toMathML(options);
};
  
Expression.prototype.toMathML = function (options) {
  throw new Error();
};
