
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
  var sign1 = "+";
  if (coefficient.isNegative()) {
    sign1 = "&minus;";
    coefficient = coefficient.negate();//?
  }
  var coefficientString = coefficient.toMathML(options);
  var precedenceOfMultiptication = new Expression.Multiplication(Expression.ZERO, Expression.ZERO).getPrecedence();
  var areBracketsRequired = coefficient.getPrecedence() < precedenceOfMultiptication; //?
  //TODO: fix
  return (isFirst && sign1 === "+" ? "" : "<mo>" + sign1 + "</mo>") +
         (coefficient.equals(Expression.ONE) ? "" : (areBracketsRequired && coefficientString !== "" ? "<mfenced open=\"(\" close=\")\"><mrow>" + coefficientString + "</mrow></mfenced>" : coefficientString) + (coefficientString !== "" ? "<mo>&times;</mo>" : "")) +
         variable.toMathML(options);
};


var decimalToMathML = function (sign, number) {
  return (sign < 0 ? "<mrow>" : "") + (sign < 0 ? "<mo>&minus;</mo>" : "") + "<mn>" + number + "</mn>" + (sign < 0 ? "</mrow>" : "");
};

var complexToMathML = function (real, imaginary, imaginarySign) {
  return real + (imaginarySign >= 0 ? "<mo>+</mo>" : "") + imaginary + "<mo>&#x2062;</mo><mi>i</mi>";
};

//TODO: move
Expression.toDecimalString = function (x, options) {
  var fractionDigits = options.fractionDigits;
  if (fractionDigits >= 0 && !Expression.has(x, Expression.Symbol) && 
                             !Expression.has(x, Expression.NonSimplifiedExpression) &&
                             !Expression.has(x, Expression.Matrix) &&
                             !Expression.has(x, Expression.Polynomial)) {
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
  var braces = options.useBraces == undefined ? ["(", ")"] : options.useBraces;
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
    result += "<mrow>";
    result += "<menclose notation=\"none\" href=\"#\" id=\"" + containerId + "\" data-matrix=\"" + Expression.escapeHTML(x.toString()) + "\" draggable=\"true\" tabindex=\"0\" contextmenu=\"matrix-menu\">";
  }

  result += "<mfenced open=\"" + braces[0] + "\" close=\"" + braces[1] + "\">";
  result += "<mrow>";
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
    result += "<mtr>";
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
        result += row;
      }
    } else {
      while (++j < cols) {
        result += "<mtd" + (cellIdGenerator != undefined ? " id=\"" + cellIdGenerator(i, j) + "\"" : "") + ">";
        if (pivotCell != undefined && i === pivotCell.i && j === pivotCell.j) {
          result += "<mstyle mathvariant=\"bold\">";
          result += "<menclose notation=\"circle\">";
        }
        if (verticalStrike !== -1 || horizontalStrike !== -1) {
          result += "<menclose notation=\"none " + (verticalStrike === j ? " " + "verticalstrike" : "") + (horizontalStrike === i ? " " + "horizontalstrike" : "") + "\">";
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
        if (verticalStrike !== -1 || horizontalStrike !== -1) {
          result += "</menclose>";
        }
        if (pivotCell != undefined && i === pivotCell.i && j === pivotCell.j) {
          result += "</menclose>";
          result += "</mstyle>";
        }
        result += "</mtd>";
      }
    }
    result += "</mtr>";
  }
  result += "</mtable>";
  result += "</mrow>";
  result += "</mfenced>";

  if (useMatrixContainer) {
    result += "</menclose>";
    result += "</mrow>";

    result += "<mrow>";
    result += "<mtext data-x=\"TODO\">";
    result += "<button type=\"button\" class=\"matrix-menu-show matrix-menu-show-new\" data-for-matrix=\"" + containerId + "\" aria-haspopup=\"true\">&#x2630;</button>";
    result += "</mtext>";
    result += "</mrow>";

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
  return "<mfenced open=\"" + "|" + "\" close=\"" + "|" + "\"><mrow>" + x.a.toMathML(options) + "</mrow></mfenced>";
};
Expression.Transpose.prototype.toMathML = function (options) {
  var x = this;
  //TODO: ^T ?
  // https://www.w3.org/TR/MathML3/chapter4.html#contm.transpose
  return "<msup><mrow>" + x.a.toMathML(options) + "</mrow><mrow><mi>T</mi></mrow></msup>";
};
Expression.SquareRoot.prototype.toMathML = function (options) {
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
  return "<msqrt><mrow>" + this.a.toMathML(Expression.setTopLevel(true, options)) + "</mrow></msqrt>";
};
Expression.CubeRoot.prototype.toMathML = function (options) {
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
  return "<mroot><mrow>" + this.a.toMathML(Expression.setTopLevel(true, options)) + "</mrow><mrow><mi>" + 3 + "</mi></mrow></mroot>";
};
Expression.NthRoot.prototype.toMathML = function (options) {
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
  return "<mroot><mrow>" + this.a.toMathML(Expression.setTopLevel(true, options)) + "</mrow><mrow><mi>" + this.n + "</mi></mrow></mroot>";
};
Expression.Function.prototype.toMathML = function (options) {
  var x = this;
  var fa = !(x.a instanceof Expression.Matrix) && !(x.a instanceof NonSimplifiedExpression && x.a.e instanceof Expression.Matrix);//?
  //TODO: fix
  return "<mrow>" +
         (typeof i18n !== "undefined" ? "<mi>" + (x.name === "rank" ? i18n.rankDenotation : (x.name === "sin" ? i18n.sinDenotation : (x.name === "tan" ? i18n.tanDenotation : x.name))) + "</mi>" : "<mi>" + x.name + "</mi>") +
         "<mo>&#x2061;</mo>" +
         (fa ? "<mfenced open=\"(\" close=\")\"><mrow>" : "") +
         x.a.toMathML(Expression.setTopLevel(true, options)) +
         (fa ? "</mrow></mfenced>" : "") +
         "</mrow>";
};
Expression.Division.prototype.toMathML = function (options) {
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
  var x = this;
  var denominator = x.getDenominator();
  var numerator = x.getNumerator();
  //???
  //if (numerator.isNegative()) {
  //  return "<mrow><mo>&minus;</mo>" + x.negate().toMathML(options) + "</mrow>";
  //}
  return "<mfrac><mrow>" + numerator.toMathML(Expression.setTopLevel(true, options)) + "</mrow><mrow>" + denominator.toMathML(Expression.setTopLevel(true, options)) + "</mrow></mfrac>";
};
Expression.Integer.prototype.toMathML = function (options) {
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
  var x = this;
  var tmp = x.toString();
  return (tmp.slice(0, 1) === "-" ? "<mo>&minus;</mo><mn>" + tmp.slice(1) + "</mn>" : "<mn>" + tmp + "</mn>");
};
Expression.BinaryOperation.prototype.toMathML = function (options) {
  var d = Expression.toDecimalString(this, options);
  if (d != undefined) {
    return d;
  }
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
      return "<msup>" + 
             "<mrow>" +
             (fa ? "<mfenced open=\"(\" close=\")\"><mrow>" : "") + a.toMathML(Expression.setTopLevel(fa || options == undefined || options.isTopLevel, options)) + (fa ? "</mrow></mfenced>" : "") +
             "</mrow>" +
             "<mrow>" +
             (fb ? "<mfenced open=\"(\" close=\")\"><mrow>" : "") + b.toMathML(Expression.setTopLevel(fb, options)) + (fb ? "</mrow></mfenced>" : "") + 
             "</mrow>" +
             "</msup>";
  }
  if (this.isNegation()) {
    // assert(fa === false);
      return "<mrow><mo>&minus;</mo>" + (fb ? "<mfenced open=\"(\" close=\")\"><mrow>" : "") + b.toMathML(Expression.setTopLevel(fb, options)) + (fb ? "</mrow></mfenced>" : "") + "</mrow>";
  }
  //TODO: fix spaces (matrix parsing)
  return "<mrow>" + 
         (fa ? "<mfenced open=\"(\" close=\")\"><mrow>" : "") + a.toMathML(Expression.setTopLevel(fa || options == undefined || options.isTopLevel, options)) + (fa ? "</mrow></mfenced>" : "") +
         "<mo>" + (s === "*" ? "&times;" : (s === "-" ? "&minus;" : s)) + "</mo>" + 
         (fb ? "<mfenced open=\"(\" close=\")\"><mrow>" : "") + b.toMathML(Expression.setTopLevel(fb, options)) + (fb ? "</mrow></mfenced>" : "") + 
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
    return "<msub>" + 
           "<mrow>" + "<mi>" + s.slice(0, i) + "</mi>" + "</mrow>" +
           "<mrow>" + indexesMathML + "</mrow>" + 
           "</msub>";
  }
  return "<mi>" + s + "</mi>";
};
Expression.Negation.prototype.toMathML = function (options) {
  var b = this.b;
  var fb = this.getPrecedence() + (this.isRightToLeftAssociative() ? -1 : 0) >= b.getPrecedence();
  fb = fb || b.isUnaryPlusMinus();
  // assert(fa === false);
  return "<mrow><mo>&minus;</mo>" + (fb ? "<mfenced open=\"(\" close=\")\"><mrow>" : "") + b.toMathML(Expression.setTopLevel(fb, options)) + (fb ? "</mrow></mfenced>" : "") + "</mrow>";
};

Condition.prototype.toMathML = function (printOptions) {
  return this._toStringInternal(function (e) {
    return e.toMathML(printOptions);
  }, "<mo>,</mo>", "<mo>&ne;</mo><mn>0</mn>", "<mo>=</mo><mn>0</mn>");
};

Expression.Complex.prototype.toMathML = function (options) {
  return "<mrow>" + this.toStringInternal(options, "<mo>&#x2062;</mo>", "<mi>i</mi>", "<mo>&minus;</mo>", "<mo>+</mo>", function (x, options) { return x.toMathML(options); }) + "</mrow>";
};

Expression.GF2.prototype.toMathML = function (options) {
  //TODO: fix
  return this.a.toMathML(options);
};
Expression.GF2Value.prototype.toMathML = function (options) {
  return "<mrow>" + "<mn>" + this.value.toString() + "</mn>" + "</mrow>";
};
Expression.Degrees.prototype.toMathML = function (options) {
  return "<mrow>" + this.value.toMathML(options) + "<mo>&#x2062;</mo><mi>&deg;</mi></mrow>";
};

NonSimplifiedExpression.prototype.toMathML = function (options) {
  //?
  //options = options.fractionDigits >= 0 ? Object.assign({}, options, {fractionDigits: -1}) : options;
  if (options.printId != undefined) {
    return "<mrow id=\"" + this.getId() + "\">" + this.e.toMathML(options) + "</mrow>";
  }
  return this.e.toMathML(options);
};
  
Expression.prototype.toMathML = function (options) {
  throw new Error();
};
