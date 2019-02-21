/*global BigInt*/

(function (global) {
  "use strict";
  if (typeof BigInt !== "undefined" && BigInt(9007199254740991) + BigInt(2) - BigInt(2) === BigInt(9007199254740991)) {
    var BigIntWrapper = function () {
    };
    BigIntWrapper.BigInt = function (x) {
      return BigInt(x);
    };
    BigIntWrapper.toNumber = function (bigint) {
      return Number(bigint);
    };
    BigIntWrapper.add = function (a, b) {
      return a + b;
    };
    BigIntWrapper.subtract = function (a, b) {
      return a - b;
    };
    BigIntWrapper.multiply = function (a, b) {
      return a * b;
    };
    BigIntWrapper.divide = function (a, b) {
      return a / b;
    };
    BigIntWrapper.remainder = function (a, b) {
      return a % b;
    };
    BigIntWrapper.unaryMinus = function (a) {
      return -a;
    };
    BigIntWrapper.lessThan = function (a, b) {
      return a < b;
    };
    global.BigIntWrapper = BigIntWrapper;
  }
}(this));

