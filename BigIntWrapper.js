/*global BigInt*/

(function (global) {
  "use strict";
  if (typeof BigInt !== "undefined" && BigInt(9007199254740991) + BigInt(2) - BigInt(2) === BigInt(9007199254740991)) {
    var BigIntWrapper = function () {
    };
    var prefix = function (radix) {
      if (radix === 2) {
        return "0b";
      }
      if (radix === 8) {
        return "0o";
      }
      if (radix === 16) {
        return "0x";
      }
      throw new RangeError();
    };
    BigIntWrapper.parseInt = function (string, radix) {
      return BigInt(radix === 10 ? string : prefix(radix) + string);
    };
    BigIntWrapper.compareTo = function (a, b) {
      return a < b ? -1 : (b < a ? +1 : 0);
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
    BigIntWrapper.negate = function (a) {
      return -a;
    };
    BigIntWrapper.fromNumber = function (number) {
      return BigInt(number);
    };
    BigIntWrapper.toNumber = function (bigint) {
      return Number(bigint);
    };
    global.BigIntWrapper = BigIntWrapper;
  }
}(this));

