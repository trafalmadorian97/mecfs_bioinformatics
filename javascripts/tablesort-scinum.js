(function () {
  // Match numbers including scientific notation:
  //  - optional sign
  //  - digits with optional decimal, or leading-decimal like ".5"
  //  - optional exponent part like "e-10"
  // Also allow commas which we'll strip before parsing.
  var SCI_NUM_RE = /^[\s]*[+-]?(?:\d{1,3}(?:,\d{3})+|\d+)?(?:\.\d+)?(?:[eE][+-]?\d+)?[\s]*$/;

  function toNumber(str) {
    // Prefer data-sort if present; Tablesort supports it.
    // (Material tables are plain HTML tables rendered from Markdown.)
    var s = String(str).trim();

    // Remove common formatting
    s = s.replace(/,/g, "");

    // If empty or not parseable, return NaN
    if (!s || !SCI_NUM_RE.test(s)) return NaN;

    return parseFloat(s);
  }

  Tablesort.extend(
    "scinum",
    function (item) {
      // Tell Tablesort when this sorter applies
      var s = String(item).trim().replace(/,/g, "");
      if (!s) return false;
      // Reject plain "e" or "+" etc — require something numeric
      return SCI_NUM_RE.test(s) && !isNaN(parseFloat(s));
    },
    function (a, b) {
      var na = toNumber(a);
      var nb = toNumber(b);

      // Push NaNs to the bottom (stable-ish behavior)
      var aNan = isNaN(na);
      var bNan = isNaN(nb);
      if (aNan && bNan) return 0;
      if (aNan) return 1;
      if (bNan) return -1;

      return na - nb;
    }
  );
})();