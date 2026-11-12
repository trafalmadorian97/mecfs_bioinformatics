// Parse TeX equations with MathJax and report any that MathJax renders as an error.
//
// Usage: node check_tex.cjs <site_config.js> <input.json> <output.json>
//
// site_config.js is the docs site's MathJax configuration script (it assigns window.MathJax).
// Its tex block is reused here, so the checker parses TeX exactly as the site does.
//
// input.json:  {"pages": [{"equations": [{"tex": "...", "display": true}, ...]}, ...]}
// output.json: {"mathjax_version": "3.2.2", "pages": [["error or null", ...], ...]}
//
// Each page's equations are parsed in document order, and the TeX state (equation labels and
// numbers) is reset between pages, as in the browser.  This catches labels defined twice on one
// page.
//
// Two kinds of failure are detected:
//   - parse errors, which MathJax renders as an merror element;
//   - undefined control sequences, which MathJax's noundefined package renders as coloured text
//     rather than as an error.  The checker sets that colour to a sentinel value so it can find
//     them without confusing them with deliberately coloured text.
//
// MathJax is loaded from the conda-forge mathjax package in the pixi environment (CONDA_PREFIX).

"use strict";

const fs = require("fs");
const vm = require("vm");

const UNDEFINED_SENTINEL_COLOR = "mjx-check-undefined-macro";

function loadSiteTexConfig(configPath) {
  const sandbox = { window: {}, document$: { subscribe: () => {} } };
  vm.runInNewContext(fs.readFileSync(configPath, "utf8"), sandbox);
  const siteConfig = sandbox.window.MathJax;
  if (!siteConfig || !siteConfig.tex) {
    throw new Error(`${configPath} does not assign window.MathJax with a tex block`);
  }
  return siteConfig.tex;
}

function unescapeXml(text) {
  return text
    .replace(/&lt;/g, "<")
    .replace(/&gt;/g, ">")
    .replace(/&quot;/g, '"')
    .replace(/&#x27;/g, "'")
    .replace(/&amp;/g, "&");
}

function findError(mathml) {
  const merror = mathml.match(/<merror data-mjx-error="([^"]*)"/);
  if (merror) {
    return unescapeXml(merror[1]);
  }
  const undefinedMacro = mathml.match(
    new RegExp(`<mtext mathcolor="${UNDEFINED_SENTINEL_COLOR}">([^<]*)</mtext>`),
  );
  if (undefinedMacro) {
    return `Undefined control sequence ${unescapeXml(undefinedMacro[1])}`;
  }
  return null;
}

async function checkEquation(MathJax, equation) {
  try {
    const mathml = await MathJax.tex2mmlPromise(equation.tex, { display: equation.display });
    return findError(mathml);
  } catch (err) {
    return `MathJax raised: ${err.message}`;
  }
}

async function main() {
  const [configPath, inputPath, outputPath] = process.argv.slice(2);
  if (!configPath || !inputPath || !outputPath) {
    throw new Error("usage: node check_tex.cjs <site_config.js> <input.json> <output.json>");
  }
  const mathjaxRoot = `${process.env.CONDA_PREFIX}/lib/mathjax/es5`;
  const siteTex = loadSiteTexConfig(configPath);
  const MathJax = await require(`${mathjaxRoot}/node-main.js`).init({
    loader: { paths: { mathjax: mathjaxRoot }, load: ["input/tex"] },
    tex: { ...siteTex, noundefined: { color: UNDEFINED_SENTINEL_COLOR } },
    startup: { typeset: false },
  });

  const input = JSON.parse(fs.readFileSync(inputPath, "utf8"));
  const pages = [];
  for (const page of input.pages) {
    MathJax.texReset();
    const errors = [];
    for (const equation of page.equations) {
      errors.push(await checkEquation(MathJax, equation));
    }
    pages.push(errors);
  }
  fs.writeFileSync(outputPath, JSON.stringify({ mathjax_version: MathJax.version, pages }));
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});
