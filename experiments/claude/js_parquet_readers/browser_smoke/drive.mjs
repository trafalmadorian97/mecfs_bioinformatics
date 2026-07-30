// Drives the Tabulator + hyparquet browser smoke test.
//
// The test page is written out at runtime rather than kept as a checked-in
// .html file, because the repo ignores *.html (and *.png, so the screenshot is
// likewise untracked). Point PARQUET at a parquet whose integer columns are
// Int32 --- INT64 decodes to BigInt and Tabulator cannot format it.
//
// Usage:
//   node drive.mjs <serve-dir> <parquet-filename>
// The serve dir must already contain the parquet and be served over HTTP, e.g.
//   python -m http.server 8765 --directory <serve-dir>

import { writeFileSync } from 'node:fs'
import { join } from 'node:path'
import { chromium } from 'playwright'

const SERVE_DIR = process.argv[2] ?? '.'
const PARQUET = process.argv[3] ?? 'h2_int32.parquet'
const ORIGIN = 'http://127.0.0.1:8765'

const PAGE = `<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <title>Tabulator + hyparquet CDN smoke test</title>
  <link href="https://unpkg.com/tabulator-tables@6.3.1/dist/css/tabulator.min.css" rel="stylesheet">
  <script src="https://unpkg.com/tabulator-tables@6.3.1/dist/js/tabulator.min.js"></script>
</head>
<body>
  <div id="status">loading…</div>
  <div id="table" style="height:600px;"></div>
  <script type="module">
    // Tabulator ships UMD (loaded above as a classic script, so it lands on
    // window); hyparquet is ESM only and must come from a CDN that serves it as
    // a browser-consumable module. Whether those two coexist is the point of
    // this test --- node alone cannot answer it.
    const report = (state, detail) => {
      document.getElementById('status').textContent = state + ': ' + detail;
      window.__result = { state, detail };
    };
    try {
      const { parquetReadObjects } = await import('https://esm.sh/hyparquet@1.26.2');
      const { compressors } = await import('https://esm.sh/hyparquet-compressors@1.1.1');

      const t0 = performance.now();
      const buffer = await (await fetch('./${PARQUET}')).arrayBuffer();
      const tFetched = performance.now();
      const rows = await parquetReadObjects({ file: buffer, compressors });
      const tDecoded = performance.now();

      // Checksum the h2 column so a silently-misdecoded BYTE_STREAM_SPLIT
      // surfaces as a wrong number rather than a plausible-looking table.
      const h2Sum = rows.reduce((acc, r) => acc + Number(r.h2), 0);

      new Tabulator('#table', {
        data: rows,
        autoColumns: true,
        height: '600px',
        layout: 'fitDataFill',
      });
      const tRendered = performance.now();

      report('OK', JSON.stringify({
        rows: rows.length,
        h2Sum,
        firstGene: rows[0].gene,
        bytes: buffer.byteLength,
        fetch_ms: Math.round(tFetched - t0),
        decode_ms: Math.round(tDecoded - tFetched),
        render_ms: Math.round(tRendered - tDecoded),
      }));
    } catch (err) {
      report('ERROR', err.name + ': ' + err.message);
    }
  </script>
</body>
</html>`

writeFileSync(join(SERVE_DIR, 'index.html'), PAGE)

const browser = await chromium.launch()
const page = await browser.newPage()
const consoleErrors = []
page.on('console', m => { if (m.type() === 'error') consoleErrors.push(m.text()) })
page.on('pageerror', e => consoleErrors.push(`pageerror: ${e.message}`))
const failedRequests = []
page.on('requestfailed', r => failedRequests.push(`${r.url()} :: ${r.failure()?.errorText}`))

await page.goto(`${ORIGIN}/index.html`, { waitUntil: 'networkidle' })
await page.waitForFunction(() => window.__result !== undefined, { timeout: 60000 })
const result = await page.evaluate(() => window.__result)
console.log('result:', result.state, result.detail)

// Row count in the DOM is the virtualisation check: Tabulator should hold only
// the visible slice, not all ~5.9k rows.
console.log('tabulator rows in DOM:', await page.locator('.tabulator-row').count())
const firstCell = await page.locator('.tabulator-row .tabulator-cell').first().textContent().catch(() => null)
console.log('first cell text:', firstCell)
console.log('console errors:', consoleErrors.length ? consoleErrors : 'none')
console.log('failed requests:', failedRequests.length ? failedRequests : 'none')
await page.screenshot({ path: join(SERVE_DIR, 'smoke.png') })
await browser.close()
