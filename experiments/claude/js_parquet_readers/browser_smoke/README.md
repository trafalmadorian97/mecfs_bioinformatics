# Browser smoke test: Tabulator + hyparquet from CDN

Verifies, in a real headless browser rather than under node, that a docs page
can load Tabulator (UMD) and hyparquet (ESM) from CDNs, fetch a
BYTE_STREAM_SPLIT + zstd parquet over HTTP, decode it, and render it.

`run_probe.py` in the parent directory establishes that the *decoding* works
under node. This test covers what node cannot: CDN ESM resolution, UMD/ESM
coexistence on one page, and whether Tabulator actually renders the decoded rows.

## Findings

- CDN ESM import of `hyparquet` + `hyparquet-compressors` via `esm.sh` works.
- Timings for the 236 KB / 5,865-row fixture: fetch 9 ms, decode 70 ms, render 4 ms.
- Tabulator virtualises: **60** row elements in the DOM for 5,865 data rows.
- **Int64 columns break Tabulator.** hyparquet correctly decodes parquet INT64
  as JavaScript `BigInt`, and Tabulator throws
  `Cannot convert a BigInt value to a number`, rendering zero rows with the
  data otherwise decoded fine. Casting the integer columns to Int32 before
  writing fixes it (`n_snps` maxes at 1,186,054, far inside Int32). The
  alternative is converting BigInt to Number in JS after decode.
- Raw float64 renders as full precision (`0.1184091243629569`), which is ugly.
  Use Tabulator column formatters to *display* rounded values while the
  underlying data — and therefore the CSV download button — stays exact.

## Running it

Needs node and a headless browser, neither of which is in the project pixi
environment. Both are obtained without touching `pixi.toml`:

```bash
# node + JS deps, in a throwaway directory
pixi exec -s nodejs -- npm install playwright tabulator-tables
pixi exec -s nodejs -- npx playwright install chromium
```

On a minimal WSL2 system the bundled chromium fails with
`error while loading shared libraries: libnspr4.so`. Rather than `sudo apt`,
pull the libraries from conda-forge into a scratch pixi project and point the
loader at them:

```bash
pixi add nspr nss dbus libdrm libxkbcommon expat alsa-lib at-spi2-core \
         cairo pango xorg-libxcomposite xorg-libxdamage xorg-libxrandr \
         xorg-libxext xorg-libxfixes xorg-libx11 xorg-libxcb
export LD_LIBRARY_PATH=<that project>/.pixi/envs/default/lib
```

Then put a parquet in a serve directory, serve it, and run the driver:

```bash
cp <some>.parquet "$SERVE_DIR/h2_int32.parquet"
python -m http.server 8765 --directory "$SERVE_DIR" &
node drive.mjs "$SERVE_DIR" h2_int32.parquet
```

`drive.mjs` writes the test page into the serve directory itself and screenshots
the result to `smoke.png` there. Neither is checked in: the repo ignores `*.html`
and `*.png`.
