// Can a browser-side JS parquet reader handle BYTE_STREAM_SPLIT + zstd?
// Tests hyparquet (pure JS) and parquet-wasm (Rust/Arrow via wasm) against the
// same fixtures, validating decoded values rather than just "did it not throw".

import { readFileSync } from 'node:fs'
import { parquetReadObjects } from 'hyparquet'
import { compressors } from 'hyparquet-compressors'

const FIXTURE_DIR = new URL('./fixtures/', import.meta.url)
const expected = JSON.parse(readFileSync(new URL('expected.json', FIXTURE_DIR)))

const FIXTURES = [
  'f64_nobss_snappy.parquet',
  'f64_nobss_zstd22.parquet',
  'f64_bss_snappy.parquet',
  'f64_bss_zstd22.parquet',
  'f32_bss_zstd22.parquet',
]

// float32 fixtures lose precision, so the checksum needs a proportional tolerance.
const tolerance = (name) => (name.startsWith('f32') ? 1e-2 : 1e-9)

function loadBuffer(name) {
  const buf = readFileSync(new URL(name, FIXTURE_DIR))
  return buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.byteLength)
}

function validate(name, rows) {
  if (rows.length !== expected.n_rows) {
    return `WRONG ROW COUNT: got ${rows.length}, want ${expected.n_rows}`
  }
  // Checksum over the whole h2 column catches silently-misdecoded floats,
  // which is exactly the failure mode a broken BYTE_STREAM_SPLIT would cause.
  const sum = rows.reduce((acc, r) => acc + Number(r.h2), 0)
  const delta = Math.abs(sum - expected.h2_sum)
  if (delta > tolerance(name)) {
    return `WRONG VALUES: h2 sum ${sum} vs ${expected.h2_sum} (delta ${delta})`
  }
  if (rows[0].gene !== expected.first.gene) {
    return `WRONG FIRST ROW: gene ${rows[0].gene} vs ${expected.first.gene}`
  }
  return `ok  (${rows.length} rows, h2 sum delta ${delta.toExponential(2)})`
}

console.log('=== hyparquet 1.26.2 (+ hyparquet-compressors 1.1.1) ===')
for (const name of FIXTURES) {
  try {
    const rows = await parquetReadObjects({ file: loadBuffer(name), compressors })
    console.log(`  ${name.padEnd(26)} ${validate(name, rows)}`)
  } catch (err) {
    console.log(`  ${name.padEnd(26)} FAILED: ${err.message}`)
  }
}

console.log('\n=== parquet-wasm 0.7.2 ===')
try {
  const wasm = await import('parquet-wasm')
  const arrow = await import('apache-arrow').catch(() => null)
  for (const name of FIXTURES) {
    try {
      const table = wasm.readParquet(new Uint8Array(loadBuffer(name)))
      let detail = `decoded, ${table.numRows} rows`
      if (arrow) {
        const t = arrow.tableFromIPC(table.intoIPCStream())
        const col = t.getChild('h2')
        let sum = 0
        for (let i = 0; i < col.length; i++) sum += Number(col.get(i))
        const delta = Math.abs(sum - expected.h2_sum)
        detail =
          delta > tolerance(name)
            ? `WRONG VALUES: h2 sum ${sum} vs ${expected.h2_sum}`
            : `ok  (${t.numRows} rows, h2 sum delta ${delta.toExponential(2)})`
      }
      console.log(`  ${name.padEnd(26)} ${detail}`)
    } catch (err) {
      console.log(`  ${name.padEnd(26)} FAILED: ${err.message}`)
    }
  }
} catch (err) {
  console.log(`  module load FAILED: ${err.message}`)
}
