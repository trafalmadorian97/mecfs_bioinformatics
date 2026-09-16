#import "@preview/cetz:0.3.4"
#import "palette.typ": ink, muted, soft

// Colour here only reinforces what the labels already say, so the diagram draws
// from the low-contrast soft palette rather than the categorical one.
#let (sand, taupe, ..) = soft

// fill: none keeps the SVG background transparent so it composites over whatever
// the page (or callout) background is, rather than baking in a white box.
#set page(width: auto, height: auto, margin: 10pt, fill: none)
#set text(size: 10pt)

// Standard normal density: liability is modelled as L = G + E with G and E
// independent normals, so L itself is normal. Units on the x-axis are standard
// deviations of the liability, which is why the density needs no scale
// parameter.
#let pdf(x) = calc.exp(-x * x / 2) / calc.sqrt(2 * calc.pi)

// The threshold, in liability standard deviations. 1.2 puts roughly 12% of the
// population above it -- a prevalence typical of a common disease, and wide
// enough that the shaded tail stays visible. Named threshold rather than tau
// because Typst math mode resolves identifiers from the surrounding scope, so
// a binding named tau would shadow the Greek letter in $tau$.
#let threshold = 1.2

// Canvas mapping: cm per liability SD, and cm per unit of density.
#let sx(x) = x * 1.15
#let sy(d) = d * 6.0

// Sample the curve finely enough that the SVG polyline reads as smooth.
#let step = 0.05
#let curve-points(from, to) = range(int((to - from) / step) + 1).map(i => {
  let x = from + i * step
  (sx(x), sy(pdf(x)))
})

#let x-min = -3.6
#let x-max = 3.6

#cetz.canvas(length: 1cm, {
  import cetz.draw: *

  // ---- Affected tail: the density to the right of the threshold ----
  // Closed back along the baseline so the region fills rather than the curve.
  line(
    ..curve-points(threshold, x-max),
    (sx(x-max), 0),
    (sx(threshold), 0),
    close: true,
    fill: sand,
    stroke: none,
  )

  // ---- Liability axis ----
  line((sx(x-min) - 0.25, 0), (sx(x-max) + 0.25, 0), stroke: 1pt + muted)

  // ---- Liability density ----
  line(..curve-points(x-min, x-max), stroke: 1.2pt + ink)

  // ---- Threshold ----
  // Drawn past the mode so the line reads as a cut through the whole
  // population, not just through the tail.
  line((sx(threshold), 0), (sx(threshold), sy(pdf(0)) + 0.3), stroke: (paint: taupe, thickness: 1.2pt, dash: "dashed"))
  content((sx(threshold), -0.32), text(size: 11pt, fill: ink)[$tau$])

  // ---- Region labels ----
  // "Unaffected" sits centred under the mode of the bulk it describes;
  // "Affected" cannot fit inside the thin tail, so it is set outside with a
  // leader line. Labels stay in ink: the soft fills are too pale for text.
  content((sx(0), sy(pdf(0)) * 0.42), text(size: 10pt, fill: ink)[Unaffected])
  content((sx(2.75), sy(pdf(0)) * 0.55), text(size: 10pt, fill: ink)[Affected])
  line(
    (sx(2.75), sy(pdf(0)) * 0.55 - 0.22),
    (sx(1.65), sy(pdf(1.65)) + 0.12),
    stroke: 0.7pt + taupe,
    mark: (end: ">", scale: 0.5, fill: taupe),
  )

  // ---- Axis label ----
  content((0, -0.85), text(size: 10pt, fill: ink)[Liability $L$])
})

#v(6pt)
#align(center, {
  set text(size: 9pt, fill: ink)
  [Liability $L = G + E$ is normally distributed across the population. \
  Individuals with $L > tau$ express the disease phenotype ($Y = 1$).]
})
