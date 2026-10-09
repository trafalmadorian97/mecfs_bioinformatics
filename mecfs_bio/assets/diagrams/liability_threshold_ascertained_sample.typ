#import "@preview/cetz:0.3.4"
#import "palette.typ": ink, muted, soft

// Colour here only reinforces what the labels already say, so the diagram draws
// from the low-contrast soft palette rather than the categorical one.
#let (sand, taupe, ..) = soft

// fill: none keeps the SVG background transparent so it composites over whatever
// the page (or callout) background is, rather than baking in a white box.
#set page(width: auto, height: auto, margin: 10pt, fill: none)
#set text(size: 10pt)

// Standard normal density of liability in the population, in liability
// standard deviations.
#let pdf(x) = calc.exp(-x * x / 2) / calc.sqrt(2 * calc.pi)

// The threshold, in liability standard deviations. Matches the population
// diagram (liability_threshold_model.typ) so the two read as a pair. Named
// threshold rather than tau because Typst math mode resolves identifiers from
// the surrounding scope, so a binding named tau would shadow the Greek letter
// in $tau$.
#let threshold = 1.2

// Population prevalence K = P(L > threshold). Typst has no erf, so integrate
// the density over the upper tail with Simpson's rule; the tail beyond 8 SD is
// negligible.
#let upper-tail-mass(from) = {
  let to = 8.0
  let n = 400
  let h = (to - from) / n
  let total = pdf(from) + pdf(to)
  for i in range(1, n) {
    total += (if calc.odd(i) { 4 } else { 2 }) * pdf(from + i * h)
  }
  total * h / 3
}
#let prevalence = upper-tail-mass(threshold)

// Fraction of cases in the ascertained sample. A balanced case/control design
// recruits equal numbers of each.
#let case-fraction = 0.5

// Within each side of the threshold, ascertainment leaves the shape of the
// liability distribution unchanged and only rescales it, so that the mass above
// the threshold becomes case-fraction and the mass below becomes its complement.
#let case-scale = case-fraction / prevalence
#let control-scale = (1 - case-fraction) / (1 - prevalence)
#let sample-pdf(x) = if x > threshold { case-scale * pdf(x) } else { control-scale * pdf(x) }

// Canvas mapping: cm per liability SD, and cm per unit of density. Same scale
// as the population diagram, so the reader can compare heights directly.
#let sx(x) = x * 1.15
#let sy(d) = d * 6.0

// Sample the curve finely enough that the SVG polyline reads as smooth.
#let step = 0.05
#let curve-points(density, from, to) = range(int(calc.round((to - from) / step)) + 1).map(i => {
  let x = from + i * step
  (sx(x), sy(density(x)))
})

#let x-min = -3.6
#let x-max = 3.6

// Height of the case density at the threshold: the peak of the whole curve.
#let case-peak = sy(case-scale * pdf(threshold))

#cetz.canvas(length: 1cm, {
  import cetz.draw: *

  // ---- Cases: the rescaled density to the right of the threshold ----
  // Closed back along the baseline so the region fills rather than the curve.
  line(
    ..curve-points(x => case-scale * pdf(x), threshold, x-max),
    (sx(x-max), 0),
    (sx(threshold), 0),
    close: true,
    fill: sand,
    stroke: none,
  )

  // ---- Liability axis ----
  line((sx(x-min) - 0.25, 0), (sx(x-max) + 0.25, 0), stroke: 1pt + muted)

  // ---- Population density, for reference ----
  line(..curve-points(pdf, x-min, x-max), stroke: (paint: muted, thickness: 0.9pt, dash: "dotted"))

  // ---- Sample density ----
  // Drawn as two pieces joined by a vertical jump at the threshold, since the
  // two sides are rescaled by different factors.
  line(..curve-points(x => control-scale * pdf(x), x-min, threshold), stroke: 1.2pt + ink)
  line(
    (sx(threshold), sy(control-scale * pdf(threshold))),
    (sx(threshold), case-peak),
    stroke: 1.2pt + ink,
  )
  line(..curve-points(x => case-scale * pdf(x), threshold, x-max), stroke: 1.2pt + ink)

  // ---- Threshold ----
  // Below the jump the dashed line is hidden by the solid density, so it
  // stays visible only where it marks the cut through the empty space above.
  line((sx(threshold), 0), (sx(threshold), case-peak + 0.3), stroke: (paint: taupe, thickness: 1.2pt, dash: "dashed"))
  content((sx(threshold), -0.32), text(size: 11pt, fill: ink)[$tau$])

  // ---- Region labels ----
  // "Cases" is set outside the tail with a leader line, as in the population
  // diagram: inside, it would collide with the dotted population curve. Labels
  // stay in ink: the soft fills are too pale for text.
  content((sx(0), sy(control-scale * pdf(0)) * 0.42), text(size: 10pt, fill: ink)[Controls])
  content((sx(2.75), case-peak * 0.6), text(size: 10pt, fill: ink)[Cases])
  line(
    (sx(2.75), case-peak * 0.6 - 0.22),
    (sx(1.5), 1.2),
    stroke: 0.7pt + taupe,
    mark: (end: ">", scale: 0.5, fill: taupe),
  )

  // ---- Population curve label ----
  // Set left of the mode, with a leader to the dotted curve, where the gap
  // between the population and sample densities is widest.
  content((sx(-2.4), sy(pdf(0)) * 1.2), text(size: 9pt, fill: muted)[Population])
  line(
    (sx(-2.15), sy(pdf(0)) * 1.2 - 0.2),
    (sx(-0.9), sy(pdf(-0.9)) + 0.06),
    stroke: 0.7pt + muted,
    mark: (end: ">", scale: 0.5, fill: muted),
  )

  // ---- Axis label ----
  content((0, -0.85), text(size: 10pt, fill: ink)[Liability $L$])
})

#v(6pt)
#align(center, {
  set text(size: 9pt, fill: ink)
  [In a sample with equal numbers of cases and controls, the density above $tau$ \
    is scaled up and the density below $tau$ is scaled down, so each side holds half the mass.]
})
