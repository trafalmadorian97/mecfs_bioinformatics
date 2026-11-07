// The coloc.abf configuration space: every pair (S_1, S_2) of causal-variant
// indices, coloured by the hypothesis (H_0 .. H_4) it belongs to.
//
// Rows index S_1 (the causal variant of trait 1), columns index S_2 (trait 2);
// index 0 means the trait has no causal variant in the locus.

// fill: none keeps the SVG background transparent so it composites over whatever
// the page (or callout) background is, rather than baking in a white box.
#set page(width: auto, height: auto, margin: 10pt, fill: none)
#set text(size: 10pt)

// Hypothesis of the configuration (S_1, S_2) = (s1, s2).
#let hypothesis(s1, s2) = {
  if s1 == 0 and s2 == 0 { 0 } else if s2 == 0 { 1 } else if s1 == 0 { 2 } else if s1 == s2 { 4 } else { 3 }
}

// Fill colour for each hypothesis H_0 .. H_4. Colour only reinforces the cell
// label, which names the hypothesis explicitly.
#let hypothesis-fill = (
  rgb("#d9d9d9"), // H_0: grey
  rgb("#f6d7a7"), // H_1: sand
  rgb("#b7dcc8"), // H_2: green
  rgb("#bcd0e8"), // H_3: blue
  rgb("#d98b8b"), // H_4: red
)

#let n = 6
#let index-labels = ([0], [1], [2], [3], [$dots$], [$Q$])

#let cells = ([$S_1 #h(2pt) \\ #h(2pt) S_2$],)
#for s2 in range(n) { cells.push(strong(index-labels.at(s2))) }
#for s1 in range(n) {
  cells.push(strong(index-labels.at(s1)))
  for s2 in range(n) {
    let h = hypothesis(s1, s2)
    cells.push(table.cell(fill: hypothesis-fill.at(h), text(size: 8.5pt, $H_#h$)))
  }
}

#table(
  columns: (46pt,) + (30pt,) * n,
  rows: 19pt,
  align: center + horizon,
  stroke: (x, y) => if x > 0 and y > 0 { 0.4pt + white } else { none },
  ..cells,
)
