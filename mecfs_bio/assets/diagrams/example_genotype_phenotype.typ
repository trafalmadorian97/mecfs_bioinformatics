#set page(width: auto, height: auto, margin: 6pt)
#set text(size: 12pt)

// A minimal node: a rounded box around its label.
#let node(label) = box(stroke: 0.6pt, inset: 8pt, radius: 3pt, label)

#grid(
  columns: 3,
  column-gutter: 10pt,
  align: horizon,
  node[Genotype], [$arrow.r$], node[Phenotype],
)
