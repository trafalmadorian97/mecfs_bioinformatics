// Shared colour palette for the project's Typst diagrams.
//
// Import what you need from a sibling diagram source, e.g.
//   #import "palette.typ": ink, muted, categorical
//   #let (primary, secondary, tertiary, ..) = categorical
// so that every diagram draws from one consistent palette.
//
// The palette separates four roles, because they have different jobs:
//
//   * chrome (ink, muted)  -- lines, labels, de-emphasis. Kept neutral so the
//     diagrams sit comfortably on the Material for MkDocs page.
//   * categorical          -- colours that ENCODE data (e.g. which branch a
//     variant arose on). Chosen for perceptual distinctness and colour-blind
//     safety (the Dark2 set), NOT to match the theme's brand colour.
//   * accent               -- the theme's brand colour, for a single emphasis
//     element (a highlight, a focus marker). Do not use it to encode data.
//   * soft                 -- low-contrast, muted colours for when colour only
//     REINFORCES meaning the diagram already carries through shape, position
//     and labels (e.g. shading a region that is also labelled). Not chosen for
//     distinctness, so never use it to encode data.

// --- Chrome / neutrals ---
#let ink = rgb("#333333") // primary line / text colour
#let muted = rgb("#8a8a8a") // secondary lines, axis labels, de-emphasis

// --- Categorical data series (Dark2) ---
// Distinguishable for the common forms of colour blindness. Index in drawing
// order, or destructure into local names in the diagram.
#let categorical = (
  rgb("#1b9e77"), // teal
  rgb("#d95f02"), // orange
  rgb("#7570b3"), // purple
  rgb("#e7298a"), // magenta
)

// --- Accent ---
// Material for MkDocs' default brand hue (indigo/deep-purple family). Reserve
// for single-accent emphasis, never for categorical encoding.
#let accent = rgb("#5e35b1")

// --- Soft / low contrast ---
// A muted, earthy set (coolors.co e2d4b7-9c9583-a1a499-b0bbbf-cadbc8). The
// light members are too pale for text on a light page: use them for fills and
// keep labels in ink. Index or destructure, as with categorical.
#let soft = (
  rgb("#e2d4b7"), // sand
  rgb("#9c9583"), // taupe
  rgb("#a1a499"), // grey-green
  rgb("#b0bbbf"), // blue-grey
  rgb("#cadbc8"), // sage
)
