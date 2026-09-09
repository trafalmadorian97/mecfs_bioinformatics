# Typst Diagram Task --- Design

Date: 2026-09-09

## Goal

Add a build-system Task that generates SVG vector diagrams from Typst
source, so that conceptual illustrations (causal DAGs, GWAS method
schematics, LD/fine-mapping cartoons, liability-threshold pictures,
etc.) can be produced reproducibly and consumed by the existing figure
system in mecfs_bio/figures for use in the MkDocs documentation.

## Background and rationale

### Why Typst rather than TikZ

- Lower build complexity. Typst is a single static binary available
  from conda-forge that emits SVG directly (source.typ -> typst ->
  output.svg). TikZ requires a multi-stage toolchain: a TeX engine
  (texlive-core, a multi-GB dependency, or tectonic) plus dvisvgm to
  convert DVI/XDV to SVG. The single-binary path is a far better fit
  for this repo's lockfile-controlled pixi environment and its
  deliberately constrained CI runner.
- Gentler learning curve. Typst + CeTZ was designed to avoid TikZ's
  well-known pain points (cryptic errors, the PGF layer leaking
  through, a bespoke coordinate dialect). It has a real scripting
  language underneath and millisecond compile times.
- Math is not a real obstacle. Typst's native math mode is
  LaTeX-adjacent (superscripts/subscripts and fraction/sum structure
  match; the main relearning is that symbols are spelled rather than
  backslashed). The equations that appear inside diagram nodes are
  short, which is where the friction is lowest. If a specific
  expression is easier in real LaTeX, the mitex Typst package renders
  literal LaTeX math strings. So "TikZ supports LaTeX equations and
  Typst does not" is not accurate.

TikZ retains an advantage only for a few specialized figure classes
(notably probabilistic-graphical-model plate notation via
tikz-bayesnet, and very mature plotting via PGFPlots). Neither is a v1
concern; actual data plots stay in Python/matplotlib. TikZ is
explicitly out of scope and not to be added preemptively.

### Decisions taken during brainstorming

- Source input: the task accepts either a path to a committed .typ
  file or an inline literal string of Typst code, behind one
  self-validating source type.
- Typst packages (CeTZ/Fletcher): fetched from the Typst online
  registry at compile time (network fetch allowed, Typst caches them).
  Exact versions are pinned in each .typ file's import lines so the
  cached fetch is deterministic. Vendoring is deferred unless CI ever
  flakes on the registry.
- V1 scope: Typst -> SVG only. No Graphviz backend, no Python drawsvg
  backend, no SVGO optimization stage. These were considered and
  deliberately deferred (YAGNI); they can be added later if a concrete
  need appears.
- Cache invalidation: NOT handled in the task. The VerifyingTraceRebuilder
  keys rebuild decisions off dependency and output asset hashes, never
  task code, and a diagram task is a leaf task with no deps. Editing
  either an inline-string source or a committed .typ file therefore
  will NOT automatically invalidate the cached SVG. This is accepted:
  the user forces diagram rebuilds manually via the
  must_rebuild_transitive argument to run(), the same convention used
  elsewhere in the project. No source-hash tracking is built.

## Components

Three new pieces plus edits to three registration points.

### New: DiagramFileMeta

Location: mecfs_bio/build_system/meta/diagram_file_meta.py

A lean FileMeta subclass (frozen attrs) describing a single SVG diagram
asset. Fields:

- id: AssetId (converter=AssetId), the unique asset id.
- sub_dir: PurePath defaulting to a diagrams subtree (see asset-store
  layout below).

The extension is fixed to .svg (a diagram is always an SVG), so it is a
constant of the meta rather than a caller-supplied field. Unlike
GWASPlotFileMeta it carries no trait or project: a conceptual diagram
is not tied to a GWAS, and omitting those fields keeps invalid states
unrepresentable. asset_id returns id.

### New: TypstDiagramTask

Location: mecfs_bio/build_system/task/typst_diagram_task.py

A frozen attrs GeneratingTask. Fields:

- meta: DiagramFileMeta
- source: TypstSource --- a small frozen wrapper holding exactly one of
  a .typ file path or a literal Typst string. Its __attrs_post_init__
  asserts that exactly one is set (not both, not neither), so the
  path-or-string choice lives behind one type with no representable
  invalid state, rather than two loose Optional parameters. It exposes
  a helper that yields a concrete .typ path when given a scratch
  directory (returning the existing path, or writing the string to a
  file under scratch).
- executor: Callable[[list[str]], str] defaulting to execute_command
  (from mecfs_bio/util/subproc/run_command.py), injected per the
  project's dependency-injection convention so tests can supply a fake
  runner --- though see Testing: the default test runs the real binary.

deps returns [] (leaf task in v1).

execute(scratch_dir, fetch, wf):
1. Resolve the source to a concrete .typ path (using scratch_dir for
   the inline-string case).
2. Choose an output path under scratch_dir, e.g. scratch_dir /
   (meta.asset_id + ".svg").
3. Run typst compile <src.typ> <out.svg> via self.executor as a
   list[str] command.
4. Return FileAsset(out.svg).

### Edit 1: register DiagramFileMeta in the Meta union

Location: mecfs_bio/build_system/meta/meta.py

Add DiagramFileMeta to the Meta union so it is a recognized asset
metadata type throughout the build system.

### Edit 2: asset-store path mapping

Location: mecfs_bio/build_system/rebuilder/metadata_to_path/simple_meta_to_path.py

Add an isinstance branch for DiagramFileMeta in
simple_meta_to_relative_path (and the corresponding import). Diagrams
are not GWAS-tied, so they get their own root-free subtree rather than
living under _GWAS. Introduce a module-level constant, e.g.

    _DIAGRAMS = PurePath("diagrams")

and return:

    _DIAGRAMS / m.sub_dir / (m.asset_id + ".svg")

(exact join to match the meta's sub_dir default). This determines where
the built SVG is cached in the asset store.

### Edit 3: figure-system integration

Location: mecfs_bio/figures/figure_exporter.py

- Add DiagramFileMeta to the ValidFigureMeta union so a diagram task is
  an accepted figure.
- Add a branch to get_figure_destination returning
  fig_dir / (meta.asset_id + ".svg") (a small get_diagram_fig_file_path
  helper mirroring the existing get_fig_file_path helpers). The default
  else-branch in export copies the single FileAsset to that
  destination, so no change to the copy logic is required.

## Source-file organization

Committed diagram sources live under a new package
mecfs_bio/diagrams/, a sibling of mecfs_bio/figures/.

Each diagram is:
- a .typ source file (omitted for inline-string diagrams), and
- a small Python module that constructs the TypstDiagramTask singleton
  (its DiagramFileMeta plus its source), following the same pattern as
  the assets/.../analysis/*.py modules that define figure task
  singletons.

Those singletons are imported into
mecfs_bio/figures/figure_task_list.py exactly like every other figure,
so diagrams flow through the existing generate/publish machinery with
no special-casing.

## Dependencies

Add typst from conda-forge to the pixi dependencies in pyproject.toml,
pinned to a minor range (for example >=0.15,<0.16). This is the only
new binary: no Node, no TeX distribution, no dvisvgm. CeTZ/Fletcher are
not pixi dependencies --- they are Typst-registry packages fetched at
compile time, with exact versions pinned inside each .typ file's import
lines.

## Data flow

    mecfs_bio/diagrams/<name>.typ   (or inline string in the task)
              |
              |  TypstDiagramTask.execute: typst compile ... <name>.svg
              v
        FileAsset(<name>.svg)  --- cached in asset store at
                                   diagrams/<sub_dir>/<id>.svg
              |
              |  FigureExporter.export copies it
              v
        <fig_dir>/<id>.svg     --- consumed by MkDocs docs

## Error handling

- Invalid source (both path and string set, or neither) fails at task
  construction via TypstSource.__attrs_post_init__ with a clear
  assertion message (shift-left / fail-fast).
- A Typst compile failure surfaces as the CalledProcessError raised by
  execute_command (non-zero return code), carrying the captured Typst
  output for diagnosis. No special handling is added; a broken diagram
  should fail the build loudly.

## Testing

Following project conventions (Task-level tests over helper-level
tests; no skipif library-presence guards since pixi guarantees typst;
inject dependencies rather than mock):

1. Render test: construct a TypstDiagramTask with a tiny real Typst
   snippet (inline string), execute it with the real typst binary in a
   scratch directory, and assert the output is a non-empty file whose
   leading bytes look like SVG/XML (starts with "<svg" or "<?xml"). A
   cheap "it actually rendered" check, not a brittle content
   assertion.
2. File-source test: same, but pointing the task at a small committed
   .typ fixture file, to exercise the path branch.
3. Construction test: TypstSource rejects the both-set and neither-set
   states.

No monkeypatching or mocking of the binary; typst is a real dependency.
Because these invoke the real compiler (and hit the package registry
only if the fixture imports a package --- the smoke fixtures should not,
to keep the test hermetic and fast), keep the smoke snippets to bare
Typst with no @preview imports.

## Out of scope (explicitly deferred)

- Graphviz and Python (drawsvg) diagram backends.
- SVGO / any SVG optimization or normalization stage.
- Vendoring CeTZ/Fletcher into a local package path.
- Automatic cache invalidation on source edits (handled manually via
  must_rebuild_transitive).
- TikZ / any TeX toolchain.
