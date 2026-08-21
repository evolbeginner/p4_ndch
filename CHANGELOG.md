# Changelog

## v0.1.0 - Initial release
### Fixed
- `--comp A,B|C,D` to select branches
- `--comp_clade A,B` to select clades
- `--refine` to enable better optimization


## v0.0.0 - Initial release

### Added
- Unrooted four-tip NDCH analysis with p4.
- IQ-TREE-style DNA model names and rate modifiers.
- Explicit branch selection with `--comp`.
- Clade-based branch selection with `--comp_clade` / `--comp-clade`.
- Comma-separated tip-pair syntax for branch and clade specs, with
  semicolons separating multiple specs.
- Configurable alignment and tree inputs via `--alignment` / `-s` and
  `--tree` / `-t`.
- One-based component labels in user-facing output.
- Optimized tree, component-labeled Newick tree, NDCH report, and JSON result
  outputs.
- Equal composition-frequency mode via `--equal-comp-freq`.
- Basic usage documentation in `README.md`.
