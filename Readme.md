# Unrooted NDCH Example

`run_unrooted_ndch.py` fits a DNA substitution model with two nucleotide
composition components on the four-tip unrooted tree `((A,B),C,D)`.

Component numbers in the output start at **1**:

- Component 1 is the background composition.
- Component 2 is assigned to selected branches.

## Run

The default input files are `alignment.nex` and `tree.nwk`:

```bash
python run_unrooted_ndch.py -m F81 --comp_clade "C,D;A,B" --pre JC
```

Use different input files with `-s/--alignment` and `-t/--tree`:

```bash
python run_unrooted_ndch.py \
  -m GTR+G \
  --alignment my_alignment.nex \
  --tree my_tree.nwk \
  --comp_clade "A,B" \
  --pre my_run
```

The script is designed for taxa named `A`, `B`, `C`, and `D`.

## Selecting branches

Use `--comp` for explicit branch specifications:

```bash
python run_unrooted_ndch.py --comp "A,B|C,D;A,B|B"
```

Use `--comp_clade` to assign Component 2 to all branches inside one or more
clades. Each clade is given as comma-separated tip names. The branch
connecting a clade to the rest of the tree is excluded:

```bash
python run_unrooted_ndch.py --comp_clade "A,B"
python run_unrooted_ndch.py --comp_clade "A,B;C,D"
```

For `--comp_clade "A,B;C,D"`, the four terminal branches receive Component 2,
while the internal branch connecting the A,B and C,D sides remains Component 1.

## Output files

For a prefix `my_run`, the script writes:

- `my_run.treefile`: optimized tree with fitted branch lengths.
- `my_run.comp.treefile`: Newick tree with component numbers as branch lengths
  (for example, `:1` or `:2`).
- `my_run.ndch`: readable report containing model parameters, all component
  frequencies and Q matrices, trees, and branch assignments.
- `my_run.result.json`: machine-readable version of the fitted results.

Add `--equal-comp-freq` to fix every component to A=C=G=T=0.25. Models use
IQ-TREE notation, including modifiers such as `+I` and `+G`.

The script requires a working p4 installation. Model parsing and
configuration are provided by the local `iqtree_models.py` module.
