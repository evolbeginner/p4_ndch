#!/usr/bin/env python3

"""
Four-tip unrooted NDCH example with command-line branch assignment.

Tree:

    ((A:1,B:1):1,C:1,D:1);

Usage with Python:

    python run_unrooted_ndch.py -m GTR+G --comp 'AB|B,AB|CD'
    python run_unrooted_ndch.py -m GTR+G --comp_clade 'AB,CD'
    python run_unrooted_ndch.py -m F81 -s alignment.nex -t tree.nwk
    python run_unrooted_ndch.py -m GTR+G --equal-comp-freq

Usage through p4:

    p4 run_unrooted_ndch.py -- -m GTR+G --comp 'AB|B,AB|CD'

Interpretation:

    comp 1 = background composition
    comp 2 = composition assigned to branches listed with --comp or
             contained within clades listed with --comp_clade

Examples:

    AB|CD
        Internal branch separating AB from CD.

    AB|A
        Branch connecting the AB ancestor to A.

    AB|B
        Branch connecting the AB ancestor to B.

    CD|C
        Branch connecting the CD side of the tree to C.

    CD|D
        Branch connecting the CD side of the tree to D.

Unspecified branches remain assigned to component 1.

The reversible DNA models and aliases from IQ-TREE's substitution-model
table are supported.  Add +I, +G, or +G<number> for rate heterogeneity.
"""

from p4 import *

import argparse
import json
import sys

from iqtree_models import configure_model, model_help_names, parse_model


# ------------------------------------------------------------
# Command-line parsing
# ------------------------------------------------------------

def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Assign component 2 to selected branches and component 1 "
            "to all other branches."
        )
    )

    parser.add_argument(
        "-m",
        "--model",
        default="HKY+G",
        help=(
            "DNA substitution model in IQ-TREE notation (default: HKY+G). "
            "Base models: %s. Modifiers: +I, +G, +G<number>."
            % model_help_names()
        )
    )

    parser.add_argument(
        "-s",
        "--alignment",
        default="alignment.nex",
        help=(
            "Input sequence alignment readable by p4 "
            "(default: alignment.nex)."
        )
    )

    parser.add_argument(
        "-t",
        "--tree",
        default="tree.nwk",
        help=(
            "Input Newick or Nexus tree readable by p4 "
            "(default: tree.nwk)."
        )
    )

    parser.add_argument(
        "--pre",
        default="iqtree",
        help=(
            "Output prefix (default: iqtree). The optimized tree and "
            "parameters are written to PREFIX.treefile and "
            "PREFIX.result.json."
        )
    )

    parser.add_argument(
        "--comp",
        default="",
        help=(
            "Comma-separated branch specifications assigned to component 2. "
            "Examples: 'AB|B', 'CD|C', 'AB|CD', "
            "or 'AB|B,CD|C'."
        )
    )

    parser.add_argument(
        "--comp_clade",
        "--comp-clade",
        default="",
        help=(
            "Comma-separated clades whose contained branches receive "
            "component 2. The branch connecting a clade to the rest of the "
            "tree is excluded. Examples: 'AB' or 'AB,CD'."
        )
    )

    parser.add_argument(
        "--equal-comp-freq",
        action="store_true",
        help=(
            "Fix every composition component to equal nucleotide "
            "frequencies (A=C=G=T=0.25) instead of estimating them."
        )
    )

    # When running through p4:
    #
    #     p4 run_unrooted_ndch.py -- --comp 'AB|B,AB|CD'
    #
    # p4 stores the arguments after the double dash in:
    #
    #     var.argvAfterDoubleDash
    #
    # When running directly with Python:
    #
    #     python run_unrooted_ndch.py --comp 'AB|B,AB|CD'
    #
    # use sys.argv[1:].
    p4_args = getattr(var, "argvAfterDoubleDash", None)

    if p4_args is not None:
        command_line_args = list(p4_args)
    else:
        command_line_args = sys.argv[1:]

    return parser.parse_args(command_line_args)


args = parse_arguments()

ALIGNMENT_FILE = args.alignment
TREE_FILE = args.tree
TREE_OUTPUT_FILE = args.pre + ".treefile"
COMPONENT_TREE_OUTPUT_FILE = args.pre + ".comp.treefile"
RESULT_FILE = args.pre + ".result.json"
REPORT_FILE = args.pre + ".ndch"


# ------------------------------------------------------------
# Read alignment and tree
# ------------------------------------------------------------

read(ALIGNMENT_FILE)
d = Data()

read(TREE_FILE)
t = var.trees[0]

t.data = d

expected_taxa = ["A", "B", "C", "D"]

if sorted(t.taxNames) != sorted(expected_taxa):
    raise RuntimeError(
        "Tree taxa do not match the expected taxa. "
        "Expected %s, found %s"
        % (expected_taxa, t.taxNames)
    )


# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

def descendant_taxa(node):
    """
    Return sorted names of all tips descended from node.
    """
    names = []

    def walk(n):
        if n.isLeaf:
            names.append(n.name)
        else:
            child = n.leftChild

            while child is not None:
                walk(child)
                child = child.sibling

    walk(node)

    return sorted(names)


def node_label(node):
    """
    Return a readable label for a node.
    """
    if node.isLeaf:
        return node.name

    return "node_%d" % node.nodeNum


def display_comp_num(internal_comp_num):
    """Convert p4's zero-based component index to a one-based output label."""
    return int(internal_comp_num) + 1


def find_tip(tree, taxon_name):
    """
    Find a tip with the specified taxon name.
    """
    for n in tree.iterNodesNoRoot():
        if n.isLeaf and n.name == taxon_name:
            return n

    raise RuntimeError(
        "Could not find tip named '%s'" % taxon_name
    )


def find_clade_node(tree, taxon_names):
    """
    Find the internal node whose descendant taxa are exactly
    taxon_names.
    """
    wanted = sorted(taxon_names)

    for n in tree.iterNodesNoRoot():
        if not n.isLeaf:
            if descendant_taxa(n) == wanted:
                return n

    raise RuntimeError(
        "Could not find an internal node subtending exactly %s"
        % wanted
    )


def split_branch_spec(specification):
    """
    Split a branch specification such as:

        AB|B
        AB|CD

    into its two sides.
    """
    specification = specification.strip()

    if not specification:
        raise ValueError(
            "An empty branch specification was supplied."
        )

    pieces = specification.split("|")

    if len(pieces) != 2:
        raise ValueError(
            "Branch specification must contain exactly one '|': %s"
            % specification
        )

    left = pieces[0].strip()
    right = pieces[1].strip()

    if not left or not right:
        raise ValueError(
            "Both sides of a branch specification are required: %s"
            % specification
        )

    return left, right


def expand_taxon_group(group, valid_taxa):
    """
    Convert a compact taxon group into a list of taxa.

    Examples:

        A    -> ["A"]
        AB   -> ["A", "B"]
        ABC  -> ["A", "B", "C"]

    This function assumes single-letter taxon names.
    """
    taxa = list(group)

    unknown_taxa = [
        taxon for taxon in taxa
        if taxon not in valid_taxa
    ]

    if unknown_taxa:
        raise ValueError(
            "Unknown taxon(s) in group '%s': %s"
            % (group, ", ".join(unknown_taxa))
        )

    if len(set(taxa)) != len(taxa):
        raise ValueError(
            "Taxon repeated in group '%s'" % group
        )

    return sorted(taxa)


def find_selected_branch(tree, specification):
    """
    Convert a user branch specification into the child node
    whose incoming branch should receive comp_1.

    Supported specifications for this tree:

        AB|CD
        AB|A
        AB|B
        CD|C
        CD|D

    The orientation is interchangeable, so these also work:

        CD|AB
        A|AB
        B|AB
        C|CD
        D|CD
    """
    left_text, right_text = split_branch_spec(specification)

    valid_taxa = set(["A", "B", "C", "D"])

    left_taxa = expand_taxon_group(
        left_text,
        valid_taxa
    )

    right_taxa = expand_taxon_group(
        right_text,
        valid_taxa
    )

    left_set = set(left_taxa)
    right_set = set(right_taxa)

    # --------------------------------------------------------
    # AB|CD:
    #
    # Select the branch leading to the internal AB node.
    # --------------------------------------------------------

    if (
        left_set == set(["A", "B"])
        and right_set == set(["C", "D"])
    ) or (
        left_set == set(["C", "D"])
        and right_set == set(["A", "B"])
    ):
        return find_clade_node(tree, ["A", "B"])


    # --------------------------------------------------------
    # AB|A:
    #
    # Select the branch leading to tip A.
    # --------------------------------------------------------

    if (
        left_set == set(["A", "B"])
        and right_set == set(["A"])
    ) or (
        left_set == set(["A"])
        and right_set == set(["A", "B"])
    ):
        return find_tip(tree, "A")


    # --------------------------------------------------------
    # AB|B:
    #
    # Select the branch leading to tip B.
    # --------------------------------------------------------

    if (
        left_set == set(["A", "B"])
        and right_set == set(["B"])
    ) or (
        left_set == set(["B"])
        and right_set == set(["A", "B"])
    ):
        return find_tip(tree, "B")


    # --------------------------------------------------------
    # CD|C and CD|D:
    #
    # Select the branch leading to the requested tip on the
    # C/D side of the internal split.
    # --------------------------------------------------------

    for taxon in ["C", "D"]:
        if (
            left_set == set(["C", "D"])
            and right_set == set([taxon])
        ) or (
            left_set == set([taxon])
            and right_set == set(["C", "D"])
        ):
            return find_tip(tree, taxon)


    raise ValueError(
        "Unsupported branch specification '%s'. "
        "For this tree, use AB|CD, AB|A, AB|B, CD|C, or CD|D."
        % specification
    )



def find_clade_branches(tree, clade_text):
    """Return branches inside a split-defined clade, excluding its boundary."""
    all_taxa = set(tree.taxNames)
    clade_taxa = set(expand_taxon_group(clade_text, all_taxa))

    if len(clade_taxa) < 2 or clade_taxa == all_taxa:
        raise ValueError(
            "Clade '%s' must contain at least two, but not all, taxa."
            % clade_text
        )

    for node in tree.iterNodesNoRoot():
        descendant_set = set(descendant_taxa(node))
        if descendant_set == clade_taxa or all_taxa - descendant_set == clade_taxa:
            break
    else:
        raise ValueError(
            "Taxa '%s' are not one side of a branch split in this tree."
            % clade_text
        )

    selected = []
    for node in tree.iterNodesNoRoot():
        descendant_set = set(descendant_taxa(node))
        opposite_set = all_taxa - descendant_set
        if descendant_set < clade_taxa or opposite_set < clade_taxa:
            selected.append(node)

    return selected

# ------------------------------------------------------------
# Parse requested branch specifications
# ------------------------------------------------------------

if args.comp.strip():
    requested_specs = [
        specification.strip()
        for specification in args.comp.split(",")
        if specification.strip()
    ]
else:
    requested_specs = []

if args.comp_clade.strip():
    requested_clades = [
        clade.strip()
        for clade in args.comp_clade.split(",")
        if clade.strip()
    ]
else:
    requested_clades = []


print("Requested component 2 branches:")

if requested_specs:
    for specification in requested_specs:
        print("  %s" % specification)
else:
    print("  none")


if len(set(requested_specs)) != len(requested_specs):
    raise ValueError(
        "The same branch specification was supplied more than once."
    )

print("Requested component 2 clades:")
if requested_clades:
    for clade in requested_clades:
        print("  %s" % clade)
else:
    print("  none")

if len(set(requested_clades)) != len(requested_clades):
    raise ValueError("The same clade was supplied more than once.")


# ------------------------------------------------------------
# Locate selected branches
# ------------------------------------------------------------

selected_nodes = {}

for specification in requested_specs:
    selected_nodes[specification] = find_selected_branch(
        t,
        specification
    )

for clade in requested_clades:
    for selected_node in find_clade_branches(t, clade):
        key = "clade %s -> %s" % (clade, node_label(selected_node))
        selected_nodes[key] = selected_node


# ------------------------------------------------------------
# Define composition components
# ------------------------------------------------------------

# comp_0 is the background composition.
#
# By default the composition values are estimated and the values below are
# starting values.  --equal-comp-freq fixes them at equal frequencies.
composition_free = 0 if args.equal_comp_freq else 1

comp_0 = t.newComp(
    partNum=0,
    free=composition_free,
    spec="specified",
    val=[0.25, 0.25, 0.25, 0.25]
)


# comp_1 is needed only for a heterogeneous run.
comp_1 = None

if selected_nodes:
    comp_1 = t.newComp(
        partNum=0,
        free=composition_free,
        spec="specified",
        val=(
            [0.25, 0.25, 0.25, 0.25]
            if args.equal_comp_freq
            else [0.40, 0.10, 0.10, 0.40]
        )
    )


# ------------------------------------------------------------
# Assign background component
# ------------------------------------------------------------

# Assign comp_0 to all branches.
#
# clade=1 applies the assignment to the whole tree.
t.setModelComponentOnNode(
    comp_0,
    node=t.root,
    clade=1
)


# ------------------------------------------------------------
# Assign comp_1 to selected branches
# ------------------------------------------------------------

for specification, selected_node in selected_nodes.items():

    print(
        "Assigning component 2 to %s -> %s; descendants = %s"
        % (
            node_label(selected_node.parent),
            node_label(selected_node),
            descendant_taxa(selected_node)
        )
    )

    # clade=0 is essential.
    #
    # It assigns comp_1 only to the branch leading to the
    # selected node, without changing descendant branches.
    t.setModelComponentOnNode(
        comp_1,
        node=selected_node,
        clade=0
    )


# ------------------------------------------------------------
# Define the nucleotide substitution model
# ------------------------------------------------------------

try:
    model_definition = parse_model(args.model)
except ValueError as error:
    raise SystemExit("error: %s" % error)

model_name = model_definition["name"]
base_model = model_definition["base_model"]
configure_model(t, model_definition)


# ------------------------------------------------------------
# Check the model
# ------------------------------------------------------------

t.modelSanityCheck()

print("\nTree with model assignments:")
t.draw(model=1)


# ------------------------------------------------------------
# Optimize likelihood
# ------------------------------------------------------------

print("\nOptimizing likelihood ...")

t.optLogLike(verbose=1)

print("\nLog likelihood: %f" % t.logLike)


# ------------------------------------------------------------
# Print substitution model, composition values, and Q matrices
# ------------------------------------------------------------

part = t.model.parts[0]
r_matrix = part.rMatrices[0]
gdasrv = part.gdasrvs[0] if part.gdasrvs else None

print("\nSubstitution model:")
print("model: %s" % model_name)
print(
    "rMatrix %d: spec=%s, free=%d, values=%s"
    % (
        r_matrix.num,
        r_matrix.spec,
        r_matrix.free,
        r_matrix.val
    )
)

if r_matrix.spec == "2p":
    print("kappa: %.10g" % r_matrix.val[0])
else:
    print("exchangeability order: AC AG AT CG CT GT")
    print("exchangeabilities: %s" % r_matrix.val)
print("IQ-TREE rate constraint code: %s" % model_definition["rate_code"])
print("p4 rate handling: %s" % model_definition["rate_mode"])
print(
    "invariant-sites proportion: %.10g (free=%d)"
    % (part.pInvar.val, part.pInvar.free)
)

if gdasrv is not None:
    print("gamma categories: %d" % part.nGammaCat)
    print("gamma alpha: %.10g (free=%d)" % (gdasrv.val[0], gdasrv.free))
    print("gamma category frequencies: %s" % gdasrv.freqs)
    print("gamma category rates: %s" % gdasrv.rates)
else:
    print("gamma categories: 1 (no among-site rate variation)")

print("\nComposition components and resulting Q matrices:")
print("Q row/column order: A C G T")

q_matrices = {}

for comp in part.comps:
    q_matrix = t.model.getBigQ(
        pNum=0,
        compNum=comp.num,
        rMatrixNum=r_matrix.num
    )
    q_matrices[comp.num] = q_matrix

    print(
        "\ncomp %d: free=%d, values=%s"
        % (
            display_comp_num(comp.num),
            comp.free,
            comp.val
        )
    )
    print("resulting Q:")
    print("             A            C            G            T")
    for state, row in zip("ACGT", q_matrix):
        print(
            "%s  %12.8f %12.8f %12.8f %12.8f"
            % ((state,) + tuple(row))
        )


# ------------------------------------------------------------
# Print branch assignments
# ------------------------------------------------------------

print("\nBranch composition assignment:")

for n in t.iterNodesNoRoot():
    parent = n.parent

    print(
        "%s -> %s : comp %d, length=%.6f, descendants=%s"
        % (
            node_label(parent),
            node_label(n),
            display_comp_num(n.parts[0].compNum),
            n.br.len,
            ",".join(descendant_taxa(n))
        )
    )


# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

t.writeNewick(TREE_OUTPUT_FILE)

# Write a second Newick tree in which each branch-length field contains its
# one-based composition component number.  Preserve the optimized lengths in
# memory; p4 itself uses zero-based component indices.
optimized_branch_lengths = {}
for n in t.iterNodesNoRoot():
    optimized_branch_lengths[n.nodeNum] = float(n.br.len)
    n.br.len = float(display_comp_num(n.parts[0].compNum))

try:
    t.writeNewick(COMPONENT_TREE_OUTPUT_FILE)
finally:
    for n in t.iterNodesNoRoot():
        n.br.len = optimized_branch_lengths[n.nodeNum]

result = {
    "logL": float(t.logLike),
    "requested_comp_2_branches": requested_specs,
    "requested_comp_2_clades": requested_clades,
    "substitution_model": {
        "name": model_name,
        "requested_name": model_definition["requested_name"],
        "base_model": base_model,
        "iqtree_rate_code": model_definition["rate_code"],
        "rate_mode": model_definition["rate_mode"],
        "equal_frequencies_in_iqtree": model_definition["equal_frequencies_in_iqtree"],
        "ndch_frequencies_optimized": not args.equal_comp_freq,
        "equal_comp_frequencies_fixed": args.equal_comp_freq,
        "invariant_sites": {
            "proportion": float(part.pInvar.val),
            "free": int(part.pInvar.free),
        },
        "r_matrix": int(r_matrix.num),
        "spec": r_matrix.spec,
        "free": int(r_matrix.free),
        "parameter_names": (
            ["kappa"] if r_matrix.spec == "2p"
            else ["AC", "AG", "AT", "CG", "CT", "GT"]
        ),
        "values": [float(value) for value in r_matrix.val],
        "gamma": (
            {
                "n_categories": int(part.nGammaCat),
                "alpha": float(gdasrv.val[0]),
                "free": int(gdasrv.free),
                "category_frequencies": [
                    float(value) for value in gdasrv.freqs
                ],
                "category_rates": [
                    float(value) for value in gdasrv.rates
                ]
            }
            if gdasrv is not None else None
        )
    },
    "q_state_order": list("ACGT"),
    "compositions": {},
    "q_matrices": {},
    "branches": {}
}


for comp in t.model.parts[0].comps:
    output_comp_num = display_comp_num(comp.num)
    result["compositions"][str(output_comp_num)] = [
        float(value) for value in comp.val
    ]
    result["q_matrices"][str(output_comp_num)] = [
        [float(value) for value in row]
        for row in q_matrices[comp.num]
    ]


for n in t.iterNodesNoRoot():
    result["branches"][node_label(n)] = {
        "parent": node_label(n.parent),
        "length": float(n.br.len),
        "comp": display_comp_num(n.parts[0].compNum),
        "descendant_taxa": descendant_taxa(n)
    }


with open(RESULT_FILE, "w") as output_file:
    json.dump(result, output_file, indent=2)


def read_text_file(file_name):
    with open(file_name) as input_file:
        return input_file.read().strip()


with open(REPORT_FILE, "w") as report:
    report.write("P4 NDCH ANALYSIS REPORT\n")
    report.write("======================\n\n")
    report.write("Input alignment: %s\n" % ALIGNMENT_FILE)
    report.write("Input tree:      %s\n" % TREE_FILE)
    report.write("Model:           %s\n" % model_name)
    report.write("Log-likelihood:  %.10f\n" % t.logLike)
    report.write("Requested branches: %s\n" % (
        ", ".join(requested_specs) if requested_specs else "none"
    ))
    report.write("Requested clades:   %s\n\n" % (
        ", ".join(requested_clades) if requested_clades else "none"
    ))

    report.write("SUBSTITUTION PROCESS\n")
    report.write("--------------------\n\n")
    report.write("Rate-matrix component: %d\n" % r_matrix.num)
    report.write("Specification: %s; free=%d\n" % (
        r_matrix.spec, r_matrix.free
    ))
    report.write("Rate order: AC AG AT CG CT GT\n")
    report.write("Rates: %s\n" % " ".join(
        "%.10g" % value for value in r_matrix.val
    ))
    report.write("Invariant-sites proportion: %.10g; free=%d\n" % (
        part.pInvar.val, part.pInvar.free
    ))
    if gdasrv is not None:
        report.write("Gamma categories: %d\n" % part.nGammaCat)
        report.write("Gamma alpha: %.10g; free=%d\n" % (
            gdasrv.val[0], gdasrv.free
        ))
        report.write("Gamma frequencies: %s\n" % " ".join(
            "%.10g" % value for value in gdasrv.freqs
        ))
        report.write("Gamma rates: %s\n" % " ".join(
            "%.10g" % value for value in gdasrv.rates
        ))
    else:
        report.write("Gamma categories: 1 (uniform rates)\n")

    report.write("\nCOMPOSITION COMPONENTS\n")
    report.write("----------------------\n")
    for comp in part.comps:
        report.write("\nComponent %d; free=%d\n" % (
            display_comp_num(comp.num), comp.free
        ))
        report.write("Frequencies (A C G T): %s\n" % " ".join(
            "%.10g" % value for value in comp.val
        ))
        report.write("Q matrix (rows/columns A C G T):\n")
        for state, row in zip("ACGT", q_matrices[comp.num]):
            report.write("  %s  %s\n" % (
                state,
                " ".join("%12.8f" % value for value in row)
            ))

    report.write("\nMAXIMUM LIKELIHOOD TREE\n")
    report.write("-----------------------\n\n")
    report.write("Optimized branch lengths:\n%s\n\n" % read_text_file(
        TREE_OUTPUT_FILE
    ))
    report.write("Component numbers encoded as branch lengths:\n%s\n\n" % (
        read_text_file(COMPONENT_TREE_OUTPUT_FILE)
    ))

    report.write("BRANCH COMPONENT ASSIGNMENTS\n")
    report.write("----------------------------\n\n")
    report.write("parent -> child  component  length  descendant_taxa\n")
    for n in t.iterNodesNoRoot():
        report.write("%s -> %s  %d  %.10g  %s\n" % (
            node_label(n.parent),
            node_label(n),
            display_comp_num(n.parts[0].compNum),
            n.br.len,
            ",".join(descendant_taxa(n))
        ))


print("\nAnalysis results written to:")
print("  Maximum-likelihood tree: %s" % TREE_OUTPUT_FILE)
print("  Component-labeled tree:  %s" % COMPONENT_TREE_OUTPUT_FILE)
print("  Model parameters:       %s" % RESULT_FILE)
print("  NDCH report:             %s" % REPORT_FILE)
