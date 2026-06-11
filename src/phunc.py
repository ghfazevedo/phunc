#!/usr/bin/env python3

import argparse
import os
import dendropy
from dendropy.simulate import treesim
from dendropy.model import reconcile
from dendropy.model import coalescent
#from dendropy.model.reconcile import monophyletic_partition_discordance
from dendropy.model.parsimony import fitch_down_pass
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
from phunc import __version__
#from scipy.stats import lognorm

# Function to parse arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculates the probability of fixation of differences in a hypothetical nuclear locus that controls phenotype under neutral divergence.")
    parser.add_argument("-t", "--tree", required=True, help="Path to the population/species tree with population size as branch annotations (Nexus format)")
    parser.add_argument("-o", "--out_dir", default="./phenofun_out", help="Output directory")
    parser.add_argument("-n", "--n_simulations", default="100", help="Number of gene trees to simulate.")
    parser.add_argument("-s", "--n_sampled_individuals", type=str, required=True, help="The number of individuals per population/species separated by comma. It should be in the same order as the populations appear in the species tree file.")
    parser.add_argument("-ts", "--target_s_statistics", type=int, required=True, help="The target s statistics as observed in the real world data to calculate the probability of generating it through the simulations.")
    parser.add_argument("-p", "--phenotype_map", type=str, required=False,
        help="Optional: Path to a tab-separated file associating species in the tree with phenotype codes. Format: one line per species, e.g. 'species1\t0', 'species2\t0', 'species3\t1', 'species4\t[0,1]'. Use [0,1] for uncertain or polymorphic states. If provided, this file will be used to assign phenotype states to taxa instead of automatically set one different state per species on tree.")
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    
    return parser.parse_args()

# Function to confirm with the user if directory exists
def confirm_proceed(message="Directory already exists. Do you want to proceed? (y/n): "):
    while True:
        response = input(message).strip().lower()
        if response == 'y':
            return True
        elif response == 'n':
            print("Exiting program.")
            return False
        else:
            print("Please enter 'y' or 'n'.")

def main():
    args = parse_arguments()

    # Convert out_dir to absolute path
    args.tree = os.path.abspath(args.tree)
    args.out_dir = os.path.abspath(args.out_dir)

    # Convert the number of individuals to a list
    n_sampled_individuals = list(map(int, args.n_sampled_individuals.split(',')))
    
    # convert string to integer
    args.n_simulations=int(args.n_simulations)

    # Check if the directory exists and confirm with the user
    if os.path.exists(args.out_dir):
        print(f"Warning: Directory '{args.out_dir}' already exists. Proceeding may erase previous outputs.")
        if not confirm_proceed():
            exit()

    # Create output directory
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    # Read tree and get taxa names
    containing_taxa = dendropy.TaxonNamespace()
    sp_tree = dendropy.Tree.get(path=args.tree, 
                                schema="nexus", 
                                preserve_underscores=True,
                                taxon_namespace=containing_taxa)

    genes_to_species = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        containing_taxon_namespace=containing_taxa,
        num_contained=n_sampled_individuals)

    # convert to containing tree
    sp_tree = reconcile.ContainingTree(sp_tree,
                contained_taxon_namespace=genes_to_species.domain_taxon_namespace,
                contained_to_containing_taxon_map=genes_to_species)

    # Simulate and save gene trees
    trees = dendropy.TreeList()
    print('Simulating trees')
    for rep in range(args.n_simulations):
        print(rep)
        gene_tree = treesim.contained_coalescent_tree(containing_tree=sp_tree, gene_to_containing_taxon_map=genes_to_species)
        trees.append(gene_tree)
    
    print('Saving newick simulated trees')
    trees.write(path= os.path.join(args.out_dir, "simulated_gene_trees.nwck"),
        schema="newick"
        )

    # Attribute code (phenotype states) to species
    taxon_state_sets_map = {}
    if args.phenotype_map:
        # Read phenotype map file
        phenotype_map = {}
        with open(args.phenotype_map, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) != 2:
                    raise ValueError(f"Invalid phenotype_map line: {line}")
                species, code_str = parts
                # Handle code as int, list, or set
                if code_str.startswith('[') and code_str.endswith(']'):
                    # Parse list, e.g. [0,1]
                    code = set(eval(code_str))
                else:
                    code = set([int(code_str)])
                phenotype_map[species] = code
        for taxon in trees.taxon_namespace:
            species = taxon.label.split()[0]
            if species not in phenotype_map:
                raise ValueError(f"Species '{species}' not found in phenotype_map file.")
            taxon_state_sets_map[taxon] = [phenotype_map[species]]
    else:
        species_to_code = {}
        current_code = 0
        for taxon in trees.taxon_namespace:
            species = taxon.label.split()[0]
            if species not in species_to_code:
                species_to_code[species] = current_code
                current_code += 1
            code = species_to_code[species]
            taxon_state_sets_map[taxon] = [set([code])]

    # Iterate over trees to calculate s
    s_count = 0 
    target_trees = dendropy.TreeList()
    s_distribution = []

    for tree in trees:
        s = fitch_down_pass(tree.postorder_node_iter(),
                            taxon_state_sets_map=taxon_state_sets_map)
        s_distribution.append(s)
        if s == args.target_s_statistics:
            s_count = s_count + 1
            target_trees.append(tree)

    #The commented code below uses the monophyletic_partition_discordance() which, contrary to statement in DendroPy manual,
    #  does not seem to be exactly the s statistics. So I chenged to the code above.
    
    ## Create the function to get species name of taxon object label of the gene trees. 
    ## This will be used for the taxa membership, since the simulated tree has taxa with names like "species1 0", "species1 1", "species2 0"
    #def mf(t):
    #    index=t.label.find(" ")
    #    return t.label[:index] 
    #for tree in trees:
    #    taxon_namespace = tree.taxon_namespace
    #    tax_parts = taxon_namespace.partition(membership_func=mf)
    #    s = monophyletic_partition_discordance(tree, taxon_namespace_partition=tax_parts)
    #    s_distribution.append(s)
    #    if s == args.target_s_statistics:
    #        s_count = s_count + 1
    #        target_trees.append(tree)


    probability = s_count/args.n_simulations

    print('Saving newick simulated trees with target s-statistics')
    target_trees.write(path= os.path.join(args.out_dir, "target_gene_trees.nwck"),
        schema="newick"
        )

    print(f"Probability of s statistics being equal to {args.target_s_statistics}: {probability}")

    results = open(os.path.join(args.out_dir, "S_statsProbs.txt"), "w")
    print(f"Probability of s statistics being equal to '{args.target_s_statistics}': '{probability}'", file=results)
    results.close()

    # Save simulated s values
    df = pd.DataFrame({'simulation': range(1, len(s_distribution) + 1), 's_statistic': s_distribution})
    df.to_csv(os.path.join(args.out_dir,"simulated_s.csv"), index=False)
    print("Simulated s values saved in", os.path.join(args.out_dir,"simulated_s.csv") )

    # Create a histogram of s values
    # Compute lower the 95% confidence interval
    lower_bound = np.percentile(s_distribution, 5.0)

    # Create histogram
    # Set number of bins to the maximum value in s_distribution
    max_s = int(max(s_distribution))
    # Create integer bins from min to max (inclusive)
    bins = np.arange(int(min(s_distribution)), max_s + 2)  # +2 to include the last value
    hist_values, bin_edges, patches = plt.hist(
        s_distribution,
        bins=bins,
        density=True,
        edgecolor='black',
        alpha=0.7,
        histtype='bar',
        rwidth=1.0  # bars touch each other
    )

    # Color bars conditionally
    for patch, left_edge in zip(patches, bin_edges[:-1]):
        if left_edge < lower_bound:
            patch.set_facecolor('red')  # Outliers in red
        else:
            patch.set_facecolor('blue')  # Normal density in blue


    # Simple histogram without 95% interval
    #plt.hist(s_distribution, bins=10, edgecolor='black', alpha=0.7)

    # Add vertical line at 'target'
    plt.axvline(x=args.target_s_statistics, color='r',
                linestyle='dashed', linewidth=2,
                label=f'Observed s = {args.target_s_statistics}')

    # Compute and plot the smooth PDF curve with gaussian
    #kde = gaussian_kde(s_distribution)  # Kernel Density Estimation
    #x_vals = np.linspace(min(s_distribution), max(s_distribution), 1000)  # X range for smooth curve
    #plt.plot(x_vals, kde(x_vals), color='black', linewidth=2, label="PDF Curve")  # Smooth PDF

    # Compute and plot the smooth PDF curve with lognormal distribution
    #shape, loc, scale = lognorm.fit(s_distribution, floc=0)  # Estimate parameters
    #x_vals = np.linspace(min(s_distribution), max(s_distribution), 1000)  # Smooth x values
    #pdf_vals = lognorm.pdf(x_vals, shape, loc, scale)  # Compute lognormal PDF
    ## Plot the lognormal PDF curve
    #plt.plot(x_vals, pdf_vals, color='black', linewidth=2, linestyle='dashed', label="Lognormal PDF")

    # Add Label at the Top, Close to the Line 
    y_top = plt.ylim()[1] * 1.01  # Position above the highest histogram bar
    plt.text(args.target_s_statistics, y_top, f"p(s={args.target_s_statistics})={probability}", color='black', ha='center', va='bottom', fontsize=12, fontweight='bold')

    # Labels and title
    plt.xlabel('Simulated s statistics')
    plt.ylabel('Probability')

    # Save as PNG and PDF
    plt.savefig(os.path.join(args.out_dir, "histogram.png"), dpi=300)
    plt.savefig(os.path.join(args.out_dir,"histogram.pdf") )

    print("Histograms saved in", os.path.join(args.out_dir,"histogram.pdf"), "and in", os.path.join(args.out_dir,"histogram.png") )

    # BAR PLOT: p(s=target | drift) vs p(s ≠ target) with evidence lines ---
    prob_target = probability
    prob_different = 1 - probability
    bar_labels = [f'p(s = {args.target_s_statistics} | drift)', f'p(s ≠ {args.target_s_statistics} | drift)']
    bar_values = [prob_target, prob_different]

    fig, ax = plt.subplots(figsize=(6, 6))
    bars = ax.bar(bar_labels, bar_values, color=['#377eb8', '#e41a1c'], edgecolor='black', alpha=0.8)

    # Evidence lines and labels (adapted from R code)
    BF_vals = np.array([3.2, 10, 100])
    prior = 0.5
    strength = BF_vals / (BF_vals + 1)
    y_lines = np.concatenate(([prior], strength))
    labels = ["no support", "weak", "substantial", "strong"]

    # Draw horizontal lines
    ax.axhline(prior, linestyle='solid', color='grey', linewidth=1)
    linestyles = ['dotted', 'dashed', (0, (5, 10))]  # solid already used for prior
    for y, ls in zip(strength, linestyles):
        ax.axhline(y, linestyle=ls, color='grey', linewidth=1)

    # Annotate evidence labels at the right edge
    y_for_labels = y_lines
    for y, label in zip(y_for_labels, labels):
        ax.text(1.05, y, label, ha='left', va='center', color='grey', fontsize=10, transform=ax.get_yaxis_transform())

    # Annotate bar values
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 5),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=12, fontweight='bold')

    ax.set_ylim(0, 1.05)
    ax.set_ylabel('Probability')
    ax.set_title('Probability of s = target and s > target')

    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "barplot.png"), dpi=300)
    plt.savefig(os.path.join(args.out_dir, "barplot.pdf"))
    print("Bar plot saved in", os.path.join(args.out_dir, "barplot.pdf"), "and in", os.path.join(args.out_dir, "barplot.png"))

if __name__ == "__main__":
    main()