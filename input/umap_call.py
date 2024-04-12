import umap
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Given a precomputed distance matrix, feed into the specified UMAP model.")
parser.add_argument("--dist", action="store", dest="distance_matrix", help="Input dataframe describing distances among all samples.")
parser.add_argument("--seed", action="store", dest="random_seed", help="Random seed for reproducibility.")
parser.add_argument("--neighbor", action="store", dest="n_neighbors", default=15, help="UMAP param: n_neighbors, default: 15.")
parser.add_argument("--mindist", action="store", dest="min_dist", default=0.1,help="UMAP param: min_dist, default: 0.1.")
parser.add_argument("--components", action="store", dest="n_components", default=2, help="UMAP param: n_components, default: 2.")
parser.add_argument("--output", action="store", dest="output", default="./umap_coords.csv", help="Output dataframe describing UMAP coordinates.")
args = parser.parse_args()

unifrac_distance = pd.read_csv(args.distance_matrix, index_col=0)
reducer = umap.UMAP(n_neighbors=int(args.n_neighbors), n_components=int(args.n_components), metric="precomputed", random_state=int(args.random_seed), min_dist=float(args.min_dist))
reducer.fit(unifrac_distance)
embedding = reducer.transform(unifrac_distance)
output_df = pd.DataFrame(embedding, index=unifrac_distance.index, columns=["UMAP_1", "UMAP_2"])
output_df.insert(0, "Sample", output_df.index)
output_df.to_csv(args.output, sep=",", header=True, index=False)
