from langchain_core.tools import tool
import pooch
import scanpy as sc
import anndata as ad

global EXAMPLE_DATA
EXAMPLE_DATA = pooch.create(
    path=pooch.os_cache("scverse_tutorials"),
    base_url="doi:10.6084/m9.figshare.22716739.v1/",
)
EXAMPLE_DATA.load_registry_from_doi()

#Tools
@tool
def load():
    """Loads sample data."""
    samples = {
    "s1d1": "s1d1_filtered_feature_bc_matrix.h5",
    "s1d3": "s1d3_filtered_feature_bc_matrix.h5",
    }
    global adatas
    adatas = {}
    for sample_id, filename in samples.items():
        path = EXAMPLE_DATA.fetch(filename)
        sample_adata = sc.read_10x_h5(path)
        sample_adata.var_names_make_unique()
        adatas[sample_id] = sample_adata
    global adata
    adata = ad.concat(adatas, label="sample")
    adata.obs_names_make_unique()
    return "Data Loaded."

@tool
def get_var_obs():
    """Get variables."""
    global adata
    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    return "Variables and observations completed."
    
@tool
def calculate_qc_metrics():
    """Calculate the quanlity control (qc) metrics."""
    global adata
    sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )
    return "Quality control completed."


@tool
def plot_violin_plot():
    """Plots a violin plot."""
    global adata
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
    )
    return "The Violin plot is drawn."

@tool
def plot_scatter():
    """Plots a scatter plot."""
    global adata
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
    return "The scatter plot is drawn."

@tool
def filter_cell_gene():
    """Filter out the cells that express to little genes and filter out genes expressed too low genes"""
    sc.pp.filter_cells(adata, min_genes=100)
    sc.pp.filter_genes(adata, min_cells=3)
    return "The genes have been filtered successfully."

@tool
def doublet_detection():
    """Identify doublets."""
    sc.external.pp.scrublet(adata, batch_key="sample")
    return "Doublet detection finished."

@tool
def normalization():
    """Normalize the data."""
    # Saving count data
    adata.layers["counts"] = adata.X.copy()
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    return "Data has been normalized."

@tool
def feature_selection():
    """Select the highly variable genes"""
    global adata
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
    sc.pl.highly_variable_genes(adata)
    return "Highly variable genes have been selected."

@tool
def pca_dimention_reduction():
    """Reduces dimentions of data using the PCA method."""
    global adata
    sc.tl.pca(adata)
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)
    sc.pl.pca(
        adata,
        color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
        dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
        ncols=2,
        size=2,
    )
    return "PCA dimention reduction complete."

@tool
def compute_neighbor_plot():
    """Compute the neighborhood graph of cells using the PCA representation of the data matrix."""
    global adata
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.pl.umap(
        adata,
        color="sample",
        # Setting a smaller point size to get prevent overlap
        size=2,
    )
    return "Neighborhood graph completed."

@tool
def clustering():
    """Use Leiden clustering directly to cluster the neighborhood graph of cell"""
    global adata
    # Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
    sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
    #sc.tl.leiden(adata, n_iterations=2)
    sc.pl.umap(adata, color=["leiden"])
    sc.pl.umap(
        adata,
        color=["leiden", "predicted_doublet", "doublet_score"],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
    )
    sc.pl.umap(
        adata,
        color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
        wspace=0.5,
        ncols=2,
    )
    return "Data has been clustered."
    
tools = [load, 
         get_var_obs, 
         calculate_qc_metrics, 
         plot_violin_plot, 
         plot_scatter, 
         filter_cell_gene, 
         doublet_detection, 
         normalization, 
         feature_selection,
         pca_dimention_reduction,
         compute_neighbor_plot,
         clustering
        ]






        