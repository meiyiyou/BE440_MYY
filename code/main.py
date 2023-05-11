from preprocess import annotate, utils, gene_sets
from analysis import gsea as gs
from analysis import network
import logging
import datetime
import click
import tomllib
import datetime

gene_set_paths = {
        "Auxin signaling": "P:auxin-activated signaling pathway", 
        "Ethylene signaling": "P:ethylene-activated signaling pathway",
        "Cytokinin signaling": "P:cytokinin-activated signaling pathway",
        "Jasmonic acid signaling": "P:jasmonic acid mediated signaling pathway",
        "Gibberellin signaling": "P:gibberellin mediated signaling pathway",
        "Abscisic acid signaling": "P:abscisic acid-activated signaling pathway",
        "Brassinosteroid signaling": "P:brassinosteroid mediated signaling pathway",
        "Salicylic acid signaling": "P:salicylic acid mediated signaling pathway"}
        
@click.group()
@click.pass_context
def cli(ctx):
    pass

@cli.command()
@click.pass_context
@click.option("--is_ctrl", default=True, help="Whether you're interested in the control or MCP-1 phenotype")
@click.option("--is_ethylene", default=True, help="Whether you're interested in looking hormone phenotype")
@click.option("--cultivar", default="golden_delicious")
@click.option("--plot/--no-plot", default=True, help="enables plotting")
def gsea(ctx, is_ctrl, is_ethylene, cultivar, plot):
    phenotype_golden_delicious_ctrl = [3, 0, 2.0349, 43.3721, 56.09, 52.0930]
    gs.gsea(
        phenotype=phenotype_golden_delicious_ctrl,
        is_ctrl=is_ctrl,
        gene_set_paths=gene_set_paths,
        title_gsea="Testing CLI",
        title_heatmap="Test heatmap",
        plot_enabled=plot
    )

@cli.command()
@click.pass_context
def preprocess(ctx):
    click.echo("Doing pre-processing stuff")

@cli.command()
@click.pass_context
@click.option("--input_xlsx", default="./data/processed/02_annotated_fpkm_v4.xlsx", help="input path for annotated FPKM")
@click.option("--output_dir", default="./data/processed/03_gene_sets")
def make_gene_sets(ctx, input_xlsx, output_dir):
    click.echo("Making new gene sets...")
    gene_sets.gene_setter(input_xlsx, output_dir, gene_set_paths.values())

@cli.command()
@click.pass_context
@click.option("--compute-corr/--no-compute-corr", default=False, help="compute all combinations of pairwise correlation")
@click.option("--compute-path", default=f"./data/processed/03_gene_sets/{datetime.datetime.now()}", help="output directory for the computed correlation matrix saved as an .npy file")
@click.option("--corr-path", default="./data/processed/03_gene_sets/corr_matrix_ctrl.npy")
@click.option("--corr-title", default="Network title")
def go_network(ctx, compute_corr, compute_path, corr_path, corr_title):
    click.echo("Generating GO network...")
    if compute_corr:
        network.compute_corr(gene_set_paths=gene_set_paths)
    network.graph_corr(corr_matrix_file=corr_path, title=corr_title)

@cli.command()
@click.pass_context
def graph_dotplot(ctx):
    gs.dotplot()

if __name__ == "__main__":
    config_project = None
    with open("config.toml", "rb") as f:
        config_project = tomllib.load(f)
    # print(data)
    cli(obj=config_project)
    # logging.basicConfig(filename="./data/test/logging.txt",
    #                     level=logging.ERROR)

    # logging.error(f"Starting session at {datetime.datetime.now()}")
    # #utils.clean_intermediates()
    # input_file = "./data/processed/01_processed_fpkm_.xlsx"
    # path_annotated_fpkm = "./data/processed/02_annotated_fpkm_v3.xlsx"
    # properties=[
    #     "keywords", "go_gene_set", "protein_existence", 
    #     "primary_accession", "uniProtkbId", "protein", "pipeline"]
    # # annotate.collect_fpkm_annotations(input_file, properties=properties, autosave_interval_count=25)
    # # annotate.preprocess_fpkm_annotations(input_file, output_xlsx=path_annotated_fpkm, properties=properties)

    # gene_ontologies = [
    #     "P:ethylene-activated signaling pathway",
    #     "P:auxin-activated signaling pathway"]
        
    # path_go_output_dir = "./data/processed/03_gene_sets"
    # gene_sets.gene_setter(input_xlsx=path_annotated_fpkm, output_dir=path_go_output_dir, gene_ontologies=gene_ontologies)