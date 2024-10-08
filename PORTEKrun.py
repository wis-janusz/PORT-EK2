import argparse
import os
import shutil
import portek

parser = argparse.ArgumentParser(
    description="Main PORT-EK v2 program. Use it to run all PORT-EK tools with the appropriate 'tool' argument."
)
parser.add_argument("tool", help="Name of the PORT-EK function you want to execute. Choose one of: new, find_k, enriched, map, classify")

parser.add_argument(
    "project_dir",
    help="path to the project directory. Must not exist for PORTEKrun.py new, must exist for all other tools",
    type=str,
)

parser.add_argument(
    "--max_k",
    help="Maximum k value to test with PORT-EK find_k. PORT-EK will test all odd k values from 5 up to and incliding max_k.",
    type=int
)

parser.add_argument(
    "--k",
    help="k value for PORT-EK enriched, map and classify",
    type=int
)

parser.add_argument(
    "--c",
    help="rarity filter threshold for PORT-EK enriched. Portek will only use k-mers that have frquency of more than c in at least one sample group",
    type=float
)

parser.add_argument(
    "--min_rmse",
    help="RMSE filter threshold for PORT-ek enriched. Portek will discard k-mers with RMSE lower than min_rmse",
    type=float
)

parser.add_argument(
    "--rare_m",
    help="Allowed mismatches when re-examining rare k-mers with PORT-ek enriched. Omit if not re-examining rare k-mers."
)


def _new_project(project_dir: str) -> None:
    if os.path.isdir(project_dir) == True:
        raise FileExistsError("This project directory already exists! PORT-EK does not allow overwriting projects. If you REALLY want to overwrite, remove the existing directory manually!")
    os.mkdir(project_dir)
    os.mkdir(f"{project_dir}/input")
    os.mkdir(f"{project_dir}/output")
    os.mkdir(f"{project_dir}/temp")
    shutil.copy2("./templates/config.yaml", project_dir)
    shutil.copy2("./templates/generate_kmers.sh", project_dir)
    print(f"New empty PORT-EK project created in {project_dir}. Please edit config.yaml and generte_kmers.sh files as required and copy input fasta files into {project_dir}/input.")


def main():
    args = parser.parse_args()
    if type(args.project_dir) != str:
        raise ValueError("Please provide a valid project directory name.")
    if args.tool == "new":
        _new_project(args.project_dir)
    elif args.tool == "find_k":
        optimal_k_finder = portek.FindOptimalKPipeline(args.project_dir, args.max_k)
        optimal_k_finder.find_optimal_k()
    elif args.tool == "enriched":
        enriched_kmers_finder = portek.EnrichedKmersPipeline(args.project_dir, args.k, args.c, args.min_rmse)
        enriched_kmers_finder.get_kmers()

    elif args.tool == "map":
        pass
    elif args.tool == "classify":
        pass
    else:
        raise ValueError("Unrecoginzed PORT-EK tool requested. Choose one of: new, find_k, enriched, map, classify.")


if __name__ == "__main__":
    main()

