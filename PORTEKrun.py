import argparse
import os
import shutil
import resource
from datetime import datetime

import portek

parser = argparse.ArgumentParser(
    description="Main PORT-EK v2 program. Use it to run all PORT-EK tools with the appropriate 'tool' argument."
)
parser.add_argument("tool", help="Name of the PORT-EK function you want to execute. Choose one of: new, find_k, enriched, map, classify")

parser.add_argument(
    "project_dir",
    help="Path to the project directory. Must not exist for PORTEKrun.py new, must exist for all other tools",
    type=str,
)

parser.add_argument(
    "--skip_kf",
    help="Skip finding k-mers in input sequences. Use with PORT-EK find_k if you already have k-mer indices",
    default=False,
    action="store_true"
)

parser.add_argument(
    "--min_k",
    help="Minimum k value to test with PORT-EK find_k (default 5)",
    type=int,
    default=5

)
parser.add_argument(
    "--max_k",
    help="Maximum k value to test with PORT-EK find_k. PORT-EK will test all odd k values from min_k up to and incliding max_k.",
    type=int
)

parser.add_argument(
    "--k",
    help="k value for PORT-EK enriched, map and classify",
    type=int
)

parser.add_argument(
    "--c",
    help="Rarity filter threshold for PORT-EK enriched. Portek will only use k-mers that have frquency of more than c in at least one sample group",
    type=float
)

parser.add_argument(
    "--min_rmse",
    help="RMSE filter threshold for PORT-ek enriched. Portek will discard k-mers with RMSE lower than min_rmse",
    type=float
)

parser.add_argument(
    "--rare_m",
    help="Allowed mismatches when re-examining rare k-mers with PORT-ek enriched. Omit if not re-examining rare k-mers.",
    type=int
)

parser.add_argument(
    "--verbose",
    help="Recieve additional output from some PORT-EK tools",
    default=False,
    action="store_true"
)

parser.add_argument(
    "--n_jobs",
    help="Number of processes used in PORT-EK find and PORT-EK enriched (default 4)",
    default=4,
    type=int
)

def _new_project(project_dir: str) -> None:
    if os.path.isdir(project_dir) == True:
        raise FileExistsError("This project directory already exists! PORT-EK does not allow overwriting projects. If you REALLY want to overwrite, remove the existing directory manually!")
    os.mkdir(project_dir)
    os.mkdir(f"{project_dir}/input")
    os.mkdir(f"{project_dir}/output")
    os.mkdir(f"{project_dir}/temp")
    shutil.copy2("./templates/config.yaml", project_dir)
    print(f"New empty PORT-EK project created in {project_dir}. Please edit config.yaml and generte_kmers.sh files as required and copy input fasta files into {project_dir}/input.")

def main():
    args = parser.parse_args()
    if type(args.project_dir) != str:
        raise ValueError("Please provide a valid project directory name.")
    
    if args.tool == "new":
        _new_project(args.project_dir)

    elif args.tool == "find_k":
        start_time = datetime.now()
        if args.skip_kf == False:
            kmer_finder = portek.KmerFinder(args.project_dir, args.min_k, args.max_k)
            kmer_finder.find_all_kmers(args.n_jobs, args.verbose)
        optimal_k_finder = portek.FindOptimalKPipeline(args.project_dir, args.min_k, args.max_k)
        optimal_k_finder.find_optimal_k(args.n_jobs, args.verbose)
        end_timeS_ARE_NOT_CANON = datetime.now()
        running_time = (end_timeS_ARE_NOT_CANON-start_time)
        print(f"Total running time: {running_time}")

    elif args.tool == "enriched":
        start_time = datetime.now()
        enriched_kmers_finder = portek.EnrichedKmersPipeline(args.project_dir, args.k, args.c, args.min_rmse)
        if args.rare_m == None:
            enriched_kmers_finder.get_kmers()
            enriched_kmers_finder.calc_kmer_stats("common")
            enriched_kmers_finder.plot_volcanos("common")       
            enriched_kmers_finder.get_enriched_kmers()
            enriched_kmers_finder.save_matrix("enriched")
        else:
            enriched_kmers_finder.get_kmers(save_rare=True)
            enriched_kmers_finder.calc_kmer_stats("common")
            enriched_kmers_finder.reexamine_rare(args.rare_m, args.n_jobs)
            enriched_kmers_finder.plot_volcanos("common")   
            enriched_kmers_finder.plot_volcanos("rare_similar")
            enriched_kmers_finder.get_enriched_kmers()
            enriched_kmers_finder.save_matrix("enriched")   

        end_timeS_ARE_NOT_CANON = datetime.now()
        running_time = (end_timeS_ARE_NOT_CANON-start_time)
        print(f"Total running time: {running_time}")


    elif args.tool == "map":
        pass

    elif args.tool == "classify":
        pass

    else:
        raise ValueError("Unrecoginzed PORT-EK tool requested. Choose one of: new, find_k, enriched, map, classify.")


if __name__ == "__main__":
    main()

