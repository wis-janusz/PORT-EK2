import argparse
import os
import shutil
from datetime import datetime

import portek

parser = argparse.ArgumentParser(
    description="Main PORT-EK v2 program. Use it to run all PORT-EK tools with the appropriate 'tool' argument."
)
parser.add_argument(
    "tool",
    help="Name of the PORT-EK function you want to execute. Choose one of: new, find_k, enriched, us, map, classify",
)

parser.add_argument(
    "project_dir",
    help="Path to the project directory. Must not exist for PORTEKrun.py new, must exist for all other tools",
    type=str,
)

parser.add_argument(
    "--max_k",
    help="Maximum k value to test with PORT-EK find_k. PORT-EK will test all odd k values from min_k up to and including max_k.",
    type=int,
    default=31,
)

parser.add_argument(
    "--k", help="k value for PORT-EK enriched, map and classify", type=int
)

parser.add_argument(
    "--verbose",
    help="Recieve additional output from some PORT-EK tools",
    default=False,
    action="store_true",
)

parser.add_argument(
    "--unsupervised",
    "-us",
    help="Performed unsupervised analysis",
    default=False,
    action="store_true",
)

parser.add_argument(
    "--n_jobs",
    help="Number of processes used in PORT-EK find and PORT-EK enriched (default 4)",
    default=4,
    type=int,
)


def _new_project(project_dir: str) -> None:
    if os.path.isdir(project_dir) == True:
        raise FileExistsError(
            "This project directory already exists! PORT-EK does not allow overwriting projects. If you REALLY want to overwrite, remove the existing directory manually!"
        )
    os.makedirs(project_dir)
    os.makedirs(f"{project_dir}/input")
    os.makedirs(f"{project_dir}/output")
    os.makedirs(f"{project_dir}/temp")
    shutil.copy2("./templates/config.yaml", project_dir)
    print(
        f"New empty PORT-EK project created in {project_dir}. Please edit config.yaml as required and copy input fasta files into {project_dir}/input."
    )


def main():
    args = parser.parse_args()
    if type(args.project_dir) != str:
        raise ValueError("Please provide a valid project directory name.")

    if args.tool == "new":
        _new_project(args.project_dir)

    elif args.tool == "find_k":
        start_time = datetime.now()
        if args.unsupervised == True:
            kmer_finder = portek.UsKmerFinder(args.project_dir, args.max_k)
            times = kmer_finder.find_all_kmers(n_jobs=args.n_jobs, verbose=args.verbose)
            optimal_k_finder = portek.UsFindOptimalKPipeline(
                args.project_dir, args.max_k, times
            )
            optimal_k_finder.find_optimal_k(n_jobs=args.n_jobs, verbose=args.verbose)
        else:
            kmer_finder = portek.KmerFinder(args.project_dir, args.max_k)
            times = kmer_finder.find_all_kmers(n_jobs=args.n_jobs, verbose=args.verbose)
            optimal_k_finder = portek.FindOptimalKPipeline(
                args.project_dir, args.max_k, times
            )
            optimal_k_finder.find_optimal_k(n_jobs=args.n_jobs, verbose=args.verbose)
        end_timeS_ARE_NOT_CANON = datetime.now()
        running_time = end_timeS_ARE_NOT_CANON - start_time
        print(f"\nTotal running time: {running_time}")

    elif args.tool == "us":
        start_time = datetime.now()
        us_kmers_finder = portek.UnsupervisedKmersPipeline(args.project_dir, args.k)
        us_kmers_finder.get_kmers(args.verbose)
        us_kmers_finder.save_matrix("common", full=True)
        us_kmers_finder.save_matrix("clusters", full=True)
        end_timeS_ARE_NOT_CANON = datetime.now()
        running_time = end_timeS_ARE_NOT_CANON - start_time
        print(f"\nTotal running time: {running_time}")

    elif args.tool == "enriched":
        start_time = datetime.now()
        enriched_kmers_finder = portek.EnrichedKmersPipeline(args.project_dir, args.k)
        enriched_kmers_finder.get_basic_kmer_stats()
        enriched_kmers_finder.calc_kmer_stats("common", verbose=args.verbose)
        enriched_kmers_finder.plot_volcanos("common")
        enriched_kmers_found = enriched_kmers_finder.get_enriched_kmers()
        if enriched_kmers_found == True:
            enriched_kmers_finder.save_counts_for_classifier()
            enriched_kmers_finder.save_matrix("common")
            enriched_kmers_finder.save_matrix("enriched")
            portek.save_kmers_fasta(
                enriched_kmers_finder.matrices["enriched"].index.to_list(),
                "enriched",
                enriched_kmers_finder.project_dir,
                enriched_kmers_finder.k,
            )
            enriched_kmers_finder.plot_PCA()
        end_timeS_ARE_NOT_CANON = datetime.now()
        running_time = end_timeS_ARE_NOT_CANON - start_time
        print(f"\nTotal running time: {running_time}")

    elif args.tool == "map":
        start_time = datetime.now()
        mapping_pipeline = portek.MappingPipeline(args.project_dir, args.k)
        mapping_pipeline.run_mapping(verbose = args.verbose)
        mapping_pipeline.analyze_mapping(verbose = args.verbose)
        mapping_pipeline.save_mappings_df()
        if args.rerun_map == True:
            mapping_pipeline.run_mapping(rerun=args.rerun, verbose = args.verbose)
            mapping_pipeline.analyze_mapping(rerun=args.rerun, verbose=args.verbose)
            

    elif args.tool == "classify":
        pass

    else:
        raise ValueError(
            "Unrecoginzed PORT-EK tool requested. Choose one of: new, find_k, enriched, map, classify."
        )


if __name__ == "__main__":
    main()
