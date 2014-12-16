import os
import Bio.SeqIO
import pandas as pd
import numpy as np
import argparse
import yaml
import datetime
import ensembler

# col_names = ['template_id', 'model', 'implicit_refined', 'explicit_refinement']

ensembler_function_types = ['build_models', 'refine_implicit_md', 'refine_explicit_md']

metadata_filenames = {
    'build_models': 'modeling-log.yaml',
    'refine_implicit_md': 'implicit-log.yaml',
    'refine_explicit_md': 'explicit-log.yaml',
}


def main():
    args, process_only_these_targets = parse_args()
    targets_dir = os.path.abspath("targets")
    models_dir = os.path.abspath("models")

    targets_fasta_filename = os.path.join(targets_dir, 'targets.fa')
    targets = list( Bio.SeqIO.parse(targets_fasta_filename, 'fasta') )

    print_col_headers(args)

    for target in targets:
        if process_only_these_targets is not None and target.id not in process_only_these_targets:
            continue

        target_models_dir = os.path.join('models', target.id)
        if not os.path.exists(target_models_dir):
            continue

        target_data = analyze_target(target_models_dir, args)
        summarized_target_data = summarize_target_data(target_data, args)
        summarized_target_data['id'] = target.id
        #print summarized_target_data
        print_target_data(summarized_target_data, args)


def parse_args():
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        '--targets', nargs='+', help='(Default: all targets) Optionally define a subset of targets to work on by providing one or more target IDs separated by spaces (e.g. "ABL1_HUMAN_D0")'
    )
    argparser.add_argument(
        '--targetsfile', type=str, help='Optionally define a filename containing a list of newline-separated target IDs. Comment targets out with "#".'
    )
    argparser.add_argument('--timings', action="store_true")
    args = argparser.parse_args()

    if args.targetsfile != None:
        with open(args.targetsfile, 'r') as targetsfile:
            process_only_these_targets = [line.strip() for line in targetsfile.readlines() if line[0] != '#']
    elif args.targets != None:
        process_only_these_targets = args.targets
    else:
        process_only_these_targets = None

    return args, process_only_these_targets


def print_col_headers(args):
    if args.timings:
        print '%15s %10s %20s %20s %20s %20s %20s' % ('template_id', 'model', 'implicit_refined', 'explicit_refined', 'modeling_timing', 'implicit_timing', 'explicit_timing')
    else:
        print '%15s %10s %20s %20s' % ('template_id', 'model', 'implicit_refined', 'explicit_refined')


def print_target_data(target_data, args):
    if args.timings:
        #print target_data.get('counts'), target_data.get('timings')
        print '%15s %10d %20d %20d %20s %20s %20s' % (target_data['id'], target_data['counts'][0], target_data['counts'][1], target_data['counts'][2], target_data['timings'][0], target_data['timings'][1], target_data['timings'][2])
    else:
        print '%15s %10d %20d %20d' % (target_data['id'], target_data['counts'][0], target_data['counts'][1], target_data['counts'][2])


def analyze_target(target_models_dir, args):
    contents = os.listdir(target_models_dir)

    target_results = {
        'data': [],
        'timings': [],
    }
    for relpath in contents:
        path = os.path.join(target_models_dir, relpath)
        if os.path.isdir(path):
            template_id = relpath
            target_results['data'].append(get_model_existances(path))
            if args.timings:
                target_results['timings'].append(get_timings(path))

    return target_results


def summarize_target_data(target_data, args):
    summarized_target_results = {}
    summarized_target_results['counts'] = calculate_target_model_counts(target_data['data'])
    if args.timings:
        summarized_target_results['timings'] = calculate_target_model_timings(target_data['timings'])
    
    return summarized_target_results


def calculate_target_model_counts(model_data):
    results = [col.sum() for col in np.array(model_data).T]
    return results


def calculate_target_model_timings(model_timings):
    timings_df = pd.DataFrame(np.array(model_timings))
    means = [timings_df[colname].mean() for colname in timings_df]
    stdevs = [timings_df[colname].std() for colname in timings_df]
    results_strs = []
    for i in range(3):
        if means[i] is np.nan:
            results_strs.append('-')
        else:
            result_str = '%s +/- %.0fs' % (seconds2timing_str(means[i]), stdevs[i])
            results_strs.append(result_str)
    return results_strs


def get_model_existances(dirpath):
    model_path = os.path.join(dirpath, 'model.pdb.gz')
    implicit_refined_path = os.path.join(dirpath, 'implicit-refined.pdb.gz')
    explicit_refined_path = os.path.join(dirpath, 'explicit-refined.pdb.gz')
    results = [
        os.path.exists(model_path),
        os.path.exists(implicit_refined_path),
        os.path.exists(explicit_refined_path),
    ]
    return results


def get_timings(dirpath):
    results = []

    for function_type in ensembler_function_types:
        metadata_filename = metadata_filenames[function_type]
        metadata_path = os.path.join(dirpath, metadata_filename)
        if os.path.exists(metadata_path):
            with open(metadata_path) as metadata_file:
                model_metadata = yaml.load(metadata_file, Loader=ensembler.core.YamlLoader)
        else:
            model_metadata = {}

        timing_str = model_metadata.get('timing')

        if timing_str is not None:
            seconds = timing_str2seconds(timing_str)
            results.append(seconds)
        else:
            results.append(np.nan)

    return results


def timing_str2seconds(timing_str):
    h, m, s = timing_str.split(':')
    seconds = float(s) + float(m)*60 + float(h)*3600
    return seconds


def seconds2timing_str(seconds):
    if seconds is np.nan:
        return None
    seconds = int(round(seconds))
    h = seconds / 3600
    remainder = seconds % 3600
    m = remainder / 60
    s = remainder % 60
    timing_str = '%d:%02d:%02d' % (h,m,s)
    return timing_str


if __name__ == '__main__':
    main()
