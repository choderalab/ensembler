import os
import Bio.SeqIO
import pandas as pd
import numpy as np
import argparse

col_names = ['template_id', 'model', 'implicit_refined', 'explicit_refinement']


def main():
    args = parse_args()
    targets_dir = os.path.abspath("targets")
    models_dir = os.path.abspath("models")

    targets_fasta_filename = os.path.join(targets_dir, 'targets.fa')
    targets = list( Bio.SeqIO.parse(targets_fasta_filename, 'fasta') )

    print '%15s %10s %20s %20s' % ('template_id', 'model', 'implicit_refined', 'explicit_refined')

    for target in targets:
        if args.targets is not None and target.id not in args.targets:
            continue

        target_models_dir = os.path.join('models', target.id)
        if not os.path.exists(target_models_dir):
            continue

        target_counts = analyze_target(target_models_dir)
        print '%15s %10d %20d %20d' % (target.id, target_counts[0], target_counts[1], target_counts[2])


def parse_args():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--targets', nargs='+')
    args = argparser.parse_args()
    return args


def analyze_target(target_models_dir):
    contents = os.listdir(target_models_dir)

    raw_data = []
    for relpath in contents:
        path = os.path.join(target_models_dir, relpath)
        if os.path.isdir(path):
            analysis = analyze_model_dir(path)
            raw_data.append(analysis)

    raw_array = np.array(raw_data)
    counts = [col.sum() for col in raw_array.T]
    return counts


def analyze_model_dir(dirpath):
    model_path = os.path.join(dirpath, 'model.pdb.gz')
    implicit_refined_path = os.path.join(dirpath, 'implicit-refined.pdb.gz')
    explicit_refined_path = os.path.join(dirpath, 'explicit-refined.pdb.gz')
    results = [
        os.path.exists(model_path),
        os.path.exists(implicit_refined_path),
        os.path.exists(explicit_refined_path),
    ]
    return results


if __name__ == '__main__':
    main()
