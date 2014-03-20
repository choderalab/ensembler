def build_models(process_only_these_targets=None, process_only_these_templates=None):
    import os
    import Bio.SeqIO

    targets_dir = os.path.abspath('targets')
    templates_dir = os.path.abspath('templates')
    models_dir = os.path.abspath('models')
    targets_fasta_filename = os.path.join(targets_dir, 'targets.fa')
    templates_fasta_filename = os.path.join(templates_dir, 'templates.fa')

    targets = Bio.SeqIO.parse(targets_fasta_filename, 'fasta')
    templates = Bio.SeqIO.parse(templates_fasta_filename, 'fasta')

    for target in targets:
        if process_only_these_targets and target.id not in process_only_these_targets: continue

        print "========================================================================="
        print "Modelling '%s'" % (target.id)
        print "========================================================================="

        models_target_dir = os.path.join(models_dir, target.id)
        if not os.path.exists(models_target_dir):
            os.mkdir(models_target_dir)

        for template in templates:
            if process_only_these_templates and template.id not in process_only_these_templates: continue
            build_model(target, template)

def build_model(target, template):
    # align target and template
    import Bio.pairwise2
    import Bio.SubsMat.MatrixInfo

    print "-------------------------------------------------------------------------"
    print "Modelling '%s' => '%s'" % (target.id, template.id)
    print "-------------------------------------------------------------------------"

    matrix = Bio.SubsMat.MatrixInfo.gonnet
    gap_open = -10
    gap_extend = -0.5

    aln = Bio.pairwise2.align.globalds(str(target.seq), str(template.seq), matrix, gap_open, gap_extend)

    print aln[0][0]
    print aln[0][1]

