nextflow.enable.dsl=2

params.top_n = 300
params.prevalence = 0.05

process PREPARE_DATA {
    output:
    path "prepare_data.done"

    script:
    """
    cd ${baseDir}
    export PYTHONPATH=src
    python src/00_prepare_data.py --top-n ${params.top_n} --prevalence ${params.prevalence}
    touch ${task.workDir}/prepare_data.done
    """
}

process DIFFERENTIAL_TAXA {
    input:
    path ready

    output:
    path "differential_taxa.done"

    script:
    """
    cd ${baseDir}
    export PYTHONPATH=src
    python src/01_differential_taxa.py
    touch ${task.workDir}/differential_taxa.done
    """
}

process TRAIN_ML {
    input:
    path ready

    output:
    path "train_ml.done"

    script:
    """
    cd ${baseDir}
    export PYTHONPATH=src
    python src/02_train_ml.py
    touch ${task.workDir}/train_ml.done
    """
}

process BIOLOGICAL_STORY {
    input:
    path diff_done
    path ml_done

    output:
    path "biological_story.done"

    script:
    """
    cd ${baseDir}
    export PYTHONPATH=src
    python src/03_gut_brain_story.py
    touch ${task.workDir}/biological_story.done
    """
}

workflow {
    prepared = PREPARE_DATA()
    diff = DIFFERENTIAL_TAXA(prepared)
    ml = TRAIN_ML(prepared)
    BIOLOGICAL_STORY(diff, ml)
}
