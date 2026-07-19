#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS/WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ISOPIPE as MAIN } from './workflows/isopipe.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATION FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def validateFullRun() {
    def problems = []

    if (!params.global_input_dir) { problems << 'missing required --global_input_dir' }
    if (!params.global_output_dir) { problems << 'missing required --global_output_dir' }
    if (!params.global_genome) { problems << 'missing required --global_genome' }
    if (!params.global_annotation) { problems << 'missing required --global_annotation' }
    if (!params.global_repeats) { problems << 'missing required --global_repeats' }
    
    if (!(params.entrypoint in ['isoseq', 'map'])) { 
      problems << 'unknown entrypoint option -> options are: isoseq, map'
    }

    if (!(params.isoseq_cluster2_mode in ['per_sample', 'multi_sample', 'both'])) { 
      problems << 'unknown isoseq_cluster2_mode option -> options are: per_sample, multi_sample, both'
    }

    if (!(params.aligner in ['minimap2', 'ultra', 'desalt', 'pbmm2'])) { 
      problems << 'unknown aligner option -> options are: minimap2, ultra, desalt, pbmm2'
    }

    if (problems) {
        error "Parameter validation failed:\n  - " + problems.join('\n  - ')
        System.exit(1)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ISOPIPE {  
    validateFullRun()

    log.info """
    isopipe v${workflow.manifest.version}
  
    Authors: ${workflow.manifest.author}
    Github:  ${workflow.manifest.homePage}

      Input     : ${params.global_input_dir}
      Output    : ${params.global_output_dir}
      Genome    : ${params.global_genome}
      Annotation: ${params.global_annotation}
      Repeats   : ${params.global_repeats}
      Aligner   : ${params.aligner}
      Database  : ${params.xorf_protein_database} (custom: ${params.xorf_custom_database})
      Profile   : ${workflow.profile}
    """.stripIndent()

    MAIN () 
}

workflow { ISOPIPE () }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
