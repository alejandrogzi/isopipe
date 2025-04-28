#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Author: Bogdan Kirilenko, 2023
// Contributor: Alejandro Gonzales-Irribarren, 2024
// Nextflow procedure to execute chain feature extraction jobs
// Joblist contains a file where each line is a separate command
// We just call these lines in parallel

// params section: basically command line arguments
params.jobs = 'NONE'  // file containing jobs
params.mem = '10.GB' // memory for each job
params.queue = 'batch' // queue name
params.time = '24.h' // time for each job
params.threads = 1 // number of threads for each job

// if still default -> nothing assigned: show usage message and quit
if (params.jobs == "NONE"){
    println("Usage: nextflow execute_jobs.nf  --jobs [jobs file] -c [config file]")
    System.exit(2);
}

// create channel lines -> we need to execute lines in parallel
Channel
    .fromPath(params.jobs)
    .splitText()
    .set{ lines }

process execute_jobs {
    executor = 'slurm'
    cpus = params.threads


    // allow each process to fail 3 times
    errorStrategy = {task.attempt <= 3 ? 'retry' : 'ignore'}
    memory = MemoryUnit.of(params.memory)

    input:
    val line

    // one line represents an independent command
    script:
    """
    ${line}
    """
}

workflow {
    execute_jobs(lines)
}
