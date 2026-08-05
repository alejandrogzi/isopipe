<p align="center">
  <p align="center">
    <img width=100 align="center" src="./assets/figures/logo.png" >
  </p>

<p align="center">
  <picture>
    <source
      media="(prefers-color-scheme: dark)"
      srcset="./assets/figures/hillerlab-dark.png"
    >
    <source
      media="(prefers-color-scheme: light)"
      srcset="./assets/figures/hillerlab-light.png"
    >
    <img
      width="200"
      alt="Hiller Lab"
      src="./assets/figures/hillerlab-light.png"
    >
  </picture>
</p>

  <span>
    <h1 align="center">
        ark
    </h1>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/ark" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/ark?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> A Reference pipeline to annotate euKaryotes at high resolution </span>
        <br>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://github.com/alejandrogzi/ark/blob/master/assets/docs/usage.md">usage</a> .
        <a href="https://github.com/hillerlab/ark/blob/main/assets/pipeline/ark.mermaid">pipeline</a> .
        <a href="https://hillerlab.com/">us</a> 
    </samp>
  </p>

</p>

---

## Usage

> [!NOTE]
> Requirements: Nextflow ≥ 25.04.6, Docker or Apptainer, Java.

```bash
git clone https://github.com/hillerlab/ark.git
cd ark
```

Edit `params.json` (set relevant options/filters), then:
```bash
# Docker
nextflow run main.nf -params-file params.json -profile docker

# Apptainer / Singularity
nextflow run main.nf -params-file params.json -profile apptainer
```

Smoke test:
```bash
nextflow run main.nf -profile test,apptainer
```

> [!NOTE]
> You can also specify these options directly in `params.json`.

A helper sh script is provided to run the pipeline on a SLURM cluster. See details below.

<details>
<summary>Click to expand</summary>


Edit the path variables at the top of `assets/hpc/ark.sh` (cache dir, container image, manifest path), then submit:

```bash
sbatch --array=1-<N> ark.sh
```

Each array task spawns one Nextflow head job that submits all compute as child SLURM jobs.

PREDICT_ORFS run as SLURM job arrays. Partition routing, array sizes, and resource tiers are documented inline in `nextflow.config` — edit there to match your cluster.

</details>

---

## Output

```
results/
├── 00_?/       *bed
└── pipeline_info/    timeline, trace, DAG
```

---

## Where to edit

| File | What |
|------|------|
| `params.json` | Genome paths, alignment settings, checkpoints — per run |
| `nextflow.config` | Compute resources, profiles, container, SLURM — rarely |
