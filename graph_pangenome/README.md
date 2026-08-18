# Pangenome Graph Construction

Four workflows for pangenome graph construction:

1. **Minigraph** — reference-guided incremental graph construction
2. **Minigraph-Cactus** — whole-genome pangenome alignment and graph construction
3. **PGGB** — reference-free all-to-all graph construction
4. **VG** — graph construction from a reference genome and an existing VCF

All scripts are run directly with `bash`; no `sbatch` wrapper is included.

> **Resource note**
>
> A Bash script itself cannot reserve memory from a cluster scheduler.
> `--threads` controls software threads. For Minigraph-Cactus, `--memory` is passed
> to Cactus/Toil as `--maxMemory`. For Minigraph, PGGB and VG, `--memory` records the
> recommended/expected memory capacity of the compute node but does not reserve RAM.
> Before running on HPC, enter or request a compute node with sufficient memory.

## Workflow overview

- **Minigraph:** whole-genome direct construction.
- **Minigraph-Cactus:** one `cactus-pangenome` workflow. Cactus retains its own internal graph splitting; the script does not manually split assemblies.
- **PGGB:** `partition-before-pggb` first, followed by sequential PGGB construction for each generated community/chromosome.
- **VG:** whole-genome `vg construct` from reference + VCF; no manual chromosome splitting.

---

## 1. Minigraph

**Script:** `01_minigraph.sh`

### Input

Prepare `assemblies.list`, with one whole-genome assembly FASTA path per line:

```text
/path/to/sample1.fa
/path/to/sample2.fa
/path/to/sample3.fa
```

The reference genome is supplied separately and should not be repeated in `assemblies.list`.

### Run

```bash
bash 01_minigraph.sh \
    reference.fa \
    assemblies.list \
    graph.gfa \
    --threads 24 \
    --memory 300G
```

### Main output

```text
graph.gfa
```

---

## 2. Minigraph-Cactus

**Script:** `02_cactus_pangenome.sh`

### Input

Prepare a two-column `seqfile.txt`:

```text
REF        /path/reference.fa
SampleA.1  /path/SampleA.hap1.fa
SampleA.2  /path/SampleA.hap2.fa
SampleB.1  /path/SampleB.hap1.fa
SampleB.2  /path/SampleB.hap2.fa
```

The reference name supplied to the script must exactly match the reference sample name in the first column.

### Run

Example for a large mammalian pangenome:

```bash
bash 02_cactus_pangenome.sh \
    seqfile.txt \
    REF \
    cattle_pg \
    cattle_pg \
    --threads 64 \
    --memory 850G \
    --workdir /path/to/large_scratch
```

### Main outputs

By default, the workflow requests:

```text
GFA
GBZ
VCF
```

To additionally generate Giraffe indexes:

```bash
bash 02_cactus_pangenome.sh \
    seqfile.txt \
    REF \
    cattle_pg \
    cattle_pg \
    --threads 64 \
    --memory 850G \
    --workdir /path/to/large_scratch \
    --giraffe
```

### Optional Cactus resource overrides

Stage-specific thread settings:

```bash
export CACTUS_MG_CORES=32
export CACTUS_MAP_CORES=8
export CACTUS_CONS_CORES=32
export CACTUS_INDEX_CORES=32
```

If a particular stage is known to underestimate memory, optional overrides can be used:

```bash
export CACTUS_MG_MEMORY=850G
export CACTUS_CONS_MEMORY=850G
export CACTUS_INDEX_MEMORY=850G
```

These stage-specific memory overrides should only be set when needed.

---

## 3. PGGB

**Script:** `03_pggb.sh`

For large mammalian whole-genome assemblies, PGGB is **partitioned first** rather than running all chromosomes together.

The workflow is:

```text
whole-genome assemblies
        ↓
partition-before-pggb
        ↓
community / chromosome FASTA files
        ↓
PGGB construction for each community
        ↓
per-community final GFA graphs
```

PanSN-style sequence names are strongly recommended.

### Run

```bash
bash 03_pggb.sh \
    all_assemblies.fa.gz \
    pggb_out \
    --threads 32 \
    --poa-threads 16 \
    --memory 256G \
    --identity 98 \
    --segment 10k \
    --min-match 47 \
    --haps 61 \
    --vcf-ref REF
```

If PanSN naming is complete and PGGB can infer haplotypes, `--haps` can be omitted.

If VCF output is not required, omit:

```text
--vcf-ref REF
```

### Main outputs

The output directory contains:

```text
partition.stdout.log
pggb_commands.sh
*.community.*.fa
*.community.*.fa.gz
*.community.*.fa.out/
*.final.gfa
```

Community/chromosome jobs are executed **sequentially** in this direct-run script to avoid multiplying peak memory usage.

---

## 4. VG

**Script:** `04_vg.sh`

VG is run directly on the whole reference genome and VCF; this workflow does **not** manually split chromosomes.

### Input

```text
reference.fa
variants.vcf.gz
```

The VCF must be bgzip-compressed. If the tabix index is absent, the script creates it automatically.

### Run

```bash
bash 04_vg.sh \
    reference.fa \
    variants.vcf.gz \
    graph.vg \
    --threads 32 \
    --memory 256G \
    --gfa
```

Structural-variant handling is enabled by default with:

```text
vg construct -S
```

If the VCF contains only SNPs/indels and SV handling is not required:

```bash
bash 04_vg.sh \
    reference.fa \
    variants.vcf.gz \
    graph.vg \
    --threads 32 \
    --memory 256G \
    --no-handle-sv
```

### Main outputs

```text
graph.vg
graph.gfa    # when --gfa is specified
```

---

## Recommended directory layout

```text
graph_pangenome/
├── 01_minigraph.sh
├── 02_cactus_pangenome.sh
├── 03_pggb_partition.sh
├── 04_vg.sh
└── README.md
```
