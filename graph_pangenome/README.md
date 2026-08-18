# Pangenome graph construction:


Important resource note
-----------------------
A normal Bash script cannot reserve RAM from a cluster scheduler. Therefore:

- `--threads` is passed directly to the graph software.
- For Minigraph-Cactus, `--memory` is passed to Cactus/Toil as `--maxMemory`.
- For Minigraph, PGGB and VG, the programs do not provide one universal
  whole-pipeline `--memory` switch. Their `--memory` argument is retained as a
  documented memory target for the current compute node / interactive job and
  is printed in the log; it does not reserve RAM by itself.

For an HPC, obtain a node with enough memory first, then run these scripts.

Only PGGB is manually partitioned
---------------------------------
- Minigraph: whole-genome direct construction.
- Minigraph-Cactus: one `cactus-pangenome` command. Cactus keeps its own default
  internal reference-chromosome splitting; this script does not manually split input.
- PGGB: `partition-before-pggb` first, then each generated community/chromosome
  command is run sequentially.
- VG: whole-genome `vg construct` directly from reference + VCF; no chromosome split.

===============================================================================
1. Minigraph
===============================================================================

Prepare `assemblies.list`, one whole-genome assembly FASTA path per line:

sample1.fa
sample2.fa
sample3.fa

Run:

bash 01_minigraph_direct.sh \
    reference.fa assemblies.list graph.gfa \
    --threads 24 --memory 300G

===============================================================================
2. Minigraph-Cactus
===============================================================================

Prepare a two-column seqfile:

REF        /path/reference.fa
SampleA.1  /path/SampleA.hap1.fa
SampleA.2  /path/SampleA.hap2.fa
SampleB.1  /path/SampleB.hap1.fa
SampleB.2  /path/SampleB.hap2.fa

Run a large mammalian example:

bash 02_cactus_pangenome_direct.sh \
    seqfile.txt REF cattle_pg cattle_pg \
    --threads 64 --memory 850G \
    --workdir /path/to/large_local_or_shared_scratch

Outputs requested by default:
- GFA
- GBZ
- VCF

Giraffe indexes are optional because they add substantial indexing work:

bash 02_cactus_pangenome_direct.sh \
    seqfile.txt REF cattle_pg cattle_pg \
    --threads 64 --memory 850G \
    --workdir /path/to/scratch \
    --giraffe

Optional stage-specific overrides can be supplied through environment variables:

export CACTUS_MG_CORES=32
export CACTUS_MAP_CORES=8
export CACTUS_CONS_CORES=32
export CACTUS_INDEX_CORES=32

If a specific stage is known to underestimate memory, optionally set:

export CACTUS_MG_MEMORY=850G
export CACTUS_CONS_MEMORY=850G
export CACTUS_INDEX_MEMORY=850G

Do not set those stage-memory values unless needed.

===============================================================================
3. PGGB
===============================================================================

For whole-genome mammalian assemblies, this script partitions first and runs the
generated community/chromosome PGGB commands sequentially.

PanSN-style sequence names are strongly recommended.

Example:

bash 03_pggb_partition_direct.sh \
    all_assemblies.fa.gz pggb_out \
    --threads 32 --poa-threads 16 --memory 256G \
    --identity 98 --segment 10k --min-match 47 \
    --haps 61 --vcf-ref REF

If PanSN naming is correct, `--haps` may be omitted.

If no VCF is required, omit `--vcf-ref`.

The script creates:
- `partition.stdout.log`
- `pggb_commands.sh`
- one `.community.*.fa[.gz]` per partition
- one `.community.*.fa[.gz].out/` directory per PGGB build

Community builds are run one by one to avoid multiplying peak memory usage.

===============================================================================
4. VG
===============================================================================

Whole-genome direct construction, no chromosome split:

bash 04_vg_from_vcf_direct.sh \
    reference.fa variants.vcf.gz graph.vg \
    --threads 32 --memory 256G --gfa

Structural-variant handling is enabled by default (`vg construct -S`).
Disable only when the input VCF contains no structural variants:

--no-handle-sv

The VCF must be bgzip-compressed. The script creates a tabix index if missing.
