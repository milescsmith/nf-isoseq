# Steps for processing this data:

**Ignore any FASTA or FASTQ files that were delivered**  These are largely useless

I had trouble remembering the steps, programs, and arguments used to
process this data.  Until I can get a Snakemake or Nextflow pipeline written,
this collection of Bash scripts will have to do.

1. Create a conda environment and install the required utilities

```
    conda create -n pb -c bioconda pbcss pbmm2 lima isoseq3 htslib
```
2. (Only if a CCS BAM file is unavailable) Circular Consensus Sequence calling

```
srun \
    --mem=128 \
    --cpus-per-task=16 \
    --partition=serial \
    ccs \
        --min-rq 0.9 \
        --log-level DEBUG \
        --log-file ccs.log \
        --num-threads {THREADS} \
        $PROJECT_FOLDER/01_raw_data/{runid}.subreads.bam $PROJECT_FOLDER/02_CSS/{runid}.ccs.bam
```
** **NOTE** Find a desktop/server/cluster with as many cores as you can get.
Even with ~50 cores, this step will take > 1 d wall time.

3. Demultiplex

You will need a FASTA file with the barcode sequences.  The barcode sequence
names need to end in `_5p` and/or `_3p`, else `lima` (or at least version 
1.10.0) run with the `--isoseq` flag will silently fail.

```
srun \
    --mem=128 \
    --cpus-per-task=16 \
    --partition=serial \
    lima \
        --isoseq \
        --dump-clips \
        --peek-guess \
        --num-threads 36 \
        --log-level INFO \
        --log-file lima.demux.log \
        $PROJECT_FOLDER/02_CSS/{runid}.ccs.bam \
        barcode.primers.fasta \
        $PROJECT_FOLDER/03_demultiplexed/demuxed.bam
```

4. Refine sequences

```
#!/bin/bash

bc_array=({barcode names matching those in the barcode.primers.fasta})

for i in {0..11}
do
    printf "${bc_array[$i]}\n"
    srun \
        --mem=128 \
        --cpus-per-task=16 \
        --partition=serial \
        isoseq3 \
            refine \
            --num-threads 16 \
            --log-file ${bc_array[$i]}.log \
            --log-level INFO \
            --verbose \
            $PROJECT_FOLDER/03_demultiplexed/demuxed.${bc_array[$i]}.bam \
            primers.fasta \
            $PROJECT_FOLDER/04_refined/refined.${bc_array[$i]}.flnc.bam &
done
```

5. Clustering

```
#!/bin/bash

bc_array=("bc1001" "bc1002" "bc1003" "bc1004" "bc1005" "bc1006" "bc1008" "bc1012" "bc1018" "bc1019" "bc1020" "bc1023")

for i in {0..11}
do
    printf "${bc_array[$i]}\n"
    srun \
        --mem=128 \
        --cpus-per-task=16 \
        --partition=serial \
        isoseq3 \
            cluster \
            --use-qvs \
            --num-threads 16 \
            --log-file ${bc_array[$i]}.log \
            --log-level INFO \
            --verbose \
            $PROJECT_FOLDER/04_refined/refined.${bc_array[$i]}.flnc.bam \
            $PROJECT_FOLDER/05_polished/polished.${bc_array[$i]}.bam \
            &
done
```

6. Mapping - use GMAP.  Minimap2 and deSalt are *suppose* to be better, but cDNA-Cupcake doesn't seem to work with their output?

```
#!/bin/bash

id_array=("590085-5-4" "590108-5-4" "541305" "541308" "541424" "541561" "541826" "510099" "550003-6-2" "500066-6-2" "500028-6-4" "541552")
bc_array=("bc1001" "bc1002" "bc1003" "bc1004" "bc1005" "bc1006" "bc1008" "bc1012" "bc1018" "bc1019" "bc1020" "bc1023")

for i in {0..11}
do
    printf "${bc_array[$i]}\n"
    srun \
        --cpus-per-task=8 \
        --mem=128 \
        --partition=highmem \
            /Volumes/guth_aci_informatics/software/gmap/bin/gmap \
                -D /Volumes/guth_aci_informatics/references/genomic/homo_sapiens/indices/gmap/gencode_v32/homo_sapiens/ \
                -d homo_sapiens \
                -f samse \
                -n 0 \
                -t 16 \
                --cross-species \
                --max-intronlength-ends 200000 \
                -z sense_force \
                $PROJECT_FOLDER/05_polished/polished.${bc_array[$i]}.hq.fasta \
                > $PROJECT_FOLDER/06_mapped/${bc_array[$i]}_hq_isoforms.fasta.sam \
                2> $PROJECT_FOLDER/06_mapped/${bc_array[$i]}_hq_isoforms.log
        &
done
```

7. Sort - While the cDNA-Cupcake documents tell you to use UNIX sort, do not do this.  It strips out the BAM header and the next step (filtering) will not work.

```
#!/bin/bash

id_array=("590085-5-4" "590108-5-4" "541305" "541308" "541424" "541561" "541826" "510099" "550003-6-2" "500066-6-2" "500028-6-4" "541552")
bc_array=("bc1001" "bc1002" "bc1003" "bc1004" "bc1005" "bc1006" "bc1008" "bc1012" "bc1018" "bc1019" "bc1020" "bc1023")

for i in {0..11}
do
    printf "${bc_array[$i]}\n"
    srun \
        --mem=128 \
        --cpus-per-task=16 \
        --partition=serial \
        samtools \
            sort \
                -O BAM \
                $PROJECT_FOLDER/06_mapped/${bc_array[$i]}_hq_isoforms.mapped.bam \
                -o $PROJECT_FOLDER/07_sorted/${bc_array[$i]}_hq_isoforms.mapped.sorted.bam
        &
done
```

8. Filter - Not sure why, but there are somehow transcripts in the mapped SAM 
file that do not exist in the pre-mapping FASTA.  Those need to be removed.
Also, we need to transcode the fasta/fastq files from gzip to bgzip.

```
id_array=("590085-5-4" "590108-5-4" "541305" "541308" "541424" "541561" "541826" "510099" "550003-6-2" "500066-6-2" "500028-6-4" "541552")
bc_array=("bc1001" "bc1002" "bc1003" "bc1004" "bc1005" "bc1006" "bc1008" "bc1012" "bc1018" "bc1019" "bc1020" "bc1023")

for i in {0..11}
do
    printf "${bc_array[$i]}\n"
    srun \
        --mem=128 \
        --cpus-per-task=16 \
        --partition=serial \
        bgzip --decompress $PROJECT_FOLDER/05_polished/${bc_array[$i]}.fastq.gz && \
        bgzip --index $PROJECT_FOLDER/05_polished/${bc_array[$i]}.fastq && \
        filter_sam \
            --fastq $PROJECT_FOLDER/05_polished/${bc_array[$i]}.fastq.gz \
            --sam $PROJECT_FOLDER/06_sorted/${bc_array[$i]}_hq_isoforms.mapped.sorted.bam \
            --prefix $PROJECT_FOLDER/07_filtered/${bc_array[$i]}_
        &
done
```

9. Collapse isoforms

```
id_array=("590085-5-4" "590108-5-4" "541305" "541308" "541424" "541561" "541826" "510099" "550003-6-2" "500066-6-2" "500028-6-4" "541552")
bc_array=("bc1001" "bc1002" "bc1003" "bc1004" "bc1005" "bc1006" "bc1008" "bc1012" "bc1018" "bc1019" "bc1020" "bc1023")

for i in {0..11}
do
    printf "${bc_array[$i]}\n"
    srun \
        --mem=128 \
        --cpus-per-task=16 \
        --partition=serial \
        collapse_isoforms_by_sam.py 
            --input $PROJECT_FOLDER/06_polished/${bc_array[$i]}.fastq \
            --sam $PROJECT_FOLDER/07_sorted/${bc_array[$i]}_filtered.sam \
            --dun-merge-5-shorter \
            --prefix $PROJECT_FOLDER/08_collapsed/${bc_array[$i]}
        &
done
```