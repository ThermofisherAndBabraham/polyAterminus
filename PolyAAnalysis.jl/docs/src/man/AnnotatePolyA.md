# AnnotatePolyA

This module annotates polyA sites if found in reference annotated regions.
If cluster option is selected - clusters terminations sites by >= k from cluster
center (either most frequent or median value) and <= m distance between clusters.

Presumptions:
1. Read name in the BAM file is tagged with a polyA length: `@XXXX:A:120:XXXX`
2. PolyA sequence is trimmed.

## Analysis

We assume that 3' end of the read matching gene sequence (`RIGHT-read`) is a correct polyA site if:

   1. polyA stretch (length is user defined) which was found in the "Right" read was not matching the
      reference sequence.

   2. Mate read is complete polyT stretch with posible 3' end reverse complement sequence matching the
      `RIGHT-read`:
      ```
      RIGHT-read - 5'-ATGCTATGCTAGTCTGATTGCTATTCGAAAAAAAAAAAAAAAAAAAA--------------------  -3'
      MATE-read  - 3'----------------------ATAAGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -5'
      ```

   3. The `RIGHT-read` might be R1 of R2, depending on the stranded library preparation kit used. If R1
   matches gene sequence, then "+" strand should be chosen, if R2 - "-".

## Result output

Results are stored in:

   1. `your_sample_detected_polyA.tsv` - Statistics of each detected polyA site without annotations.

   2. `your_sample_annotated_polyA.bed` - Annotated polyA sites in `bed` format where columns:
      ```
         Chr	Start	End	Name	Counts	Strand	Feature	Median	Min	Max	Biotype
      ```
      If site does not match any features in the reference gff3 file `NA` is returned.
      One gene might have several polyadelynation sites, therefore to separate such site number is
      added to the gene name: `NAME.1, NAME.2 ...`. Same applies to the `NA` values.

   3. `your_sample_detected_clusters_polyA.tsv` - only if cluster option is selected.
      ```
         Chrmosome	Start	End	Strand	ClusterCenter	ClusterSize	ClusterMedian	ClusterMean	ClusterMin	ClusterMax	Cluster1stQuartile	Cluster3rdQuartile	TSMedian	TSMean	TSMin	TSMax	TS1stQuartile	TS3rdQuartile	Smaple	TotalReads	PassedReads	PolyAReads
      ```
   4. `your_sample_annotated_polyA_clusters.bed` - only if cluster option is selected.
      ```
         Chr	Start	End	GeneName	ClusterSize	Strand	Feature	ClusterCenter	Biotype	ClusterMedian	ClusterMean	ClusterMin	ClusterMax	Cluster1stQuartile	Cluster3rdQuartile	TSMedian	TSMean	TSMin	TSMax	TS1stQuartile	TS3rdQuartile
      ```

`BED` files can be easily loaded on `IGV` for reviewing and comparison between
different outputs.
