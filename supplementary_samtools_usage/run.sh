#!/bin/bash


######################################################################
## 1. Testing aln.sam file type
######################################################################
echo "[ Type of file:aln.sam ]" >&2 
file aln.sam


######################################################################
## 2 .sam to bam file conversion (bam files compress more efficiently)
######################################################################
qual=20
echo "[ sam => bam (minimum mapping quality $qual) ]" >&2 
samtools view -q $qual -bS -o aln.bam aln.sam


######################################################################
## 3. Testing aln.bam file type
######################################################################
echo "[ Type of file:aln.bam ]" >&2 
file aln.bam


######################################################################
## 4. Printing header information
######################################################################
echo "[ Printing header of bam file only ]" >&2
samtools view -H aln.bam


######################################################################
## 5. Example extracting column 2 from body of sam file
######################################################################
echo "[ Examining sam flag column using perl (0 indexed )]" >&2 
samtools view -S aln.sam  | perl -lane 'print @F[1]' | head -n 10
echo "[ Examining sam flag column using awk (1 indexed )]" >&2 
samtools view -S aln.sam  | awk '{print $2}' | head -n 10 


######################################################################
## 6. Example extracting column 2 from body of bam file
######################################################################
echo "[ Examining bam flag column using perl (0 indexed )]" >&2 
samtools view aln.bam  | perl -lane 'print @F[1]' | head -n 10
echo "[ Examining bam flag column using awk (1 indexed )]" >&2 
samtools view aln.sam  | awk '{print $2}' | head -n 10


######################################################################
## 7. Code snippet to explore other columns
######################################################################
for i in {1..11}; 
do
  echo "[ Command to look at Column $i ]" >&2 
  echo "samtools view -S aln.sam  | awk '{print \$$i}' | head -n 10" >&2 
done


######################################################################
## 8. Sorting the bam file
######################################################################
echo "[ sorting bam file ]" >&2 
samtools sort aln.bam aln_sorted


######################################################################
## 9. There currently is no index file
######################################################################
echo "[ list files ending in '.bai' extension ]" >&2 
ls -l | grep "\.bai$"


######################################################################
## 10. Create an index 
######################################################################
echo "[ create an index file ]" >&2
samtools index aln_sorted.bam


######################################################################
## 11. index file created
######################################################################
echo "[ list files ending in '.bai' extension]" >&2
ls -l | grep "\.bai$"


######################################################################
## 12. alignment statistics (idxstats)
######################################################################
echo "[ bam index stats on a sorted bam file (idxstats) ]" >&2
echo "refname                 length  mapped" >&2
samtools idxstats aln_sorted.bam


######################################################################
## 13. alignment statistics (idxstats)
######################################################################
echo "[ simple stats on a bam file (flagstat) ]" >&2
samtools flagstat aln_sorted.bam


######################################################################
## 14. Counting matches in a bam file
######################################################################
echo "[ Count matches in aln_sored.bam ]" >&2
samtools view -c aln_sorted.bam


######################################################################
## 15. samtools text alignment viewer (requires sorted indexed bam file)
######################################################################
echo "[ Running text alignment viewer in 5 seconds, type ? in samtools tview for help ]" >&2
sleep 5
samtools tview aln_sorted.bam ../05.reference/pool1.fasta


######################################################################
## 16. bam to sam conversion (some tools require sam format)
######################################################################
echo "[ bam => sam ]" >&2
samtools view -o aln_sorted.sam aln_sorted.bam


######################################################################
## 17. merging several bam files (454 example)
######################################################################
echo "[ merging multiple bam files into one bam file ]" >&2
samtools merge pooled1.bam ../11.Samtools_Filter_Pool1/Pool1_Pop[1-8].filtered.bam 


######################################################################
## 18. Filter on specific reference region
######################################################################
echo "[ Printing header of bam file only ]" >&2
samtools view -H aln.bam


######################################################################
## 19. Filter on reference CO_Pool1_contig00004, region
######################################################################
echo "[ Filter on reference:CO_Pool1_contig00004, region 500-1000 ]" >&2
samtools view aln_sorted.bam CO_Pool1_contig00004:500-1000 | head -n 10


######################################################################
## 20. How many records are in the region 500-1000?
######################################################################
echo "[ How many records in region? ]" >&2
samtools view aln_sorted.bam CO_Pool1_contig00004:500-1000 | wc -l


######################################################################
## 21. Piping samtools commands (doesnt store intermediate files)
######################################################################
echo "[ SAM directly to a sorted BAM file ]" >&2
samtools view -bS aln.sam | samtools sort - aln_sorted


######################################################################
## 22. Multiple piping
######################################################################
echo "[ SAM directly to a sorted BAM file with no duplicates ]" >&2

#samtools view -bS aln.sam | samtools sort - aln_sorted
#samtools rmdup -S aln_sorted.bam aln_sorted_unique.bam

#samtools view -bS aln.sam | 
#samtools sort - - |
#samtools rmdup -S - aln_sorted_rmdup.bam

