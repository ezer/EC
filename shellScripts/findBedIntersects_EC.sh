#first, merge and sort the two replicates for LUX17
cat LUX_ZT10_17C_rep.bed LUX_ZT10_17C.bed >LUX_ZT10_17C_merged.bed

#then sort these
sortBed -i LUX_ZT10_17C_merged.bed >LUX_ZT10_17C_merged_sort.bed

#then merge it
bedtools merge -i LUX_ZT10_17C_merged_sort.bed >LUX_ZT10_17C_merged_sort_mergeSegments.bed

####Now I need a huge file, but first I need to merge the 3 other EC components
cat LUX_ZT10_22C.bed ELF3_ZT10_22C.bed ELF4_ZT10_22C.bed >EC_22C.bed
sortBed -i EC_22C.bed >EC_22C_sorted.bed
bedtools merge -i EC_22C_sorted.bed > EC_22C_sorted_mergeSegments.bed

cat EC_22C_sorted_mergeSegments.bed LUX_ZT10_17C_merged_sort_mergeSegments.bed >EC_ALL.bed
sortBed -i EC_ALL.bed >EC_ALL_sort.bed
bedtools merge -i EC_ALL_sort.bed >EC_ALL_sort_mergeSegments.bed

#Now I need to intersect:
bedtools intersect -loj -names LUX17 LUX22 ELF4 ELF3 -a EC_ALL_sort_mergeSegments.bed -b LUX_ZT10_17C_merged_sort_mergeSegments.bed LUX_ZT10_22C.bed ELF4_ZT10_22C.bed ELF3_ZT10_22C.bed >EC_ALL_intersect_TP3_ZT10.bed
