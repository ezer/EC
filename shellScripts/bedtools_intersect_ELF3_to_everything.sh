#bedtools intersect -loj -a phyb_diffpeak_c3.0_cond2.bed  -b ELF3pos.bed >phyb_diffpeak_c3.0_cond2_vs_ELF3_loj_intersect.bed
#bedtools intersect -loj -a phyb_diffpeak_c3.0_cond2.bed  -b ELF4pos.bed >phyb_diffpeak_c3.0_cond2_vs_ELF4_loj_intersect.bed
#bedtools intersect -loj -a phyb_diffpeak_c3.0_cond2.bed  -b LUXpos.bed >phyb_diffpeak_c3.0_cond2_vs_LUX_loj_intersect.bed
#bedtools intersect -loj -a phyb_diffpeak_c3.0_cond2.bed  -b PIF5pos.bed >phyb_diffpeak_c3.0_cond2_vs_PIF5_loj_intersect.bed




bedtools slop -i ELF3pos2.bed -g /home/Reference_genomes/Arabidopsis_thaliana_Ensembl_TAIR10/Ensembl/TAIR10/GenomeStudio/Arabidopsis_thaliana/Ensembl-TAIR10/ChromInfo.txt  -b 100 >ELF3_slop100bp.bed

#other EC components
bedtools intersect -loj -a ELF3pos2.bed -b LUXpos2.bed >ELF3_bLUX_loj_intersect.bed
bedtools intersect -loj -a ELF3_slop100bp.bed -b LUXpos2.bed >ELF3_bLUX_slop100_loj_intersect.bed

bedtools intersect -loj -a ELF3pos2.bed -b ELF4pos2.bed >ELF3_bELF4_loj_intersect.bed
bedtools intersect -loj -a ELF3_slop100bp.bed -b ELF4pos2.bed >ELF3_bELF4_slop100_loj_intersect.bed

#LUX17C-fc5.narrowPeak.curated.bed
bedtools intersect -loj -a ELF3pos2.bed -b LUX17C-fc5.narrowPeak.curated.bed >ELF3_bLUX17_loj_intersect.bed
bedtools intersect -loj -a ELF3_slop100bp.bed -b LUX17C-fc5.narrowPeak.curated.bed >ELF3_bLUX17_slop100_loj_intersect.bed


#phyB peaks
bedtools intersect -loj -a ELF3pos2.bed -b phyB_17C.bed >ELF3_bphyB_17_loj_intersect.bed
bedtools intersect -loj -a ELF3_slop100bp.bed -b phyB_17C.bed >ELF3_bphyB_17_slop100_loj_intersect.bed

bedtools intersect -loj -a ELF3pos2.bed -b phyB_27C.bed >ELF3_bphyB_27_loj_intersect.bed
bedtools intersect -loj -a ELF3_slop100bp.bed -b phyB_27C.bed >ELF3_bphyB_27_slop100_loj_intersect.bed

#each PIF
#bedtools intersect -loj -a ELF3pos.bed -b PIF1pos.bed >ELF3_PIF1_loj_intersect.bed
#bedtools intersect -loj -a ELF3_slop100bp.bed -b PIF1pos.bed >ELF3_PIF1_slop100_loj_intersect.bed

#bedtools intersect -loj -a ELF3pos.bed -b PIF3pos.bed >ELF3_PIF3_loj_intersect.bed
#bedtools intersect -loj -a ELF3_slop100bp.bed -b PIF3pos.bed >ELF3_PIF3_slop100_loj_intersect.bed

#bedtools intersect -loj -a ELF3pos.bed -b PIF4pos.bed >ELF3_PIF4_loj_intersect.bed
#bedtools intersect -loj -a ELF3_slop100bp.bed -b PIF4pos.bed >ELF3_PIF4_slop100_loj_intersect.bed

#bedtools intersect -loj -a ELF3pos.bed -b PIF5pos.bed >ELF3_PIF5_loj_intersect.bed
#bedtools intersect -loj -a ELF3_slop100bp.bed -b PIF5pos.bed >ELF3_PIF5_slop100_loj_intersect.bed

#bedtools intersect -loj -a phyB_17C_slop_b_100bp.bed  -b phyB_27C_slop_b_100bp.bed >phyb_a17C_b27C_slop_b_100bp_loj_intersect.bed
#bedtools intersect -loj -a phyB_27C_slop_b_100bp.bed  -b phyB_17C_slop_b_100bp.bed >phyb_a27C_b17C_slop_b_100bp_loj_intersect.bed

#bedtools intersect -wa -wb \
# -a phyb_diffpeak_c3.0_cond2.bed \
# -b "ELF3pos.bed" "ELF4pos.bed" "LUXpos.bed" "PIF5pos.bed" \
# -names "ELF3" "ELF4" "LUX" "PIF5" >phyb_diffpeak_c3.0_cond2_vs_ELF3_3_4_5_wb_intersect.bed

#bedtools intersect -loj -a phyb_diffpeak_c3.0_cond2.bed -b ELF3pos.bed ELF4pos.bed LUXpos.bed PIF5pos.bed -names ELF3 ELF4 LUX PIF5 >phyb_diffpeak_c3.0_cond2_vs_ELF3_3_4_5_loj_intersect.bed
