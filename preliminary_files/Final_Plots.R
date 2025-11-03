# Plotting examples 
# Full Analysis----
# Substrate degradation traits----

Figure_2 = ggarrange(Figure_test_1, Figure_test_2,Figure_test_1.ISO,Figure_test_2.ISO, 
          labels = c("A","B","C","D"),
          ncol = 2, nrow = 2) + theme(panel.background = element_blank())
Figure_2
png("Output_Data/Figures/Figure_2.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(Figure_2)
dev.off()

Figure_2.new = ggarrange(Figure_test_1.new, Figure_test_2.new,Figure_test_1.ISO.new,Figure_test_2.ISO.new, 
                     labels = c("A","B","C","D"),
                     ncol = 2, nrow = 2) + theme(panel.background = element_blank())
Figure_2.new
png("Output_Data/Figures/Figure_2.new.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(Figure_2.new)
dev.off()

pdf("Output_Data/Figures/Figure_2_overleaf.pdf",
    width=4.81*2.75,height=2.21*2.75)
print(Figure_2.new)
dev.off()

# Monomer Uptake----

Figure_3 = ggarrange(Figure_test_3, Figure_test_4,Figure_test_5,Figure_test_3.ISO,
                     Figure_test_4.ISO,Figure_test_5.ISO,
                     labels = c("A","B","C","D","E","F"),
                     ncol = 3, nrow = 2) + theme(panel.background = element_blank())
Figure_3
png("Output_Data/Figures/Figure_3.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(Figure_3)
dev.off()

Figure_3.new = ggarrange(Figure_test_3.new, Figure_test_4.new,Figure_test_5.new,Figure_test_3.ISO.new,
                     Figure_test_4.ISO.new,Figure_test_5.ISO.new,
                     labels = c("A","B","C","D","E","F"),
                     ncol = 3, nrow = 2) + theme(panel.background = element_blank())
Figure_3.new
png("Output_Data/Figures/Figure_3.new.png",
    width=3500*1.5,height=1969*1.5,res=300)
print(Figure_3.new)
dev.off()

pdf("Output_Data/Figures/Figure_3_overleaf.pdf",
    width=4.81*2.75,height=2.21*2.75)
print(Figure_3.new)
dev.off()

# Stress Tolerance----

Figure_4 = ggarrange(Figure_test_6,Figure_test_7,Figure_test_8,Figure_test_9,
                     Figure_test_6.ISO,Figure_test_7.ISO,Figure_test_8.ISO,
                     Figure_test_9.ISO,
                     labels = c("A","B","C","D","E","F","G","H"),
                     ncol = 4, nrow = 2) + theme(panel.background = element_blank())
Figure_4
png("Output_Data/Figures/Figure_4.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(Figure_4)
dev.off()

Figure_4.new = ggarrange(Figure_test_6.new,Figure_test_7.new,Figure_test_8.new,Figure_test_9.new,
                     Figure_test_6.ISO.new,Figure_test_7.ISO.new,Figure_test_8.ISO.new,
                     Figure_test_9.ISO.new,
                     labels = c("A","B","C","D","E","F","G","H"),
                     ncol = 4, nrow = 2) + theme(panel.background = element_blank())
Figure_4.new

png("Output_Data/Figures/Figure_4.new.png",
    width=3500*1.5,height=1969*1.5,res=300)
print(Figure_4.new)
dev.off()

pdf("Output_Data/Figures/Figure_4_overleaf.pdf",
    width=4.81*3,height=2.21*3)
print(Figure_4.new)
dev.off()

# Cross-correlations----

pdf("Output_Data/Figures/Figure_5.pdf",
    width=4.81*3,height=2.21*3)
print(Figure_5)
dev.off()

pdf("Output_Data/Figures/Figure_5.a.pdf",
    width=4.81*3,height=2.21*3)
print(Figure_5.a)
dev.off()

pdf("Output_Data/Figures/Figure_6.pdf",
    width=4.81*3,height=2.21*3)
print(Figure_6)
dev.off()

pdf("Output_Data/Figures/Figure_6.a.pdf",
    width=4.81*3,height=2.21*3)
print(Figure_6.a)
dev.off()

# Minimum Generation Time----

Figure_10 = ggarrange(Figure_test_10a, Figure_test_10b,Figure_test_11, 
                      labels = c("A","B","C"),
                      ncol = 3, nrow = 1) + theme(panel.background = element_blank())
Figure_10
png("Output_Data/Figures/Figure_10.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(Figure_10)
dev.off()

Figure_10.new = ggarrange(Figure_test_10a.new, Figure_test_10.ISO.new,Figure_test_11.ISO.new, 
                      labels = c("A","B","C"),
                      ncol = 3, nrow = 1) + theme(panel.background = element_blank())
Figure_10.new
png("Output_Data/Figures/Figure_10.new.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(Figure_10.new)
dev.off()

pdf("Output_Data/Figures/Figure_10_overleaf.pdf",
    width=4.81*2.75,height=2.21*2.75)
print(Figure_10.new)
dev.off()

# Optimal Growth Temperature----

Figure_10 = ggarrange(Figure_test_10a, Figure_test_10b,Figure_test_11, 
                      labels = c("A","B","C"),
                      ncol = 3, nrow = 1) + theme(panel.background = element_blank())
Figure_10
png("Output_Data/Figures/Figure_10.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(Figure_10)
dev.off()

Figure_11.new = ggarrange(Figure_test_11.new, Figure_test_12.ISO.new, 
                      labels = c("A","B"),
                      ncol = 2, nrow = 1) + theme(panel.background = element_blank())
Figure_11.new
png("Output_Data/Figures/Figure_11.new.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(Figure_11.new)
dev.off()

pdf("Output_Data/Figures/Figure_11_overleaf.pdf",
    width=4.81*2.75,height=2.21*2.75)
print(Figure_11.new)
dev.off()

# CUE plots----

Figure_11 = ggarrange(Figure_test_1.phylum,Figure_test_1.Class,
                      Figure_test_1.Order,Figure_test_1.Family,Figure_test_1.Genus, 
                     labels = c("A","B","C","D","E"),
                     ncol = 3, nrow = 2) + theme(panel.background = element_blank())
Figure_11

pdf("Output_Data/Figures/Figure_7_power_law_new.pdf",
    width=14,height=14*3/5)
print(Figure_11)
dev.off()

Figure_12.new = ggarrange(Figure_test_1.phylum.new,Figure_test_1.Class.new,
                      Figure_test_1.Order.new,Figure_test_1.Family.new,Figure_test_1.Genus.new, 
                      labels = c("A","B","C","D","E"),
                      ncol = 3, nrow = 2) + theme(panel.background = element_blank())
Figure_12.new
png("Output_Data/Figures/Figure_12.new.png",
    width=3500*1.5,height=1969*1.5,res=300)
print(Figure_12.new)
dev.off()

pdf("Output_Data/Figures/Figure_12.new.pdf",
    width=4.81*2.75,height=2.21*2.75)
print(Figure_12.new)
dev.off()

#

Figure_12 = ggarrange(Figure_test_10b, Figure_test_11,Figure_test_1.CUE, 
                      labels = c("A","B","C"),
                      ncol = 3, nrow = 1) + theme(panel.background = element_blank())
Figure_12
png("Output_Data/Figures/Figure_12.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(Figure_12)
dev.off()

# Summary statistics ----

Figure_S5 = ggarrange(Figure_S5.A,Figure_S5.B,
                         labels = c("A","B"),
                         ncol = 2, nrow = 1) + theme(panel.background = element_blank())
Figure_S5

png("Output_Data/Figures/Figure_S5.png",
    width=3500*0.5,height=1969*0.5,res=300)
print(Figure_S5)
dev.off()

pdf("Output_Data/Figures/Figure_S5_overleaf.pdf",
    width=4.81*1.5,height=2.21*1.5)
print(Figure_S5)
dev.off()

# Enzyme encoding genes ----

Figure_2_power_law = ggarrange(Figure_CAZy, Figure_Protein,Figure_CAZy.total.ISO,
                               Figure_Protein.enzyme.total.ISO, 
                               labels = c("A","B","C","D"),
                               ncol = 2, nrow = 2) + theme(panel.background = element_blank())
Figure_2_power_law
pdf("Output_Data/Figures/Figure_2_power_law_new.pdf",
    width=14,height=14*3/5)
print(Figure_2_power_law)
dev.off()

# Transporters encoding genes ----

Figure_3_power_law = ggarrange(Figure_tranport.total, Figure_GH.trasnport,
                               Figure_amino.transport,Figure_tranport.total_ISO,
                               Figure_GH.transporter,Figure_Amino.transporter,
                               labels = c("A","B","C","D","E","F"),
                               ncol = 3, nrow = 2) + theme(panel.background = element_blank())
Figure_3_power_law
pdf("Output_Data/Figures/Figure_3_power_law_new.pdf",
    width=15,height=15*3/5)
print(Figure_3_power_law)
dev.off()

# Stress Tolerance encoding genes ----

Figure_4_power_law = ggarrange(Figure_osmolyte,Figure_biofilm,Figure_temp_fast,
                               Figure_pH_fast,Figure_osmolyte_ISO,Figure_Biofilm_ISO,
                               Figure_Temp.Tol_ISO,Figure_pH.Tol_ISO,
                               labels = c("A","B","C","D","E","F","G","H"),
                               ncol = 4, nrow = 2) + theme(panel.background = element_blank())
Figure_4_power_law

pdf("Output_Data/Figures/Figure_4_power_law_new.pdf",
    width=16,height=16*3/5)
print(Figure_4_power_law)
dev.off()

# Life History encoding genes (Yield) ----

Figure_5_power_law = ggarrange(Figure_mgr_fast, Figure_mgr_ISO, Figure_Yield_ISO, 
                               labels = c("A","B","C"),
                               ncol = 3, nrow = 1) + theme(panel.background = element_blank())

Figure_5_power_law

pdf("Output_Data/Figures/Figure_5_power_law_new.pdf",
    width=12,height=12*3/5)
print(Figure_5_power_law)
dev.off()

# Life History encoding genes (Temperature) ----

Figure_6_power_law = ggarrange(Figure_ogt_fast, Figure_OGT_ISO, 
                               labels = c("A","B"),
                               ncol = 2, nrow = 1) + theme(panel.background = element_blank())

Figure_6_power_law

pdf("Output_Data/Figures/Figure_6_power_law_new.pdf",
    width=12,height=12*3/5)
print(Figure_6_power_law)
dev.off()

# A/S ratio vs genome size ----

Figure_7_A_S = ggarrange(Figure_genome_A_S_MAGs, Figure_genome_A_S_ISO, 
                               labels = c("A","B"),
                               ncol = 2, nrow = 1) + theme(panel.background = element_blank())

Figure_7_A_S

pdf("Output_Data/Figures/Figure_7_A_S.pdf",
    width=12,height=12*3/5)
print(Figure_7_A_S)
dev.off()

# Specific tradeoffs ----

Figure_8_AS = ggarrange(Figure_Total_pH_MAGs,Figure_Total_pH_ISO,
                        Figure_amino_pH_MAGs,Figure_amino_pH_ISO, 
                        Figure_amino_GH_MAGs,Figure_amino_GH_ISO,
                        Figure_protein_pH_MAGs,Figure_protein_pH_ISO,
                        Figure_protein_CAZy_MAGs,Figure_CAZy_pH_ISO,
                               labels = c("A","B","C","D", "E", "F", "G", "H","I", "J"),
                               ncol = 2, nrow = 5) + theme(panel.background = element_blank())
Figure_8_AS
pdf("Output_Data/Figures/Figure_8_AS_pH.pdf",
    width=12,height=12*3/5)
print(Figure_8_AS)
dev.off()

Figure_9_AS = ggarrange(Figure_Total_temp_MAGs,Figure_Total_temp_ISO,
                        Figure_amino_temp_MAGs,Figure_amino_temp_ISO, 
                        Figure_GH_temp_MAGs,Figure_GH_temp_ISO,
                        Figure_protein_temp_MAGs,Figure_protein_temp_ISO,
                        Figure_CAZy_temp_MAGs,Figure_CAZy_temp_ISO,
                        labels = c("A","B","C","D", "E", "F", "G", "H","I", "J"),
                        ncol = 2, nrow = 5) + theme(panel.background = element_blank())
Figure_9_AS
pdf("Output_Data/Figures/Figure_9_AS_temp.pdf",
    width=12,height=12*3/5)
print(Figure_9_AS)
dev.off()

Figure_10_AS = ggarrange(Figure_Total_biofilm_MAGs,Figure_Total_biofilm_ISO,
                        Figure_amino_biofilm_MAGs,Figure_amino_biofilm_ISO, 
                        Figure_GH_biofilm_MAGs,Figure_GH_biofilm_ISO,
                        Figure_protein_biofilm_MAGs,Figure_protein_biofilm_ISO,
                        Figure_CAZy_biofilm_MAGs,Figure_CAZy_biofilm_ISO,
                        labels = c("A","B","C","D", "E", "F", "G", "H","I", "J"),
                        ncol = 2, nrow = 5) + theme(panel.background = element_blank())
Figure_10_AS
pdf("Output_Data/Figures/Figure_10_AS_biofilm.pdf",
    width=12,height=12*3/5)
print(Figure_10_AS)
dev.off()

Figure_11_AS = ggarrange(Figure_Total_osmolyte_MAGs,Figure_Total_osmolyte_ISO,
                         Figure_amino_osmolyte_MAGs,Figure_amino_osmolyte_ISO, 
                         Figure_GH_osmolyte_MAGs,Figure_GH_osmolyte_ISO,
                         Figure_protein_osmolyte_MAGs,Figure_protein_osmolyte_ISO,
                         Figure_CAZy_osmolyte_MAGs,Figure_CAZy_osmolyte_ISO,
                         labels = c("A","B","C","D", "E", "F", "G", "H","I", "J"),
                         ncol = 2, nrow = 5) + theme(panel.background = element_blank())
Figure_11_AS
pdf("Output_Data/Figures/Figure_11_AS_osmolyte.pdf",
    width=12,height=12*3/5)
print(Figure_11_AS)
dev.off()

