# Plotting examples 

# Substrate degradation traits----

Figure_2 = ggarrange(Figure_test_1, Figure_test_2,Figure_test_1.ISO,Figure_test_2.ISO, 
          labels = c("A","B","C","D"),
          ncol = 2, nrow = 2) + theme(panel.background = element_blank())
Figure_2
png("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Figures/Figure_2.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(Figure_2)
dev.off()

# Monomer Uptake----

Figure_3 = ggarrange(Figure_test_3, Figure_test_4,Figure_test_5,Figure_test_3.ISO,
                     Figure_test_4.ISO,Figure_test_5.ISO,
                     labels = c("A","B*","C*","D","E","F"),
                     ncol = 3, nrow = 2) + theme(panel.background = element_blank())
Figure_3
png("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Figures/Figure_3.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(Figure_3)
dev.off()

# Stress Tolerance----

Figure_4 = ggarrange(Figure_test_6,Figure_test_7,Figure_test_8,Figure_test_9,
                     Figure_test_6.ISO,Figure_test_7.ISO,Figure_test_8.ISO,
                     Figure_test_9.ISO,
                     labels = c("A","B","C","D","E","F**","G**","H"),
                     ncol = 4, nrow = 2) + theme(panel.background = element_blank())
Figure_4
png("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Figures/Figure_4.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(Figure_4)
dev.off()

# Cross-correlations----

png("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Figures/Figure_5.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(Figure_5)
dev.off()

png("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Figures/Figure_6.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(Figure_6)
dev.off()

