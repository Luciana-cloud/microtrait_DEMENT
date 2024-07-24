# Plotting examples 

Figure_1 = ggarrange(Figure_test_1, Figure_test_2, Figure_test_3,Figure_test_6,
          Figure_test_1.ISO,Figure_test_2.ISO,Figure_test_3.ISO,Figure_test_6.ISO, 
          labels = c("A","B","C","D","E","F","G","H"),
          ncol = 4, nrow = 2) + theme(panel.background = element_blank())
Figure_1
png("C:/luciana_datos/UCI/Project_2 (microtrait-dement)/MICROTRAIT_DEMENT/Figures/Figure_1.png",
    width=3500*1.35,height=1969*1.35,res=300)
print(Figure_1)
dev.off()
