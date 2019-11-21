# 
# ggsave(fig1,file="img/fig1.png",width=300,height=400,units="mm",dpi=300)
# ggsave(sfig1x,file="img/sfig1x.png",width=300,height=600,units="mm",dpi=100)
# ggsave(sfig1y,file="img/sfig1y.png",width=300,height=600,units="mm",dpi=100)
# 
# # source('code_analysis/03_burst_detection_figures.R')
# ggsave("img/fig2.png",fig2,height=16,width=12)
# ggsave("img/sfig2x1.png",sfig2x1,height=450,width=1000,units = "mm")
# ggsave("img/sfig2x2.png",sfig2x2,height=450,width=1000,units = "mm")
# ggsave("img/sfig2x3.png",sfig2x3,height=450,width=1000,units = "mm")
# ggsave("img/sfig2x4.png",sfig2x4,height=225.5,width=1000,units = "mm")
# 
# 
# ggsave(fig3,file="img/fig3.png",width=350,height=185,units="mm",dpi=600)
# 
# ggsave(fig4,file="img/fig4.png",width = 400,height = 195,units = "mm",dpi=600)  
# 
# ggsave(fig5,file="img/fig5.png",width=400,height=195,units="mm",dpi=600)
# 
# ggsave(fig7,file="img/fig7.png",width=288,height=216,units="mm",dpi=600)

## Make it in TIFF format...
# 
# ggsave(fig1,file="img/fig1.tiff",width=300,height=400,units="mm",dpi=120)
# 
# # source('code_analysis/03_burst_detection_figures.R')
# ggsave("img/fig2.tiff",fig2,height=16,width=12,dpi=120)
# 
# 
# ggsave(fig3,file="img/fig3.tiff",width=350,height=185,units="mm",dpi=150)
# 
# ggsave(fig4,file="img/fig4.tiff",width = 400,height = 195,units = "mm",dpi=150)  
# 
# ggsave(fig5,file="img/fig5.png",width=400,height=195,units="mm",dpi=600)
# 
# ggsave(fig7,file="img/fig7.png",width=288,height=216,units="mm",dpi=600)



# Desperately trying to make it better resolution - then use PACE (https://pacev2.apexcovantage.com/Upload)

ggsave(fig1,file="img/fig1x.tiff",width=300,height=300,units="mm",dpi=300)

# source('code_analysis/03_burst_detection_figures.R')
ggsave("img/fig2x.tiff",fig2,height=16,width=12,dpi=300)


ggsave(fig3,file="img/fig3x.tiff",width=400,height=195,units="mm",dpi=300)

ggsave(fig4,file="img/fig4x.tiff",width = 400,height = 195,units = "mm",dpi=300)  

ggsave(fig5,file="img/fig5x.tiff",width=400,height=195,units="mm",dpi=300)

ggsave(fig7,file="img/fig7xx.tiff",width=130,height=150,units="mm",dpi=300)
