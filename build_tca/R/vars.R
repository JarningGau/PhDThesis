htca.markers <- c(
  "VIM","CD14","VWF","ACTA2","DLK1","INSL3","AMH", "SOX9","CCL21","CLDN3",
  "DAZL","UTF1","ID4","KIT","STRA8","DMRT1",
  "PRDM9","ZCWPW1","DMC1","RAD51AP2","SPO11","MLH3",
  "PIWIL1","POU5F2","CCNA1","SUN5",
  "RUNX2","TNP1","PRM1","PRM2"
)

mtca.markers <- c(
  "Vim","Cd14","Vwf","Acta2","Cd34","Sox9",
  "Dazl","Pou5f1","Mki67","Utf1","Dnmt3l","Id4","Gfra1","Nanos3","Esx1","Kit","Stra8","Dmrt1",
  "Prdm9","Zcwpw1","Dmc1","Rad51ap2","Sycp3","Mlh3",
  "Piwil1","Pou5f2","Ccna1","Sun5",
  "Tex36","Tnp1","Prm1","Prm2")

corrected.mmu.age <- c(
  "E6.5" = "6.5 dpc",
  "E7.5" = "7.5 dpc",
  "E8.5" = "8.5 dpc",
  "E9.5" = "9.5 dpc",
  "E10.5" = "10.5 dpc",
  "E11.5" = "11.5 dpc",
  "E12.5" = "12.5 dpc",
  "E13.5" = "13.5 dpc",
  "E14.5" = "14.5 dpc",
  "E15.5" = "15.5 dpc",
  "E16.5" = "16.5 dpc",
  "E17.5" = "17.5 dpc",
  "E18.5" = "18.5 dpc",
  "P0" = "0 dpp",
  "P2" = "2 dpp",
  "P3" = "3 dpp",
  "P6" = "6 dpp",
  "P7" = "7 dpp",
  "PND0" = "0 dpp",
  "PND1" = "1 dpp",
  "PND2" = "2 dpp",
  "PND3" = "3 dpp",
  "PND4" = "4 dpp",
  "PND5" = "5 dpp",
  "PND6" = "6 dpp",
  "PND7" = "7 dpp",
  "PND8" = "8 dpp",
  "PND10" = "10 dpp",
  "PND12" = "12 dpp",
  "PND14" = "14 dpp",
  "PND30" = "30 dpp",
  "PND35" = "35 dpp",
  "7 wk" = "Adult",
  "7.5 wk" = "Adult",
  "8 wk" = "Adult",
  "9 wk" = "Adult",
  "10 wk" = "Adult",
  "11 wk" = "Adult",
  "14 wk" = "Adult",
  "15 wk" = "Adult",
  "16 wk" = "Adult",
  "18 wk" = "Adult",
  "20 wk" = "Adult",
  "7 and 11 wk" = "Adult",
  "67 dpp" = "Adult"
)

mmu.age.levels <- c(paste(6:18+0.5, "dpc"), paste(c(0:8,10,12,14,15,18,20,25,30,35), "dpp"), "Adult")

