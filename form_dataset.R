library(data.table)

## 1-0
load("~/genome/step1_omni/snp_info.RData");
load("~/genome/step2_data0/genotype_with_grp.RData");

IMP <- function(x, bank) {
  pm1 <- nrow(x)
  pm2 <- ncol(x)
  out <- x;
  for (i in 1:pm2) {
    out[is.na(out[ ,i]),i] <- bank[i]
  };
  return(out)
};



# 2-1 SNP from qc_HWE
Raw.HWE <- fread("~/genome/step3_data2/train2.hwe", header=TRUE) # @____@
Raw.HWE <- as.data.frame(Raw.HWE)
AFF.HWE <- Raw.HWE[Raw.HWE$TEST=="AFF", ]

Subset1_Xrs <- as.character(AFF.HWE[AFF.HWE$P<0.0001,2]) # total = 102
rm(AFF.HWE)
rm(Raw.HWE)

# 2-2 SNP from papers
Subset2_Xrs <- read.table("~/Dropbox/shared/dataset/snplist_external_manual review.txt")[ ,1]
Subset2_Xrs <- as.character(unique(Subset2_Xrs)) # total = 618
sum(Subset2_Xrs %in% SNP_Info2$Xrs) # total left = 41

# 2-3 SNP from standard analysis
Std.Ass <- as.data.frame(fread("~/genome/standard_test/trial_1.assoc", header=T))
Subset3_Xrs <- Std.Ass[Std.Ass$P<0.0001,2] # 34

# 2-4 combine two sources

# Dataset1
SNP.Keep1 <- c(Subset2_Xrs)
SNP.Keep1 <- SNP.Keep1[SNP.Keep1 %in% SNP_Info2$Xrs] # total = 41

# Dataset2
SNP.Keep2 <- unique(c(Subset2_Xrs, Subset3_Xrs))
SNP.Keep2 <- SNP.Keep2[SNP.Keep2 %in% SNP_Info2$Xrs] # total = 73

# Dataset3
SNP.Keep3 <- unique(c(Subset1_Xrs, Subset2_Xrs, Subset3_Xrs))
SNP.Keep3 <- SNP.Keep3[SNP.Keep3 %in% SNP_Info2$Xrs] # total = 169


## 3-0
Orig_Name <- colnames(Raw_Geno)
Orig_Name <- Orig_Name[-c(1:6)]
colnames(Raw_Geno)[-c(1:6)] <- SNP_Info2$Xrs

Var_Keep <- SNP_Info2$Xrs
Var_Keep <- Var_Keep %in% SNP.Keep3
Var_Keep <- c(rep(TRUE, 6), Var_Keep)

Data_Grp2 <- Raw_Geno[ ,Var_Keep]

# 3-1 imputation
Data.Cont <- Data_Grp2[Data_Grp2$PHENOTYPE==1, -c(1:6)]
Data.Case <- Data_Grp2[Data_Grp2$PHENOTYPE==2, -c(1:6)]

Mean.Cont <- as.vector(colMeans(as.matrix(Data.Cont), na.rm = TRUE))
Mean.Case <- as.vector(colMeans(as.matrix(Data.Case), na.rm = TRUE))

IMP.Cont <- IMP(x=Data.Cont, bank=Mean.Cont)
sum(is.na(IMP.Cont))
IMP.Case <- IMP(x=Data.Case, bank=Mean.Case)
sum(is.na(IMP.Case))

Data_Impu <- rbind(IMP.Case, IMP.Cont);
Data_Impu <- cbind(Data_Grp2[ ,1:6], Data_Impu);

# 4-0 Finalize the training and the test data
Data_Train <- Data_Impu[CV_Fold$GROUP!=2, ]
Data_Train$PHENOTYPE <- Data_Train$PHENOTYPE - 1

Data_Test <- Data_Impu[CV_Fold$GROUP==2, ]
Data_Test$PHENOTYPE <- Data_Test$PHENOTYPE - 1

drop <- ls()
drop <- drop[!(drop %in% c("Data_Grp2",
                           "Data_Impu",
                           "Data_Train",
                           "Data_Test",
                           "SNP.Keep1",
                           "SNP.Keep2",
                           "SNP.Keep3",
                           "Subset1_Xrs",
                           "Subset2_Xrs",
                           "Subset3_Xrs"))]
rm(list=drop)
rm(drop)
save.image("data_ready_2.RData")



