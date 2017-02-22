library(data.table);

omni <- fread("omni.raw", header=TRUE);
hwe <- fread("omni.hwe", header=TRUE);
map <- fread("~/genome/original/adgwas.map", header=FALSE);
(ncol(omni)-6)==nrow(hwe)/3;
nrow(omni)==364;

omni <- as.data.frame(omni);
hwe <- as.data.frame(hwe);
map <- as.data.frame(map);
hwe.unaff <- hwe[hwe$TEST=="UNAFF", ];
hwe.aff <- hwe[hwe$TEST=="AFF", ];

Name1 <- hwe.unaff$SNP; # no minor allele;
Name2 <- colnames(omni)[7:length(omni)]; # with minor allele;

SNP_Info <- data.frame(RS=Name1, RS_MA=Name2,
                           CHR=map$V1, DIST=map$V3, POS=map$V4);
SNP_Info$RS <- as.character(SNP_Info$RS)
SNP_Info$RS_MA <- as.character(SNP_Info$RS_MA)
drop <- ls()
drop <- drop[drop!="SNP_Info"]
rm(list=drop)
rm(drop)
save.image("snp_info.RData")



