works_with_R("3.2.2",
             ggplot2="1.0.1",
             "tdhock/namedCapture@a31a38be12d4dad4aec6257a47c0e2307679589c",
             data.table="1.9.6")

all.txt.vec <- c(oncoscan.txt.vec, Sys.glob("aCGH/*/*"))

chr.pos.pattern <- paste0(
  "(?<chrom>chr.*?)",
  ":",
  "(?<chromStart>[0-9]+)",
  "-",
  "(?<chromEnd>[0-9]*)")

profile.list <- list()

oncoscan.txt.vec <- Sys.glob("OncoScan/*")
for(oncoscan.txt in oncoscan.txt.vec){
  oncoscan.dt <- fread(oncoscan.txt)
  oncoscan.dt[, chrom := {
    factor(Chromosome, 1:24, paste0("chr", c(1:22, "X", "Y")))
  }]
  ## ASSUME chrom 23 = X and 24 = Y.
  oncoscan.dt[, chromStart := Position]
  oncoscan.dt[, chromEnd := Position+1]
  oncoscan.dt[, logratio := logRatio]
  profile.list[[oncoscan.txt]] <- oncoscan.dt
}

system("head -30 aCGH/Agilent_v2/08_318.txt|nl")
aCGH.dt <- fread("aCGH/Agilent_v1/03_141.txt", skip=9)

aCGH.txt.vec <- Sys.glob("aCGH/Agilent_v1/*")
for(aCGH.txt in aCGH.txt.vec){
  aCGH.dt <- fread(aCGH.txt, skip=9)
  stopifnot(aCGH.dt$FeatureNum[1]==1)
  match.df <- str_match_named(aCGH.dt$SystematicName, chr.pos.pattern, list(
    chromStart=as.integer, chromEnd=as.integer))
  ## ASSUME SystematicName which does not match a genome position is
  ## not important.
  ignore <- is.na(match.df$chrom)
  na.names <- aCGH.dt$SystematicName[ignore]
  probes <- data.table(match.df, aCGH.dt)[!ignore, ]
  probes[, logratio := LogRatio]
  profile.list[[aCGH.txt]] <- probes
}

aCGH.txt <- "aCGH/Agilent_v2/08_318.txt"
cmd <- paste("cat", aCGH.txt, "| sed 's/\\s*$//'")
aCGH.dt <- fread(cmd)

aCGH.txt.vec <- Sys.glob("aCGH/Agilent_v2/*")
for(aCGH.txt in aCGH.txt.vec){
  cmd <- paste("cat", aCGH.txt, "| sed 's/\\s*$//'")
  aCGH.dt <- fread(cmd)
  logratio.name <- grep("log ratio", names(aCGH.dt), value=TRUE)
  aCGH.dt$logratio <- aCGH.dt[[logratio.name]]
  aCGH.dt[, chrom := paste0("chr", ChrName)]
  aCGH.dt[, chromStart := Start]
  aCGH.dt[, chromEnd := Stop]
  profile.list[[aCGH.txt]] <- aCGH.dt
}

## Both Nimblegen files types can be read using the following code.
aCGH.dt <- fread("aCGH/Nimblegen_72K/01_003.txt")
aCGH.dt <- fread("aCGH/Nimblegen_720K/01_300.txt")

ggplot()+
  geom_point(aes(POSITION, RATIO, color="RATIO"),
             pch=1,
             data=aCGH.dt)+
  geom_point(aes(POSITION, RATIO_CORRECTED, color="RATIO_CORRECTED"),
             pch=1,
             data=aCGH.dt)+
  facet_grid(CHROMOSOME ~ .)+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))

aCGH.txt.vec <- Sys.glob("aCGH/Nimblegen_*/*")
for(aCGH.txt in aCGH.txt.vec){
  aCGH.dt <- fread(aCGH.txt)
  ## ASSUME RATIO_CORRECTED is logratio
  aCGH.dt[, logratio := RATIO_CORRECTED]
  aCGH.dt[, chrom := CHROMOSOME]
  aCGH.dt[, chromStart := POSITION]
  aCGH.dt[, chromEnd := POSITION+1]
  profile.list[[aCGH.txt]] <- aCGH.dt
}

not.processed <- all.txt.vec[!all.txt.vec %in% names(profile.list)]
stopifnot(length(not.processed) == 0)

just.sample <- gsub("logratio.|.txt", "", basename(names(profile.list)))
just.dir <- sub("aCGH/", "", dirname(names(profile.list)))
sample.id.vec <- gsub("_", "-", paste0(just.dir, "_", just.sample))
stopifnot(length(unique(sample.id.vec)) == length(all.txt.vec))

header.tmp <-
  paste('track',
        'type=bedGraph',
        'db=hg19',
        'export=yes',
        'visibility=full',
        'maxSegments=20',
        'alwaysZero=on',
        'share=curie.fr',
        'graphType=points',
        'yLineMark=0',
        'yLineOnOff=on',
        'name=%s',
        'description="%s"')

dir.create("bedGraph", showWarnings=FALSE)
file.i.vec <- seq_along(profile.list)
valid.chroms <- paste0("chr", c(1:22, "X", "Y"))

for(file.i in file.i.vec){
  sample.id <- sample.id.vec[[file.i]]
  file.txt <- names(profile.list)[[file.i]]
  header <- sprintf(header.tmp, sample.id, file.txt)
  profile <- profile.list[[file.i]]
  maybe.na <- profile[chrom %in% valid.chroms, data.table(
    chrom,
    chromStart=sprintf("%d", as.integer(chromStart)),
    chromEnd=sprintf("%d", as.integer(chromEnd)),
    logratio)]
  some.na <- apply(is.na(maybe.na), 1, any)
  maybe.na[some.na,]
  not.na <- maybe.na[!some.na, ]
  stopifnot(!is.na(not.na))
  cat(sprintf("%4d / %4d writing %d probes to %s\n",
              file.i, length(profile.list),
              nrow(not.na), bedGraph.gz))
  bedGraph.gz <- file.path("bedGraph", paste0(sample.id, ".bedGraph.gz"))
  con <- gzfile(bedGraph.gz, "w")
  writeLines(header, con)
  write.table(not.na, con, quote=FALSE, row.names=FALSE, col.names=FALSE)
  close(con)
}

## PROBLEM: chromStart=59160767 > max(chr19)=59128983 on Nimblegen 72K
## array. Probably not hg19. sample 01_005 uploaded.

## Same for 720K 10_296 chromStart=59160767 > max(chr19)=59128983.

