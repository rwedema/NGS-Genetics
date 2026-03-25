library(ggplot2)
library(DNAcopy)

df <- read.table("/students/2025-2026/Thema05/eyelid_carcinoma/data/varscan_output/subset58_57_varscan/output.copynumber",
                 header=T, stringsAsFactors = F)

df$pos <- df$chr_start
CNA.object <- CNA(genomdat=df$log2_ratio, chrom=df$chrom, maploc=df$pos, data.type="logratio")

sm <- smooth.CNA(CNA.object)
seg <- segment(sm, verbose=1, alpha=0.01)
segments <- seg$output

head(segments)

chrlen <- aggregate(chr_stop ~ chrom, data=df, FUN=max)
chrlen$cumstart <- c(0, cumsum(as.numeric(head(chrlen$chr_stop, -1))))
chrlen$cumend <- chrlen$cumstart + chrlen$length
chrlen$midpoint <- (chrlen$cumstart + chrlen$cumend) / 2

names(chrlen)[2] <- "length"

df <- merge(df, chrlen[,c("chrom", "cumstart")], by="chrom")
df$cumpos <- df$pos + df$cumstart
segments <- merge(segments, chrlen[,c("chrom", "cumstart")], by.x="chrom", by.y="chrom")
segments$cumstart <- segments$loc.start + segments$cumstart
segments$cumend <- segments$loc.end + segments$cumstart

p <- ggplot() +
  geom_point(data=df, aes(x=cumpos, y=log2_ratio), size=0.4, alpha=0.6) +
  geom_segment(data=segments, aes(x=cumstart, xend=cumend, y=seg.mean, yend=seg.mean), color="red", linewidth=0.6) +
  scale_x_continuous(
    breaks = chrlen$midpoint,
    labels = chrlen$chrom,
    expand = c(0.01, 0.01)
  ) +
  labs(x="", y="log2(tumor/normal)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  

print(p)

## Chromo Plot

The `vcfR` package contains functions to create a 'chromo plot' that shows an overview of both the variants and statistics regarding mapping quality. Creating this plot requires extra data besides the variant file, namely:
  * the reference sequence (genome or gene) in FASTA format and
* a file describing *features* (genes, transcripts, CDSs, introns, exons, etc.) on this reference genome (in GFF format).

For this project we can use a subset of these files that only span our region of interest (a single gene for instance). 

**Note:** this is experimental and not recommended right now; subsetting all required data for a single gene works, but currently is not displayed properly.

```{r, include=FALSE}
vcf_data <- read.delim("solution/variants_cardio.vcf", header = FALSE, comment.char = "#")
freq <- strsplit(vcf_data$V8, split = ";")
freq <- as.numeric(gsub(x = do.call(rbind, freq)[,2], pattern = "AF=", replacement = ""))
vcf_GRanges <- GRanges(seqnames = vcf_data$V1, 
                       ranges = IRanges(start = vcf_data$V2, 
                                        end = vcf_data$V2))
mcols(vcf_GRanges) <- DataFrame(AF = freq)

# Determine overlap with BED data
overlaps <- findOverlaps(vcf_GRanges, bed_data)

queryHits(overlaps)[vcf_GRanges[queryHits(overlaps)]$AF < 0.8]

overlap_df <- DataFrame(overlaps)
overlap_df$Gene <- names(bed_data)[overlap_df$subjectHits]
# plot variants per gene
queryHits_tab_df <- as.data.frame(table(overlap_df$subjectHits))
queryHits_tab_df$Var1 <- names(bed_data)[queryHits_tab_df$Var1]
#barplot(queryHits_tab_df$Freq, names.arg = queryHits_tab_df$Var1, cex.names = 0.7, las = 2)
```



```{r echo=FALSE, message=FALSE, warning=FALSE, include=FALSE}
# Get the ranges (including introns) for all genes:
panel_ranges <- unlist(range(bed_data))

# Get sequences for all genes
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# Only run to install once
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)

Hsapiens <- BSgenome.Hsapiens.UCSC.hg38

## random tx subset
tx <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)

## extract sequence
seq <- getSeq(Hsapiens, "chr6")#panel_ranges)

## add names
names(seq) <- names(panel_ranges)
```

```{r, include=FALSE, echo=FALSE}
library(vcfR)
library(ape)
library(rtracklayer)

vcf_file <- "solution/variants.vcf"
dna <- as.DNAbin(seq$DSP)
names(dna) <- "chr6"
vcf <- read.vcfR(vcf_file)
vcf_chr6 <- vcf[getCHROM(vcf) == "chr6"]

export.gff(bed_data$DSP, con = "solution/bed_data.gff")
gff <- read.table("solution/bed_data.gff", sep="\t", quote="")

chrom <- create.chromR(name='DSP', vcf=vcf_chr6, seq=dna, ann=gff)
plot(chrom)
chromoqc(chrom, xlim=c(min(gff$V4), max(gff$V5)))
```


