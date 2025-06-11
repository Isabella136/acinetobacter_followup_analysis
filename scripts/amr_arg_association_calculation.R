# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Missing argument(s)")
} else if (length(args) > 4) {
    stop("Extra argument(s)")
}
blast.matrix <- args[1]
disk.diffusion <- args[2]
phenotype.clsi <- args[3]
output <- args[4]

# Open blast matrix csv file as data frame
blast.matrix.df <- read.csv(blast.matrix, header = TRUE, row.names = 1)

# Open disk diffusion tsv file as data frame
disk.diffusion.df <- read.csv(disk.diffusion, header = TRUE, sep = '\t', row.names = 1, na.strings = "-1")

# Open phenotype clsi csv file as data frame
phenotype.clsi.df <- read.csv(phenotype.clsi, header = TRUE, row.names = 1)

# Turn into a binary matrix
antimicrobials.with.known.stds <- as.character(rownames(phenotype.clsi.df))
antimicrobials.tested <- as.character(colnames(disk.diffusion.df))
for (antimicrobial in antimicrobials.with.known.stds)
{
    if (antimicrobial %in% antimicrobials.tested) {
        disk.diffusion.df[[antimicrobial]][disk.diffusion.df[[antimicrobial]] <= phenotype.clsi.df[antimicrobial,"R"] & !is.na(disk.diffusion.df[[antimicrobial]])] <- 1
        disk.diffusion.df[[antimicrobial]][disk.diffusion.df[[antimicrobial]] > phenotype.clsi.df[antimicrobial,"R"] & !is.na(disk.diffusion.df[[antimicrobial]])] <- 0
    }
}

# Remove antimicrobials that don't have a clsi standard
for (antimicrobial in antimicrobials.tested)
{
    if (! antimicrobial %in% antimicrobials.with.known.stds ) {
        disk.diffusion.df <- disk.diffusion.df[, colnames(disk.diffusion.df) != antimicrobial]
    }
}

# Order row names for blast.matrix.df and retrieve needed rows for disk.diffusion.df
blast.matrix.df <- blast.matrix.df[order(rownames(blast.matrix.df)),]

if ('-' %in% rownames(blast.matrix.df)[1]) {
    disk.diffusion.df <- disk.diffusion.df[rownames(blast.matrix.df),]
} else {
    rownames(disk.diffusion.df) <- sapply(strsplit(rownames(disk.diffusion.df), split = ' ', fixed = TRUE), "[" , 3)
    print(rownames(disk.diffusion.df))
    disk.diffusion.df <- disk.diffusion.df[rownames(blast.matrix.df),]
}

# Calculate Fisher
fisher.matrix <- matrix(
    data = NA, 
    nrow = length(colnames(blast.matrix.df)), 
    ncol = length(colnames(disk.diffusion.df)),
    dimnames= list(colnames(blast.matrix.df), colnames(disk.diffusion.df)))

for (row.gene in colnames(blast.matrix.df)) {
    for (col.antimicrobial in colnames(disk.diffusion.df)) {
        yes.res.yes.gene <- length(disk.diffusion.df[[col.antimicrobial]][!is.na(disk.diffusion.df[[col.antimicrobial]]) & disk.diffusion.df[[col.antimicrobial]] == 1 & blast.matrix.df[[row.gene]] == 1])
        yes.res.no.gene <- length(disk.diffusion.df[[col.antimicrobial]][!is.na(disk.diffusion.df[[col.antimicrobial]]) & disk.diffusion.df[[col.antimicrobial]] == 1 & blast.matrix.df[[row.gene]] == 0])
        no.res.yes.gene <- length(disk.diffusion.df[[col.antimicrobial]][!is.na(disk.diffusion.df[[col.antimicrobial]]) & disk.diffusion.df[[col.antimicrobial]] == 0 & blast.matrix.df[[row.gene]] == 1])
        no.res.no.gene <- length(disk.diffusion.df[[col.antimicrobial]][!is.na(disk.diffusion.df[[col.antimicrobial]]) & disk.diffusion.df[[col.antimicrobial]] == 0 & blast.matrix.df[[row.gene]] == 0])

        # If gene found in 10 or more susceptible isolates, skip
        if (no.res.yes.gene > 10) {
            next
        }

        contingency.table <- matrix(
            data = c(yes.res.yes.gene, yes.res.no.gene, no.res.yes.gene, no.res.no.gene), 
            nrow = 2, 
            ncol = 2)
        fisher.results <- fisher.test(contingency.table)

        # If p-value <= 0.01, add to matrix
        if (fisher.results$p.value <= 0.01) {
            fisher.matrix[row.gene, col.antimicrobial] <- fisher.results$estimate
            print(row.gene)
            print(col.antimicrobial)  
            print(yes.res.yes.gene)
            print(yes.res.no.gene)  
            print(no.res.yes.gene)
            print(no.res.no.gene)  
            print(fisher.matrix[row.gene, col.antimicrobial])
        }
    }
}

write.csv(data.frame(fisher.matrix), output)