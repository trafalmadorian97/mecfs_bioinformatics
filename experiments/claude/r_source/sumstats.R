function (files, ref, trait.names = NULL, se.logit, OLS = NULL, 
    linprob = NULL, N = NULL, betas = NULL, info.filter = 0.6, 
    maf.filter = 0.01, keep.indel = FALSE, parallel = FALSE, 
    cores = NULL, ambig = FALSE, direct.filter = FALSE) 
{
    if (is.list(files)) {
        wrn <- paste0("DeprecationWarning: In future versions a list of filenames will no longer be accepted.\n", 
            "                    Please change files to a vector to ensure future compatibility.")
        warning(wrn)
        files_ <- c()
        for (i in 1:length(files)) {
            files_ <- c(files_, files[[i]])
        }
        files <- files_
    }
    len <- length(files)
    if (is.null(N)) 
        N <- rep(NA, len)
    if (is.null(OLS)) {
        OLS <- rep(FALSE, len)
    }
    if (is.null(linprob)) {
        linprob <- rep(FALSE, len)
    }
    if (is.null(betas)) {
        betas <- rep(FALSE, len)
    }
    .check_file_exists(files)
    .check_file_exists(ref)
    .check_equal_length(files, trait.names)
    .check_equal_length(files, OLS)
    .check_equal_length(files, linprob)
    .check_equal_length(files, N)
    .check_equal_length(files, betas)
    .check_range(info.filter)
    .check_range(maf.filter)
    .check_range(N, min = 0, max = Inf, allowNA = TRUE)
    .check_boolean(keep.indel)
    .check_boolean(parallel)
    if (!is.null(cores)) 
        .check_range(cores, min = 0, max = Inf)
    begin.time <- Sys.time()
    filenames <- files
    ref2 <- ref
    if (is.null(trait.names)) {
        names.beta <- paste0("beta.", 1:len)
        names.se <- paste0("se.", 1:len)
    }
    else {
        names.beta <- paste0("beta.", trait.names)
        names.se <- paste0("se.", trait.names)
    }
    log2 <- paste(trait.names, collapse = "_")
    log2 <- str_remove_all(log2, "/")
    if (object.size(log2) > 200) {
        log2 <- substr(log2, 1, 100)
    }
    log.file <- file(paste0(log2, "_sumstats.log"), open = "wt")
    .LOG("The preparation of ", length(trait.names), " summary statistics for use in Genomic SEM began at: ", 
        begin.time, file = log.file)
    .LOG("Please note that the files should be in the same order that they were listed for the ldsc function", 
        file = log.file)
    .LOG("Reading in reference file", file = log.file)
    ref <- fread(ref, header = T, data.table = F)
    .LOG("Applying MAF filer of ", maf.filter, " to the reference file.", 
        file = log.file)
    ref <- subset(ref, ref$MAF >= maf.filter)
    if (ambig) {
        .LOG("Removing strand ambiguous SNPs.", file = log.file)
        ref <- subset(ref, (ref$A1 != "T" | ref$A2 != "A"))
        ref <- subset(ref, (ref$A1 != "A" | ref$A2 != "T"))
        ref <- subset(ref, (ref$A1 != "C" | ref$A2 != "G"))
        ref <- subset(ref, (ref$A1 != "G" | ref$A2 != "C"))
    }
    data.frame.out <- ref
    if (!parallel) {
        files <- lapply(files, read.table, header = T, quote = "\"", 
            fill = T, colClasses = c(P = "character"), na.string = c(".", 
                NA, "NA", ""))
        .LOG("All files loaded into R!", file = log.file)
        Output <- list()
        for (i in 1:len) {
            Output[[i]] <- .sumstats_main(i, utilfuncs = NULL, 
                filenames[i], trait.names[i], N[i], keep.indel, 
                OLS[i], betas[i], info.filter, linprob[i], se.logit[i], 
                names.beta[i], names.se[i], ref, ref2, files[[i]], 
                log.file, direct.filter)
        }
    }
    else {
        if (is.null(cores)) {
            int <- detectCores() - 1
        }
        else {
            int <- cores
        }
        if (int > length(filenames)) {
            .LOG("Number of requested cores(", int, ") greater than the number of files (", 
                length(filenames), "). Deferring to the lowest number", 
                file = log.file)
            int <- length(filenames)
        }
        cl <- makeCluster(int, type = "PSOCK")
        registerDoParallel(cl)
        on.exit(stopCluster(cl))
        utilfuncs <- list()
        utilfuncs[[".get_renamed_colnames"]] <- .get_renamed_colnames
        utilfuncs[[".LOG"]] <- .LOG
        utilfuncs[["inner_join"]] <- inner_join
        .LOG("As parallel sumstats was requested, logs of each file will be saved separately", 
            file = log.file)
        Output <- foreach(i = 1:length(filenames), .export = c(".sumstats_main"), 
            .packages = c("stringr")) %dopar% {
            .sumstats_main(i, utilfuncs, filenames[i], trait.names[i], 
                N[i], keep.indel, OLS[i], betas[i], info.filter, 
                linprob[i], se.logit[i], names.beta[i], names.se[i], 
                ref, ref2, NULL, NULL, direct.filter)
        }
    }
    for (i in 1:len) {
        data.frame.out <- suppressWarnings(inner_join(data.frame.out, 
            Output[[i]], by = "SNP"))
    }
    end.time <- Sys.time()
    total.time <- difftime(time1 = end.time, time2 = begin.time, 
        units = "sec")
    mins <- floor(floor(total.time)/60)
    secs <- total.time - mins * 60
    .LOG("     ", file = log.file, print = FALSE)
    .LOG("After merging across all summary statistics using listwise deletion, performing QC, and merging with the reference file, there are ", 
        nrow(data.frame.out), " SNPs left in the final multivariate summary statistics file", 
        file = log.file)
    .LOG("Sumstats finished running at ", end.time, file = log.file)
    .LOG("Running sumstats for all files took ", mins, " minutes and ", 
        secs, " seconds", file = log.file)
    .LOG("Please check the log file ", paste0(log2, "_sumstats.log"), 
        " to ensure that all columns were interpreted correctly and no warnings were issued for any of the summary statistics files.", 
        file = log.file)
    flush(log.file)
    close(log.file)
    return(data.frame.out)
}
