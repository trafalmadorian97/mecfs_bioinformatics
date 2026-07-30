function (SE_SNP, I_LD, varSNP, GC, coords, k, i) 
{
    V_SNP <- diag(k)
    if (GC == "conserv") {
        for (p in 1:nrow(coords)) {
            x <- coords[p, 1]
            y <- coords[p, 2]
            if (x != y) {
                V_SNP[x, y] <- (SE_SNP[i, y] * SE_SNP[i, x] * 
                  I_LD[x, y] * I_LD[x, x] * I_LD[y, y] * varSNP[i]^2)
            }
            if (x == y) {
                V_SNP[x, x] <- (SE_SNP[i, x] * I_LD[x, x] * varSNP[i])^2
            }
        }
    }
    if (GC == "standard") {
        for (p in 1:nrow(coords)) {
            x <- coords[p, 1]
            y <- coords[p, 2]
            if (x != y) {
                V_SNP[x, y] <- (SE_SNP[i, y] * SE_SNP[i, x] * 
                  I_LD[x, y] * sqrt(I_LD[x, x]) * sqrt(I_LD[y, 
                  y]) * varSNP[i]^2)
            }
            if (x == y) {
                V_SNP[x, x] <- (SE_SNP[i, x] * sqrt(I_LD[x, x]) * 
                  varSNP[i])^2
            }
        }
    }
    if (GC == "none") {
        for (p in 1:nrow(coords)) {
            x <- coords[p, 1]
            y <- coords[p, 2]
            if (x != y) {
                V_SNP[x, y] <- (SE_SNP[i, y] * SE_SNP[i, x] * 
                  I_LD[x, y] * varSNP[i]^2)
            }
            if (x == y) {
                V_SNP[x, x] <- (SE_SNP[i, x] * varSNP[i])^2
            }
        }
    }
    return(V_SNP)
}
