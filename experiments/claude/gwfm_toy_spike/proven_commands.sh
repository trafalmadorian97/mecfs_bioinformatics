# Proven GWFM toy pipeline (all exit 0, causals rs11/rs171/rs331 PIP=1.0)
gctb --bfile toy --make-block-ldm --block-info blocks.txt --out ldm --thread 1
gctb --ldm ldm --gwas-summary toy.ma --make-ldm-eigen --ldm-eigen-cutoff 0.995 --thread 1 --out matched
gctb --gwfm RC --ldm-eigen matched --gwas-summary toy.ma --annot annot.txt --gene-map genemap.txt --thread 1 --out gwfm --chain-length 3000 --burn-in 1000
gctb --cs --pwld-file ldm/rsq0.5.pwld --pip 0.9 --pep 0.7 --gene-map genemap.txt --mcmc-samples gwfm --out gwfm
