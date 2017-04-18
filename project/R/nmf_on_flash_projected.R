library(NMF)

data <- get(load("flash_projected.rda"))

out <- nmf(t(data), rank=20)
save(out, file="nmf_on_flash_projected_gtex_K_20.rda")



