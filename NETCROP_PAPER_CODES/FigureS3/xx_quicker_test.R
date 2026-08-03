setwd(file.path(here::here("NETCROP_PAPER_CODES", "FigureS3")))

HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))

source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "SBM_DCBM_helpers.R"))

RUN_NAME <- "xx_quicker_FigureS3_small_networks"
run.paths <- netcrop_output_action(file.path("output", RUN_NAME),
                                   file.path("logs", RUN_NAME))
OUTPUT_DIR <- run.paths$output_dir
LOG_DIR <- run.paths$log_dir
OUTPUT_ACTION <- run.paths$action

################################################################################

detected.cores <- parallel::detectCores()
ncore <- if (is.na(detected.cores)) 1L else max(1L, floor(detected.cores/2)) # it can be set to anything else
nsim <- 10L

n.all <- 100*(1:2)
K.all <- c(3)
beta.all <- 0.3
alpha <- 1
# pi.all <- c("balanced", "unbalanced")
pi.all <- c("balanced")

edge.density.all <- 0.25
p.test.all <- 0.1
o.range.all <- 0

mod.all <- c('SBM', 'DCBM')

big.loop <- as.data.frame(expand.grid(nn = n.all, KK = K.all, bbb = beta.all,
                                      pipi = pi.all, ee = edge.density.all, mod = mod.all))

output.file <- file.path(OUTPUT_DIR, "FigureS3_small_networks.csv")
loopers <- 1L:nrow(big.loop)
if (identical(OUTPUT_ACTION, "resume") && file.exists(output.file)) {
  existing <- readr::read_csv(output.file, show_col_types = FALSE)
  keys <- c("looper", "nsim", "algorithm", "R", "p_test", "o_range")
  existing <- existing[!duplicated(existing[keys]), , drop = FALSE]
  expected <- nsim * (3L + 2L + 2L)
  counts <- table(existing$looper)
  complete <- as.integer(names(counts)[counts >= expected])
  existing <- existing[existing$looper %in% complete, , drop = FALSE]
  readr::write_csv(existing, output.file)
  loopers <- loopers[!loopers %in% complete]
  if (length(loopers) == 0L) message("All Figure S3 configurations are complete.")
}

for(looper in loopers){
  bloop <- big.loop[looper, ]
  nn <- bloop$nn
  KK <- bloop$KK
  bbb <- bloop$bbb
  pipi <- bloop$pipi
  ee <- bloop$ee
  mod <- bloop$mod

  if(pipi == "balanced"){
    PI <- rep(1/KK, KK)
  }else{
    PI <- c(0.6, rep(0.4/(KK-1), KK-1))
  }

  avg.deg <- nn * ee
  generator <- if (mod == "SBM") SBM.gen else DCBM.gen
  BB <- alpha * diag(1-bbb, KK) + bbb

  cat("\n--------------------------------------------")
  cat("\n---", paste(looper, nn, KK, bbb, pipi, ee, mod, sep = ":"), "Started---\n")
  all.out <- parallel::mclapply(1:nsim, function(ii){
    cat("--", ii, ": gen --\n")
    # on.exit(gc())
    t0 <- proc.time()

    net <- generator(n = nn, K = KK, g = NULL, B = BB,
                     # beta = bbb,
                     PI = PI,
                     avg.deg = avg.deg, ncore = 1, seed = 100*ii)

    lambda <- mean(rowSums(net$A))

    nc.all <- tibble::tibble()

    for(tete in p.test.all){
      for(ooo in o.range.all){
        param.out <- netcrop_param(p.test = tete, n = nn, o.range = ooo)

        ss <- param.out$s
        oo <- param.out$o

        for(rr in c(1, 3, 5)){
          cat("--", ii, ": netcrop --", rr, "\n")
          err <- 0
          err <- tryCatch({
            ram.nc <- peakRAM::peakRAM({time.nc <- system.time({
              out.nc <- netcrop_blockmodel(
                A = net$A, K.CAND = 5, s = ss, o = oo, R = rr,
                laplace = F, dc.est = 2,
                loss = c("l2"),
                mod.cand = c("SBM", "DCBM"), ncore = 1, seed = 500 + 20*rr*ii,
                rngR = TRUE)
            })[3]})$Peak_RAM_Used_MiB},
            error = function(e){-1})

          if(err == -1){
            out.nc <- list(l2.model = NA_character_)
            time.nc <- -1
            ram.nc <- -1
          }

          netcrop_status(ii, nsim, "NETCROP", time.nc,
                         out.nc$l2.model, rr)

          nc.all <- dplyr::bind_rows(nc.all, tibble::tibble(
            looper = looper,
            nsim = ii, model = mod, n = nn, K = KK, PI = pipi,
            beta = bbb, rho = ee, p_test = tete, o_range = ooo,
            algorithm = "NETCROP",
            s = ss, o = oo, R = rr, cv = -1,
            best.l2 = out.nc$l2.model,
            run_time = time.nc, ram_MiB = ram.nc
          ))
        } #rr
      } #ooo
    } #tete

    for(stst in c(1, 20)){
      cat("--", ii, ": NCV --", stst, "\n")

      err <- 0
      err <- tryCatch({
        ram.ncv <- peakRAM::peakRAM({time.ncv <- system.time({
          out.ncv <- NCV.stability.BM(
            A = net$A, max.K = 5, cv = 3, R = stst, tau = 0,
            laplace = F, dc.est = 2,
            loss = c("l2"), ncore = 1, seed = 500 + 20*stst*ii)
        })[3]})$Peak_RAM_Used_MiB},
        error = function(e){-1})

      if(err == -1){
        out.ncv <- list(best.l2.stable = NA_character_)
        time.ncv <- -1
        ram.ncv <- -1
      }

      netcrop_status(ii, nsim, "NCV", time.ncv,
                     out.ncv$best.l2.stable, stst)

      nc.all <- dplyr::bind_rows(nc.all, tibble::tibble(
        looper = looper,
        nsim = ii, model = mod, n = nn, K = KK, PI = pipi,
        beta = bbb, rho = ee, p_test = -1, o_range = -1,
        algorithm = "NCV",
        s = -1, o = -1, R = stst, cv = 3,
        best.l2 = out.ncv$best.l2.stable,
        run_time = time.ncv, ram_MiB = ram.ncv
      ))
      cat("--", ii, ": ECV --", stst, "\n")

      err <- 0
      err <- tryCatch({
        ram.ecv <- peakRAM::peakRAM({time.ecv <- system.time({
          out.ecv <- ECV.stability.BM(
            A = net$A, max.K = 5, train.p = 0.9,
            cv = 3, R = stst, tau = 0, dc.est = 2,
            loss = c("l2"), ncore = 1, seed = 500 + 20*stst*ii)
        })[3]})$Peak_RAM_Used_MiB},
        error = function(e){-1})

      if(err == -1){
        out.ecv <- list(best.l2.stable = NA_character_)
        time.ecv <- -1
        ram.ecv <- -1
      }

      netcrop_status(ii, nsim, "ECV", time.ecv,
                     out.ecv$best.l2.stable, stst)

      nc.all <- dplyr::bind_rows(nc.all, tibble::tibble(
        looper = looper,
        nsim = ii, model = mod, n = nn, K = KK, PI = pipi,
        beta = bbb, rho = ee, p_test = -1, o_range = -1,
        algorithm = "ECV",
        s = -1, o = -1, R = stst, cv = 3,
        best.l2 = out.ecv$best.l2.stable,
        run_time = time.ecv, ram_MiB = ram.ecv
      ))

      # gc()
    }

    t1 <- proc.time()
    list(data = nc.all, elapsed = (t1 - t0)[3])
  }, mc.cores = ncore)

  configuration.output <- dplyr::bind_rows(lapply(all.out, `[[`, "data"))
  readr::write_csv(configuration.output, output.file,
                   append = file.exists(output.file))
  elapsed <- vapply(all.out, function(x) as.numeric(x$elapsed), numeric(1))
  # gc()
  cat("\n---", paste(looper, nn, KK, bbb, pipi, ee, mod, sep = ":"), "Ended---Time: ",
      round(max(elapsed, na.rm = T), 2), "---")
  cat("\n--------------------------------------------")
  cat("\n--- ", looper, " / ", nrow(big.loop), " = ",
      round(100 * looper/nrow(big.loop)), "% done ---" )
  cat("\n--------------------------------------------")

} # looper

################################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

setwd(file.path(here::here("NETCROP_PAPER_CODES", "FigureS3")))
all.out <- readr::read_csv(output.file, show_col_types = FALSE)

all.sum <- all.out |>
  tidyr::complete(R = tidyr::full_seq(R, 1)) |>
  mutate(true_model = paste0(model, '-', K),
         algorithm = factor(algorithm, levels = c("NCV", "ECV", "NETCROP")),
         full.algorithm = factor(paste0(algorithm, "+", "R=",R),
                                 levels = c(
                                   "NCV+R=1", "NCV+R=20",
                                   "ECV+R=1", "ECV+R=20",
                                   paste0("NETCROP+R=", c(1, 3, 5))
                                 )),
         model = factor(model, levels = c('SBM', 'DCBM'))) |>
  group_by(model, n, algorithm, full.algorithm, R) |>
  summarize(
    nsim = n(),
    accuracy = 100*mean(best.l2 == true_model),
    mean_time = log(mean(run_time)),
    mean_ram = mean(ram_MiB),
    .groups = "drop"
  )

all.sum.long <- all.sum |>
  pivot_longer(cols = c(accuracy, mean_time, mean_ram),
               names_to = "Measure", values_to = "Value") |>
  mutate(Measure = factor(case_when(
    Measure == "accuracy" ~ "Accuracy (%)",
    Measure == "mean_time" ~ "log(Mean runtime (sec.))",
    Measure == "mean_ram" ~ "Mean RAM usage (MiB)"
  ), levels = c("Accuracy (%)", "log(Mean runtime (sec.))", "Mean RAM usage (MiB)")))

S3_small_SBM <- all.sum.long |>
  filter(model == "SBM") |>
  mutate(R = factor(R)) |>
  ggplot(aes(x = n, y = Value,
             color = algorithm, fill = algorithm,
             linetype = R,
             shape = algorithm)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ Measure, scales = "free") +
  labs(x = "Number of nodes (n)",
       color = "Algorithm", fill = "Algorithm",
       shape = "Algorithm"
  ) +
  cowplot::theme_minimal_grid() +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 10),
    strip.background = element_blank()
  )

S3_small_DCBM <- all.sum.long |>
  filter(model == "DCBM") |>
  mutate(R = factor(R)) |>
  ggplot(aes(x = n, y = Value,
             color = algorithm, fill = algorithm,
             linetype = R,
             shape = algorithm)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ Measure, scales = "free") +
  labs(x = "Number of nodes (n)",
       color = "Algorithm", fill = "Algorithm",
       shape = "Algorithm"
  ) +
  cowplot::theme_minimal_grid() +
  theme(
    legend.position = "bottom",
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank(),
    strip.text = element_text(size = 10),
    strip.background = element_blank()
  )


S3_small_plot <- cowplot::plot_grid(
  S3_small_SBM, S3_small_DCBM,
  align = "hv",
  nrow = 2,
  # common.legend = T,
  labels = c("SBM-3", "DCBM-3"),
  label_size = 12
  # font.label = list(size = 12)
)

print(S3_small_plot)

ggsave(file.path(OUTPUT_DIR, "Figure3_small_plot.png"), plot = S3_small_plot,
       device = "png", height = 4.5, width = 8, units = "in", dpi = 600,
       bg = "white")
