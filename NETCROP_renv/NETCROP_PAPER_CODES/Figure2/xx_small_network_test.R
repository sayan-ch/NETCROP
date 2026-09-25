setwd(file.path(here::here("NETCROP_PAPER_CODES", "Figure2")))

HELPERS_DIR <- file.path(here::here("NETCROP_PAPER_CODES", "helpers"))
source(file.path(HELPERS_DIR, "General_helpers.R"))
source(file.path(HELPERS_DIR, "PARTUNE_RSC_helpers.R"))

RUN_NAME <- "xx_small_network_test"
run.paths <- netcrop_output_action(file.path("output", RUN_NAME),
                                   file.path("logs", RUN_NAME))
OUTPUT_DIR <- run.paths$output_dir
LOG_DIR <- run.paths$log_dir
OUTPUT_ACTION <- run.paths$action

detected.cores <- parallel::detectCores()
ncore <- if (is.na(detected.cores) | .Platform$OS.type == "windows" ) 1L else max(1L, floor(detected.cores/2))
nsim <- 2L
n <- 500L
K <- 3L
K.max <- 5L
beta <- 1 / 3
rho <- 0.3

p.test <- 0.1
param.out <- netcrop_param(p.test = p.test, n = n, o.range = 0)
s <- param.out$s
o <- param.out$o

results.file <- file.path(OUTPUT_DIR, "small_figure2_results.rds")
resume.state <- netcrop_resume_rds(results.file, nsim, OUTPUT_ACTION)
results <- resume.state$results

for (sim in resume.state$simulations) {
  net <- DCBM.gen(n = n, K = K, beta = beta, rho = rho,
                  ncore = ncore, seed = 100 + sim)

  netcrop.time <- system.time({
    netcrop.result <- netcrop.tune.regsp(
      A = net$A, K = K, tau.cand = seq(0, 2, 0.1), DCBM = TRUE,
      s = s, o = o, R = 5, laplace = TRUE, dc.est = 2,
      loss = "l2", true.g = net$member,
      ncore = ncore, seed = 200 + 10 * sim
    )
  })[3]
  netcrop_status(sim, nsim, "NETCROP", netcrop.time,
                 netcrop.result$croissant.all.accu["l2"])

  dk.time <- system.time({
    dk.result <- DKest(
      A = net$A, K = K, true.g = net$member,
      tau.cand = seq(0, 0.1, 0.01), laplace = TRUE,
      DCBM = TRUE, DC.est = 2, ncore = ncore
    )
  })[3]
  best.dk <- dk.result[which.min(dk.result[, "DK.stat"]), "tau.cand"]
  netcrop_status(sim, nsim, "Davis-Kahan", dk.time, best.dk)

  results[[sim]] <- list(
    net = net, nc.out = netcrop.result, nc.time = netcrop.time,
    dk.out = dk.result, dk.time = dk.time,
    n = n, K = K, K.max = K.max, p.test = p.test
  )
  saveRDS(results[[sim]], file.path(LOG_DIR,
    paste0("small_figure2_sim", sim, ".rds")))
  saveRDS(results, results.file)
}

message("Small Figure 2 test completed. Results: ", normalizePath(OUTPUT_DIR))

################################################################################
list.all <- readRDS(results.file)
## plotting
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(RColorBrewer)

showtext::showtext_auto()

final.out <- lapply(list.all, `[[`, "nc.out")
dk.out <- lapply(list.all, `[[`, "dk.out")

mat.out <- list()
for(ii in 1:length(final.out)){
  mat.out[[ii]] <- 100*c(
    final.out[[ii]]$all.accu[1],
    max(final.out[[ii]]$all.accu),
    dk.out[[ii]][,'accu'][which.min(dk.out[[ii]][,'DK.stat'])],
    # final.out[[ii]]$all.accu[which.min(dk.out[[ii]][,'DK.stat'])],
    final.out[[ii]]$croissant.all.accu["l2"],
    final.out[[ii]]$croissant.all.accu["l2.mean"],
    final.out[[ii]]$croissant.all.accu["l2.mode"],
    final.out[[ii]]$croissant.all.accu["bin.dev"],
    final.out[[ii]]$croissant.all.accu["bin.dev.mean"],
    final.out[[ii]]$croissant.all.accu["bin.dev.mode"],
    final.out[[ii]]$croissant.all.accu["AUC"],
    final.out[[ii]]$croissant.all.accu["AUC.mean"],
    final.out[[ii]]$croissant.all.accu["AUC.mode"],
    final.out[[ii]]$croissant.all.accu["pair.NMI.loss"],
    final.out[[ii]]$croissant.all.accu["pair.NMI.loss.mean"],
    final.out[[ii]]$croissant.all.accu["pair.NMI.loss.mode"],
    final.out[[ii]]$croissant.all.accu["pair.hemming.loss"],
    final.out[[ii]]$croissant.all.accu["pair.hemming.loss.mean"],
    final.out[[ii]]$croissant.all.accu["pair.hemming.loss.mode"]
  )
  
  names(mat.out[[ii]]) <- c("0", "Oracle",
                            "Davis-Kahan Estimator",
                            "NETCROP(l2)", "NETCROP(l2)-Mean", "NETCROP(l2)-Mode",
                            "NETCROP(bd)", "NETCROP(bd)-Mean", "NETCROP(bd)-Mode",
                            "NETCROP(AUC)", "NETCROP(AUC)-Mean", "NETCROP(AUC)-Mode",
                            "NETCROP(NMI)", "NETCROP(NMI)-Mean", "NETCROP(NMI)-Mode",
                            "NETCROP(Hemming)", "NETCROP(Hemming)-Mean", "NETCROP(Hemming)-Mode"
  )
  
}


mat.all <- do.call(cbind, mat.out)

mat.median <- apply(mat.all, 1, mean)
mat.med.sd <- apply(mat.all, 1, sd)

plot.out <- tibble(
  tau = factor(rownames(mat.all),
               levels = c("0", "Oracle",
                          "Davis-Kahan Estimator",
                          "NETCROP(l2)", "NETCROP(l2)-Mean", "NETCROP(l2)-Mode",
                          "NETCROP(bd)", "NETCROP(bd)-Mean",
                          "NETCROP(bd)-Mode",
                          "NETCROP(AUC)", "NETCROP(AUC)-Mean", "NETCROP(AUC)-Mode",
                          "NETCROP(NMI)", "NETCROP(NMI)-Mean",
                          "NETCROP(NMI)-Mode", "NETCROP(Hemming)",
                          "NETCROP(Hemming)-Mean", "NETCROP(Hemming)-Mode")),
  accuracy = mat.median,
  sd = mat.med.sd,
  xxmin = mat.median - mat.med.sd,
  xxmax = mat.median + mat.med.sd
)


################################################################################
## only 0, oracle, l2
plot.l2 <- plot.out |> filter(
  tau %notin% tau[grepl("NMI|Hemming|bd|AUC", tau)]
) |>
  mutate(
    tau = factor(tau,
                 levels = c("0", "Oracle",
                            "Davis-Kahan Estimator",
                            "NETCROP(l2)", "NETCROP(l2)-Mean", "NETCROP(l2)-Mode"
                 )
    ))

plot.l2$group <- factor(
  if_else(
    plot.l2$tau %in% c("0", "Oracle", "Davis-Kahan Estimator"),
    "Group 1", "Group 2"
  )
)

color.l2 <- plot.l2 |>
  ggplot(aes(x = accuracy, y = tau,
             color = group
  )) +  # Map color to the 'group' variable
  geom_point() +
  geom_errorbar(aes(xmin = xxmin, xmax = xxmax), orientation = "y") +
  geom_vline(xintercept = mat.median["Oracle"], linetype = "dashed") +
  # scale_x_continuous(limits = c(70, 100)) +
  scale_y_discrete(
    limits = rev(levels(plot.l2$tau)),
    labels = c(
      "0" = "0",
      "Oracle" = "Oracle",
      "Davis-Kahan Estimator" = "Davis-Kahan Estimator",
      "NETCROP(l2)" = expression(paste("NETCROP(", italic(l[2]), ")")),
      "NETCROP(l2)-Mean" = expression(paste("NETCROP(", italic(l[2]), ")-Mean")),
      "NETCROP(l2)-Mode" = expression(paste("NETCROP(", italic(l[2]), ")-Mode"))
      # "NETCROP(bd)" = "NETCROP(bd)",
      # "NETCROP(bd)-Mean" = "NETCROP(bd)-Mean",
      # "NETCROP(bd)-Mode" = "NETCROP(bd)-Mode",
      # "NETCROP(AUC)" = "NETCROP(AUC)",
      # "NETCROP(AUC)-Mean" = "NETCROP(AUC)-Mean",
      # "NETCROP(AUC)-Mode" = "NETCROP(AUC)-Mode"
    )) +
  scale_color_manual(
    values = c("Group 1" = "blue",  # Color for Group 1 (0, Oracle, Davis-Kahan)
               "Group 2" = "red"  # Color for Group 2 (l2)
    )) +
  labs(
    # title = "Clustering Accuracy of Regularized Spectral Clustering",
    x = expression("Clustering Accuracy (%) [ Mean \u00b1 SD]"),
    y = expression(tau)
    # y = quote(tau)
  ) +
  cowplot::theme_minimal_hgrid() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5, size = 16, color = "black"),
    axis.text.y = element_text(family = "sans", size = 13,
                               face = "bold",
                               color = c(
                                 rep("darkred", 3),
                                 rep("darkblue", 3)
                               )
    ),  # Monospace font with clean styling
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black"),
    # axis.title.y = element_text(),
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
    legend.title = element_blank(),  # Remove legend title for a cleaner look
    legend.position = "none",  # Place the legend on top
    # panel.grid.major = element_blank(),  # Remove major gridlines for a cleaner look
    panel.grid.minor = element_blank(),  # Remove minor gridlines for a cleaner look
    plot.margin = margin(10, 10, 10, 10)  # Adjust plot margin for a better layout
  )

print(color.l2)