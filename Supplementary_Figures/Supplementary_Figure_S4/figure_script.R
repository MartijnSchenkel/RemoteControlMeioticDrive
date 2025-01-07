library(tidyverse)
library(viridis)
library(ggpattern)
library(cowplot)

# Administrative setup
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

# Read data
dA <- read.table("data_S4A.txt", T) %>% as_tibble()
dB <- read.table("data_S4B.txt", T) %>% as_tibble()
dC <- read.table("data_S4C.txt", T) %>% as_tibble()
dD <- read.table("data_S4D.txt", T) %>% as_tibble()
dE <- read.table("data_S4E.txt", T) %>% as_tibble()

# Plot A
dA2 <- dA %>% gather(X, Y, A, D, value = "Frequency", key = "Allele")

dA2end <- dA2 %>% filter(G == 500)
dA2end
delta <- 0.055
dtA <- tibble(text = c('X', 'Y',
                       'italic(A)[2]', 'italic(D)[2]'),
              x = c(450, 450, 350, 450), y = c(0.329 - delta, 0.671 + delta, 
                                               0.476 - delta, 0.476 - delta),
              Allele = unique(dA2$Allele))

dtA


pA <- dA2 %>% ggplot(aes(G, Frequency, color = Allele)) + geom_line(linewidth = 0.8) +
  scale_color_viridis(discrete = T, end = 0.8, 
                      labels = c(bquote(italic(A)[2]), 
                                 bquote(italic(D)[2]), 
                                 "X", "Y")) + 
  geom_hline(yintercept = 0.5, color = "black", linetype = 2, linewidth = 0.8) + 
  scale_x_continuous(limits = c(0, 500), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1), expand = c(0,0)) + 
  labs(x = "Generations") +
  geom_text(data = dtA, aes(x, y, color = Allele, label = text), parse = T, inherit.aes = F) +
  
  theme(panel.background = element_rect(fill = "white", color = "black"),
        axis.line = element_line(color = "black"))
pA

# Plot B


dB2 <- dB %>% gather(DA2D2, DA2Y, DD2Y, value = "LD", key = "Haplotype")


dB2end <- dB2 %>% filter(G == 500)

delta <- 0.04

dtB <- tibble(text = c('italic(A)[2]~italic(D)[2]',
                       'italic(A)[2]~Y',
                       'italic(D)[2]~Y'),
              x = c(950, 850, 950), y = c(0.977 - delta, 0.0147 + delta, 0.0147 + delta),
              Haplotype = unique(dB2$Haplotype))

dtB
dB2 <- dB %>% filter(G <= 500) %>% gather(DA2D2, DA2Y, DD2Y, value = "LD", key = "Haplotype")


dB2end <- dB2 %>% filter(G == max(G))
delta <- 0.06
dtB <- tibble(text = c(bquote(italic(Assister-Distorter)),
                       bquote(italic(Assister)*"-XY"),
                       bquote(italic(Distorter)*"-XY")),
              x = c(390,420,420), y = c(0.979- delta, 0.342 + delta, 0.337 - delta),
              Haplotype = unique(dB2$Haplotype))

dtB
pB <- dB2 %>% filter(G <= 500) %>% ggplot(aes(G, LD, color = Haplotype)) + geom_line(linewidth = 0.8) +
  scale_color_viridis(discrete = T, end = 0.8, 
                      labels = c(bquote(italic(A)[2]~italic(D)[2]), 
                                 bquote(italic(A)[2]~"Y"),
                                 bquote(italic(D)[2]~"Y"))) + 
  geom_hline(yintercept = 0, color = "black", linetype = 2, linewidth = 0.8) + 
  scale_x_continuous(limits = c(0, 500), expand = c(0,0)) + 
  scale_y_continuous(limits = c(-0.2, 1), 
                     breaks = seq(-0.2, 1, 0.2), 
                     expand = c(0,0)) + 
  labs(x = "Generations", y = "Linkage Disequilibrium") +
  geom_text(data = dtB, aes(x, y, color = Haplotype, label = text), parse = T, inherit.aes = F) +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        axis.line = element_line(color = "black"))
pB

# Plot C
dC2 <- dC %>% mutate(tot = A1D1X + A2D1X + A1D1Y + A2D1Y + 
                       A1D2X + A2D2X + A1D2Y + A2D2Y)
dC2 <- dC2 %>% mutate(A1D1X = A1D1X / tot,
                      A2D1X = A2D1X / tot,
                      A1D1Y = A1D1Y / tot,
                      A2D1Y = A2D1Y / tot,
                      
                      A1D2X = A1D2X / tot,
                      A2D2X = A2D2X / tot,
                      A1D2Y = A1D2Y / tot,
                      A2D2Y = A2D2Y / tot)

dC3 <- dC2 %>% gather(A1D1X, A2D1X, A1D1Y, A2D1Y,
                      A1D2X, A2D2X, A1D2Y, A2D2Y,
                      value = "Frequency", key = "Haplotype")

dC4 <- dC3 %>% mutate(Haplotype = factor(Haplotype, c("A1D1X", "A2D1X", "A1D1Y", "A2D1Y",
                                                      "A1D1O", "A2D1O", 
                                                      "A1D2X", "A2D2X", "A1D2Y", "A2D2Y",
                                                      "A1D2O", "A2D2O")),
                      Hatch = ifelse(Haplotype %in% c("A1D2X", "A2D2X", "A1D2Y", "A2D2Y",
                                                      "A1D2O", "A2D2O"), 1, 0))
max(dC3$tot)


text <- c(bquote(italic(A)[1]~italic(D)[1]~"X"), 
          bquote(italic(A)[1]~italic(D)[1]~"Y"),
          # bquote(italic(A)[1]~italic(D)[2]~"X"),
          # bquote(italic(A)[1]~italic(D)[2]~"Y"),
          bquote(italic(A)[2]~italic(D)[1]~"X"),
          bquote(italic(A)[2]~italic(D)[1]~"Y"),
          bquote(italic(A)[2]~italic(D)[2]~"X"),
          bquote(italic(A)[2]~italic(D)[2]~"Y"))
yvals <- c(0.85, 0.6, 0.6, 0.20, 0.435, 0.215)
xvals <- c(250, 25, 250, 25, 250, 250)
angle <- c(0, 1, 0, 1, 0, 0)
dtC <- tibble(text, yvals, xvals)

vircol <- viridis_pal()(8)

pC <- dC4 %>% filter(Hatch == 1,G <= 500) %>% ggplot(aes(G, Frequency, fill = Haplotype)) + 
  geom_area(data = dC4 %>% filter(G <= 500), aes(G, Frequency, fill = Haplotype),
            color = "black", alpha = 0.8) +
  geom_area_pattern(aes(pattern = Hatch), 
                    color = 'black', alpha = 0, pattern = "stripe",
                    pattern_angle = 45, pattern_color = "black", 
                    pattern_density = 0.001,
                    pattern_spacing = 0.03) +
  geom_rect(aes(xmin = 210, xmax = 290, ymin = 0.2, ymax = 0.23), fill = "white") +
  geom_rect(aes(xmin = 210, xmax = 290, ymin = 0.42, ymax = 0.45), fill = "white") +
  
  annotate('rect', xmin = 210, xmax = 290, ymin = 0.2, ymax = 0.23, 
           alpha = 0.8, fill = vircol[8]) +
  annotate('rect', xmin = 210, xmax = 290, ymin = 0.42, ymax = 0.45, 
           alpha = 0.8, fill = vircol[6]) +
  
  # geom_rect(aes(xmin = 420, xmax = 500, ymin = 0.113, ymax = 0.115), fill = "white") +
  # annotate('rect', xmin = 420, xmax = 500, ymin = 0.113, ymax = 0.115, 
  #          alpha = 0.8, fill = vircol[6]) +
  # 
  
  geom_area(data = dC4, aes(G, Frequency, fill = Haplotype),
            color = "black", alpha = 0) +
  
  theme(plot.background = element_blank(),
        panel.background = element_blank()) + 
  scale_fill_viridis(discrete = T, 
                     labels = c(bquote(italic(A)[1]~italic(D)[1]~"X"), 
                                bquote(italic(A)[1]~italic(D)[1]~"Y"),
                                bquote(italic(A)[1]~italic(D)[2]~"X"),
                                bquote(italic(A)[1]~italic(D)[2]~"Y"),
                                bquote(italic(A)[2]~italic(D)[1]~"X"),
                                bquote(italic(A)[2]~italic(D)[1]~"Y"),
                                bquote(italic(A)[2]~italic(D)[2]~"X"),
                                bquote(italic(A)[2]~italic(D)[2]~"Y"))) +
  geom_text(data = dtC, aes(x = xvals, y = yvals, label = text, angle = angle*90), 
            inherit.aes = F, parse = T) +
  scale_x_continuous(limits = c(0, 501), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0, 1.001), expand =c(0.0, 0),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
  labs(x = "Generations", y = "Frequency") +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        axis.line = element_line(color = "black"))
# pC

# Plot D
pD <- dD %>% ggplot(aes(k, c, z = DF)) + geom_contour_filled(color = "black", bins=10) + 
  scale_fill_viridis(option = "A", discrete = T, 
                     
                     labels = seq(0, 1, 0.1)) +
  scale_x_continuous(limits = c(0, 0.5), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0.5, 1), expand = c(0, 0),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
  labs(x = bquote(italic(trans)~"drive ("*italic(k)[D]*")"), 
       y = bquote("Cosegregation rate ("*italic('c')*")"),
       fill = bquote("["*italic(D)[2]~"]")) +
  theme(panel.background = element_rect(fill = "white", color = "black"))
pD

# Plot E
pE <- dE %>% ggplot(aes(k, c, z = SR)) + geom_contour_filled(color = "black", bins=10) + 
  scale_fill_viridis(option = "mako", discrete = T, 
                     labels = seq(0.5, 1, 0.05)) +
  scale_x_continuous(limits = c(0, 0.5), expand = c(0,0)) + 
  scale_y_continuous(limits = c(0.5, 1), expand = c(0, 0),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
  labs(x = bquote(italic(trans)~"drive ("*italic(k)[D]*")"), 
       y = bquote("Cosegregation rate ("*italic('c')*")"),
       fill = "Sex ratio") +
  theme(panel.background = element_rect(fill = "white", color = "black"))
pE


legend_D <- get_legend(
  pD + 
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(legend.key.size = unit(0.5, "cm")) +
    theme(legend.key.height = unit(0.35, "cm")) # Adjust size here
)

legend_E <- get_legend(
  pE + 
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(legend.key.size = unit(0.5, "cm"))+
    theme(legend.key.height = unit(0.35, "cm")) # Adjust size here
)

topleft <- plot_grid(pA + theme(legend.position = "none"),
                     pB + theme(legend.position = "none"),
                     ncol = 1, nrow = 2, labels = c("A", "B"))
topright  <-  plot_grid(pC + theme(legend.position = "none"), labels = "C")


top <- plot_grid(topleft, NULL, topright,NULL, ncol = 4, 
                 rel_widths = c(0.425, 0.075, 0.45, 0.075))

bottom <- plot_grid(pD+ theme(legend.position = "none"), legend_D, 
                    pE + theme(legend.position = "none"), legend_E, 
                    ncol = 4, 
                    labels = c("D", "", "E", ""), 
                    rel_widths = c(0.425, 0.075, 0.45, 0.075))
plot <-  plot_grid(top, bottom, ncol = 1, nrow = 2, rel_heights = c(0.6, 0.4))
# plot

# pdf(file = "2024_10_22_Figure_C.pdf", width = 8, height = 10)
# plot
# dev.off()

pdf(file = "2025_01_07_Figure_S4.pdf", width = 10, height = 10)
plot
dev.off()

png(file = "2025_01_07_Figure_S4.png", width = 10, height = 10, units = "in", res = 1000)
plot
dev.off()
