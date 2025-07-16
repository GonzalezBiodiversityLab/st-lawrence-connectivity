library(ggplot2)
library(dplyr)
library(readr)
library(RColorBrewer)
library(ggpubr)
library(tidyr)

# Read percent change summary (from your generated data)
summary_stats_con <- read_csv("./Circuitscape/outputs/percent_change_summary.csv") %>%
  mutate(Metric = "Connectivity")
summary_stats_hab <- read_csv("./habitat/percent_change_summary.csv") %>%
  mutate(Metric = "Habitat")


# Add Metric and 2010 baseline value
summary_stats <- summary_stats_con %>%
  rbind(summary_stats_hab) %>%
  mutate(
#    Metric = "Connectivity",
    `2010` = 0  # Always zero for baseline
  ) %>%
  rename(
    Scenario = scenario,
    Species = species,
    `2110` = mean_percent_change,
    SE = se_percent_change
  ) %>%
  select(Metric, Scenario, Species, `2010`, `2110`, SE)


# Set plotting theme
theme_set(
  theme_bw() +
    theme(
      aspect.ratio = 1,
      panel.spacing = unit(4, "lines"),
      panel.border = element_rect(linewidth = 2),
      strip.text.x.top = element_blank(),
      axis.title = element_text(size = 50),
      axis.text = element_text(size = 50),
      legend.position = "none",
      legend.text = element_text(size = 50),
      axis.ticks = element_line(linewidth = 1.5),
      axis.ticks.length = unit(.75, "cm"),
      plot.margin = unit(c(0, 0, 0, 0), "lines")
    )
)


# Reshape into long format for plotting
perc_chng_hab <- summary_stats %>%
  filter(Metric == 'Habitat') %>%
  pivot_longer(cols = c(`2010`, `2110`), names_to = "Year", values_to = "PctChng") %>%
  mutate(
    Scenario = factor(Scenario, levels = c("NCNC", "NC45", "NC85", "BAUNC", "BAU45", "BAU85")),
    Year = factor(Year, levels = c("2010", "2110")),
    SE = ifelse(Year == "2110", SE, 0)  # No error bar at baseline
  )


#scn_colours <- c(brewer.pal(n = 6, name = 'RdBu')[4:6], brewer.pal(n = 6, name = 'RdBu')[3:1])


# Plot percent habitat change
# Define nudges for points and error bars by scenario
nudge_points <- tibble::tribble(
  ~Scenario, ~Species, ~point_nudge,
  "NCNC",  "BLBR",  0.00,
  "NC45",  "BLBR",  0.00,
  "NC85",  "BLBR",  0.00,
  "BAUNC", "BLBR",  0.00,
  "BAU45", "BLBR",  0.00,
  "BAU85", "BLBR",  0.00,
  
  "NCNC",  "MAAM",  0.00,
  "NC45",  "MAAM",  0.00,
  "NC85",  "MAAM",  0.00,
  "BAUNC", "MAAM",  0.00,
  "BAU45", "MAAM",  0.00,
  "BAU85", "MAAM",  0.00,
  
  "NCNC",  "PLCI",  0.00,
  "NC45",  "PLCI",  0.00,
  "NC85",  "PLCI",  0.00,
  "BAUNC", "PLCI",  0.00,
  "BAU45", "PLCI",  0.00,
  "BAU85", "PLCI",  0.00,
  
  "NCNC",  "RANA",  0.00,
  "NC45",  "RANA",  0.03,
  "NC85",  "RANA",  0.06,
  "BAUNC", "RANA",  0.00,
  "BAU45", "RANA",  0.03,
  "BAU85", "RANA",  0.06,
  
  "NCNC",  "URAM",  0.00,
  "NC45",  "URAM",  0.03,
  "NC85",  "URAM",  0.06,
  "BAUNC", "URAM",  0.00,
  "BAU45", "URAM",  0.03,
  "BAU85", "URAM",  0.06
)

nudge_ci <- tibble::tribble(
  ~Scenario, ~Species, ~ci_nudge,
  "NCNC",  "BLBR",  0.00,
  "NC45",  "BLBR",  0.06,
  "NC85",  "BLBR",  0.00,
  "BAUNC", "BLBR",  0.00,
  "BAU45", "BLBR",  0.06,
  "BAU85", "BLBR",  0.00,
  
  "NCNC",  "MAAM",  0.00,
  "NC45",  "MAAM",  0.00,
  "NC85",  "MAAM",  0.00,
  "BAUNC", "MAAM",  0.00,
  "BAU45", "MAAM",  0.00,
  "BAU85", "MAAM",  0.00,
  
  "NCNC",  "PLCI",  0.00,
  "NC45",  "PLCI",  0.06,
  "NC85",  "PLCI",  0.00,
  "BAUNC", "PLCI",  0.00,
  "BAU45", "PLCI",  0.06,
  "BAU85", "PLCI",  0.00,
  
  "NCNC",  "RANA",  0.00,
  "NC45",  "RANA",  0.03,
  "NC85",  "RANA",  0.06,
  "BAUNC", "RANA",  0.00,
  "BAU45", "RANA",  0.03,
  "BAU85", "RANA",  0.06,
  
  "NCNC",  "URAM",  0.00,
  "NC45",  "URAM",  0.03,
  "NC85",  "URAM",  0.06,
  "BAUNC", "URAM",  0.00,
  "BAU45", "URAM",  0.03,
  "BAU85", "URAM",  0.06
)


# Create a column with the x-axis nudge values
perc_chng_hab <- perc_chng_hab %>%
  left_join(nudge_points, by = c("Scenario", "Species")) %>%
  left_join(nudge_ci, by = c("Scenario", "Species")) %>%
  mutate(
    point_nudge = ifelse(is.na(point_nudge), 0, point_nudge),
    ci_nudge = ifelse(is.na(ci_nudge), 0, ci_nudge)
  ) %>%
  mutate(
    x_pos_point = as.numeric(factor(Year, levels = c("2010", "2110"))) + ifelse(Year == "2110", point_nudge, 0),
    x_pos_ci = as.numeric(factor(Year, levels = c("2010", "2110"))) + ifelse(Year == "2110", ci_nudge, 0)
  )


# Ensure Scenario is a factor with consistent levels and correct order
scenario_levels <- c("NCNC", "NC45", "NC85", "BAUNC", "BAU45", "BAU85")

perc_chng_hab <- perc_chng_hab %>%
  mutate(
    Scenario = factor(Scenario, levels = scenario_levels)
  )

# Make sure your colours align with these levels
scn_colours <- c(
  brewer.pal(n = 6, name = 'RdBu')[4:6],  # NC scenarios
  brewer.pal(n = 6, name = 'RdBu')[3:1]   # BAU scenarios
)
names(scn_colours) <- scenario_levels

# Plot
p1 <- ggplot(data = perc_chng_hab, aes(x = Year, y = PctChng, color = Scenario, group = Scenario)) +
  geom_line(linewidth = 3) +
  geom_point(aes(x = x_pos_point), size = 7) +
  geom_errorbar(
    data = filter(perc_chng_hab, Year == "2110"),
    aes(
      x = x_pos_ci,
      ymin = PctChng - 1.96 * SE,
      ymax = PctChng + 1.96 * SE,
      color = Scenario
    ),
    inherit.aes = FALSE,
    width = 0.05,
    linewidth = 1.5
  ) +
  scale_color_manual(values = scn_colours) +
  facet_grid(. ~ Species) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(x = "Year", y = "Change in habitat area (%)") +
  scale_x_discrete(expand = expansion(mult = 0.25)) +
  theme(
    axis.title = element_text(size = 50),
    axis.text = element_text(size = 50),
    legend.text = element_text(size = 50),
    axis.ticks = element_line(linewidth = 1.5),
    axis.ticks.length = unit(.75, "cm"),
    plot.margin = unit(c(0, 0, 0, 0), "lines")
  )







# Create plot

# Reshape into long format for plotting
perc_chng_connect <- summary_stats %>%
  filter(Metric == 'Connectivity') %>%
  pivot_longer(cols = c(`2010`, `2110`), names_to = "Year", values_to = "PctChng") %>%
  mutate(
    Scenario = factor(Scenario, levels = c("NCNC", "NC45", "NC85", "BAUNC", "BAU45", "BAU85")),
    Year = factor(Year, levels = c("2010", "2110")),
    SE = ifelse(Year == "2110", SE, 0)  # No error bar at baseline
  )


# Define nudges for points and error bars by scenario
nudge_points <- tibble::tribble(
  ~Scenario, ~Species, ~point_nudge,
  "NCNC",  "BLBR",  0.00,
  "NC45",  "BLBR",  0.00,
  "NC85",  "BLBR",  0.00,
  "BAUNC", "BLBR",  0.00,
  "BAU45", "BLBR",  0.00,
  "BAU85", "BLBR",  0.00,
  
  "NCNC",  "MAAM",  0.00,
  "NC45",  "MAAM",  0.00,
  "NC85",  "MAAM",  0.00,
  "BAUNC", "MAAM",  0.00,
  "BAU45", "MAAM",  0.00,
  "BAU85", "MAAM",  0.00,
  
  "NCNC",  "PLCI",  0.00,
  "NC45",  "PLCI",  0.00,
  "NC85",  "PLCI",  0.00,
  "BAUNC", "PLCI",  0.00,
  "BAU45", "PLCI",  0.00,
  "BAU85", "PLCI",  0.00,
  
  "NCNC",  "RANA",  0.00,
  "NC45",  "RANA",  0.03,
  "NC85",  "RANA",  0.06,
  "BAUNC", "RANA",  0.00,
  "BAU45", "RANA",  0.03,
  "BAU85", "RANA",  0.06,
  
  "NCNC",  "URAM",  0.00,
  "NC45",  "URAM",  0.03,
  "NC85",  "URAM",  0.06,
  "BAUNC", "URAM",  0.00,
  "BAU45", "URAM",  0.03,
  "BAU85", "URAM",  0.06
)

nudge_ci <- tibble::tribble(
  ~Scenario, ~Species, ~ci_nudge,
  "NCNC",  "BLBR",  0.00,
  "NC45",  "BLBR",  0.06,
  "NC85",  "BLBR",  0.00,
  "BAUNC", "BLBR",  0.00,
  "BAU45", "BLBR",  0.06,
  "BAU85", "BLBR",  0.00,
  
  "NCNC",  "MAAM",  0.00,
  "NC45",  "MAAM",  0.06,
  "NC85",  "MAAM",  0.00,
  "BAUNC", "MAAM",  0.00,
  "BAU45", "MAAM",  0.06,
  "BAU85", "MAAM",  0.00,
  
  "NCNC",  "PLCI",  0.00,
  "NC45",  "PLCI",  0.06,
  "NC85",  "PLCI",  0.00,
  "BAUNC", "PLCI",  0.00,
  "BAU45", "PLCI",  0.06,
  "BAU85", "PLCI",  0.00,
  
  "NCNC",  "RANA",  0.00,
  "NC45",  "RANA",  0.03,
  "NC85",  "RANA",  0.06,
  "BAUNC", "RANA",  0.00,
  "BAU45", "RANA",  0.03,
  "BAU85", "RANA",  0.06,
  
  "NCNC",  "URAM",  0.00,
  "NC45",  "URAM",  0.03,
  "NC85",  "URAM",  0.06,
  "BAUNC", "URAM",  0.00,
  "BAU45", "URAM",  0.03,
  "BAU85", "URAM",  0.06
)


# Create a column with the x-axis nudge values
# perc_chng_connect <- perc_chng_connect %>%
#   left_join(nudge_table, by = c("Scenario", "Species")) %>%
#   mutate(
#     x_nudge = ifelse(Year == "2110", x_nudge, 0),
#     x_pos = as.numeric(factor(Year, levels = c("2010", "2110"))) + x_nudge
#   )

perc_chng_connect <- perc_chng_connect %>%
  left_join(nudge_points, by = c("Scenario", "Species")) %>%
  left_join(nudge_ci, by = c("Scenario", "Species")) %>%
  mutate(
    point_nudge = ifelse(is.na(point_nudge), 0, point_nudge),
    ci_nudge = ifelse(is.na(ci_nudge), 0, ci_nudge)
  ) %>%
  mutate(
    x_pos_point = as.numeric(factor(Year, levels = c("2010", "2110"))) + ifelse(Year == "2110", point_nudge, 0),
    x_pos_ci = as.numeric(factor(Year, levels = c("2010", "2110"))) + ifelse(Year == "2110", ci_nudge, 0)
  )


# Ensure Scenario is a factor with consistent levels and correct order
scenario_levels <- c("NCNC", "NC45", "NC85", "BAUNC", "BAU45", "BAU85")

perc_chng_connect <- perc_chng_connect %>%
  mutate(
    Scenario = factor(Scenario, levels = scenario_levels)
  )

# Make sure your colours align with these levels
scn_colours <- c(
  brewer.pal(n = 6, name = 'RdBu')[4:6],  # NC scenarios
  brewer.pal(n = 6, name = 'RdBu')[3:1]   # BAU scenarios
)
names(scn_colours) <- scenario_levels

# Plot
p2 <- ggplot(data = perc_chng_connect, aes(x = Year, y = PctChng, color = Scenario, group = Scenario)) +
  geom_line(linewidth = 3) +
  geom_point(aes(x = x_pos_point), size = 7) +
  geom_errorbar(
    data = filter(perc_chng_connect, Year == "2110"),
    aes(
      x = x_pos_ci,
      ymin = PctChng - 1.96 * SE,
      ymax = PctChng + 1.96 * SE,
      color = Scenario
    ),
    inherit.aes = FALSE,
    width = 0.05,
    linewidth = 1.5
  ) +
  scale_color_manual(values = scn_colours) +
  facet_grid(. ~ Species) +
  guides(colour = guide_legend(nrow = 1)) +
  labs(x = "Year", y = "Change in habitat connectivity (%)") +
  scale_x_discrete(expand = expansion(mult = 0.25)) +
  theme(
    axis.title = element_text(size = 50),
    axis.text = element_text(size = 50),
    legend.text = element_text(size = 50),
    axis.ticks = element_line(linewidth = 1.5),
    axis.ticks.length = unit(.75, "cm"),
    plot.margin = unit(c(0, 0, 0, 0), "lines")
  )


figure <- ggarrange(p1, p2, nrow = 2,
                    common.legend = TRUE,
                    legend = "bottom",
                    align = "v"
)


#figure

ggexport(figure, 
         filename = "ggfigure3.png",
         width = 3600,
         height = 1800
)
