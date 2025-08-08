library(ggplot2)
library(ggplot2)
library(forcats)
library(patchwork)

# Forest plot function for PRS analysis (used in my thesis)
# 
# Creates a two-panel forest plot:
# - Upper panel: OR for standardized PRS (per SE)
# - Lower panel: OR for quintile comparison (e.g., top 20% vs rest)
# 
# data: Data frame must contain the following columns: term, var, OR, CI_l, CI_u, p
# xmin: Minimum x-axis limit
# xmax: Maximum x-axis limit  
# title: Plot title
# std_prs: Term name for standardized PRS (default: "std_prs")
# ntile_comparison: Term name for quintile comparison (default: "q5_vs_rest")


plot_ics <- function(data, xmin, xmax, title, 
                     std_prs = "std_prs", 
                     ntile_comparison = "q5_vs_rest") {

  # Upper panel: Standardized PRS results
  data_stprs <- subset(data, term == std_prs)
  
  p1 <- ggplot(data_stprs, aes(y = fct_rev(factor(var)))) +
    geom_point(aes(x = OR), shape = 15, size = 3, color = "black") +
    geom_text(aes(x = OR, label = round(OR, 2)), 
              hjust = 0.5, vjust = -0.7, size = 3.5) +
    geom_linerange(aes(xmin = CI_l, xmax = CI_u), size = 0.5) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(x = "OR (95% IC) por cada SE", y = "") +
    theme_classic() +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11)
    ) +
    coord_cartesian(xlim = c(xmin, xmax))
  
  # Variable names panel
  p2 <- ggplot(data_stprs, aes(y = fct_rev(factor(var)))) +
    ggtitle(title) +
    geom_text(aes(x = 0, label = var), 
              hjust = 0, fontface = "bold", family = "sans", size = 4) +
    labs(y = "") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0, family = "sans", color = "black", 
                               size = 12, margin = margin(b = 10))
    )
  
  # P-values panel
  p3 <- ggplot(data_stprs, aes(y = fct_rev(factor(var)))) +
    geom_text(aes(x = 0, label = p), 
              hjust = 0, family = "sans", size = 3.5) +
    theme_void()
  
  # Layout for upper panel
  layout1 <- c(
    area(t = 0, l = 0, b = 30, r = 3),   # Variable names
    area(t = 0, l = 4, b = 30, r = 10),  # Forest plot
    area(t = 0, l = 11, b = 30, r = 12)  # P-values
  )
  
  plot_upper <- p2 + p1 + p3 + plot_layout(design = layout1) +
    plot_annotation(theme = theme(
      panel.background = element_rect(fill = 'transparent', colour = NA),
      plot.background = element_rect(fill = 'transparent', colour = NA)
    ))
  
  # Lower panel: Quintile comparison results
  data_quintile <- subset(data, term == ntile_comparison)
  
  p11 <- ggplot(data_quintile, aes(y = fct_rev(factor(var)))) +
    geom_point(aes(x = OR), shape = 15, size = 3, color = "black") +
    geom_text(aes(x = OR, label = round(OR, 2)), 
              hjust = 0.5, vjust = -0.7, size = 3.5) +
    geom_linerange(aes(xmin = CI_l, xmax = CI_u), size = 0.5) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", alpha = 0.7) +
    labs(x = "OR (95% IC) Q5 vs resto", y = "") +
    theme_classic() +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 10),
      axis.title.x = element_text(size = 11)
    ) +
    coord_cartesian(xlim = c(xmin, xmax))
  
  # Variable names for lower panel
  p22 <- ggplot(data_quintile, aes(y = fct_rev(factor(var)))) +
    geom_text(aes(x = 0, label = var), 
              hjust = 0, fontface = "bold", family = "sans", size = 4) +
    labs(y = "") +
    theme_void()
  
  # P-values for lower panel
  p33 <- ggplot(data_quintile, aes(y = fct_rev(factor(var)))) +
    geom_text(aes(x = 0, label = p), 
              hjust = 0, family = "sans", size = 3.5) +
    theme_void()
  
  # Layout for lower panel
  layout2 <- c(
    area(t = 0, l = 0, b = 30, r = 3),   # Variable names
    area(t = 0, l = 4, b = 30, r = 10),  # Forest plot
    area(t = 0, l = 11, b = 30, r = 12)  # P-values
  )
  
  plot_lower <- p22 + p11 + p33 + plot_layout(design = layout2) +
    plot_annotation(theme = theme(
      panel.background = element_rect(fill = 'transparent', colour = NA),
      plot.background = element_rect(fill = 'transparent', colour = NA)
    ))
  
  # Combine both panels vertically
  final_plot <- plot_upper / plot_lower + 
    plot_layout(heights = c(1, 1)) +
    plot_annotation(theme = theme(
      panel.background = element_rect(fill = 'transparent', colour = NA),
      plot.background = element_rect(fill = 'transparent', colour = NA)
    ))
  
  return(final_plot)
}

# Example:
# plot_ics(data = my_data, xmin = 0.5, xmax = 2.0, title = "PRS associations for the Spanish cohort", std_prs = "std_prs", ntile_comparison = "q5_vs_rest")
