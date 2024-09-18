#' Forest plot
#' 
#' Create table with forest plot embedded; used only to plot HR for pollutants
#'
#' @author Ensor Palacios, \email{ensorrafael.palacios@@bristol.ac.uk}
#'
#' ----------------------------------------------------------------------------

# Load packages
library(tidyverse)
library(patchwork)
library(RColorBrewer)

plot_forest <- function(df_data, save_name) {
    # From https://www.khstats.com/blog/forest-plots/#just-the-code
    # Order pollutants
    level_order <- fct_rev(df_data[['explanatory']])
    df_data[['explanatory']] <- factor(df_data[['explanatory']], levels = level_order)

    # Prepare data
    df_data <- df_data %>%
        mutate( # Create color groups
            colormap = case_when(grepl('pm25', explanatory) ~ 1,
                grepl('pmcoarse', explanatory) ~ 2,
                grepl('pmabs', explanatory) ~ 3,
                grepl('pm10', explanatory) ~ 4,
                grepl('no2', explanatory) ~ 5,
                grepl('no', explanatory) ~ 6,
                grepl('appc', explanatory) ~ 7,
                grepl('all', explanatory) ~ 8),
            colormap = map(colormap, function (x) brewer.pal(8, name = 'RdYlBu')[x]) %>% unlist,
            # Create background rectangles (not used anymore)
            position = explanatory %>% as.numeric,
            top = position + 0.5,
            bottom = position - 0.5,
            # Recode p-values
            p = case_when(
                p < 0.001 ~ '<0.001',
                round(p, 2) == .05 ~ as.character(round(p, 3)),
                p < .01 ~ str_pad(
                    as.character(round(p, 3)),
                    width = 4,
                    pad = '0',
                    side = 'right'
                ),
                round(p, 2) == 1 ~ as.character('1.00'),
                TRUE ~ str_pad(
                    as.character(round(p, 2)),
                    width = 4,
                    pad = '0',
                    side = 'right')
            ),
            .keep = 'all')

    # Rename air pollutants
    new_names <- c('"PM"[25]*" iqr"', '"PM"[25]*" 2q" ', '"PM"[25]*" 3q"', '"PM"[25]*" 4q"',
        '"PM"[coarse]*" iqr"', '"PM"[coarse]*" 2q"', '"PM"[coarse]*" 3q"', '"PM"[coarse]*" 4q"',
        '"PM"[abs]*" iqr"', '"PM"[abs]*" 2q"', '"PM"[abs]*" 3q"', '"PM"[abs]*" 4q"',
        '"PM"["10"]*" iqr"', '"PM"[10]*" 2q"', '"PM"[10]*" 3q"', '"PM"[10]*" 4q"',
        '"NO"[2]*" iqr"', '"NO"[2]*" 2q"', '"NO"[2]*" 3q"', '"NO"[2]*" 4q"',
        '"NO iqr"', '"NO 2q"', '"NO 3q"', '"NO 4q"',
        '"pollution score"'
    ) %>% as.expression
    levels(df_data$explanatory) <- new_names

    # Get height plot
    height <- df_data %>% dim %>% .[1] + 1
    
    # Plot pollutants with corresponding HR (CI)
    plot_left <- df_data %>% 
        ggplot(aes(y = fct_rev(explanatory))) +
#        geom_rect(aes(xmin = -Inf, 
#                      xmax = +Inf, 
#                      ymin = bottom, 
#                      ymax = top, 
#                      fill = rev(colormap)), 
#                  alpha = 0.2,
#                  show.legend = F) +
        geom_text(aes(x = 0, label = explanatory), hjust = 0, parse=T) +
        geom_text(aes(x = 1, label = str_glue('{round(HR, 2)} ({round(L95, 2)}-{round(U95, 2)})')), hjust = 0) +
        annotate('text', x = 0.2, y =  height , label = 'Pollutant', fontface = 'bold') + 
        annotate('text', x = 1.3, y =  height, label = 'Hazard ratio (CI)', fontface = 'bold') + 
        coord_cartesian(xlim = c(0, 3),
                        ylim = c(0, height + 1)) +
        theme_void()
 
    # Forest plot
    xlim <- c(min(0.9, min(df_data$L95)), max(1.4, max(df_data$U95)))
    plot_mid <- df_data %>% 
        ggplot(aes(y = fct_rev(explanatory))) +
#        geom_rect(aes(xmin = -Inf, 
#                      xmax = +Inf, 
#                      ymin = bottom, 
#                      ymax = top, 
#                      fill = rev(colormap)), 
#                  alpha = 0.00) +
        geom_text(aes(x= 0, label =  explanatory), hjust = 0, fontface = 'bold') +
        theme_classic() +
        geom_point(aes(x=HR, color = colormap), shape=16, size = 3) +
        geom_linerange(aes(xmin=L95, xmax=U95, color = colormap), linewidth = 1) +
        geom_vline(xintercept = 1, linetype="dashed") +
        coord_cartesian(xlim = c(xlim[1], xlim[2]),
                        ylim = c(0, height + 1)) +
        theme(axis.line.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none") +
        labs(x = "HR")

    # Plot p-values 
    plot_right <- df_data %>% 
        ggplot() + 
#        geom_rect(aes(xmin = -Inf, 
#                      xmax = +Inf, 
#                      ymin = bottom, 
#                      ymax = top, 
#                      fill = rev(colormap)), 
#                  alpha = 0.2,
#                  show.legend = F) +
        coord_cartesian(ylim = c(0, height + 1)) +
        geom_text(aes(x = 0, 
                      y = fct_rev(explanatory), 
                      label = p)) +
        annotate('text', x = 0, y = height, label = 'p-value', fontface = 'bold') + 
        theme_void()

    # Create composite plot 
    layout <- c(
                area(t = 0, l = 0, b =30, r = 3),
                area(t = 0, l = 3, b =30, r = 4),
                area(t = 0, l = 5, b =30, r = 5)
    )
    plot_left + plot_mid + plot_right + plot_layout(design = layout)

    # Export plot
    ggsave(str_glue('{save_name}_forest.eps'), width=9, height=6)

#    cairo_ps(str_glue('{save_name}_forest.eps'), width=7, heigh=7)
#    plot_left + plot_mid + plot_right + plot_layout(design = layout)
#    dev.off()
}
