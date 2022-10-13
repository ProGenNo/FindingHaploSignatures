library(conflicted)
library(glue)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scico)
library(grid)

conflict_prefer("filter", "dplyr")

raw_data_frame <- read.table(
    file = "220824_all_multivar.txt",
    sep = "\t",
    header = T,
    stringsAsFactors = F
)

parse_array <- function(x) {
    return(
        strsplit(
            x = substr(
                x = x, 
                start = 2, 
                stop = nchar(x = x) - 1
            ), 
            split = ", "
        )[[1]]
    )
}

for(i in 1:nrow(raw_data_frame)) {
    raw_data <- raw_data_frame[i, ]

    sequence <- raw_data$peptide[1]
    PSMId <- raw_data$PSMId[1]

    observed_spectrum_df <- data.frame(
        observed_mz = as.numeric(parse_array(raw_data$observed_mz)), 
        observed_intensity = as.numeric(parse_array(raw_data$observed_intensity)), 
        observed_annotation = parse_array(raw_data$observed_annotation),  
        stringsAsFactors = F
    ) %>% 
        mutate(
            observed_annotation = ifelse(
                test = observed_annotation == "None",
                yes = NA,
                no = observed_annotation
            )
        )

    predicted_spectrum_df <- data.frame(
        predicted_mz = as.numeric(parse_array(raw_data$predicted_mz)), 
        predicted_intensity = as.numeric(parse_array(raw_data$predicted_intensity)), 
        predicted_annotation = parse_array(raw_data$predicted_annotation), 
        stringsAsFactors = F
    )

    annotation_threshold_ppm <- 10

    peak_matching_df <- observed_spectrum_df %>% 
        mutate(
            dummy = 0 # Trick: merging on a dummy variable produces all combinations of rows
        ) %>% 
        left_join(
            predicted_spectrum_df %>% 
                mutate(
                    dummy = 0
                ),
            by = "dummy"
        ) %>% 
        mutate(
            error_ppm = 1000000 * (observed_mz - predicted_mz) / predicted_mz
        ) %>% 
        filter(
            abs(error_ppm) <= annotation_threshold_ppm
        ) %>% 
        select(
            -dummy
        )

    median_observed <- median(peak_matching_df$observed_intensity)
    median_predicted <- median(peak_matching_df$predicted_intensity)

    observed_spectrum_df <- observed_spectrum_df %>% 
        mutate(
            scaled_observed_intensity = observed_intensity / median_observed
        )

    predicted_spectrum_df <- predicted_spectrum_df %>% 
        mutate(
            scaled_predicted_intensity = predicted_intensity / median_predicted
        )

    observed_plot_df <- observed_spectrum_df %>% 
    rename(
        mz = observed_mz,
        intensity = observed_intensity,
        annotation = observed_annotation,
        scaled_intensity = scaled_observed_intensity
    ) %>% 
    mutate(
        category = ifelse(
            test = mz %in% peak_matching_df$observed_mz,
            yes = "observed_predicted",
            no = "observed"
        )
    )

    predicted_plot_df <- predicted_spectrum_df %>% 
        rename(
            mz = predicted_mz,
            intensity = predicted_intensity,
            annotation = predicted_annotation,
            scaled_intensity = scaled_predicted_intensity
        ) %>% 
        mutate(
            category = "predicted",
            scaled_intensity = -scaled_intensity
        )

    plot_df <- rbind(observed_plot_df, predicted_plot_df)

    plot_df$category <- factor(plot_df$category, c("observed", "observed_predicted", "predicted"))

    # categoryColors <- c("grey80", scico(n = 2, palette = "cork", begin = 0.2, end = 0.8))
    categoryColors <- c("grey80", "blue3", "forestgreen", "black")

    plot_df <- plot_df %>% 
        arrange(
            category, intensity
        )

    plot <- ggplot() +
        theme_bw(
            base_size = 48
        ) +
        geom_segment(
            data = plot_df,
            mapping = aes(
                x = mz,
                xend = mz,
                y = 0,
                yend = scaled_intensity,
                col = category
            ),
            size = 3.5
        ) +
        geom_hline(
            yintercept = 0,
            col = "black"
        ) +
        geom_text(
            data = predicted_plot_df,
            mapping = aes(
                x = mz,
                y = scaled_intensity - 0.1,
                label = annotation
            ),
            col = categoryColors[4],
            vjust = 1,
            size = 10.5
        ) +
        scale_color_manual(
            values = categoryColors
        ) + 
        scale_x_continuous(
            name = "m/z"
        ) + 
        scale_y_continuous(
            name = "Scaled Intensity"
        ) +
        theme(
            legend.position = "none",
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank()
        )

    png(
        filename = paste(i, sequence, PSMId, ".png", sep="_"),
        width = 1800,
        height = 1200   
    )
    grid.draw(plot)
    dev.off()
}
