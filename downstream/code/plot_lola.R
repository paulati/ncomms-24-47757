# Author: Anabella Trigila
# ------


# plot_lola.R

# Install and load ggthemr package
# install.packages("rphylopic")
library(ggthemr)
library(rphylopic)

# Define custom color palette
random_colours <- c("#111111", "#65ADC2", "#233B43", "#E84646", "#C29365", "#362C21", "#316675", "#168E7F", "#109B37")
my_palette <- define_palette(
  swatch = random_colours,
  gradient = c(lower = random_colours[1L], upper = random_colours[2L]),
  gridline = "#111111"
)

# Set ggthemr theme
ggthemr(my_palette)

# Load necessary libraries
library(ggplot2)

# Function to create a horizontal bar graph with LOLA results
library(ggplot2)
library(dplyr)


plot_lola <- function(lola_output, value_col = "pValueLog", img_path = NULL, char_limit = 25) {

  # Convert value_col to a symbol for tidy evaluation
  value_sym <- sym(value_col)

  # Truncate descriptions to the specified character limit
  lola_output <- lola_output %>%
    mutate(description = ifelse(nchar(description) > char_limit,
                                paste0(substr(description, 1, char_limit - 3), "..."),
                                description))

plot <- lola_output %>%
    head(n = 10) %>%
    arrange(!!value_sym) %>%
    ggplot(aes(y = !!value_sym, x = reorder(description, !!value_sym))) +
    geom_col(width = 0.5) +
    coord_flip() +
    xlab("") +
    ylab(paste("", value_col, "", sep = "")) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),  # Remove major gridlines
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      axis.title.y = element_blank(),      # Remove y-axis label
      axis.title.x = element_text(size = 12),  # Change x-axis label
      axis.text.x = element_text(size = 10)  # Adjust x-axis text size
    )

  # If an image path is provided, add the PhyloPic image
  if (!is.null(img_path)) {
    img <- png::readPNG(img_path)
    plot <- ggpubr::ggarrange(plot) + add_phylopic(img,x = 0.90, y = 0.20, height = 0.10)
  }

  return(plot)
}

# Example usage
# plot_lola_mammals(lola.output)

