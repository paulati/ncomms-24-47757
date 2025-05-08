# plot_bargraph.R

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



# Function to create a horizontal bar graph
plot_horizontal_bar <- function(data, x_var, y_var, x_label, png_image) {

  img <- png::readPNG(png_image)

plot_1 <-   ggplot(data, aes(x = {{x_var}}, y = reorder({{y_var}}, {{x_var}}))) +
    labs(x = x_label) +
    geom_bar(stat = "identity") +
    scale_x_continuous(labels = scales::number_format(scale = 1, accuracy = 1)) +  # Format x-axis labels
    theme(
      panel.grid.major = element_blank(),  # Remove major gridlines
      panel.grid.minor = element_blank(),  # Remove minor gridlines
      axis.title.y = element_blank(),      # Remove y-axis label
      axis.title.x = element_text(size = 12),  # Change x-axis label
      axis.text.x = element_text(size = 10)  # Adjust x-axis text size
    )

ggpubr::ggarrange(plot_1) + add_phylopic(img,
                                         alpha = 1,
                                         x = 0.90,
                                         y = 0.13,
                                         height = 0.05)


}

