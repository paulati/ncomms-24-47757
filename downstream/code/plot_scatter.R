# Load necessary libraries
library(ggthemr)
library(ggpubr)
library(rphylopic)
library(png)

# Define custom color palette
random_colours <- c("#111111", "#65ADC2", "#233B43", "#E84646", "#C29365", "#362C21", "#316675", "#168E7F", "#109B37")
my_palette <- ggthemr::define_palette(
  swatch = random_colours,
  gradient = c(lower = random_colours[1L], upper = random_colours[2L]),
  gridline = "#111111"
)

# Set ggthemr theme
ggthemr::ggthemr(my_palette)

# Function to create a scatter plot with optional PNG image
plot_scatter <- function(data, x_var, y_var, name_var, x_label, y_label, png_image = NULL, number_to_show = 25) {
  # Check if the PNG image path is provided
  if (!is.null(png_image)) {
    img <- png::readPNG(png_image)
    # Create the scatter plot without the PNG image
    scatter_plot <- ggpubr::ggscatter(
      data,
      x = x_var,
      y = y_var,
      add = "reg.line",
      add.params = list(color = my_palette$swatch[1], fill = my_palette$swatch[2]),
      conf.int = TRUE,
      cor.coef = TRUE,
      cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
      label = name_var,
      label.select = list(criteria = paste(y_var, "> ", number_to_show)),
      font.label = list(size = 9, face = "italic"),
      repel = TRUE
    ) +
      ggplot2::ylab(y_label) +
      ggplot2::xlab(x_label)

    # Use ggarrange to arrange the scatter plot and add_phylopic
    scatter_plot <- ggpubr::ggarrange(scatter_plot) +
      rphylopic::add_phylopic(img, alpha = 1, x = 0.90, y = 0.18, height = 0.1)

    return(scatter_plot)
  } else {
    # Create the scatter plot without the PNG image
    scatter_plot <- ggpubr::ggscatter(
      data,
      x = x_var,
      y = y_var,
      add = "reg.line",
      add.params = list(color = my_palette$swatch[1], fill = my_palette$swatch[2]),
      conf.int = TRUE,
      cor.coef = TRUE,
      cor.coeff.args = list(method = "pearson", label.x = 3, label.sep = "\n"),
      label = name_var,
      label.select = list(criteria = paste(y_var, "> ", number_to_show)),
      font.label = list(size = 9, face = "italic"),
      repel = TRUE
    ) +
      ggplot2::ylab(y_label) +
      ggplot2::xlab(x_label)

    return(scatter_plot +  my_palette$theme)
  }
}
