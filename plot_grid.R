require(readr)
require(ggplot)
to_plot <- read_tsv("output.txt", col_names = F)
colnames(to_plot) <- c("Birth", "Death", "Type")
ggplot(to_plot, aes(x = Birth, y = Death)) +
  geom_point(aes(col = Type)) +
  scale_color_manual(values = c("blue", "green"))
