require(readr)
require(ggplot)
# if only a single output file is present:
single_output_file <- TRUE
if(single_output_file) {
  to_plot <- read_tsv("output.txt", col_names = F)
  colnames(to_plot) <- c("Birth", "Death", "Type")
}

# if multiple output files in many directories are present:
dir_path <- "/home/p123456/simulation/"
f <- list.files(dir_path, pattern = "output.txt", recursive = TRUE)
for(x in f) {
  file_name <- paste0(dir_path, x)
  if(file.exists(file_name) {
    vx <- read_tsv(file_name, col_names = F)
    colnames(vx) <- c("Birth", "Death", "Type")
    to_plot <- rbind(to_plot, vx)
  }
}
     
ggplot(to_plot, aes(x = Birth, y = Death)) +
  geom_point(aes(col = Type)) +
  scale_color_manual(values = c("blue", "green"))
