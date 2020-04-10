require(tidyverse)

vx <- read_tsv("output.txt")

p1 <- ggplot(vx, aes(x = birth_virus, y = death_virus, fill = num_cancer_cells)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal()

vx %>%
  gather(key = "statistic", value = "value", c(num_normal_cells, num_cancer_cells, num_infected_cells)) %>%
  ggplot(aes(x = birth_virus, y = death_virus, fill = value)) +
    geom_tile()  +
  scale_fill_viridis_c() +
  theme_minimal() +
  facet_wrap(~statistic)

vx %>%
  ggplot(aes(x = birth_virus, y = death_virus, fill = outcome)) +
  geom_tile()  +
  scale_fill_viridis_d() +
  theme_minimal()


