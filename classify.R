library(tidyverse)
results <- read_tsv("3d.txt")

classify_result <- function(num_normal,
                            num_cancer,
                            num_infected,
                            num_resistant,
                            num_empty) {

  # translate to frequencies:
  total_number_cells <- num_normal + num_cancer +
                        num_infected + num_resistant + num_empty
  freq_normal <- num_normal / total_number_cells
  freq_cancer <- (num_cancer + num_resistant) / total_number_cells
  freq_infected <- num_infected / total_number_cells
  freq_resistant <- num_resistant / total_number_cells

  lower_limit <- 1 / total_number_cells

  if (freq_normal == 0 && freq_infected == 0) {
    return("Cancer_wins")
  }
  if (freq_cancer == 0) {
    return("Cancer_eradicated")
  }

  if (freq_cancer > freq_normal && freq_cancer > freq_infected &&
      freq_infected > 0) {
    return("Cancer_winning")
  }
  if (freq_infected > freq_normal && freq_infected > freq_cancer &&
      freq_cancer > 0) {
    return("Virus_winning")
  }
  if (freq_normal > freq_infected && freq_normal > freq_cancer &&
      freq_cancer > 0 && freq_infected > 0) {
    return("Normal_winning")
  }

  # could not match
  return("No_winner")
}

results$outcome_2 <- rep(NA, length(results$outcome))
for(i in seq_along(results$outcome_2)) {
  results$outcome_2[i] <- classify_result(results$num_normal_cells[i],
                                          results$num_cancer_cells[i],
                                          results$num_infected_cells[i],
                                          results$num_resistant_cells[i],
                                          results$num_empty_cells[i])
}
table(results$outcome_2)

# if you want, you can actually replace with labels for a more formal figure:
results$outcome_2[ results$outcome_2 == "Cancer_wins"]       <- "A"
results$outcome_2[ results$outcome_2 == "Cancer_eradicated"] <- "B"
results$outcome_2[ results$outcome_2 == "Virus_winning"]     <- "B-"
results$outcome_2[ results$outcome_2 == "Cancer_winning"]    <- "A-"
results$outcome_2[ results$outcome_2 == "Normal_winning"]    <- "B+"
results$outcome_2[ results$outcome_2 == "No_winner"]         <- "C"


