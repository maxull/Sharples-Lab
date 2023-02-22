dxa %>%
        ungroup() %>%
        pivot_longer(cols = age:LBM, names_to = "variable", values_to = "value") %>%
        group_by(sex, include, variable) %>%
        summarise(m = mean(value), 
                  s = sd(value)) %>%
        ungroup() %>%
        mutate(summary = paste0(sprintf("%.1f", m),
                                " (",
                                sprintf("%.1f", s), 
                                ")")) %>%
        select(sex, include, variable, summary) %>%
        pivot_wider(id_cols = variable, names_from = c(sex, include), values_from = summary ) %>%
        # sort the rows -- create a factor level
        mutate(variable = factor(variable, levels = c("age", "height", "weight", "LBM"))) %>%
        
        mutate(Variable = c("Age (yrs)", "Stature (cm)", "Body mass (kg)", "Lean body mass (%)")) %>%
        
        arrange(variable) %>%
        
        select(Variable, female_excl:male_incl) %>%
        kable(format = "html", col.names = c(" ", 
                                             "Excluded",
                                             "Included", 
                                             "Excluded", 
                                             "Included"), 
              caption = "Table 1. Participant characteristics") %>%
        add_header_above(c(" " = 1, "Female" = 2, "Male" = 2)) %>% 
        kable_classic(full_width = F, html_font = "Cambria")
