library(tidyverse)
library(knitr)
library(kableExtra)

download.file("https://www.dropbox.com/s/g2t97j8edqvvktn/tr003_dxa.csv?raw=1", 
              destfile = "dxa_data.csv")
dxa_data <- read_csv("dxa_data.csv")

#---------------

#We have some DXA data, here we wish to present lean body mass as a percentage of the whole body mass
# LBM = (LEAN/WHOLE)*100 =>

dxa <- dxa_data %>%
        rowwise()%>%
        mutate(LBM = (lean.whole/(fat.whole+BMD.whole+lean.whole))*100)%>%
        select(subject, age, height, weight, LBM, sex, include)%>%
                print()

#When making a summary table, we first need to format out data suitable for tables:
#therefore we need summary data such as means and SD:

dxa %>%
        ungroup()%>%
        dplyr::group_by(m.age = mean(age), 
                 s.age = sd(age),
                 m.height = mean(height),
                 s.height = sd(height),
                 m.weight = mean(weight),
                 s.weight = sd(weight),
                 m.lbm = mean(LBM), 
                 s.lbm = sd(LBM))%>%
        print()

#---- The above script is somewhat long and untidy. Below is a more efficient way to summarize:

summary_table <- dxa %>%
        ungroup()%>%
                        pivot_longer(cols = age:LBM, names_to = "variable", values_to = "value")%>%
                        group_by(sex, include, variable)%>%
        summarise(m = mean(value), 
                  s = sd(value))%>%
        print()

# We essentially just take our main variables and allocate them to have values in their own column for easier calculation


summary_table <- dxa %>%
        ungroup()%>%
        pivot_longer(cols = age:LBM, names_to = "variable", values_to = "value")%>%
        group_by(sex, include, variable)%>%
        summarise(m = mean(value), 
                  s = sd(value))%>%
        #next we want our numbers to be represented by ONE digit [sprintf] with both m and s in the same cell [paste0]:
        ungroup()%>%
        mutate(summary = paste0(sprintf("%.1f", m), 
                                "(", 
                                sprintf("%.1f", s),
                                ")")) %>%
        
       #next we tidy so that the variables: sex and include are numerical in relation to the variable column:
        
       select(sex, include, variable, summary)%>%
       pivot_wider(id_cols = variable, names_from = c(sex, include), values_from = summary)%>%
               print()
        
        # This gave us a dataset whereby the variable column remain with sex and include merged into their own variables
        # With the values containing the mean and SD as such: "mean(SD)"  
        # we did not, however, get numbers represented by one digit, rather they have one decimal..

#--------


```{r my_table, results ="asis"}
summary_table %>%
        kable(col.names = c("Variable", 
                                                "Female excluded",
                                                "Female included", 
                                                "Male excluded", 
                                                "Male included"), 
                                        caption = "Participant characteristics")
```
#This formats so that we get more than pure html code and actually a table

#----

summary_table %>%
        kable(format = "html", col.names = c(" ", 
                                             "Excluded",
                                             "Included", 
                                             "Excluded", 
                                             "Included"), 
              caption = "Participant characteristics")%>%
        add_header_above(c(" " = 1, "Female" = 2, "Male" = 2))

#Now we can format this table to look more professional:
#First we manage the date similar as done previously:

dxa %>%
        ungroup() %>%
        pivot_longer(cols = age:LBM, names_to = "variable", values_to = "value") %>%
        group_by(sex, include, variable) %>%
        summarise(m = mean(value), 
                  s = sd(value)) %>%
        ungroup() %>%
        mutate(summary = paste0(sprintf("%.1f", m),
                                " \u00b1 ",
                                                        #\u00b1 is the code for +/- symbol
                                sprintf("%.1f", s), 
                                "")) %>%
        select(sex, include, variable, summary) %>%
        pivot_wider(id_cols = variable, names_from = c(sex, include), values_from = summary )%>%
        
        # Then we sort the rows by using a factor level:
        
        mutate(variable = factor(variable, levels = c("age", "height", "weight", "LBM"))) %>%
        arrange(variable) %>%
        mutate(Variable = c("Age (yrs)", "Stature (cm)", "Body mass (kg)", "Lean body mass (%)")) %>%
        select(Variable, female_excl:male_incl) %>%
        
        #Then we make our table anew:
        kable(format = "html", col.names = c(" ", 
                                             "Excluded",
                                             "Included", 
                                             "Excluded", 
                                             "Included"), 
              caption = "Table 1. Participant characteristics") %>%
        add_header_above(c(" " = 1, "Female" = 2, "Male" = 2)) %>%
        footnote(general = "Values are Mean (SD)")%>%
        #Finally we make the format into a classical paper-like table:
        
        kable_classic(full_width = F, html_font = "Cambria")


#------------------------

library(flextable)










