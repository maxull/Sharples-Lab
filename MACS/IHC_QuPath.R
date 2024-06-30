###
###
###             QuPath analysis if nuclei
###
###

# load data

MACS001_1.1 <- read_csv('/Users/maxullrich/OneDrive - UGent/Work/Projects/MACS/DATA/IHC/MACS_001_1.1.csv')

colnames(MACS001_1.1) <- gsub(" ", "_", colnames(MACS001_1.1))
colnames(MACS001_1.1) <- gsub(":", "", colnames(MACS001_1.1))
colnames(MACS001_1.1) <- gsub("\\^", "", colnames(MACS001_1.1))

MACS001_1.1 %>% 
        ggplot(aes(x = Nucleus_Area))+
        geom_histogram()

# check distribution of nuclei, myonuclei and satellite cell signal

MACS001_1.1 %>% 
        ggplot(aes(x = Nucleus_Blue_max))+
        geom_histogram()

MACS001_1.1 %>% 
        ggplot(aes(x = Nucleus_Green_max))+
        geom_histogram()

MACS001_1.1 %>% 
        ggplot(aes(x = Nucleus_Red_max))+
        geom_histogram()

# focus on the nuclei with the stronges signal first

# filter blue > 250

MACS001_1.1 %>% 
        filter(Nucleus_Blue_max > 250) %>% 
        ggplot(aes(x = Nucleus_Red_max, y = Nucleus_Green_max))+
        geom_point()

# threshold for green / satellite cell stain 240 
# threshold for red / myonuclei stain 120

MACS001_1.1 %>% 
        filter(Nucleus_Blue_max > 250) %>% 
        mutate(Myonuclei = ifelse(Nucleus_Red_max > 120, 'Myonuclei', ''))
