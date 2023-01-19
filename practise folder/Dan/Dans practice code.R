# Dan's coding journey
#
# good luck mate :)


## PACKAGE INSTALLATION
# input the following code and type the package of interest:

        install.packages("...")
        
# N.B. it is possible to download packages directly from CRAN repository via ´packages > install´ in the bottom-right window.
        
        

## RUNNING PACKAGES
# important to first run packages (cmd+enter) before proceeding with data import and creating figures. 
# N.B. it is possible to run multiple packages simultaneously by highlight all packages of interest if on separate lines like below 
        
        library(readxl)
        library(tidyverse)
        library(ggplot2)

# N.B.it is also possible to run multiple packages simultaneoulsy by using ´;´ symbol to separate packages when on the same line like below.

        library(readxl);library(tidyverse);library(ggplot2)
        
        
# to identify the current working environment use code below:
        
        getwd()
        
        
## IMPORTING DATA FROM EXCEL 
# N.B. it is important to ensure data in raw excel files starts from row 1, column 1.
# refer to https://dhammarstrom.github.io/IDR4000/lesson_5_import_data.html 
# important to first run ´readxl´package. This will enable function to read excel files.
# after running package, insert the code below to access local(?) or github files.
                
        RMA_MVC<-read_excel("/Users/danieltu/OneDrive - nih.no/NIH/Study 1 - Human Repeated Muscle Atrophy (RMA)/Data (RMA)/RMA_MVC.xlsx")
      
# Will automatically open the first tab within excel file.
# select the ´tab´ button within quotation marks for dropdown menu of folders. Simply repeat until you reach the xlsx./csv. file of interest.
# it is also possible to open excel files directly from the bottom-right window. Select ´Home > ... > Import Dataset > Import´. The dataset should now appear in a new tab within script window. However, much more effecient to access data via code rather than manually importing. 
# if end folder within path is not desired, manually input ´/users´code to create new path

