#### DOCUMENTATION ####

##* Contact ####
# Author : Elika Garg <elika.garg@mail.mcgill.ca>
# OK reasons to contact : trouble with script, customization to the script, understanding its working (after you have read through the annotations and comments and referred to RSelenium documentation as well), other feedback, showering praises
# NOT OK reasons to contact : trouble with Selenium (please utilize the power of Internet for troubleshooting that)

##* Installation ####
# Libraries to install : tidyverse, magrittr, stringr, data.table, RSelenium
# Install selenium in your machine as per instructions from RSelenium - this can be a long and cumbersome process, often requiring administrative permissions, stay strong and may the force be with you

##* Troubleshoot ####
# Note : If you see this error - "docker: Cannot connect to the Docker daemon at unix:///var/run/docker.sock. Is the docker daemon running?."- go to terminal, login as 'su' and run "systemctl start docker" - this happens when docker is inactive and needs to be activated manually

##* Runtime ####
# Time taken for 1 gene is around 30s, apart from the overheads (~1min per part)

##* Caution ####
# This script only works for the current webpage structure, if it changes in any manner (e.g., menu options), the script will have to be edited
# In case of trouble, please monitor the screenshot displays to ensure that the webpage is being navigated correctly
# Always do a random manual check at the end to compare machine output with human output
# This is a hard-coded script, so it is not very flexible and breaks easily - everytime something breaks, that is patched up but new bugs may emerge ..

##* Usage ####
# What this script CAN do : extract information from a webpage
# What this script CANNOT do : download information from a webpage

#### CODE ####

##* Libraries ####

library(tidyverse); library(magrittr); library(stringr); library(data.table)

require(RSelenium)

##* Wrapper function ####

# the wrapper was written to clean the temporary files and enable the use of a master list
# in case of interruptions, re-running the wrapper also takes care of salvaging the data from the temporary files
# although, you may have to run the wrapper a few times if the list is huge, as it sometimes gets interrupted for reasons unknown to the cosmos

##> Arguments ####

# input_list = character vector of gene names
# name_file = name of the file which would be used to save logs, temporaty files, and final results; this can include the full path but should exclude file extensions
# screenshot = if TRUE, displays screenshots at every step, default is FALSE (note : since it saves a screenshot for every step for every gene, R soon runs out of memory if the input list it too long)
# use_master = name of the master file (with path if not in the working directory, and .csv extenstion) to use (if it exists, or create if it doesn't) and add to; this master file will comprise all the unique genese from runs so far, default is NA (in which case it is not used)
# increase_wait = number of seconds that will be added to the waiting time at each step (one number is added to all); use this if the internet speed is slow

##> Results ####

# if the script runs smoothly, you would get one final results file (.csv), and the temporary files (.log and .txt) will be removed
# the results file will consist of all unique and non-NA genes from the input list 
# all files created with this script have file names that end with "_REGN" before the file extension

## Log file
# you can check the log file dynamically (e.g. by using the terminal watch command) to see the progress of the script
# it also displays an extra line whenever the given gene name does not match the found name, keep in mind that the gene name found on the website will be used to get the locations
# this can be very useful if there's an interruption

## Temporary file
# this file is written line-by-line as the for loop runs through the gene list, so that if there is an interruption you can salvage this data

## Results file
# this collates the temporary file results into a table and processes them a little for a cleaner look
# the number of genes in the file is added to the name at the end
# Given_Name = gene name from the input list (with whitespaces removed)
# Gene = the first gene name found in the search box output; the results are presented for this gene
# Chromosome = chromosome of the gene
# First_Variant_BP = base pair location of the first Variant ID
# Last_Variant_BP = base pair location of the last Variant ID
# Gene_Full_Location = full location of the gene as-is from the search box output

## Master file
# if the master file option is selected, collects and collates all output files

##> Example ####

# system(command = "mkdir REGN_Test") # this is to create the testing directory
# run_ncbivv(input_list = c("NUP210", "c11orf10", "RASSF1"), name_file = "REGN_Test/ncbivv_test1", use_master = "REGN_Test/master_test.csv") # this generates 1 ncbivv_test_3_output_REGN.csv [all 3 genes]) and 1 master_test.csv [all 3 genes]
# run_ncbivv(input_list = c("NUP210", "SOX21", "Aph1c"), name_file = "REGN_Test/ncbivv_test2", use_master = "REGN_Test/master_test.csv") # this generates 1 ncbivv_test2_3_output_REGN.csv [all 3 genes]) and updates the master_test.csv [all 5 genes]
# system(command = "rm -r REGN_Test") # this is to remove the testing directory


run_ncbivv <- function(input_list, name_file, screenshot = FALSE, use_master = NA, increase_wait = 0){
while(length(str_subset(list.files(path = dirname(name_file), pattern = basename(name_file)), "_output_REGN.csv"))<1){ # while loop runs until the final output file has been created
  closeAllConnections() # closes any open connections, to deal with temporary files
  Sys.sleep(5+increase_wait)
  input_list = str_trim(input_list) %>% unique() %>% str_replace_na() %>% grep("^NA$|^$", ., value = T, ignore.case = T, invert = T) # removes whitespaces and only retains unique, valid entries
  
  deal_dups <- function(dff) {
    dff_mixed <- dff %>% 
      distinct() %>% 
      group_by(Given_Name) %>% 
      dplyr::filter((Last_Variant_BP-First_Variant_BP)==max(Last_Variant_BP-First_Variant_BP) | is.na(Last_Variant_BP-First_Variant_BP))
    dff_singles <- dff_mixed %>% 
      add_count() %>% 
      dplyr::filter(n==1) %>% 
      dplyr::select(-n)
    dff_doubles <- dff_mixed %>% 
      add_count() %>% 
      dplyr::filter(n>1) %>% 
      dplyr::select(-n)
    if(nrow(dff_doubles)>0) dff_doubles %>% 
      summarise_all(funs(first(na.omit(.)))) %>% 
      full_join(dff_singles) %>% 
      return() else return(dff_singles)
  }
  
  if(!is.na(use_master)) {
  if(!file.exists(use_master) & length(list.files(path = dirname(name_file), pattern = "_output_REGN.csv"))>0) { # if specified master file does not exist, then creates one using all the output REGN files
    message(" # # # # # # # # # # # # # # # # # # # # # ")
    message("Creating master file : ", use_master)
    message("# # # # # # # # # # # # # # # # # # # # # #")
    list.files(path = dirname(name_file), pattern = "_output_REGN.csv", full.names = T) %>%
      lapply(., read_csv) %>%
      rbindlist(.) %>%
      as.data.frame() %>%
      deal_dups() %>%
      write_csv(use_master)
  } 
  if(file.exists(use_master)) {# then uses the newly created or existing master file
    message(" # # # # # # # # # # # # # # # # # # # # # ")
    message("Using master file : ", use_master)
    message("# # # # # # # # # # # # # # # # # # # # # #")
    master_sofar <- read_csv(use_master) %>% 
      inner_join(data.frame(Given_Name = input_list)) %>% 
      deal_dups()
    write_csv(master_sofar, paste0(name_file, "_", nrow(master_sofar), "_frommaster_REGN.csv"))
  }
  }# creates a subset of the master file that is part of the input list, so that the genes that have previously been run, need not be run again
  # however, there might exist some updates to some genes, therefore, the largest difference in the first and last positions is used
  
  temp_sofar <- str_subset(list.files(path = dirname(name_file), pattern = paste0(basename(name_file), "_"), full.names = T), "_REGN.txt")
  if(length(temp_sofar)>0) { # finds and uses temporary text files
    tryCatch(lapply(temp_sofar, read_tsv) %>%
      rbindlist(.) %>% 
      as.data.frame() %>% 
      mutate(Chromosome = str_match(Gene_Full_Location, "Chr(.*?):")[,2]) %>%
      mutate_at(c("First_Variant_BP", "Last_Variant_BP"), funs(str_replace_all(., ",", ""))) %>%
      dplyr::select(Given_Name, Gene, Chromosome, First_Variant_BP, Last_Variant_BP, Gene_Full_Location) %>% 
      write_csv(paste0(name_file,"_rescue_REGN.csv")), error = function(e) e)
  } # creates a rescue file from the temporary text files
  results_sofar <- NULL
    check_sofar <- str_subset(list.files(path = dirname(name_file), pattern = paste0(basename(name_file), "_"), full.names = T), "_REGN.csv")
    if(length(check_sofar)>0) results_sofar <- lapply(check_sofar, read_csv) %>% rbindlist(.) # finds and uses incomplete output files
  
  temp_to_remove <- str_subset(list.files(path = dirname(name_file), pattern = paste0(basename(name_file), "_"), full.names = T), "_REGN")
  
  message(" # # # # # # # # # # # # # # # # # # # # # ")
  message("Using and removing temporary files :")
  print(temp_to_remove)
  message("# # # # # # # # # # # # # # # # # # # # # #")
  
  input_bad <- temp_to_remove %>% 
    str_subset(".log") %>% # finds and uses temporary log files
    sapply(., read_lines) %>% 
    unlist() %>% 
    str_subset("Skipping") %>% 
    str_remove("Skipping : No results found - ") # finds the genes that were skipped, and adds them in the output with NA values, so that these are not re-run several times
  
  tryCatch(results_sofar %<>% 
             as.data.frame %>% 
             full_join(data.frame(Given_Name = input_bad)) %>% 
             inner_join(data.frame(Given_Name = input_list)) %>% 
             distinct(), error = function(e) e)
  
  walk(temp_to_remove, file.remove)
  
  write_csv(as.data.frame(results_sofar), paste0(name_file, "_", nrow(results_sofar), "_incomplete_REGN.csv"))
  # creates an incomplete file using all the temporary log and text files, as well as the rescued and master subset files
  
  if(setequal(input_list, results_sofar[[1]])) { # if the input list is complete, produces final output file
    message(" # # # # # # # # # # # # # # # # # # # # # ")
    message("Completed !!!", "\nFinished number of genes = ", length(input_list))
    message("# # # # # # # # # # # # # # # # # # # # # #")
    write_csv(results_sofar, paste0(name_file, "_", nrow(results_sofar), "_output_REGN.csv"))
    temp_to_remove <- grep("_output_REGN.csv$", str_subset(list.files(path = dirname(name_file), pattern = basename(name_file), full.names = T), "_REGN"), invert = TRUE, value = TRUE)
    
    message(" # # # # # # # # # # # # # # # # # # # # # ")
    message("Using and removing temporary files :")
    print(temp_to_remove)
    message("# # # # # # # # # # # # # # # # # # # # # #")
    
    walk(temp_to_remove, file.remove) 
    if(!is.na(use_master)) { # creates or updates the master file
      message(" # # # # # # # # # # # # # # # # # # # # # ")
      message("Creating/Updating master file :", use_master)
      message("# # # # # # # # # # # # # # # # # # # # # #")
      list.files(path = dirname(name_file), pattern = "_output_REGN.csv", full.names = T) %>%
      lapply(., read_csv) %>%
      rbindlist(.) %>%
      as.data.frame() %>%
      deal_dups() %>%
      write_csv(use_master)}
    } else {
  
  input_list_trunc = setdiff(input_list, results_sofar[[1]]) # if there is remnant input list, runs the remaining genes
  message(" # # # # # # # # # # # # # # # # # # # # # ")
  message("Length of input list total = ", length(input_list))
  message("Length of input list done = ", length(results_sofar[[1]]))
  message("Length of input list remaining = ", length(input_list_trunc))
  message("# # # # # # # # # # # # # # # # # # # # # #")
  
  ncbivv(input_list = input_list_trunc, name_file = name_file, screenshot = screenshot, increase_wait = increase_wait)
  Sys.sleep(10+increase_wait)
    }
}
}

##* Main function ####

ncbivv <- function(input_list, name_file, screenshot = FALSE, increase_wait = 0) {
  
  #### Setup webpage in docker browser ####
  
  message("Starting docker by running this in the terminal : \n 
          sudo docker run -d -p 4445:4444 selenium/standalone-firefox:2.53.1")
  
  system(command = "sudo docker run -d -p 4445:4444 selenium/standalone-firefox:2.53.1")
  
  Sys.sleep(10+increase_wait)
  
  remDr <- remoteDriver(port = 4445) ## starts server ##
  remDr$open()
  
  if(screenshot) message("Caution : displaying screenshots, if the list is very long (>100), might run out of memory")
  
  remDr$navigate("https://www.ncbi.nlm.nih.gov/variation/view/") ## opens NCBI Variation Viewer webpage ##
  Sys.sleep(5+increase_wait); if(screenshot) remDr$screenshot(display = TRUE) # this command is used variably throughout, it put the system on sleep for the number of seconds indicated so that the browser has enough time to load the webpage according to the preceding command, then it displays a screenshot of the current webpage
  
  ####| "Pick Assembly" ####
  ## "Pick Assembly" section starts ##
  ## it selects the second option that is the "GRCh37.p13: Annotation Release 105" assembly (if the options change, the script will need to be edited) ## 
  
  assembly_click <- remDr$findElement(using = "css selector", ".gb-region-navigator .gb-section:nth-child(1) span:nth-child(1)")
  assembly_click$clickElement()
  Sys.sleep(5+increase_wait); if(screenshot) remDr$screenshot(display = TRUE)
  
  assembly_menu <- remDr$findElement(using = "css", ".gb-ds-item")
  assembly_menu$clickElement()
  Sys.sleep(5+increase_wait); if(screenshot) remDr$screenshot(display = TRUE)
  
  assembly_default <- remDr$findElement(using = "css", ".gb-selected")
  Sys.sleep(5+increase_wait); if(screenshot) remDr$screenshot(display = TRUE)
  
  remDr$mouseMoveToLocation(webElement = assembly_default)
  remDr$sendKeysToActiveElement(list(key = "down_arrow", key = "enter"))# this part will need editing if "Pick Assembly" menu options change
  Sys.sleep(5+increase_wait); if(screenshot) remDr$screenshot(display = TRUE)
  
  ## "Pick Assembly" section ends ##
  
  ncbivv_log <- file(paste0(name_file, "_REGN.log"), "w+") # logs script messages to a text file, use this file to investigate sources of error
  
  sink(ncbivv_log, type = "message")
  
  message("Processing : ")
  
  ncbivv_con <- file(paste0(name_file, "_REGN.txt"), "w+") # writes temporary results in a text file, use this file to recover output if the script is interrupted and does not produce neatly collated results
  
  output_matrix <- matrix(data = NA, nrow = length(input_list), ncol = 5, dimnames = list(input_list, c("Given_Name", "Gene", "First_Variant_BP", "Last_Variant_BP", "Gene_Full_Location")))
  
  write_lines(str_c(c("Given_Name", "Gene", "First_Variant_BP", "Last_Variant_BP", "Gene_Full_Location"), collapse = "\t"), ncbivv_con, append = TRUE)
  
  #### Gene list web-scraping ####
  
  for(listed_genes in seq_along(input_list)){
    
    message(listed_genes, " ", input_list[listed_genes])
    
    output_matrix[listed_genes, "Given_Name"] <- input_list[listed_genes]
    
    ####| "Search" ####
    ## gene search and selection ##
    
    search_click <- remDr$findElement(using = "css", "#loc-search")
    search_click$clearElement() ## clears search box
    if(screenshot) remDr$screenshot(display = TRUE)
    search_click$sendKeysToElement(list(input_list[listed_genes], key = "enter"))
    Sys.sleep(10+increase_wait)
    try(remDr$screenshot(display = FALSE), TRUE)
    
    if(remDr$status > 0) {
      message("Skipping : No results found - ", input_list[listed_genes])
      message(remDr$statusCodes[which(remDr$statusCodes$Code == remDr$status),])
      remDr$refresh()
      Sys.sleep(5+increase_wait)
      next
    }
    
    gene_name <- remDr$findElement(using = "css", ".grid-data tr:nth-child(1) .gbs-name")
    
    ## this subsection deals with the cases when the gene name in the input list is not found and NCBI suggests an alternative ##
    ## the first alternative is automatically selected and the information for this alternative is extracted ##
    ## (but the "given name" is saved as well) ##
    
    if(gene_name$getElementText()[[1]] != str_to_upper(input_list[listed_genes])) {
      message("Different names : Given = ", input_list[listed_genes], " , Found = ", gene_name$getElementText()[[1]])
      gene_select <- remDr$findElement(using = "css", "tr:nth-child(1) .ctrl span")
      gene_select$clickElement()
      Sys.sleep(10+increase_wait); if(screenshot) remDr$screenshot(display = TRUE)
      }
    output_matrix[listed_genes, "Gene"] <- gene_name$getElementText()[[1]]
    
    gene_position <- remDr$findElement(using = "css", ".grid-data tr:nth-child(1) .gbs-location span") ## gets full gene location/position from the search box result##
    
    output_matrix[listed_genes, "Gene_Full_Location"] <- gene_position$getElementText()[[1]]
    
    ####| "Source database" ####
    ## "dbSNP" ##
    
    dbSNP_click <- remDr$findElement(using = "css", "#SourceDB_ID_value_0") ## checks the "dbSNP" box in "Source database"
    if(!dbSNP_click$isElementSelected()[[1]]) {
      dbSNP_click$clickElement()
      Sys.sleep(10+increase_wait); if(screenshot) remDr$screenshot(display = TRUE)
    }
    
    ####| "Variant type" ####
    ## "single nucleotide variant" ##
    
    varSNP_click <- remDr$findElement(using = "css", "#VarType_ID_value_1483")
    if(!varSNP_click$isElementSelected()[[1]]) {
      varSNP_click$clickElement()
      Sys.sleep(10+increase_wait); if(screenshot) remDr$screenshot(display = TRUE)
    }
    
    #### Variation data table ####
    ## gene locations ##
    
    ## the ultimate goal of the script is to get the 'Location' of the first and the last 'Variant ID' from the first and the last pages, respectively
    
    dbSNP_start <- remDr$findElement(using = "css", ".vdt-location")
    output_matrix[listed_genes, "First_Variant_BP"] <- dbSNP_start$getElementText()[[1]] # first variant location
    
    ## this subsection gets the 'Location' of the last 'Variant ID' from the variation data table
    ## it needed to be convoluted because if the number of pages is 1, 2, or more the navigation options are different
    
    numpages <- remDr$findElement(using = "css", ".nav-controls .ui-ncbigrid-paged-endPage") ## gets the number of total pages in the variation data table
    numpages = as.numeric(gsub(",", "", numpages$getElementText()[[1]]))
    
    if(numpages > 2){ ## if the number of pages is more than 2, it goes to the 'last' page
      dbSNP_endclick <- remDr$findElement(using = "css", ".nav-controls .ui-ncbigrid-paged-pageControl-last")
      dbSNP_endclick$clickElement()
      Sys.sleep(5+increase_wait); if(screenshot) remDr$screenshot(display = TRUE)
    } else if(numpages == 2){ ## if the number of pages is more than 2, it goes to the 'next' page
      dbSNP_endclick <- remDr$findElement(using = "css", ".nav-controls .next")
      dbSNP_endclick$clickElement()
      Sys.sleep(5+increase_wait); if(screenshot) remDr$screenshot(display = TRUE)
    }  ## if the number of pages is 1, it stays on the 'first' page
    
    dbSNP_end <- remDr$findElements(using = "css", ".vdt-location")
    output_matrix[listed_genes, "Last_Variant_BP"] <- dbSNP_end[[length(dbSNP_end)]]$getElementText()[[1]] # last variant location
    
    write_lines(str_c(output_matrix[listed_genes,], collapse = "\t"), ncbivv_con, append = TRUE)
  }
  
  on.exit(sink(NULL, type = "message"), add = TRUE)
  
  on.exit(close(ncbivv_con), add = TRUE)
  
  # this chunk cleans up the data generated from the for loop above to write a csv file
  # in case there were interruptions and all you have is the temporary text file, then you can supply that information and use this code to get the clean output
  
  output_matrix %<>% 
    as.data.frame() %>% 
    mutate(Chromosome = str_match(Gene_Full_Location, "Chr(.*?):")[,2]) %>%
    mutate_at(c("First_Variant_BP", "Last_Variant_BP"), funs(str_replace_all(., ",", ""))) %>%
    dplyr::select(Given_Name, Gene, Chromosome, First_Variant_BP, Last_Variant_BP, Gene_Full_Location)
  
  write_csv(output_matrix, paste0(name_file, "_", nrow(output_matrix), "_REGN.csv")) # add the total number of entries to the name of the file
  
  message("Stopping docker by running this in the terminal : \n 
        sudo docker stop  $(sudo docker ps -q)")
  
  on.exit(system(command = "sudo docker stop  $(sudo docker ps -q)"), add = TRUE)
  
  Sys.sleep(10+increase_wait)
  
  on.exit(close(ncbivv_log), add = TRUE)
}
