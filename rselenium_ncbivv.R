#### DOCUMENTATION ####

##* Contact ####
# Author : Elika Garg <elika.garg@mail.mcgill.ca>
# OK reasons to contact : trouble with script, customization to the script, understanding its working (after you have read through the annotations and comments and referred to RSelenium documentation as well), other feedback, showering praises
# NOT OK reasons to contact : trouble with Selenium (please utilize the power of Internet for troubleshooting that)

##* Installation ####
# Libraries to install : tidyverse, magrittr, stringr, RSelenium
# Install selenium in your machine as per instructions from RSelenium - this can be a long and cumbersome process, often requiring administrative permissions, stay strong and may the force be with you

##* Troubleshoot ####
# Note : If you see this error - "docker: Cannot connect to the Docker daemon at unix:///var/run/docker.sock. Is the docker daemon running?."- go to terminal, login as 'su' and run "systemctl start docker" - this happens when docker is inactive and needs to be activated manually

##* Runtime ####
# Time taken for 1 gene is around 30s, apart from the overheads (~1min per part)

##* Caution ####
# This script only works for the current webpage structure, if it changes in any manner (e.g., menu options), the script will have to be edited
# Please monitor the screenshot displays to ensure that the webpage is being navigated correctly
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

# the wrapper was written to deal with huge lists which run into memory issues (the docker browser cache builds up and the cleaning commands were not helpful)
# therefore, the wrapper function divided the original list into chunks (size defined by 'splitter') and runs each chunk in a separate block, then collates all results

##> Arguments ####

# input_list = character vector of gene names
# name_file = name of the file which would be used to save logs, temporaty files, and final results; this can include the full path but should exclude file extensions
# splitter = chunk sizes for big lists, default is 100

##> Results ####

# if the script runs smoothly, you would get a log file (.log), a temporary file (.txt), and a results file (.csv) for every part, plus a final collated results file (.csv)

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

##> Example ####

# run_ncbivv(input_list = c("NUP210", "c11orf10", "RASSF1"), name_file = "ncbivv_test", splitter = 2)
# this generates 2 log files (ncbivv_test_3_part1.log, ncbivv_test_3_part2.log), 2 temporary files (ncbivv_test_3_part1.txt, ncbivv_test_3_part2.txt), and 3 results files (ncbivv_test_3_part1_2.csv [2 genes because splitter was 2], ncbivv_test_3_part2_1.csv [1 gene because total input length was 3], ncbivv_test_3.csv [all 3 genes])


run_ncbivv <- function(input_list, name_file, splitter = 100){
  
  try(system(command = "sudo docker ps -q"), system(command = "sudo docker stop  $(sudo docker ps -q)")) # if the docker is already running, stops it first (this can happen when there is an interruption in the script and it exits before stopping the docker)
  
  input_list = str_trim(input_list) %>% unique() # removes whitespaces and only retains unique entries
  
  # if the splitter value is more than the list, then the list is divided up into chunks and the results are collated at the end
  if(length(input_list) > splitter){
    split_input_list <- split(input_list, ceiling(seq_along(input_list)/splitter))
    ncbivv_result <- mapply(ncbivv, 
                            input_list = split_input_list, 
                            name_file = paste0(name_file, "_", length(input_list), "_part", seq_along(split_input_list)), # each chunk is also saved with the length of the total input (unique entries) and the part/chunk number added to the name
                            SIMPLIFY = FALSE) %>% 
      data.table::rbindlist(.)
    write_csv(ncbivv_result, paste0(name_file, "_", length(input_list), ".csv")) # collated result has the length of the total input (unique entries) added to the name
  } else 
    ncbivv_result <- ncbivv(input_list = input_list, name_file = name_file)
  return(ncbivv_result)
}

##* Main function ####

ncbivv <- function(input_list, name_file) {
  
  #### Setup webpage in docker browser ####
  
  message("Starting docker by running this in the terminal : \n 
          sudo docker run -d -p 4445:4444 selenium/standalone-firefox:2.53.1")
  
  system(command = "sudo docker run -d -p 4445:4444 selenium/standalone-firefox:2.53.1")
  
  Sys.sleep(10)
  
  remDr <- remoteDriver(port = 4445) ## starts server ##
  remDr$open()
  
  remDr$navigate("https://www.ncbi.nlm.nih.gov/variation/view/") ## opens NCBI Variation Viewer webpage ##
  Sys.sleep(5); remDr$screenshot(display = TRUE) # this command is used variably throughout, it put the system on sleep for the number of seconds indicated so that the browser has enough time to load the webpage according to the preceding command, then it displays a screenshot of the current webpage
  
  ####| "Pick Assembly" ####
  ## "Pick Assembly" section starts ##
  ## it selects the second option that is the "GRCh37.p13: Annotation Release 105" assembly (if the options change, the script will need to be edited) ## 
  
  assembly_click <- remDr$findElement(using = "css selector", ".gb-region-navigator .gb-section:nth-child(1) span:nth-child(1)")
  assembly_click$clickElement()
  Sys.sleep(5); remDr$screenshot(display = TRUE)
  
  assembly_menu <- remDr$findElement(using = "css", ".gb-ds-item")
  assembly_menu$clickElement()
  Sys.sleep(5); remDr$screenshot(display = TRUE)
  
  assembly_default <- remDr$findElement(using = "css", ".gb-selected")
  Sys.sleep(5); remDr$screenshot(display = TRUE)
  
  remDr$mouseMoveToLocation(webElement = assembly_default)
  remDr$sendKeysToActiveElement(list(key = "down_arrow", key = "enter"))# this part will need editing if "Pick Assembly" menu options change
  Sys.sleep(5); remDr$screenshot(display = TRUE)
  
  ## "Pick Assembly" section ends ##
  
  ncbivv_log <- file(paste0(name_file, ".log"), "w+") # logs script messages to a text file, use this file to investigate sources of error
  
  sink(ncbivv_log, type = "message")
  
  message("Processing : ")
  
  ncbivv_con <- file(paste0(name_file, ".txt"), "w+") # writes temporary results in a text file, use this file to recover output if the script is interrupted and does not produce neatly collated results
  
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
    remDr$screenshot(display = TRUE)
    search_click$sendKeysToElement(list(input_list[listed_genes], key = "enter"))
    Sys.sleep(10)
    try(remDr$screenshot(display = TRUE), TRUE)
    if(remDr$status > 0) {
      message("Skipping : No results found")
      message(remDr$statusCodes[which(remDr$statusCodes$Code == remDr$status),])
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
      Sys.sleep(10); remDr$screenshot(display = TRUE)}
    output_matrix[listed_genes, "Gene"] <- gene_name$getElementText()[[1]]
    
    gene_position <- remDr$findElement(using = "css", ".grid-data tr:nth-child(1) .gbs-location span") ## gets full gene location/position from the search box result##

    output_matrix[listed_genes, "Gene_Full_Location"] <- gene_position$getElementText()[[1]]
    
    ####| "Source database" ####
    ## "dbSNP" ##
    
    dbSNP_click <- remDr$findElement(using = "css", "#SourceDB_ID_value_0") ## checks the "dbSNP" box in "Source database"
    if(!dbSNP_click$isElementSelected()[[1]]) {
      dbSNP_click$clickElement()
      Sys.sleep(10); remDr$screenshot(display = TRUE)
    }
    
    ####| "Variant type" ####
    ## "single nucleotide variant" ##
    
    varSNP_click <- remDr$findElement(using = "css", "#VarType_ID_value_1483")
    if(!varSNP_click$isElementSelected()[[1]]) {
      varSNP_click$clickElement()
      Sys.sleep(10); remDr$screenshot(display = TRUE)
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
      Sys.sleep(5); remDr$screenshot(display = TRUE)
    } else if(numpages == 2){ ## if the number of pages is more than 2, it goes to the 'next' page
      dbSNP_endclick <- remDr$findElement(using = "css", ".nav-controls .next")
      dbSNP_endclick$clickElement()
      Sys.sleep(5); remDr$screenshot(display = TRUE)
    }  ## if the number of pages is 1, it stays on the 'first' page
    
    dbSNP_end <- remDr$findElements(using = "css", ".vdt-location")
    output_matrix[listed_genes, "Last_Variant_BP"] <- dbSNP_end[[length(dbSNP_end)]]$getElementText()[[1]] # last variant location
    
    write_lines(str_c(output_matrix[listed_genes,], collapse = "\t"), ncbivv_con, append = TRUE)
  }
  
  sink(NULL, type = "message")
  
  close(ncbivv_con)  
  
  # this chunk cleans up the data generated from the for loop above to write a csv file
  # in case there were interruptions and all you have is the temporary text file, then you can supply that information and use this code to get the clean output
  
  output_matrix %<>% 
    as.data.frame() %>% 
    mutate(Chromosome = str_match(Gene_Full_Location, "Chr(.*?):")[,2]) %>%
    mutate_at(c("First_Variant_BP", "Last_Variant_BP"), funs(str_replace_all(., ",", ""))) %>%
    dplyr::select(Given_Name, Gene, Chromosome, First_Variant_BP, Last_Variant_BP, Gene_Full_Location)
  
  write_csv(output_matrix, paste0(name_file, "_", length(input_list), ".csv")) # add the total number of entries to the name of the file
  
  remDr$close()
  remDr$closeServer()
  remDr$quit()
  
  message("Stopping docker by running this in the terminal : \n 
        sudo docker stop  $(sudo docker ps -q)")
  
  system(command = "sudo docker stop  $(sudo docker ps -q)")
  
  Sys.sleep(10)
  
  close(ncbivv_log)
  
  return(output_matrix)
}