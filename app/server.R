
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
# setwd("D:/One Drive Jithin/OneDrive - St John's National Academy of Health Sciences/PsychoPy/R-EEG Analysis/BoubaKikiApp")
library(shiny)
library(edf)
library(dplyr)
library(ggplot2)
library(reshape2)
library(signal)
shinyServer(function(input, output) {
  observeEvent(input$goButton,{
  
  edf_server <- input$edf_file
  dataset <- read.edf(edf_server$datapath)
  # dataset <- read.edf("20170417114942_S7-MMM2-20170417.edf")
  # str(dataset)
  

  #Classifying trigger as start of shape or start of sound
  trigger <- as.data.frame(dataset$events)
  trigger$type <- with(trigger,ifelse(substr(annotation,8,8)=="#","Sound","-"))
  
  yes <- c("1")
  no <- c("2")
  
  trigger$match <- with(trigger,ifelse(substr(annotation,9,9) %in% yes,"Yes",ifelse(substr(annotation,9,9) %in% no,"No","-")))
  
  sound <- trigger[trigger$type=="Sound",]
  
  i.pre <- 400
  i.pre <- input$averaging
  
  i.stim_start <- 150
  i.stim_start <- input$stim_start
  
  i.stim_stop <- 1000
  i.stim_stop <- input$stim_stop
  
  sound$baseline_start <- sound$onset-(i.pre/1000) #400ms pre-audio
  sound$stimulus_start <- sound$onset-i.stim_start/1000 #150ms pre-audio stimulus interval
  sound$stimulus_end <- sound$onset+ i.stim_stop/1000 #1000ms post-audio stimulus interval
  
  
  #Change the $__ and =c("__",)
  ch1 <- as.data.frame(dataset$signal$F3,col.names=c("F3","t"))
  ch2 <- as.data.frame(dataset$signal$Fz,col.names=c("Fz","t"))
  ch3 <- as.data.frame(dataset$signal$F4,col.names=c("F4","t"))
  ch4 <- as.data.frame(dataset$signal$C3,col.names=c("C3","t"))
  ch5 <- as.data.frame(dataset$signal$Cz,col.names=c("Cz","t"))
  ch6 <- as.data.frame(dataset$signal$C4,col.names=c("C4","t"))
  ch7 <- as.data.frame(dataset$signal$P3,col.names=c("P3","t"))
  ch8 <- as.data.frame(dataset$signal$P4,col.names=c("P4","t"))
  
  ch_all <- merge(ch1,ch2,by="t")
  ch_all <- merge(ch_all,ch3,by="t")
  ch_all <- merge(ch_all,ch4,by="t")
  ch_all <- merge(ch_all,ch5,by="t")
  ch_all <- merge(ch_all,ch6,by="t")
  ch_all <- merge(ch_all,ch7,by="t")
  ch_all <- merge(ch_all,ch8,by="t")
  
  rm(ch1,ch2,ch3,ch4,ch5,ch6,ch7,ch8)
  

  
  bf <- signal::butter(2,W=c(0.3/250,30/250),type="pass",plane="z")
  for (i in 2:9) {
    ch_all[,i] <- signal::filtfilt(bf,ch_all[,i])
    
  }
  
  #Function to calculate epoch average for each sound trigger
  epoch_baseline <- function(channel_data,base_start,stim_onset,stim_start,stim_end,stim_record,stim_annotation,stim_match){
    #Stimulus data
    stim_data <- channel_data[channel_data[,1]>=stim_start&channel_data[,1]<=stim_end,]
    stim_data$type = "Stimulus"
    
    
    #Pre data
    base_data <- channel_data[channel_data[,1]>=base_start&channel_data[,1]<=stim_onset,]
    base_data$type = "Baseline"
    
    #Averages
    base_avg <- base_data[,2:10] %>% group_by(type) %>% summarise_all(.,funs(mean(.,na.rm=TRUE)))
    
    colnames(base_avg) <- c("type",paste("base_",colnames(ch_all)[2:9]))
    
    base_avg <- base_avg[rep(1,each=length(stim_data[,1])),]
    stim_data[,2:9] <- stim_data[,2:9] - base_avg[,2:9]
    stim_data$record <- stim_record
    stim_data$annotation <- stim_annotation
    stim_data$match <- stim_match
    stim_data$index <- seq(1:length(stim_data[,"match"]))
    return(stim_data)
  }
  
  # function(channel_data,base_start,stim_onset,stim_start,stim_end,stim_record)
  
  consolidated <- epoch_baseline(ch_all,sound$baseline_start[1],sound$onset[1],sound$stimulus_start[1],sound$stimulus_end[1],sound$record[1],sound$annotation[1],sound$match[1])
  
  for (i in 2: length(sound$onset)){
    consolidated <- rbind(consolidated,epoch_baseline(ch_all,sound$baseline_start[i],sound$onset[i],sound$stimulus_start[i],sound$stimulus_end[i],sound$record[i],sound$annotation[i],sound$match[i]))
  }
  
  for (i in 2:9) {
    ch_all[,i] <- ifelse(ch_all[,i]>-200000&ch_all[,i]<200000,ch_all[,i],NA)
  }
  
  
  consolidated_waves <- consolidated[,c(2:9,12:14)] %>% group_by(annotation,match,index) %>% summarise_all(.,funs(mean(.,na.rm=TRUE)))
  
  summary <- consolidated_waves %>% group_by(annotation) %>% summarize(index = max(index))
  # rm(ch_all,ch1,ch2,ch3,ch4,ch5,ch6,ch7,ch8)

 consolidated_waves[consolidated_waves$annotation=="Trigger#1",] -> k1y
 consolidated_waves[consolidated_waves$annotation=="Trigger#2",] -> k1n

  # consolidated_waves[consolidated_waves$annotation=="Trigger#3",]
  # consolidated_waves[consolidated_waves$annotation=="Trigger#4",]

  # b1mm <- cbind(rep("Difference",length(b1y$match)),b1y$index,b1n[,c(4:11)]-b1y[,c(4:11)])
  # colnames(b1mm) <- colnames(b1y)[2:11]
  # b1 <- rbind(cbind(b1y[,2:11]),cbind(b1n[,2:11]),b1mm[,1:10])

  
  k1mm <- cbind(rep("Difference",length(k1y$match)),k1y$index,k1n[,c(4:11)]-k1y[,c(4:11)])
  colnames(k1mm) <- colnames(k1y)[2:11]
  k1 <- rbind(cbind(k1y[,2:11]),cbind(k1n[,2:11]),k1mm[,1:10])
  
  rm(b1y,b1n,k1y,k1n,b2y,b2n,k2y,k2n)
  rm(b1mm,k1mm,b2mm,k2mm)
  
  
  # b1$ms <- seq(from=-150,to=1000,by=2)
  k1$ms <- seq(from=-150,to=1000,by=2)
  # b2$ms <- seq(from=-150,to=1000,by=2)
  # k2$ms <- seq(from=-150,to=1000,by=2)
  
  # b1$id <- "Bouba1"
  k1$id <- "Kiki1"
  # b2$id <- "Bouba2"
  # k2$id <- "Kiki2"
  
  # bk <- rbind(b1,k1,b2,k2)
  # bk <- rbind(b1,k1)
  bk <- k1
    
  
  output$channel_plot <- renderPlot({
    channel = "Cz"
    channel = input$channel
    
    x.min = 300
    x.min <-  input$range[1]
    
    x.max = 500
    x.max <- input$range[2]
    ggplot(bk,aes_string(x="ms",y=channel,col="match")) + geom_line() + facet_grid(~id) + geom_vline(xintercept =x.min,colour="orange") + geom_vline(xintercept = x.max,colour="black")
    })
  
  
  output$statement <- renderText({paste0("This table can be downloaded by clicking the 'Download Result' button ")})
  
  output$downloadresult <- downloadHandler(
    filename = function(){
      paste0("Result_",edf_server$name,"_",Sys.Date(),".csv")
    },
    content = function(file) {
      write.csv(bk,file) 
    }
  )
  
  output$result_table <- renderTable({
    
    channel = "Cz"
    channel = input$channel
    
    x.min = 300
    x.min <-  input$range[1]
    
    x.max = 500
    x.max <- input$range[2]
    
    columns <- c(channel,"id","match")
    mean_min <- bk[bk$ms>=x.min & bk$ms <= x.max,columns] %>% group_by(id,match) %>% summarize_all(funs(mean,min))
    
    mean_bk <- reshape2::dcast(data=mean_min,id~match,value.var="mean")
    colnames(mean_bk) <- c("ID","Mean of Difference","Mean of Mismatch","Mean of Match")
    min_bk <- reshape2::dcast(data=mean_min,id~match,value.var="min")
    colnames(min_bk) <- c("ID","Min of Difference","Min of Mismatch","Min of Match")
    mean_min <- cbind(mean_bk,min_bk[,2:4])
    min_index <- bk[bk$ms>=x.min & bk$ms <= x.max & bk$match=="Difference",columns] %>% group_by(id,match) %>% summarize_all(funs(which.min))
    min_index$ms <- x.min + (min_index[,3]-1)*2
    mean_min$Latency = min_index[,4]
    mean_min
  })
  
})
})