library(shinyalert) # pop-up error messages
library(shinycssloaders) # spinning plot loading icon
library(icesTAF) # contains mkdir
library(tidyverse)
library(shiny)

# ui <- basicPage(
#     plotOutput("plot1",
#                click = "plot_click",
#                dblclick = "plot_dblclick",
#                brush = "plot_brush"
#     ),
#     verbatimTextOutput("data"),
#     tableOutput("data_table")
# )
# 
# server <- function(input, output) {
#     values <- reactiveValues()
#     #values$df_all <- mtcars
# 
#     output$plot1 <- renderPlot({
#         mtcars %>%
#             ggplot(aes(x = wt, y = mpg)) +
#             geom_point()+
#             facet_wrap(~cyl)
#     })
# 
# 
#   values$new_vals <- tibble(dye = mtcars$cyl, assignment = rep('none', times = length(mtcars$cyl))) %>% 
#                           distinct(dye, .keep_all = TRUE)
# 
#     observeEvent({input$plot_dblclick$panelvar1}, {
#         values$new_vals <- values$new_vals %>%
#                             mutate( assignment = replace(assignment, dye == input$plot_dblclick$panelvar1, "hit")) 
#         })
# 
#         observeEvent({input$plot_click$panelvar1}, {
#           values$new_vals <- values$new_vals %>%
#             mutate(assignment = replace(assignment, dye == input$plot_click$panelvar1, "sensitive")) %>%
#             mutate( assignment = replace(assignment, dye == input$plot_brush$panelvar1, "none")) 
#             })
# 
#         observeEvent(input$plot_brush$panelvar1, {
#           values$new_vals <- values$new_vals %>%
#             mutate( assignment = replace(assignment, dye == input$plot_brush$panelvar1, "none")) 
#                })
# 
#        output$data_table <-  renderTable({  values$new_vals  })
# 
# }
# 
# shinyApp(ui, server)
# 

facet_no_y_theme <- theme( # for the first hit-calling plot, the most austere
    text = element_text(size = 4),
    axis.title = element_blank(), # don't label the axes
    axis.text.x = element_text(), # don't label the numbers on the axes
    axis.ticks = element_blank(), # dont have ticks
    legend.position = "right", # put a legent on the right of the plot
    plot.title = element_text(lineheight=.8, face="bold", size = 12), # have a title
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_blank(),
    panel.spacing.x = unit(0.1, "lines"),
    panel.spacing.y = unit(0.1, "lines"),
    aspect.ratio = (1/1.618)
)

facet_wrap_linetyped2 <- function(df_melt, title, facets_wide) {
    p <- ggplot(df_melt, aes(x = Temperature, # temperature on X
                             y = value, # RFU on y
                             color = channel_f, # colored by the state
                             linetype = type,
                             group = dye_conc_type_channel # group means series, as in, this defines the unique data sets
    )) +
        geom_line(size = 0.3, alpha = 0.8) + # change the line type depending on the dye concentration # linetype = df_melt$conc #
        facet_wrap(~dye, scales = "free", ncol = facets_wide) +
        labs(title = title, color = "Channel") +
        theme_bw() +
        scale_color_manual( values = c("Cy5.5" = "#67000d", "Cy5" = "#a50f15", "ROX" = "#ef3b2c", "TAMRA" = "#f16913", "JOE" = "#74c476", "FAM" = "#2171b5")) +
        scale_linetype_manual(values = c("dashed", "solid")) +
        facet_no_y_theme
    
    p # return the plot
}

save_screen_meta_data <- function (path_in,
                                   protein,
                                   exp_num,
                                   daughter_num,
                                   script_version,
                                   instrument,
                                   plate_type,
                                   plate_lot,
                                   additional_notes,
                                   save_outputs,
                                   buffer_type,
                                   facets_wide
) {
    
    # define and/or create directory logic
    path_int <- paste0(path_in, "/intermediate/") # holds intermediate data files
    path_fin <- paste0(path_in, "/final/") # holds all final materials
    path_fin_facets <- paste0(path_fin, "/hits/") # holds all final materials
    
    #mkdir(path_raw) # if exists, fails silently ("FALSE")
    mkdir(path_int)
    mkdir(path_fin)
    mkdir(path_fin_facets)
    
    # write the read_me files
    # generate a read me  and session Info fils e for this experiment with relevant specs immediately
    fileConn <- file(paste0(path_fin, exp_num, "_", as.character(base::Sys.Date()),"_",  "readme.txt"))
    writeLines(c(paste("Experiment", paste0(exp_num, protein), sep = ": "),
                 paste("Dye daughter plates:", daughter_num, sep = ": "),
                 paste("Buffer used for screen:", buffer_type, sep = ":"),
                 paste("Hit calling script version", script_version, sep = ": "),
                 paste("Instrument used", instrument, sep = ": "),
                 paste("Plate type", plate_type, sep = ": "),
                 paste("Plate lot", plate_lot, sep = ": "),
                 paste("Analysis completed on", as.character(base::Sys.Date()), sep = ": "),
                 additional_notes), fileConn)
    close(fileConn)
    
    writeLines(capture.output(sessionInfo()), paste0(path_fin, exp_num, as.character(base::Sys.Date()),"_","sessionInfo.txt")) # save the session info into the final materials folder
}

make_channel_vec <- function( df ) { # make the vector which specifies channel for each reading
    channels <- df %>%
        group_by(`0`) %>%
        filter(n() == 1) %>%
        select(`0`) %>%
        as_vector()
    
    n_meas <- df %>%  # the number of wells measured per channel (always the same for all channels )
        group_by(`0`) %>%
        filter(n() > 1) %>%
        tally() %>%
        nrow()
    
    rep(channels , each = n_meas + 1) # add one, for the row which will still contain the channel itself
}

read_qTower <- function( file_path ) {
    
    df_raw <- read_csv(file_path,
                       col_names = c(0:500) %>% as.character()) %>% # read in 500 columns; this should exceed any actual run, and fill in columsn to right as NAs
        select_if(~sum(!is.na(.)) > 0) #%>% # remove all columns which are all NAs
    
    
    df <- df_raw %>%
        drop_na( tail(names(.), 1) %>% var() ) %>% # drop the header, which is empty in the tailing cols
        mutate( channel = make_channel_vec(.)) %>% # add channel as a column
        filter(!`0` %in% .$channel) %>%
        rename(well = `0`) %>%
        mutate_at(vars(-well, -channel), as.numeric) %>%
        pivot_longer(-c("well", "channel"), names_to = "Temperature", values_to = "value") %>%
        mutate_at(vars(well, channel), as.character) %>%
        mutate_at(vars(Temperature, value), as.numeric) %>%
        mutate(channel_f = factor(.$channel, levels = c("FAM", "JOE", "TAMRA", "ROX", "Cy5", "Cy5.5", "SyproOrange")))
}

# new daughter layout function
df_to_layout <- function(df, layout_type) {
    df_m <-   set_names( df ,  c("type","row",as.numeric( df [1,-c(1,2)]))) %>%
        . [ -1 , -1] %>%
        reshape2::melt( . ,id.vars = "row") %>%
        mutate( . , well = as_vector(map2( . $row,  . $variable, paste0)) ) %>%
        set_names( . , c("row", "column", layout_type, "well"))
    df_m
}

make_layout <- function( filename ) { # from path to raw layout to a final fomatted layout file
    # read the layout file, and split each layout into an individual
    layout_list <- data.table::fread( filename, header = TRUE) %>%
        as_tibble() %>%
        split( . ,  . $Type)
    
    # put into a merge-able form
    layout <- df_to_layout(layout_list[[1]], names(layout_list)[[1]])[c(1,2,4)] # initialize the list
    for (i in c(1:length(layout_list))) {
        layout <- layout %>%
            mutate("var" =  as_vector(df_to_layout(layout_list[[i]], layout_type = names(layout_list)[[i]])[3] )) %>% # append the column of interest
            set_names(c(names(layout), names(layout_list)[[i]])) # rename based on the column of interest
    }
    layout <- layout %>%
        unite("condition", c(4:ncol(.)), remove = FALSE) %>% # create a unique column, used to define groups after averaging
        mutate_if(is.factor, as.character)
    
    layout
}



# Define server logic required to draw a histogram
server <- function(input, output) {
    
    values <- reactiveValues()
    
    observeEvent(input$analyze_button, {
        protein_file_path <- input$protein_file
        
        # ##### make necessary directories and save the README file
        # save_screen_meta_data(input$path_in,
        #                       input$protein,
        #                       input$exp_num,
        #                       input$daughter_num,
        #                       input$script_version,
        #                       input$instrument,
        #                       input$plate_type,
        #                       input$plate_lot,
        #                       input$additional_notes,
        #                       input$save_outputs,
        #                       input$buffer_type,
        #                       input$facets_wide)
        
        
        ##### read in the layout files
        tryCatch({
            if (input$use_layout_protein == "layout_821") {
                values$protein_layout <- readRDS("layout_821.rds")
                
            } else if (input$use_layout_protein == "layout_715") {
                values$protein_layout <- readRDS("layout_715.rds")
                
            } else if (input$use_layout_protein == "layout_606") {
                values$protein_layout <- readRDS("layout_606.rds")
                
            } else {
                protein_layout_path <- input$protein_layout$datapath
                print(protein_layout_path)
                values$protein_layout <-  make_layout(input$protein_layout$datapath)
            }
            values$protein_layout %>% head()
        }, # end tryCatch function
        error = function(e){
            shinyalert("Protein layout file failed to upload!", "Please check file and try again.")
            values$protein_layout <- NULL
        }
        )
        
        tryCatch({
            if (input$use_layout_buffer == "layout_821") {
                print("buffer layout 821")
                values$buffer_layout <- readRDS("layout_821.rds")
                
            } else if (input$use_layout_buffer == "layout_715") {
                print("buffer layout 715")
                values$buffer_layout <- readRDS("layout_715.rds")
                
            } else if (input$use_layout_buffer == "layout_606") {
                print("buffer layout 606")
                values$buffer_layout <- readRDS("layout_606.rds")
                
            } else {
                print("buffer layout custom")
                buffer_layout_path <- input$buffer_layout$datapath
                print(buffer_layout_path)
                values$buffer_layout <-  make_layout(input$buffer_layout$datapath)
            }
            values$buffer_layout %>% head()
        }, # end tryCatch function
        error = function(e){
            shinyalert("Buffer layout file failed to upload!", "Please check file and try again.")
            values$buffer_layout <- NULL
        }
        )
        
        ##### read in the raw data files
        tryCatch({
            protein_path <- input$protein_file$datapath
            print(protein_path)
            values$protein_file <-  read_qTower(protein_path) %>%
                mutate(type = rep('protein', times = nrow(.)),
                       screen = rep(input$exp_num, times = nrow(.)),
                       Temperature = .$Temperature + (input$start_T-1))
            
            values$protein_file %>% head()
            
        }, # end tryCatch function
        error = function(e){
            shinyalert("Protein data file failed to upload!", "Please check file and try again.")
            values$protein_file <- NULL
        }
        )
        
        tryCatch({
            buffer_path <- input$buffer_file$datapath
            print(buffer_path)
            values$buffer_file <-  read_qTower(buffer_path) %>%
                mutate(type = rep('buffer', times = nrow(.)),
                       screen =rep(input$exp_num, times = nrow(.)),
                       Temperature = .$Temperature + (input$start_T-1))
            
            values$buffer_file %>% head()
        }, # end tryCatch function
        error = function(e){
            shinyalert("Buffer data file failed to upload!", "Please check file and try again.")
            values$buffer_file <- NULL
        }
        )
        
        # join layouts with data files, and join buffer and protein files
        tryCatch({
            values$buffer_file_layout <- full_join(values$buffer_file, values$buffer_layout, by = "well")
            values$protein_file_layout <- full_join(values$protein_file, values$protein_layout, by = "well")
            
            values$df_all <- bind_rows(values$protein_file_layout, values$buffer_file_layout) %>%
                unite(dye_conc_type_channel, c(dye, conc, type, channel), sep = "-", remove = FALSE)
            
            # if (input$save_outputs == TRUE) {
            # 
            # filename = function() {
            #   int_out_path <- paste0(input$path_in, "intermediate/", input$exp_num, "_", as.character(base::Sys.Date()))
            #   df_save_name <- paste0(int_out_path, "_df_all.rds")
            #   df_save_name
            #   #paste(input$dataset, ".csv", sep = "")
            # },
            # content = function(file) {
            #   write.csv(values$df_all, file, row.names = FALSE)
            # }
            # 
            #     # int_out_path <- paste0(input$path_in, "intermediate/", input$exp_num, "_", as.character(base::Sys.Date()))
            #     # df_save_name <- paste0(int_out_path, "_df_all.rds")
            #     #write_rds(values$df_all, df_save_name)
            # }
            
        }, # end tryCatch function
        error = function(e){
            shinyalert("Protein, buffer, and layout files don't match!", "Please check files and try again.")
            values$protein_file_layout <- NULL
            values$buffer_file_layout <- NULL
            values$df_all <- NULL
        }
        )
        
        output$download_files <- downloadHandler(
            filename = function() {
                
                paste0(input$exp_num, "_", as.character(base::Sys.Date()), "_df_all.rds")
            },
            content = function(file) {
                # # # write the read_me
                # fileConn <- file(paste0(input$exp_num, "_", as.character(base::Sys.Date()),"_",  "readme.txt"))
                # writeLines(c(paste("Experiment", paste0(input$exp_num, input$protein), sep = ": "),
                #              paste("Dye daughter plates:", input$daughter_num, sep = ": "),
                #              paste("Buffer used for screen:", input$buffer_type, sep = ":"),
                #              paste("Hit calling script version", input$script_version, sep = ": "),
                #              paste("Instrument used", input$instrument, sep = ": "),
                #              paste("Plate type", input$plate_type, sep = ": "),
                #              paste("Plate lot", input$plate_lot, sep = ": "),
                #              paste("Analysis completed on", as.character(base::Sys.Date()), sep = ": "),
                #              additional_notes), fileConn)
                # close(fileConn)
                
                # write the session info
                # writeLines(capture.output(sessionInfo()), paste0(path_fin, exp_num, as.character(base::Sys.Date()),"_","sessionInfo.txt")) # save the session info into the final materials folder
                
                # write the RDS
                write_rds(values$df_all, file)
            }
        )
        
        output$download_summary_files <- downloadHandler(
            filename = function() {
                
                paste0(input$exp_num, "_", as.character(base::Sys.Date()), "readme.txt")
            },
            content = function(file) {
                writeLines(c(paste("Experiment", paste0(input$exp_num, input$protein), sep = ": "),
                             paste("Dye daughter plates", input$daughter_num, sep = ": "),
                             paste("Buffer used for screen", input$buffer_type, sep = ":"),
                             paste("Hit calling script version", input$script_version, sep = ": "),
                             paste("Instrument used", input$instrument, sep = ": "),
                             paste("Plate type", input$plate_type, sep = ": "),
                             paste("Plate lot", input$plate_lot, sep = ": "),
                             paste("Analysis completed on", as.character(base::Sys.Date()), sep = ": "),
                             paste("Additional notes", input$additional_notes, sep = ": "),
                             paste("Session info", capture.output(sessionInfo()))
                ), 
                file)
            }
        )
        
        output$download_plot <- downloadHandler(
            filename = function() {
                
                paste0(input$exp_num, "_", as.character(base::Sys.Date()), "all_screen_plot.pdf")
            },
            content = function(file) {
                ### make and save the primary plot
                values$plot_title_save <- paste0(input$exp_num, "--", as.character(base::Sys.Date()), "_full_screen_data")
                values$p_save <- facet_wrap_linetyped2(values$df_all, 
                                                       values$plot_title_save, 
                                                       facets_wide = 14)
                
                values$plot_ratio <- values$df_all %>%
                    select(dye) %>%
                    distinct() %>%
                    nrow()
                
                ggsave(file, values$p_save, width = 10, height =  (1/1.618 * values$plot_ratio/14 + 1) )
                
            }
        )
        
        output$save_dt_button <- downloadHandler(
            filename = function() {
                paste0(input$exp_num, "_", as.character(base::Sys.Date()), "_manually_called_hits.csv")
            },
            content = function(file) {
                write.csv(values$new_vals %>% arrange(assignment), file)
                # write_rds(values$df_all, file, row.names = FALSE)
            }
        )
        
        # downloadButton("downloadData", "Download")
        
        
        
    }) # end observe event analyze_button
    
    # observeEvent(values$df_all, {
    #     req(values$df_all)
    #     req(input$save_outputs == TRUE)
    #     ### make and save the primary plot
    #     values$plot_title <- paste0(input$exp_num, "--", as.character(base::Sys.Date()), "_full_screen_data")
    #     values$plot_save_name <- paste0(input$path_in, "final/", values$plot_title, ".pdf")
    #     
    #     values$p_save <- facet_wrap_linetyped2(values$df_all, values$plot_title, facets_wide = 14)
    #     values$plot_ratio <- values$df_all %>%
    #         select(dye) %>%
    #         distinct() %>%
    #         nrow()
    #     
    #     #fin_out_path <- paste0(input$path_in, "final/", input$exp_num, "_", as.character(base::Sys.Date()))
    #     ggsave(values$plot_save_name, values$p_save, width = 10, height =  (1/1.618 * values$plot_ratio/14 + 1) )
    # })
    
    
    #### create the interactive hit-calling plot 
    output$plot <- renderPlot(
        # facet_wrap_linetyped2(values$df_all, title = values$plot_title, facets_wide = 6 )
        facet_wrap_linetyped2(values$df_all %>% filter(dye %in% c("A009", "A010", "A011")), title = values$plot_title, facets_wide = 6 ) # testing, only deal with three wells of data
    ) 
    
    #### plot interaction hit calls
    observeEvent(values$df_all, {
        values$new_vals <- tibble(dye = values$df_all$dye, assignment = rep('none', times = length(values$df_all$dye))) %>%
            distinct(dye, .keep_all = TRUE)
    })
    
    observeEvent({input$plot_dblclick$panelvar1}, {
        values$new_vals <- values$new_vals %>%
            mutate( assignment = replace(assignment, dye == input$plot_dblclick$panelvar1, "hit"))
    })
    
    observeEvent({input$plot_click$panelvar1}, {
        values$new_vals <- values$new_vals %>%
            mutate(assignment = replace(assignment, dye == input$plot_click$panelvar1, "sensitive")) %>%
            mutate( assignment = replace(assignment, dye == input$plot_brush$panelvar1, "none"))
    })
    
    observeEvent(input$plot_brush$panelvar1, {
        values$new_vals <- values$new_vals %>%
            mutate( assignment = replace(assignment, dye == input$plot_brush$panelvar1, "none"))
    })
    
    output$data_table <-  DT::renderDataTable({  values$new_vals %>%
            filter(assignment != "none") %>%
            arrange(assignment)})
    
    
    #### for the validation plot
    output$plot_hit <- renderPlot({
        values$hit_dyes <- values$new_vals %>%
            filter(assignment != "none") %>%
            select(dye) %>%
            as_vector()
        
        facet_wrap_linetyped2(values$df_all %>%
                                  filter( dye  %in%  values$hit_dyes ),
                              title = paste0(values$plot_title, "hits"),
                              facets_wide = 6 )
    })
    
    ### the hit determination
    observeEvent(values$df_all, {
        values$new_vals_hit <- tibble(dye = values$df_all$dye, assignment = rep('none', times = length(values$df_all$dye))) %>%
            distinct(dye, .keep_all = TRUE)
    })
    
    observeEvent({input$plot_hit_dblclick$panelvar1}, {
        values$new_vals_hit <- values$new_vals_hit %>%
            mutate( assignment = replace(assignment, dye == input$plot_hit_dblclick$panelvar1, "validate"))
    })
    
    observeEvent({input$plot_hit_click$panelvar1}, {
        values$new_vals_hit <- values$new_vals_hit %>%
            mutate(assignment = replace(assignment, dye == input$plot_hit_click$panelvar1, "validate_maybe")) %>%
            mutate( assignment = replace(assignment, dye == input$plot_hit_brush$panelvar1, "none"))
    })
    
    observeEvent(input$plot_hit_brush$panelvar1, {
        values$new_vals_hit <- values$new_vals_hit %>%
            mutate( assignment = replace(assignment, dye == input$plot_hit_brush$panelvar1, "none"))
    })
    
    output$data_table_hits <-  DT::renderDataTable({  values$new_vals_hit %>%
            filter(assignment != "none") %>%
            arrange(assignment)})
    
    output$download_files_hit <- downloadHandler(
        filename = function() {
            #int_out_path <- paste0(input$path_in, "intermediate/", input$exp_num, "_", as.character(base::Sys.Date()))
            paste0(input$exp_num, "_", as.character(base::Sys.Date()), "_validation_list.csv")
        },
        content = function(file) {
            
            write_csv(values$new_vals_hit %>%
                          arrange(desc(assignment)), file)
        }
    )
    
    
    # observeEvent(input$save_dt_button, { # write the hit-list to the final folder when this button is pressed
    #   print("saving hit df")
    #   write.csv(values$new_vals %>% arrange(assignment), file = paste0(input$path_in, "final/",input$exp_num, "--", as.character(base::Sys.Date()), "_", input$protein, "_dye_screen", input$daughter_num, "_", "manually_called_hits.csv")) # the reformatted and renamed buffer+dye data
    #   saveRDS(values$new_vals %>% arrange(assignment), file = paste0(input$path_in, "final/" ,input$exp_num, "--", as.character(base::Sys.Date()), "_", input$protein, "_dye_screen", input$daughter_num, "_", "manually_called_hits.rds")) # the reformatted and renamed buffer+dye data
    # })
    
    #### have the option to use a pre-stored background file for various buffers
    
    #### have the option to use pre-stored layout files for various mothers
    
}
#
# Define UI for application that draws a histogram
ui <- navbarPage(useShinyalert(),
                 tabPanel("Instructions",
                          fluidRow(
                              # left panel: interactive sliders
                              column(3),
                              column(6,
                                     tags$h1("Instructions for hit calling"), 
                                     tags$p("Here's the new hit-calling app for the dye screens! It's not in it's final form, but it's useable for now."),
                                     tags$br(),
                                     tags$li(tags$strong("Tab 1: Entering screening information."),
                                             "In the first panel, copy information from the notebook entry for the scren into the requested windows. Having this information associated with the processed files is critical for doing retrospective analyses for things like buffer compatibility, so it's really important that these entires are as accurate as possible. Whatever path is set as the 'Path to Data' is where the outputs will be saved. You can either upload a layout file, or if you used a specific daughter plate, you can select those from the radio buttons!"),
                                     tags$br(),
                                     tags$li(tags$strong("Tab 2: Pick hits."),
                                             "Once you click 'Analyze' you can proceed to the second step, which is visualizing and selecing hits. To select something as a DSF hit, double-click that panel. To select something as protein-sensitive, but not a DSF hit, click once. To remove a dye from the hit list, click and drag anywhere in the window. The most recent version of the selected hits will be shown in a table to the right of the plot. Once you've finished calling hits, you can download your assignments by pressing the 'Save hit list to final folder' button."),
                                     tags$br(),
                                     tags$li(tags$strong("Tab 2: Downloading results."),
                                             "This part is super clunky right now, sorry! Do the following: in the folder where the raw data is kept, make two new folders, one titled 'intermediate' and one titled 'final'. Download the read me, plot, and hit definitions into the final folder. Downlaod the raw data into the intermediate folder. You should be able to leave the file names as-is! When you click the 'download plot' button, it will look like it's not doing anything, but just wait a minute. Not sure why this is so slow, but it does work....."),
                                     tags$br(),
                                     tags$p("Ultimately, this website will do hit calling for you using the hit-calling neural net Zach and Kangway wrote. For now, we're still doing it manually. Good luck!")
                              ),
                              
                              column(3)
                              
                              
                          )),
                 #"UCSF Dye Screen Processing",
                 tabPanel("Enter screening information",
                          column(3,
                                 p("Basic screen information", style = "font-size: 16px;",align = "center"),
                                 textInput("protein", "Protein", "Cedilla_buffer"),
                                 textInput("exp_num", "Experiment number", "Ds0014"),
                                 textInput("buffer", "Buffer used", "TEAD/BCL-6: 50mM Tris pH8.0, 150mM NaCl, 1mM TCEP"),
                                 textAreaInput("additional_notes", "Additional notes","Making sure everything looks ok with Alya's buffer screen in the 0062000 Axygen plates before proceeding to buffer screening.", height = 255)
                          ),
                          
                          column(3,
                                 p("Technical notes", style = "font-size: 16px;",align = "center"),
                                 numericInput("start_T", "Starting temperature (C)", 25),
                                 numericInput("end_T", "Ending temperature (C)", 94),
                                 textInput("daughter_num", "Dye daughter number", "715"),
                                 textInput("script_version", "Analysis script version", "hit_calling_v13.Rmd"),
                                 textInput("instrument", "qPCR instrument used", "qTower384g"),
                                 textInput("plate_type", "Type of PCR plate used", "Axygen PCR-284-LC480WNFBC"),
                                 textInput("plate_lot", "PCR plate lot number", "lot 00620000")
                                 # ,
                                 # checkboxInput("save_outputs", "Save outputs", value = TRUE, width = NULL)
                          ),
                          
                          column(3,
                                 p("Screening files", style = "font-size: 16px;",align = "center"),
                                 #textAreaInput("path_in", "Path to data", "~/Box Sync/data/dye_screens/dye_screening_code/mock_screen/"),
                                 
                                 tags$hr(),
                                 p("Select protein screening data"),
                                 fileInput("protein_file", "Upload protein screen (csv)",
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),
                                 fileInput("protein_layout", "Upload protein layout (csv)",
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),
                                 p("Or:"),
                                 radioButtons("use_layout_protein", "Use daughter layout for protein plate:",
                                              c("Custom layout" = "custom", #selected by default
                                                "Daughter 715" = "layout_715",
                                                "Daughter 821" = "layout_821")),
                                 
                                 tags$hr(),
                                 p("Select buffer screening data"),
                                 # radioButtons("use_buffer_data", "Use pre-existing buffer control data:",
                                 #              c("Custom-matched buffer" = "custom", # selected by default
                                 #                "Buffer P1 with 715" = "layout_715",
                                 #                "Buffer P1 with 821" = "layout_821")),
                                 #
                                 # p("Or"),
                                 fileInput("buffer_file", "Upload buffer screen (csv)",
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),
                                 
                                 fileInput("buffer_layout", "Upload buffer layout (csv)",
                                           accept = c(
                                               "text/csv",
                                               "text/comma-separated-values,text/plain",
                                               ".csv")
                                 ),
                                 p("Or"),
                                 radioButtons("use_layout_buffer", "Use daughter layout for buffer plate:",
                                              c("Custom layout" = "custom", #selected by default
                                                "Daughter 715" = "layout_715",
                                                "Daughter 821" = "layout_821"
                                              )),
                                 
                                 actionButton("analyze_button", "Analyze")
                          )),
                 tabPanel("Pick hits",
                          fluidRow(
                              # left panel: interactive sliders
                              column(9,
                                     wellPanel(id = "facet_plot",
                                               #plotOutput("plot1", brush = "plot_brush", height = "5000px"),
                                               plotOutput("plot",
                                                          click = "plot_click",
                                                          dblclick = "plot_dblclick",
                                                          brush = "plot_brush",
                                                          height = "5000px"
                                               ) %>% withSpinner(color="#525252"),
                                               style = "overflow-y:scroll; max-height: 600px")
                              ),
                              column(3,
                                     downloadButton("download_files", "Download processed data (save in 'intermediate' folder)"), 
                                     downloadButton("download_summary_files", "Download read-me (save in 'final' folder)"), 
                                     downloadButton("download_plot", "Download plot (takes a minute or so, save in 'final' folder)"),
                                     downloadButton("save_dt_button", "Save hit list to final folder"), 
                                     DT::dataTableOutput("data_table")))
                 ),
                 tabPanel("Select dyes to validate",
                          fluidRow(
                              # left panel: interactive sliders
                              column(9,
                                     wellPanel(id = "facet_plot_hit",
                                               #plotOutput("plot1", brush = "plot_brush", height = "5000px"),
                                               plotOutput("plot_hit",
                                                          click = "plot_hit_click",
                                                          dblclick = "plot_hit_dblclick",
                                                          brush = "plot_hit_brush",
                                                          height = "5000px"
                                               ) %>% withSpinner(color="#525252"),
                                               style = "overflow-y:scroll; max-height: 600px")
                              ),
                              column(3,
                                     downloadButton("download_files_hit", "Download list of dyes to validate"), 
                                     downloadButton("download_layout_hits", "Download layout validate plate layout"), 
                                     DT::dataTableOutput("data_table_hits")))
                 )
                 
)

shinyApp(ui, server)


