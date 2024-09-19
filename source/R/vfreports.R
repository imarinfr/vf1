#' @rdname sfa
#' @title Single Field Reporting
#' @description Generates of one-page reports of single field analyses
#' @details
#' \itemize{
#'   \item\code{vfsfa} saves a pdf with one-page reports of single field analyses
#'   \item\code{vfsfashiny} generates interactive one-page reports of single field
#'     analyses based on Shiny
#' }
#' @param vf visual field data
#' @param td the total deviation values. If \code{NULL} (default) then use
#'           visualFields normative values
#' @param tdp the total deviation probability values. If \code{NULL} (default)
#'            then use visualFields normative values
#' @param pd the pattern deviation values. If \code{NULL} (default) then use
#'           visualFields normative values
#' @param pdp the pattern deviation probability values. If \code{NULL} (default)
#'            then use visualFields normative values
#' @param file The pdf file name where to save the one-page reports of single field analysis
#' @param ... other graphical arguments
#' @return No return value
#' @export
vfsfa <- function(vf, td = NULL, tdp = NULL, pd = NULL, pdp = NULL, file, ...) {
  # always sort by ID, eye, date, and time
  vf <- vfsort(vf)
  defmar <- par("mar") # read default par
  defps <- par("ps")
  on.exit(par(mar = defmar, ps = defps)) # reset default par on exit, even if the code crashes
  pdf(file, width = 8.27, height = 11.69)
  par(mar = c(0, 0, 0, 0))
  for(i in 1:nrow(vf)) {
    par(ps = 10)
    scrlist <- mountlayoutsfa()
    vfiter <- vfselect(vf, i)
    screen(scrlist$title)
    filltitle("Single Field Analysis")
    screen(scrlist$info)
    fillinfosfa(vfiter)
    screen(scrlist$vftxt)
    text(0.50, 1, "Sensitivities", adj = c(0.5, 1), font = 2)
    screen(scrlist$tdtxt)
    text(0.50, 1, "Total Deviation", adj = c(0.5, 1), font = 2)
    screen(scrlist$pdtxt)
    text(0.50, 1, "Pattern Deviation", adj = c(0.5, 1), font = 2)
    screen(scrlist$coltxt)
    text(0.50, 1, "Color Scale", adj = c(0.5, 1), font = 2)
    screen(scrlist$vf)
    vfplot(vfiter, td = td, tdp = tdp, pd = pd, pdp = pdp, type = "s", mar = c(0, 0, 0, 0), ps = 8)
    screen(scrlist$td)
    vfplot(vfiter, td = td, tdp = tdp, pd = pd, pdp = pdp, type = "tds", mar = c(0, 0, 0, 0), ps = 8)
    screen(scrlist$pd)
    vfplot(vfiter, td = td, tdp = tdp, pd = pd, pdp = pdp, type = "pd", mar = c(0, 0, 0, 0), ps = 8)
    screen(scrlist$col)
    drawcolscalesfa(getgpar()$colmap$map$probs, getgpar()$colmap$map$cols, ps = 6, ...)
    screen(scrlist$foot)
    fillfoot()
    close.screen(all.screens = TRUE)
  }
  invisible(dev.off())
}

#' @rdname spa
#' @title Series Progession Analysis
#' @description Generation of one-page reports of series progression analyses
#' \itemize{
#'   \item\code{vfspa} saves a pdf with one-page reports of series progression analyses
#'   \item\code{vfspashiny} generates interactive one-page reports of series progression
#'     analyses based on Shiny
#' }
#' @param vf visual field data
#' @param file The pdf file name where to save the one-page reports of single field analysis
#' @param type Type of data to use. It can be `\code{s}`, `\code{td}`, or
#' `\code{pd}`.
#' @param nperm Number of permutations. Default is 7!
#' @param trunc value for the Truncated Product Method (see reference).
#' Default is 1
#' @param testSlope slope, or slopes, to test as null hypothesis. Default is 0.
#' if a single value, then the same null hypothesis is used for all locations.
#' If a vector of values, then (for \code{plr} and \code{poplr}) each
#' location of the visual field will have a different null hypothesis. The length
#' of testSlope must be 1 or equal to the number of locations to be used in the PLR
#' or PoPLR analysis
#' @param ... other graphical arguments
#' @return No return value
#' @references
#' N. O'Leary, B. C. Chauhan, and P. H. Artes. \emph{Visual field progression in
#' glaucoma: estimating the overall significance of deterioration with permutation
#' analyses of pointwise linear regression (PoPLR)}. Investigative Ophthalmology
#' and Visual Science, 53, 2012
#' @export
vfspa <- function(vf, file, type = "td", nperm = factorial(7),
                  trunc = 1, testSlope = 0, ...) {
  # sort
  vf <- vfsort(vf)
  # run regression analyses
  res <- runregressions(vf, type, nperm, trunc, testSlope)
  defmar <- par("mar") # read default par
  defps  <- par("ps")
  on.exit(par(mar = defmar, ps = defps)) # reset default par on exit, even if the code crashes
  pdf(file, width = 8.27, height = 11.69)
  par(mar = c(0, 0, 0, 0))
  par(ps = 10)
  for(i in 1:length(res)) {
    # for each subject/eye
    scrlist <- mountlayoutspa()
    vfeye <- vffilter(vf, !!sym("id") == res[[i]]$id, !!sym("eye") == res[[i]]$eye)
    screen(scrlist$title)
    filltitle("Series Progression Analysis")
    screen(scrlist$info)
    fillinfospa(res[[i]], type, nperm, trunc, testSlope)
    screen(scrlist$sparktxt)
    text(0.50, 1, "Sparklines", adj = c(0.5, 1), font = 2)
    screen(scrlist$inttxt)
    text(0.50, 1, "Baseline", adj = c(0.5, 1), font = 2)
    screen(scrlist$sltxt)
    text(0.50, 1, "Slopes", adj = c(0.5, 1), font = 2)
    screen(scrlist$poplrtxt)
    text(0.50, 1, "PoPLR histograms", adj = c(0.5, 1), font = 2)
    screen(scrlist$coltxt)
    text(0.50, 1, "Color Scale", adj = c(0.5, 1), font = 2)
    screen(scrlist$mstxt)
    text(0.50, 1, "Mean sensitivity (MS)", adj = c(0.5, 1), font = 2)
    screen(scrlist$mdtxt)
    text(0.50, 1, "Mean Deviation (MD)", adj = c(0.5, 1), font = 2)
    screen(scrlist$ghtxt)
    text(0.50, 1, "General Height (GH)", adj = c(0.5, 1), font = 2)
    screen(scrlist$int)
    vfint <- vfselect(vfeye, sel = 1) # get first
    vfint[,getvfcols()] <- plr(vfeye)$int
    vfplot(vfint, type = "tds", mar = c(0, 0, 0, 0), ps = 8)
    screen(scrlist$sl)
    vals <- switch(type, "s" = vfeye, "td" = gettd(vfeye), "pd" = getpd(gettd(vfeye)))
    vfplotplr(vals, mar = c(0, 0, 0, 0), ps = 8)
    screen(scrlist$col)
    drawcolscalespa(getgpar()$progcolmap$b$map$probs, getgpar()$progcolmap$b$map$cols, ps = 6, ...)
    screen(scrlist$spark)
    vfsparklines(vals, mar = c(0, 0, 0, 0), ps = 8)
    screen(scrlist$poplrb)
    drawhist(res[[i]]$poplr, alternative = "GT")
    screen(scrlist$poplrw)
    drawhist(res[[i]]$poplr, alternative = "LT")
    screen(scrlist$ms)
    drawgi(res[[i]]$ms, ylab = "Mean Sensitivity", ps = 8)
    screen(scrlist$md)
    drawgi(res[[i]]$md, ylab = "Mean Deviation", ps = 8)
    screen(scrlist$gh)
    drawgi(res[[i]]$gh, ylab = "General Height", ps = 8)
    screen(scrlist$foot)
    fillfoot()
    close.screen(all.screens = TRUE)
  }
  invisible(dev.off())
}

#' @rdname sfa
#' @export
vfsfashiny <- function(vf, ...) {
  # sort
  vf <- vfsort(vf)
  # user interface
  ui <- fluidPage(
    useShinyjs(),
    titlePanel("Single Field Analysis"),
    sidebarLayout(
      sidebarPanel(
        column(7, selectInput("id",  label = "Patient ID", choices = unique(vf$id), selected = vf$id[1])),
        column(5, selectInput("eye", label = "Eye", choices = NULL)),
        div(dataTableOutput("vfdata"), style = "font-size: 85%")
      ),
      mainPanel(
        column(12, align = "center", htmlOutput("info")),
        br(),
        tabsetPanel(type = "pills",
          tabPanel("Sensitivity",       plotOutput("s"),  plotOutput("dum", height = "60px")),
          tabPanel("Total Deviation",   plotOutput("td"), plotOutput("ctd", height = "60px")),
          tabPanel("Pattern Deviation", plotOutput("pd"), plotOutput("cpd", height = "60px"))
        ),
        column(6, align = "center", actionButton("prevbtn", icon = icon("arrow-left"), "")),
        column(6, align = "center", actionButton("nextbtn", icon = icon("arrow-right"), ""))
      )
    )
  )
  # server
  server <- function(input, output, session) {
    vfdata <- vffilter(vf, !!sym("id") == vf$id[1] & !!sym("eye") == vf$eye[1])
    vfsel  <- reactiveVal(vfselect(vfdata, 1))
    ####################
    # EVENTS
    ####################
    # select Patient ID, update eye options
    # if a new Patient ID is selected
    observeEvent(input$id, {
      vfdata <<- vffilter(vf, !!sym("id") == input$id) # get data from subject
      eyes <- sort(unique(vfdata$eye)) # get eyes for this subject
      vfdata <<- vffilter(vfdata, !!sym("eye") == eyes[1]) # refine selection get data only for first eye
      vfsel(vfselect(vfdata, 1)) # choose first record
      output$vfdata <- rendervftable(vfdata, 1)
      updateSelectInput(session, "eye", choices = eyes)
    })
    # if a different eye is selected
    observeEvent(input$eye, {
      vfdata <<- vffilter(vf, !!sym("id") == input$id & !!sym("eye") == input$eye) # get data for the selected eye
      vfsel(vfselect(vfdata, 1)) # choose first record
      output$vfdata <- rendervftable(vfdata, 1)
    }, ignoreInit = TRUE)
    # update data from selected table row
    observeEvent(input$vfdata_rows_selected, vfsel(vfselect(vfdata, input$vfdata_rows_selected)), ignoreInit = TRUE)
    # previous record
    observeEvent(input$prevbtn, {
      selected <- input$vfdata_rows_selected
      if(selected > 1) selected <- selected - 1
      vfsel(vfselect(vfdata, selected))
      output$vfdata <- rendervftable(vfdata, selected)
    })
    # next record
    observeEvent(input$nextbtn, {
      selected <- input$vfdata_rows_selected
      if(selected < nrow(vfdata)) selected <- selected + 1
      vfsel(vfselect(vfdata, selected))
      output$vfdata <- rendervftable(vfdata, selected)
    })
    ####################
    # OUTPUT
    ####################
    # show patient's info and global indices
    output$info <- renderText(fillinfosfashiny(vfsel()))
    # show the sensitivity plots (if that tab is selected in the UI)
    output$s    <- renderPlot(vfplot(vfsel(), type = "s"))
    # show the total-deviation maps (if that tab is selected in the UI)
    output$td   <- renderPlot(vfplot(vfsel(), type = "tds"))
    # show the pattern-deviation maps (if that tab is selected in the UI)
    output$pd   <- renderPlot(vfplot(vfsel(), type = "pd"))
    # show col scales for TD or PD probability maps
    output$ctd  <- renderPlot(drawcolscalesfa(getgpar()$colmap$map$probs, getgpar()$colmap$map$cols, ps = 8, ...)) # col scale for TD probability maps
    output$cpd  <- renderPlot(drawcolscalesfa(getgpar()$colmap$map$probs, getgpar()$colmap$map$cols, ps = 8, ...)) # col scale for TD probability maps. Redundant but necessary
    output$dum  <- renderPlot({}) # for senstivity values, we need not plot col scale, but need this for consistent formatting
  }
  shinyApp(ui, server)
}

#' @rdname spa
#' @references
#' N. O'Leary, B. C. Chauhan, and P. H. Artes. \emph{Visual field progression in
#' glaucoma: estimating the overall significance of deterioration with permutation
#' analyses of pointwise linear regression (PoPLR)}. Investigative Ophthalmology
#' and Visual Science, 53, 2012
#' @export
vfspashiny <- function(vf, type = "td", nperm = factorial(7),
                       trunc = 1, testSlope = 0, ...) {
  vf <- vfsort(vf)
  # get eyes analyzed 
  vfeye <- unique(data.frame(id = vf$id, eye = vf$eye, stringsAsFactors = FALSE))
  # run regression analyses
  res <- runregressions(vf, type, nperm, trunc, testSlope)
  ui <- fluidPage(
    titlePanel("Series Progression Analysis"),
    sidebarLayout(
      sidebarPanel(
        column(7, selectInput("id",  label = "Patient ID", choices = unique(vfeye$id), selected = unique(vfeye$id)[1])),
        column(5, selectInput("eye", label = "Eye", choices = NULL)),
        br(), br(), br()
      ),
      mainPanel(
        column(12, align = "center", htmlOutput("info")),
        tabsetPanel(type = "pills",
          tabPanel("Baseline", plotOutput("bas"), plotOutput("cbas", height = "60px")),
          tabPanel("PLR", plotOutput("plr"), plotOutput("cplr", height = "60px")),
          tabPanel("Sparklines", plotOutput("skl")),
          tabPanel("PoPLR",
            column(width = 12, align = "center",
              plotOutput("poplrr", width = "250px", height = "250px"),
              plotOutput("poplrl", width = "250px", height = "250px")
            )
          ),
          tabPanel("Global trends",
            column(width = 4, plotOutput("ms", height = "250px")),
            column(width = 4, plotOutput("md", height = "250px")),
            column(width = 4, plotOutput("gh", height = "250px"))
          )
        )
      )
    )
  )

  server <- function(input, output, session) {
    ressel   <- reactiveVal(res[[1]]) # select first subject/eye to show results
    vfseries <- reactiveVal(vf[vf$id == res[[1]]$id & vf$eye == res[[1]]$eye,])
    ####################
    # EVENTS
    ####################
    # select Patient ID, update eye options
    # if a new Patient ID is selected
    # update available eyes for analysis
    observeEvent(input$id, {
      ressel(res[[which(vfeye$id == input$id & vfeye$eye == vfeye$eye[vfeye$id == input$id][1])]])
      updateSelectInput(session, "eye", choices = vfeye$eye[vfeye$id == input$id])
    })
    # update selected eye for which to show the analysis
    observeEvent(input$eye, ressel(res[[which(vfeye$id == input$id & vfeye$eye == input$eye)]]), ignoreInit = TRUE)
    # update vfseries
    observeEvent(ressel(), vfseries(vf[vf$id == ressel()$id & vf$eye == ressel()$eye,]))
    ####################
    # OUTPUT
    ####################
    # show patient's info and analysis
    output$info <- renderText(fillinfospashiny(ressel()))
    # show the baseline plot
    output$bas  <- renderPlot({
      vfint <- vfselect(vfseries(), sel = 1) # get first
      vfint[,getvfcols()] <- plr(vfseries())$int
      vfplot(vfint, type = "tds", mar = c(0, 0, 0, 0))
    })
    # show the pointwise linear regression plot
    output$plr  <- renderPlot(vfplotplr(vfseries(), type, mar = c(0, 0, 0, 0)))
    # show the sparklines plot
    output$skl  <- renderPlot(vfsparklines(vfseries(), type, mar = c(0, 0, 0, 0)))
    # show the probability scales
    output$cbas <- renderPlot(drawcolscalesfa(getgpar()$colmap$map$probs, getgpar()$colmap$map$cols, ps = 8, ...))
    output$cplr <- renderPlot(drawcolscalesfa(getgpar()$progcolmap$b$map$probs, getgpar()$progcolmap$b$map$cols, ps = 8, ...))
    # PoPLR histograms
    output$poplrl <- renderPlot(drawhist(ressel()$poplr, "LT"))
    output$poplrr <- renderPlot(drawhist(ressel()$poplr, "GT"))
    # global indices progression analysis
    output$ms <- renderPlot(drawgi(ressel()$ms, "Mean Sensitivity"))
    output$md <- renderPlot(drawgi(ressel()$md, "Mean Deviation"))
    output$gh <- renderPlot(drawgi(ressel()$gh, "General Height"))
  }
  shinyApp(ui, server)
}

########################################################
# internal helper functions for pdf and shiny reports
########################################################
#' @noRd
mountlayoutsfa <- function() {
  # all the boxes are defined in mm divided by the width and height
  # of the page, also in mm
  boxtitle  <- c(10 / 210, 200 / 210, 267 / 297, 287 / 297)
  boxinfo   <- c(10 / 210,  65 / 210, 130 / 297, 265 / 297)
  boxsenst  <- c(65 / 210, 200 / 210, 260 / 297, 265 / 297)
  boxsens   <- c(65 / 210, 200 / 210, 190 / 297, 263 / 297)
  boxtdt    <- c(65 / 210, 200 / 210, 182 / 297, 185 / 297)
  boxtd     <- c(65 / 210, 200 / 210, 110 / 297, 183 / 297)
  boxpdt    <- c(65 / 210, 200 / 210, 102 / 297, 105 / 297)
  boxpd     <- c(65 / 210, 200 / 210,  30 / 297, 103 / 297)
  boxcolort <- c(65 / 210, 200 / 210,  25 / 297,  28 / 297)
  boxcolor  <- c(65 / 210, 200 / 210,  15 / 297,  25 / 297)
  boxfoot   <- c(10 / 210, 200 / 210,  10 / 297,  15 / 297)
  scr <- split.screen(rbind(boxtitle, boxinfo, boxsenst, boxsens,
                            boxtdt, boxtd, boxpdt, boxpd,
                            boxcolort, boxcolor, boxfoot))
  return(list(title  = scr[1], info = scr[2],
              vftxt  = scr[3], vf   = scr[4],
              tdtxt  = scr[5], td   = scr[6],
              pdtxt  = scr[7], pd   = scr[8],
              coltxt = scr[9], col  = scr[10],
              foot   = scr[11]))
}

#' @noRd
mountlayoutspa <- function() {
  # all the boxes are defined in mm divided by the width and height
  # of the page, also in mm
  boxtitle    <- c( 10 / 210, 200 / 210, 267 / 297, 287 / 297)
  boxinfo     <- c( 20 / 210, 100 / 210, 140 / 297, 265 / 297)
  boxinttxt   <- c(100 / 210, 200 / 210, 260 / 297, 265 / 297)
  boxint      <- c(100 / 210, 200 / 210, 182 / 297, 262 / 297)
  boxsltxt    <- c(100 / 210, 200 / 210, 170 / 297, 175 / 297)
  boxsl       <- c(100 / 210, 200 / 210,  92 / 297, 172 / 297)
  boxcoltxt   <- c(100 / 210, 200 / 210,  84 / 297,  87 / 297)
  boxcol      <- c(100 / 210, 200 / 210,  74 / 297,  82 / 297)
  boxsparktxt <- c( 10 / 210,  60 / 210, 130 / 297, 135 / 297)
  boxSpark    <- c( 10 / 210,  60 / 210,  80 / 297, 125 / 297)
  boxpoplrtxt <- c( 60 / 210, 100 / 210, 130 / 297, 135 / 297)
  boxpoplrb   <- c( 60 / 210, 100 / 210, 103 / 297, 130 / 297)
  boxpoplrw   <- c( 60 / 210, 100 / 210,  75 / 297, 102 / 297)
  boxmstxt    <- c( 10 / 210,  70 / 210,  63 / 297,  68 / 297)
  boxms       <- c( 10 / 210,  70 / 210,  16 / 297,  63 / 297)
  boxmdtxt    <- c( 75 / 210, 135 / 210,  63 / 297,  68 / 297)
  boxmd       <- c( 75 / 210, 135 / 210,  16 / 297,  63 / 297)
  boxghtxt    <- c(140 / 210, 200 / 210,  63 / 297,  68 / 297)
  boxgh       <- c(140 / 210, 200 / 210,  16 / 297,  63 / 297)
  boxfoot     <- c( 10 / 210, 200 / 210,  10 / 297,  15 / 297)
  scr <- split.screen(rbind(boxtitle, boxinfo, boxinttxt, boxint, boxsltxt, boxsl,
                            boxcoltxt, boxcol, boxsparktxt, boxSpark,
                            boxpoplrtxt, boxpoplrb, boxpoplrw,
                            boxmstxt, boxms, boxmdtxt, boxmd, boxghtxt, boxgh,
                            boxfoot))
  return(list(title    = scr[1],  info   = scr[2],
              inttxt   = scr[3],  int    = scr[4],
              sltxt    = scr[5],  sl     = scr[6],
              coltxt   = scr[7],  col    = scr[8],
              sparktxt = scr[9],  spark  = scr[10],
              poplrtxt = scr[11], poplrb = scr[12], poplrw = scr[13],
              mstxt    = scr[14], ms     = scr[15],
              mdtxt    = scr[16], md     = scr[17],
              ghtxt    = scr[18], gh     = scr[19],
              foot     = scr[20]))
}

#' @noRd
filltitle <- function(txt) {
  nvinfo <- getnv()$info
  text(0.50, 1.00, txt, adj = c(0.5, 1), cex = 1.5, font = 2)
  if(isnotempty(nvinfo$perimetry)) perim_txt <- nvinfo$perimetry
  if(isnotempty(nvinfo$perimetry)) perim_txt <- paste(perim_txt, "with the", nvinfo$strategy, "strategy")
  if(isnotempty(perim_txt)) text(0.01, 0.60, toTitleCase(perim_txt), adj = c(0, 1))
  if(isnotempty(nvinfo$size)) { # if there is info about size
    stim_txt <- paste("Stimulus size:", nvinfo$size)
    if(is.numeric(nvinfo$size)) # if numeric assume that value is in degrees
      stim_txt <- paste(stim_txt, "deg.")
    text(0.01, 0.10, stim_txt, adj = c(0, 0))
  }
  if(isnotempty(nvinfo$name)) text(0.99, 0.10, nvinfo$name, adj = c(1, 0))
  segments(0, 0, 1, 0)
}

#' @noRd
fillinfosfa <- function(vf) {
  g  <- getgl(vf)
  gp <- getglp(g)
  # format global indices
  ms  <- ifelse(isnotempty(g$msens), paste(round(g$msens, 1), "dB"), "")
  msp <- ifelse(isnotempty(gp$msens), paste0("(p < ", gp$msens, ")"), "")
  md <- ifelse(isnotempty(g$tmd), paste(round(g$tmd, 1), "dB"),   "")
  mdp <- ifelse(isnotempty(gp$tmd), paste0("(p < ", gp$tmd, ")"),    "")
  psd <- ifelse(isnotempty(g$psd), paste(round(g$psd, 1), "dB"),   "")
  psdp <- ifelse(isnotempty(gp$psd), paste0("(p < ", gp$psd, ")"),    "")
  vfi <- ifelse(isnotempty(g$vfi), paste(round(g$vfi, 1), "%"),   "")
  vfip <- ifelse(isnotempty(gp$vfi), paste0("(p < ", gp$vfi, ")"),    "")
  # format reliability indices
  fpr <- ifelse(isnotempty(vf$fpr), paste(round(100 * vf$fpr), "%"), "")
  fnr <- ifelse(isnotempty(vf$fnr), paste(round(100 * vf$fnr), "%"), "")
  fl <- ifelse(isnotempty(vf$fl), paste(round(100 * vf$fl), "%"),  "")
  # print in the pdf file
  rect(0, 0, 1, 1, col = "#F6F6F6", border = NA)
  text(0.50, 0.98, "Patient Information", adj = c(0.5, 1), cex = 1.2, font = 2)
  text(0.05, 0.90, "Patient ID:\nEye:\nAge:\n\nDate:\nTime:",   adj = c(0, 1), font = 2)
  text(0.95, 0.90, paste0(vf$id, "\n", vf$eye, "\n", vf$age, "\n\n", vf$date, "\n", vf$time), adj = c(1, 1))
  segments(0.05, 0.66, 0.95, 0.66, col = "#BBBBBB")
  text(0.50, 0.64, "Global indices", adj = c(0.5, 1), cex = 1.2, font = 2)
  text(0.05, 0.56, "MS:\n\nMD:\n\nPSD:\n\nVFI:",   adj = c(0, 1))
  text(0.52, 0.56, paste0(ms,  "\n\n", md,  "\n\n", psd,  "\n\n", vfi), adj = c(1, 1))
  text(0.95, 0.56, paste0(msp, "\n\n", mdp, "\n\n", psdp, "\n\n", vfip), adj = c(1, 1))
  segments(0.05, 0.30, 0.95, 0.30, col = "#BBBBBB")
  text(0.50, 0.28, "Reliability indices", adj = c(0.5, 1), cex = 1.2, font = 2)
  text(0.05, 0.20, "False Positives:\n\nFalse Negatives:\n\nFixation Losses:", adj = c(0, 1))
  text(0.80, 0.20, paste0(fpr, "\n\n", fnr, "\n\n", fl),   adj = c(1, 1))
}

#' @noRd
fillinfospa <- function(res, type, nperm, trunc, testSlope) {
  # format global regression results
  # mean sensitivity
  msint <- format(round(res$ms$int, 1),  nsmall = 1)
  mssl <- format(round(res$ms$sl, 2),   nsmall = 2)
  msp <- format(round(res$ms$pval, 3), nsmall = 1)
  idxp <- which(mssl >= 0)
  idxn <- which(mssl < 0)
  mssl[idxp] <- paste("+", mssl[idxp])
  mssl[idxn] <- paste("-", substr(mssl[idxn], 2, nchar(mssl[idxn])))
  # mean deviation
  mdint <- format(round(res$md$int, 1),  nsmall = 1)
  mdsl <- format(round(res$md$sl, 2),   nsmall = 2)
  mdp <- format(round(res$md$pval, 3), nsmall = 1)
  idxp <- which(mdsl >= 0)
  idxn <- which(mdsl < 0)
  mdsl[idxp] <- paste("+", mdsl[idxp])
  mdsl[idxn] <- paste("-", substr(mdsl[idxn], 2, nchar(mdsl[idxn])))
  # general height
  ghint <- format(round(res$gh$int, 1),  nsmall = 1)
  ghsl <- format(round(res$gh$sl, 2),   nsmall = 2)
  ghp <- format(round(res$gh$pval, 3), nsmall = 1)
  idxp <- which(ghsl >= 0)
  idxn <- which(ghsl < 0)
  ghsl[idxp] <- paste("+", ghsl[idxp])
  ghsl[idxn] <- paste("-", substr(ghsl[idxn], 2, nchar(ghsl[idxn])))
  # format PoPLR results
  if(type == "s")  poplrtype <- "Sensitivity"
  if(type == "td") poplrtype <- "Total Deviation"
  if(type == "pd") poplrtype <- "Pattern Deviation"
  sl <- format(round(res$poplr$csl,  1), nsmall = 1)
  pl <- format(round(res$poplr$cslp, 3), nsmall = 1)
  sr <- format(round(res$poplr$csr,  1), nsmall = 1)
  pr <- format(round(res$poplr$csrp, 3), nsmall = 1)
  # print in the pdf file
  rect(0, 0, 1, 1, col = "#F6F6F6", border = NA)
  text(0.50, 0.96, "Patient Information", adj = c(0.5, 1), cex = 1.2, font = 2)
  text(0.05, 0.88, "Patient ID:\nEye:\nDate from:\nDate to:\nAge:", adj = c(0, 1), font = 2)
  text(0.95, 0.88, paste0(res$id, "\n", res$eye, "\n", res$dateStart, "\n", res$dateEnd,
                          "\nfrom ", res$ageStart, " to ", res$ageEnd), adj = c(1, 1))
  segments(0.02, 0.68, 0.98, 0.68, col = "#BBBBBB")
  text(0.50, 0.64, "Global indices", adj = c(0.5, 1), cex = 1.2, font = 2)
  text(0.05, 0.56, "MS:\nMD:\nGH:", adj = c(0, 1), font = 2)
  text(0.95, 0.56, paste0(msint, " ", mssl, " y (p < ", msp, ")\n",
                          mdint, " ", mdsl, " y (p < ", mdp, ")\n",
                          ghint, " ", ghsl, " y (p < ", ghp, ")"), adj = c(1, 1))
  segments(0.02, 0.42, 0.98, 0.42, col = "#BBBBBB")
  text(0.50, 0.38, "PoPLR slope analysis", adj = c(0.5, 1), cex = 1.2, font = 2)
  text(0.05, 0.30, paste0("data type:\npermutations:\ntruncation:\ntest slope:\n\nS (right):\nS (left):"),
       adj = c(0, 1), font = 2)
  text(0.95, 0.30, paste0(poplrtype, "\n", nperm, "\n", trunc, "\n", testSlope, "\n\n",
                          sr, " (p < ", pr, ")\n", sl, " (p < ", pl, ")"), adj = c(1, 1))
}

#' @noRd
drawcolscalesfa <- function(probs, cols, ...) {
  if(!(0 %in% probs)) {
    probs <- c(0, probs)
    cols  <- c("#000000", cols)
  }
  colrgb <- col2rgb(cols) / 255
  txtcol <- rep("#000000", length(probs))
  txtcol[(0.2126 * colrgb[1,]
          + 0.7152 * colrgb[2,]
          + 0.0722 * colrgb[3,]) < 0.4] <- "#FFFFFF"
  pol <- NULL
  y <- c(0.5, 0.5, -0.5, -0.5)
  xini <- (26 - length(probs)) / 2
  xend <- 25 - xini
  pol[1] <- list(data.frame(x = c(xini, xini + 1, xini + 1, xini), y = y))
  for(i in 2:length(probs)) {
    xl <- pol[[i-1]]$x[2]
    xu <- xl + 1
    pol[i] <- list(data.frame(x = c(xl, xu, xu, xl), y = y))
  }
  x <- xini + 1:length(probs)
  y <- rep(0, length(probs))
  par(mar = c(0, 0, 0, 0), ...)
  plot(x, y, typ = "n", ann = FALSE, axes = FALSE,
       xlim = c(1, 25), ylim = c(-0.25, 0.25), asp = 1)
  for(i in 1:length(x)) polygon(pol[[i]], border = NA, col = cols[i])
  text(x - diff(x)[1] / 2, y, probs, col = txtcol)
}

#' @noRd
drawcolscalespa <- function(probs, cols, ...) {
  if(!(0 %in% probs)) {
    probs <- c(0, probs)
    cols  <- c("#000000", cols)
  }
  colrgb <- col2rgb(cols) / 255
  txtcol <- rep("#000000", length(probs))
  txtcol[(0.2126 * colrgb[1,]
          + 0.7152 * colrgb[2,]
          + 0.0722 * colrgb[3,]) < 0.4] <- "#FFFFFF"
  pol <- NULL
  y <- c(0.5, 0.5, -0.5, -0.5)
  xini <- (17 - length(probs)) / 2
  xend <- 16 - xini
  pol[1] <- list(data.frame(x = c(xini, xini + 1, xini + 1, xini), y = y))
  for(i in 2:length(probs)) {
    xl <- pol[[i-1]]$x[2]
    xu <- xl + 1
    pol[i] <- list(data.frame(x = c(xl, xu, xu, xl), y = y))
  }
  x <- xini + 1:length(probs)
  y <- rep(0, length(probs))
  plot(x, y, typ = "n", ann = FALSE, axes = FALSE,
       xlim = c(1, 17), ylim = c(-0.25, 0.25), asp = 1)
  for(i in 1:length(x)) polygon(pol[[i]], border = NA, col = cols[i])
  text(x - diff(x)[1] / 2, y, probs, col = txtcol)
}

#' @noRd
fillfoot <- function() {
  box()
  vf_txt <- paste0("visualFields version ", packageVersion("visualFields"), 
                   " (", packageDate("visualFields"), ")")
  r_txt  <- version$version.string
  text(0.01, 0.50, vf_txt, adj = c(0, 0.5))
  text(0.99, 0.50, r_txt,  adj = c(1, 0.5))
}

#' @noRd
isnotempty <- function(field) {
  if(is.null(field) || is.na(field) || length(field) != 1 || field == "") return(FALSE)
  return(TRUE)
}

# Scripts for series progression analysis
#' @noRd
runregressions <- function(vf, type, nperm, trunc, testSlope) {
  # run all global and PoPLR analyses a priori so that it does not take too long
  # when we move to see the analyses for different eyes
  uid <- unique(data.frame(id = vf$id, eye = vf$eye))
  print("running PoPLR analyses. Please be patient, it is a slow process")
  pb <- txtProgressBar(min = 0, max = nrow(uid), initial = 0, style = 3) # show text progress bar
  res <- lapply(1:nrow(uid), function(i) {
    vfiter <- vffilter(vf, !!sym("id") == uid$id[i], !!sym("eye") == uid$eye[i])
    g  <- getgl(vfiter)
    ms <- glr(g[,c("date", "msens")], testSlope = testSlope)
    md <- glr(g[,c("date", "tmd")], testSlope = testSlope)
    gh <- glr(g[,c("date", "gh")], testSlope = testSlope)
    vals <- switch(type, "s" = vfiter, "td" = gettd(vfiter), "pd" = getpd(gettd(vfiter)))
    res <- poplr(vals, nperm = nperm, trunc = trunc, testSlope = testSlope)
    setTxtProgressBar(pb, i)
    return(list(id = uid$id[i], eye = as.character(uid$eye[i]),
                dateStart = vfiter$date[1], dateEnd = vfiter$date[res$nvisits],
                ageStart = vfiter$age[1], ageEnd  = vfiter$age[res$nvisits],
                poplr = list(int = res$int, sl = res$sl, pval = res$pval,
                             csl = res$csl, cslp = res$cslp, cslall = res$cstats$cslall,
                             csr = res$csr, csrp = res$csrp, csrall = res$cstats$csrall),
                ms = ms, md = md, gh = gh))
  })
  close(pb)
  return(res)
}

#' @noRd
drawgi <- function(res, ylab, ...) {
  xlim <- c(res$years[1], res$years[length(res$years)])
  if(res$sl <= 0) {
    ylim <- c(max(res$pred) - 5, max(res$pred) + 3)
  } else {
    ylim <- c(min(res$pred) - 3, min(res$pred) + 5)
  }
  firstTick <- ceiling(ylim[1])
  tickMarks <- c(firstTick, firstTick + 2, firstTick + 4, firstTick + 6)
  tooLarge  <- res$years[res$data > ylim[2]]
  tooSmall  <- res$years[res$data < ylim[1]]
  defpar <- par(no.readonly = TRUE) # read default par
  on.exit(par(defpar))              # reset default par on exit, even if the code crashes
  par(plt = c( 0.3, 1, 0.3, 1 ), mgp = c( 1.5, 0.5, 0 ), ...)
  plot(res$years, res$data, type = "n", axes = FALSE, ann = FALSE, xlim = xlim, ylim = ylim, ...)
  axis(1, las = 1, tcl = -0.3, lwd = 0.5, lwd.ticks = 0.5)
  axis(2, las = 1, tcl = -0.3, lwd = 0.5, lwd.ticks = 0.5, at = tickMarks)
  years <- res$years
  dat <- res$data
  years[dat < ylim[1] | dat > ylim[2]] <- NA
  dat[dat < ylim[1] | dat > ylim[2]]   <- NA
  points(years, dat, pch = 19, col = "gray50", ...)
  lines(res$years, res$pred, lwd = 3, ...)
  if(length(tooSmall) > 0) points(tooSmall, ylim[1] * rep(1, length(tooSmall)), pch = 10, col = "red", ...)
  if(length(tooLarge) > 0) points(tooLarge, ylim[2] * rep(1, length(tooLarge)), pch = 10, col = "red", ...)
  title(xlab = "years")
  title(ylab = ylab)
}

#' @noRd
drawhist <- function(res, alternative, ...) {
  maxsp <- 5
  sep   <- 1 / 10 # separation between bins in a tenth so we have 100 bins
  if(alternative == "LT") {
    coltxt  <- "#FF0000"
    colhist <- "#FF000080"
    s       <- res$csl
    sp      <- res$cslall
  } else {
    coltxt  <- "#008000"
    colhist <- "#00800040"
    s       <- res$csr
    sp      <- res$csrall
  }
  s  <- s  / length(res$pval) # average by the number of locations here
  sp <- sp / length(res$pval)
  sp[sp > maxsp] <- maxsp
  if(s > maxsp) s <- maxsp
  breaks <- seq(0, maxsp, by = sep)
  ymax <- max(hist(sp, breaks = breaks, plot = FALSE)$density)
  xlim <- c(0, maxsp)
  ylim <- c(0, 1.1 * ymax)
  defpar <- par(no.readonly = TRUE) # read default par
  on.exit(par(defpar))              # reset default par on exit, even if the code crashes
  par(plt = c(0, 1, 0.4, 1), mgp = c(0.6, 0.1, 0), ...)
  hist(sp, breaks = breaks, freq = FALSE, main = "", xlim = xlim, ylim = ylim,
       xlab = "", ylab = "", lty = 0, col = colhist, axes = FALSE, ann = FALSE)
  axis(1, at = 0:maxsp, tcl = -0.2, lwd = 0.5, lwd.ticks = 0.5)
  title(xlab = "S / n")
  lines(c(s, s), 0.8 * c(0, ymax), col = coltxt)
  points(s, 0.8 * ymax, pch = 21, col = coltxt, bg = coltxt)
  if(alternative == "LT") text(x = maxsp, y = ymax, label = "S (left)", adj = c(1, 1), font = 2, ...)
  if(alternative == "GT") text(x = maxsp, y = ymax, label = "S (right)", adj = c(1, 1), font = 2, ...)
}

#' @noRd
fillinfosfashiny <- function(vfsel) {
  # fill information table
  infoTab <- matrix("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;", 4, 13) # leave large space between blocks
  cssTab  <- matrix("", 4, 13)
  infoTab[,1] <- c("<b>Patient ID:</b>", "<b>Eye:</b>", "<b>Date:</b>", "<b>Time:</b>")
  infoTab[,2] <- "&nbsp;" # leave small space between columns
  infoTab[,3] <- c(vfsel$id, vfsel$eye, format(vfsel$date), vfsel$time)
  # global indices
  infoTab[,5] <- c("<b>MS:</b>", "<b>MD:</b>", "<b>PSD:</b>", "<b>VFI:</b>")
  infoTab[,6] <- "&nbsp;" # leave small space between columns
  g    <- getgl(vfsel)
  gp   <- getglp(g)
  ms   <- paste(round(g$msens),     " dB")
  md   <- paste(round(g$tmd),       " dB")
  psd  <- paste(round(g$psd),       " dB")
  vfi  <- paste(round(g$vfi),       "  %")
  msp  <- paste0("(p < ", gp$msens, ")")
  mdp  <- paste0("(p < ", gp$tmd,   ")")
  psdp <- paste0("(p < ", gp$psd,   ")")
  vfip <- paste0("(p < ", gp$vfi,   ")")
  infoTab[,7] <- c(ms, md, psd, vfi)
  infoTab[,8] <- "&nbsp;" # leave small space between columns
  infoTab[,9] <- c(msp, mdp, psdp, vfip)
  # check if text color needs to change (if p-value < 0.05)
  cssTab[c(gp$msens, gp$tmd, gp$psd, gp$vfi) < 0.05, 7:9] <- "color: #fb0000; font-weight: 700;"
  # false positive and negative rates and fixation losses
  fpr <- paste(round(100 * vfsel$fpr), "%")
  fnr <- paste(round(100 * vfsel$fnr), "%")
  fl  <- paste(round(100 * vfsel$fl),  "%")
  infoTab[,11] <- c("<b>Test duration:</b>", "<b>False positives:</b>", "<b>False negatives:</b>", "<b>Fixation losses:</b>")
  infoTab[,12] <- "&nbsp;" # leave small space between columns
  infoTab[,13] <- c(vfsel$duration, fpr, fnr, fl)
  # check if text color needs to change (if FPR, FNR, FL are greater than 15%)
  cssTab[c(0, vfsel$fpr, vfsel$fnr, vfsel$fl) > 0.15, 13] <- "color: #fb0000; font-weight: 700;"
  return(htmlTable(infoTab, css.cell = cssTab,
                   css.table = "margin-top: 0em; margin-bottom: 1em;",
                   align = c("l", "c", "r", "c", "l", "c", "r", "c", "r", "c", "l", "c", "r"),
                   cgroup = c("Patient information", "", "Global indices", "", "Reliability indices"),
                   n.cgroup = c(3, 1, 5, 1, 3),
                   col.rgroup = c("none", "#F7F7F7")))
}

#' @noRd
rendervftable <- function(vfdata, selected)
  renderDataTable(vfdata[,c("date", "time")], server = FALSE,
                  rownames = FALSE, colnames = c("Test date", "Test time"),
                  selection = list(mode = "single", selected = selected),
                  options = list(paging = FALSE, searching = FALSE))

#' @noRd
fillinfospashiny <- function(selres) {
  # fill information table
  infoTab <- matrix("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;", 4, 15) # leave large space between blocks
  cssTab  <- matrix("", 4, 15)
  infoTab[,1] <- c("<b>Patient ID:</b>", "<b>Eye:</b>", "<b>From:</b>", "<b>To:</b>")
  infoTab[,2] <- "&nbsp;" # leave small space between columns
  infoTab[,3] <- c(selres$id, selres$eye, format(selres$dateStart), format(selres$dateEnd))
  # progression global indices
  infoTab[,5] <- c("<b>MS:</b>", "<b>MD:</b>", "<b>PSD:</b>", "")
  infoTab[,6] <- "&nbsp;" # leave small space between columns
  ms  <- paste(round(selres$ms$sl, 2), "dB / y")
  md  <- paste(round(selres$md$sl, 2), "dB / y")
  gh  <- paste(round(selres$gh$sl, 2), "dB / y")
  msp <- paste0("(p < ", round(selres$ms$pval, 3), ")")
  mdp <- paste0("(p < ", round(selres$md$pval, 3), ")")
  ghp <- paste0("(p < ", round(selres$gh$pval, 3), ")")
  infoTab[,7] <- c(ms, md, gh, "")
  infoTab[,8] <- "&nbsp;" # leave small space between columns
  infoTab[,9] <- c(msp, mdp, ghp, "")
  cssTab[c(selres$ms$pval, selres$md$pval, selres$gh$pval, 1) < 0.05, 7:9] <- "color: #fb0000; font-weight: 700;"
  # progression PoPLR
  csl  <- round(selres$poplr$csl)
  csr  <- round(selres$poplr$csr)
  cslp <- paste0("(p < ", round(selres$poplr$cslp, 3), ")")
  csrp <- paste0("(p < ", round(selres$poplr$csrp, 3), ")")
  infoTab[,11] <- c("<b>S left:</b>", "<b>S right:</b>", "", "")
  infoTab[,12] <- "&nbsp;" # leave small space between columns
  infoTab[,13] <- c(csl, csr, "", "")
  infoTab[,14] <- "&nbsp;" # leave small space between columns
  infoTab[,15] <- c(cslp, csrp, "", "")
  cssTab[c(selres$poplr$cslp, 1, 1, 1) < 0.05, 12:15] <- "color: #fb0000; font-weight: 700;"
  cssTab[c(1, selres$poplr$csrp, 1, 1) < 0.05, 12:15] <- "color: #228b22; font-weight: 700;"
  return(htmlTable(infoTab, css.cell = cssTab,
                   css.table = "margin-top: 0em; margin-bottom: 1em;",
                   align = c("l", "c", "r", "c", "l", "c", "r", "c", "r", "c", "l", "c", "r", "c", "r"),
                   cgroup = c("Patient information", "", "Global slopes", "", "PoPLR"),
                   n.cgroup = c(3, 1, 5, 1, 5),
                   col.rgroup = c("none", "#F7F7F7")))
}