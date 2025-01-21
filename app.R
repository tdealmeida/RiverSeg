## packages necessaries ----
library(shiny)
library(shinyBS)
library(changepoint)
# devtools::install_github("lvaudor/hubr")
library(hubr)
library(dplyr)
library(RPostgres)
library(ecp)
library(cpm)
library(reticulate)
library(zoo)
library(Rbeast)
library(readr)
library(sf)
library(leaflet)
library(ggplot2)
library(EnvCpt)
library(kcpRS)
library(cumSeg)
library(trend)
library(bfast)



##### fonction #####
# Fonction récursive pour PETTITT
pettitt_test <- function(data, index_offset, alpha, min_seg) {
  result <- list()
  
  # Vérifier que le segment est assez long
  if (length(data) < min_seg) {
    result$change_points <- NULL
    return(result)
  }
  
  # Appliquer le test de Pettitt
  test <- pettitt.test(data)
  
  # Si le test est significatif, un point de rupture est détecté
  if (test$p.value < alpha) {
    cp <- test$estimate
    cp_global <- cp + index_offset
    
    # Appliquer récursivement la fonction aux segments gauche et droit
    left <- pettitt_test(data[1:cp], index_offset = index_offset, alpha = alpha, min_seg = min_seg)
    right <- pettitt_test(data[(cp + 1):length(data)], index_offset = cp_global, alpha = alpha, min_seg = min_seg)
    
    # Combiner les points de rupture trouvés
    result$change_points <- c(left$change_points, cp_global, right$change_points)
  } else {
    result$change_points <- NULL
  }
  
  return(result)
}

# Fonction récursive pour BUISHANDS U test
buishand_U_test <- function(data, index_offset, alpha, min_seg) {
  result <- list()
  
  # Vérifier que le segment est assez long
  n <- length(data)
  if (n < min_seg) {
    result$change_points <- NULL
    return(result)
  }
  
  # Appliquer le test de Buishand
  test <- bu.test(data)
  
  # Si le test est significatif, un point de rupture est détecté
  if (test$p.value < alpha) {
    cp <- test$estimate
    cp_global <- cp + index_offset - 1  # Ajuster pour l'indice global
    
    # Éviter de détecter des points de rupture aux extrémités
    if (cp <= min_seg || (n - cp) < min_seg) {
      result$change_points <- NULL
      return(result)
    }
    
    # Appliquer récursivement la fonction aux segments gauche et droit
    left <- buishand_U_test(data[1:cp], index_offset = index_offset, alpha = alpha, min_seg = min_seg)
    right <- buishand_U_test(data[(cp + 1):n], index_offset = cp_global + 1, alpha = alpha, min_seg = min_seg)
    
    # Combiner les points de rupture trouvés
    result$change_points <- c(left$change_points, cp_global, right$change_points)
  } else {
    result$change_points <- NULL
  }
  
  return(result)
}

# Fonction récursive pour Buishands Range test
buishand_range_test <- function(data, index_offset = 1, alpha = 0.05, min_seg = 10) {
  result <- list()
  
  # Vérifier que le segment est assez long
  n <- length(data)
  if (n < min_seg) {
    result$change_points <- NULL
    return(result)
  }
  
  # Appliquer le test de Buishand avec 'br.test'
  test <- br.test(data)
  
  # Si le test est significatif, un point de rupture est détecté
  if (test$p.value < alpha) {
    cp <- test$estimate
    cp_global <- cp + index_offset - 1  # Ajuster pour l'indice global
    
    # Éviter de détecter des points de rupture aux extrémités
    if (cp <= min_seg || (n - cp) < min_seg) {
      result$change_points <- NULL
      return(result)
    }
    
    # Appliquer récursivement la fonction aux segments gauche et droit
    left <- buishand_range_test(data[1:cp], index_offset = index_offset, alpha = alpha, min_seg = min_seg)
    right <- buishand_range_test(data[(cp + 1):n], index_offset = cp_global + 1, alpha = alpha, min_seg = min_seg)
    
    # Combiner les points de rupture trouvés
    result$change_points <- c(left$change_points, cp_global, right$change_points)
  } else {
    result$change_points <- NULL
  }
  
  return(result)
}

# Fonction récursive pour SNHT
compute_snht <- function(data) {
  n <- length(data)
  mean_total <- mean(data)
  sigma2 <- var(data)
  T_values <- numeric(n)
  
  for (k in 2:(n - 2)) {
    n1 <- k
    n2 <- n - k
    mean1 <- mean(data[1:k])
    mean2 <- mean(data[(k + 1):n])
    T_values[k] <- n1 * (mean1 - mean_total)^2 / sigma2 + n2 * (mean2 - mean_total)^2 / sigma2
  }
  
  return(T_values)
}
snht_test <- function(data, index_offset, alpha, min_seg) {
  result <- list()
  
  # Vérifier que le segment est assez long
  n <- length(data)
  if (n < (2 * min_seg)) {  # Besoin d'au moins 2*min_seg pour diviser le segment
    result$change_points <- NULL
    return(result)
  }
  
  # Calculer les statistiques SNHT
  T_values <- compute_snht(data)
  
  # Trouver la position du maximum de la statistique de test
  T_max <- max(T_values, na.rm = TRUE)
  cp <- which.max(T_values)
  
  # Calculer la valeur critique pour le niveau de signification
  # La distribution de T sous H0 est approximativement une distribution du chi-deux à 1 degré de liberté
  critical_value <- qchisq(1 - alpha, df = 1)
  
  # Vérifier si la statistique dépasse la valeur critique
  if (T_max > critical_value) {
    cp_global <- cp + index_offset - 1  # Ajuster pour l'indice global
    
    # Éviter de détecter des points de rupture aux extrémités
    if (cp <= min_seg || (n - cp) < min_seg) {
      result$change_points <- NULL
      return(result)
    }
    
    # Appliquer récursivement la fonction aux segments gauche et droit
    left <- snht_test(data[1:cp], index_offset = index_offset, alpha = alpha, min_seg = min_seg)
    right <- snht_test(data[(cp + 1):n], index_offset = cp_global + 1, alpha = alpha, min_seg = min_seg)
    
    # Combiner les points de rupture trouvés
    result$change_points <- c(left$change_points, cp_global, right$change_points)
  } else {
    result$change_points <- NULL
  }
  
  return(result)
}




##### appli #####

# reticulate::py_config()
# reticulate::virtualenv_create("r-reticulate")
# reticulate::virtualenv_install("r-reticulate", packages = c("numpy", "pandas", "ruptures"))
# 
# ruptures <- import("ruptures")



con <- DBI::dbConnect(RPostgres::Postgres(),
                      host = Sys.getenv("DBMAPDO_HOST"),
                      port = Sys.getenv("DBMAPDO_PORT"),
                      dbname = Sys.getenv("DBMAPDO_NAME"),
                      user      = Sys.getenv("DBMAPDO_USER"),
                      password  = Sys.getenv("DBMAPDO_PASS"))


query <- "SELECT * FROM network_metrics"

Data <- sf::st_read(dsn = con, query = query)

dbDisconnect(con)




str(Data)
Data <- arrange(Data,measure)
# varnum=colnames(Data)[which(purrr::map_df(Data, class)=="numeric")]
varnum=colnames(Data)

# Data <- conn()  # obtention des données 
# varnum=colnames(Data) #récupération des colonnes 


## définir l'interface utilisateur Shiny
ui <- fluidPage(
  titlePanel("Segmentation & détection de rupture"),
  
  sidebarLayout(
    sidebarPanel(
      bsCollapse(
        id = "collapseExample",
        multiple = TRUE,  # Permettre plusieurs panneaux ouverts
        
        # Panel Import des Données
        bsCollapsePanel(
          title = "Import des Données",
          radioButtons("Data", "Type de données",
                       choices = c("Exemple", "Simulations", "Import de données")),
          
          checkboxInput("log", "Transformation logarithmique", FALSE),
          
          # Si 'Exemple' est sélectionné
          conditionalPanel(
            condition = "input.Data == 'Exemple'",
            selectInput("riviere", "Sélectionnez la rivière", 
                        choices = sort(unique(Data$toponyme)), 
                        selected = "le Drac"),
            selectInput("metrique", "Sélectionnez la variable", 
                        choices = varnum,
                        selected = "active_channel_width")
          ),
          
          # Si 'Simulations' est sélectionné
          conditionalPanel(
            condition = "input.Data == 'Simulations'",
            numericInput("N", "Longueur de la série (N)", value = 1000),
            numericInput("K", "Nombre de segments (K)", value = 5),
            numericInput("Theta", "Theta (θ)", value = 3),
            numericInput("Ratio", "Ratio de Sigma / Theta (r)", value = 1),
            numericInput("heterosced", "Hétéroscédasticité (H)", value = 1),
            sliderInput("Delta", "Distance inter-segments x N (δ)", 
                        min = 0, max = 0.05, value = 0, step = 0.01)
          ),
          
          # Si 'Import de données' est sélectionné
          conditionalPanel(
            condition = "input.Data == 'Import de données'",
            fileInput("file1", "Choisir un fichier CSV", 
                      accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
            checkboxInput("header", "Avec en-têtes", TRUE),
            textInput("na.strings", "Données manquantes", value = NA),
            selectInput("sep", "Séparation de colonnes", 
                        choices = c(";", ",", ".", "\t"), selected = ","),
            selectInput("dec", "Séparation des décimales", 
                        choices = c(".", ","), selected = "."),
            # selectInput("selectedCol", "Choisir une variable", choices = NULL),
            uiOutput("colonneUI")
          )
        ),
        
        # Panel Paramètres de Segmentation
        bsCollapsePanel(
          title = "Paramètres de Segmentation",
          selectInput("methode_segmentation", "Méthode de segmentation",
                      choices = c("Changepoint", "Hubert", "E.divisive", "E.agglo", 
                                  "CPM", "GGS", "EnvCpt", "Rbeast", "Cumseg", "KcpRS", 
                                  "Pettitt", "Buishands U test", "Buishands Range test",
                                  "Bfast", "SNHT", "Ruptures (python)"),
                      selected = "Changepoint"),
          
          # Si 'Changepoint' est sélectionné
          conditionalPanel(
            condition = "input.methode_segmentation == 'Changepoint'",
            selectInput("methode", "Algorithme", 
                        choices = c("PELT", "BinSeg", "AMOC", "CTM", "SegNeigh"), 
                        selected = "PELT"),
            radioButtons("type_changement", "Type de changement", 
                         choices = c("mean", "var", "meanvar"), selected = "mean"),
            selectInput("penalty", "Pénalité", 
                        choices = c("None", "SIC", "BIC", "MBIC", "AIC", 
                                    "Hannan-Quinn", "Asymptotic", "Manual", "CROPS"),
                        selected = "Manual")
          ),
          
          # Si 'SegNeigh' est sélectionné
          conditionalPanel(
            condition = "input.methode == 'SegNeigh'",
            numericInput("Q", "Paramètre Q", value = 1)
          ),
          
          # Si 'Manual' est sélectionné comme pénalité
          conditionalPanel(
            condition = "input.penalty == 'Manual'",
            selectInput("penalty_manual", "Choisir la pénalité manuelle :", 
                        choices = list("Pénalité par défaut" = "default", 
                                       "Pénalité log(n)" = "log (n)", 
                                       "Autre" = "other"))
          ),
          
          # Si 'log (n)' est sélectionné
          conditionalPanel(
            condition = "input.penalty_manual == 'log (n)'",
            sliderInput("custom_penalty_log", "Pénalité personnalisée log(n) x :", 
                        min = 1, max = 50, value = 5, step = 1)
          ),
          
          # Si 'other' est sélectionné
          conditionalPanel(
            condition = "input.penalty_manual == 'other'",
            numericInput("custom_penalty_other", "Pénalité personnalisée :", 
                         min = 0, value = 5, step = 1)
          ),
          
          # Si 'CROPS' est sélectionné
          conditionalPanel(
            condition = "input.penalty == 'CROPS'",
            sliderInput("penalty_crops", "Valeur de pénalité (pen.value)", 
                        min = 0, max = 10000, value = c(1, 1000)),
            numericInput("ncpts", "Nombre de segments optimal", value = NULL)
          ),
          
          # Si 'Asymptotic' est sélectionné
          conditionalPanel(
            condition = "input.penalty == 'Asymptotic'",
            sliderInput("penalty_asymptotic", "Valeur de pénalité (pen.value)", 
                        min = 0.01, max = 1, value = 0.05)
          ),
          
          # Si 'Hubert' est sélectionné
          conditionalPanel(
            condition = "input.methode_segmentation == 'Hubert'",
            numericInput("alpha", "Alpha", value = 0.05)
          ),
          
          # Si 'E.divisive' est sélectionné
          conditionalPanel(
            condition = "input.methode_segmentation == 'E.divisive'",
            numericInput("sig.lvl", "Niveau de signification (sig.lvl)", value = 0.05),
            numericInput("R", "Nombre de permutations (R)", value = 199),
            numericInput("k", "Nombre de ruptures (k)", value = NULL),
            numericInput("min.size", "Taille minimale (min.size)", value = 30),
            numericInput("alpha", "Alpha", value = 1)
          ),
          
          # Si 'E.agglo' est sélectionné
          conditionalPanel(
            condition = "input.methode_segmentation == 'E.agglo'",
            numericInput("alpha_agglo", "Alpha", value = 1)
          ),
          
          # Si 'CPM' est sélectionné
          conditionalPanel(
            condition = "input.methode_segmentation == 'CPM'",
            selectInput("cpmType", "Type de CPM (cpmType)", 
                        choices = c("Student", "Barlett", "GLR", "Exponential", 
                                    "GLRAdjusted", "FET", "Mann-Whitney", 
                                    "Mood", "Lepage", "Kolmogorov-Smirnov", 
                                    "Cramer-Von-Mises"), selected = "Student"),
            sliderInput("ARL0", "ARL0", min = 400, max = 50000, value = 3000, step = 100),
            sliderInput("startup", "Startup", min = 20, max = 1000, value = 20)
          ),
          
          # Si 'KcpRS' est sélectionné
          conditionalPanel(
            condition = "input.methode_segmentation == 'KcpRS'",
            numericInput("wsize", "Taille de la fenêtre (wsize)", value = 5),
            numericInput("nperm", "Nombre de permutations (nperm)", value = 199),
            numericInput("kmax_kcp", "Kmax", value = 30),
            numericInput("alpha_kcp", "Alpha", value = 0.05)
          ),
          
          # Si 'Pettitt' est sélectionné
          conditionalPanel(
            condition = "input.methode_segmentation == 'Pettitt'",
            numericInput("alpha_pettitt", "Alpha", value = 0.05),
            numericInput("min_seg_pettitt", "Segment minimum", value = 5)
          ),
          
          # Si 'Buishands U test' est sélectionné
          conditionalPanel(
            condition = "input.methode_segmentation == 'Buishands U test'",
            numericInput("alpha_buishand_U", "Alpha", value = 0.05),
            numericInput("min_seg_buishand_U", "Segment minimum", value = 5)
          ),
          
          # Si 'Buishands Range test' est sélectionné
          conditionalPanel(
            condition = "input.methode_segmentation == 'Buishands Range test'",
            numericInput("alpha_buishand_range", "Alpha", value = 0.05),
            numericInput("min_seg_buishand_range", "Segment minimum", value = 5)
          ),
          
          # Si 'SHNT' est sélectionné
          conditionalPanel(
            condition = "input.methode_segmentation == 'SNHT'",
            numericInput("alpha_snht", "Alpha", value = 0.05),
            numericInput("min_seg_snht", "Segment minimum", value = 5)
          ),
          
          # Si 'Bfast' est sélectionné
          conditionalPanel(
            condition = "input.methode_segmentation == 'Bfast'",
            numericInput("h_bfast", "Valeur de h", value = 0.05)
          )
          
          # # Si 'Ruptures (Python)' est sélectionné
          # conditionalPanel(
          #   condition = "input.methode_segmentation == 'Ruptures (python)'",
          #   selectInput("algo_ruptures", "Algorithmes",
          #               choices = c("Pelt",
          #                           "Binary Segmentation",
          #                           "Bottom-Up",
          #                           "Window-Based"),
          #               selected = "Pelt"),
          #   selectInput("model_ruptures", "Model",
          #               choices = c("l1","l2", "normal", "rbf", "linear","cosine","clinear","rank","mahalanobis","ar"),
          #               selected = "l2"),
          #   numericInput("pen_ruptures", "Penalty", value = 3)
          # )
          
        )
      ),
      
      # Afficher les cartes géologiques et aériennes si Exemple est sélectionné
      conditionalPanel(
        condition = "input.Data == 'Exemple'",
        checkboxInput("show_geology", "Afficher la carte géologique", FALSE),
        checkboxInput("show_satellite", "Afficher la carte aérienne", FALSE)
      ),
      
      hr(),
      checkboxInput("segment_button", HTML("<b>Activer la segmentation</b>"), FALSE),
      
      hr(),
      h4("Code de segmentation :"),
      verbatimTextOutput("code_display"),
      
      actionButton("copy_button", "Copier le code"),
      
      tags$script(HTML("
    Shiny.addCustomMessageHandler('copyText', function(code) {
      var textArea = document.createElement('textarea');
      textArea.value = code;
      document.body.appendChild(textArea);
      textArea.select();
      document.execCommand('copy');
      document.body.removeChild(textArea);
    });
  ")),
      
    ),
    
    #graphique principal à droite
    mainPanel(
      tabsetPanel( 
        tabPanel("Segmentation",
                 plotOutput("plot"),
                 
                 conditionalPanel(
                   condition = "input.Data == 'Simulations'",
                   uiOutput("tableaux")
                 ),
                 conditionalPanel(
                   condition = "input.Data == 'Import de données'",
                   tableOutput("dataPreview")
                 ),
                 conditionalPanel(
                   condition = "input.Data == 'Exemple'",
                   leafletOutput("map")
                 )),
        
        tabPanel("doc")
        
      )
    )
    
  ) # End of sidebarLayout
) # End of fluidPage



### définir le serveur Shiny 
server <- function(input, output,session) {
  
  real_data = reactive({  # reactive pour obtenir les Exemple
    
    data <- Data[Data$toponyme == input$riviere, ]
    while (anyNA(data[1,])) {      # Supprimer les lignes au début et à la fin contenant des NA
      data <- data[-1,]
    }
    while (anyNA(data[nrow(data),])) {
      data <- data[-nrow(data),]
    }
    data <- data[rowSums(!is.na(data)) >= 2, ]      # Supprimer les lignes qui n'ont pas au moins deux éléments valides pour l'interpolation
    data <- data %>% mutate_if(is.numeric, na.approx)      # Appliquer l'interpolation
  })
  
  simulate_data = reactive({  # reactive pour obtenir les Simulations
    
    simulate <- function(N, K, ratio,theta,delta,heterosced) {
      
      margin = ceiling(N / 40) # marge de distance entre les segments
      if (K - 1 > (N - margin - 1)) { # Vérification que le nombre de segments est possible
        k = seq_len(N - margin - 1) # Si ce n'est pas possible, on met un segment par point
      } else { 
        k = numeric(K - 1) # Sinon on génère les segments
        while (TRUE) { # Boucle pour générer les segments
          k <- sort(sample(setdiff(seq_len(N - margin), (N - margin):N), K - 1, replace = FALSE)) # Génération des segments
          if (all(diff(k) >= margin+1)) break # Vérification que l'espacement minimal est respecté
        }
      }
      
      
      # création des segments
      tib_seg <- tibble::tibble( # Création des segments
        segment = 1:K, # Numéro du segment
        mu = rnorm(K, 0, theta),  # Choix de la moyenne de chaque segment en fonction de theta
        sigma = ratio * theta * (1 - heterosced) + ratio * abs(mu) / sqrt(2/pi) * heterosced) # Choix de l'écart type de chaque segment 
      
      # obtention des numeros de segments
      series <- rep(0, N) # Création d'une série de 0
      series[k] <- 1 # Ajout des segments
      series <- cumsum(series) + 1 # cumsum pour avoir les numéros de segments
      
      # création de série entière
      tib <- tibble::tibble( # Création de la série entière
        segment = series, # Numéro du segment
        l = 1:N # Numéro de l'observation
      ) %>% 
        left_join(tib_seg, by = "segment") %>%  # Jointure avec les segments
        mutate(
          x = mu + rnorm(N, 0, sd = sigma) # création de la série entière en ajoutant du bruit (sigma)
        )
      
      # Ajout d'une distance entre les segments en remplaçant par des NA en fonction de N
      last_segment <- tib$segment[1] # Initialisation du dernier segment
      for (i in 2:nrow(tib)) { # Boucle pour ajouter des NA
        if (!is.na(tib$segment[i]) && tib$segment[i] != last_segment 
            && tib$segment[i] != 1 && !is.na(last_segment)) { # Si le segment est différent du dernier segment
          if (delta > 0) { # Si delta est supérieur à 0
            end_row <- min(i + delta * N - 1, nrow(tib)) # Calcul de la fin de la zone à remplacer par des NA
            tib[i:end_row, ] <- NA # Remplacement par des NA
          }
        }
        last_segment <- tib$segment[i] # Mise à jour du dernier segment
      }
      na_indices <- which(is.na(tib$x)) # selection des zones à interpoler
      tib$x <- na.approx(tib$x) # interpolation linéaire des NA
      tib$x[na_indices] <- tib$x[na_indices] + # ajout du bruit dans les zons interpolés
        rnorm(length(na_indices), 0, mean(tib_seg$sigma))
      
      
      # # Calcule de theta chapeau
      # mean_by_segment <- tib %>%
      #   group_by(segment) %>%
      #   summarise(mean_segment = mean(x))
      # global_mean <- mean(mean_by_segment$mean_segment)
      # theta_hat <- sd(mean_by_segment$mean_segment)
      
      # Calcul de la moyenne et de l'écart type global + ratio
      Emu_k <- mean(tib$mu, na.rm = TRUE)
      mu_k <- tib_seg$mu
      Esigma_ktheta <- mean(tib$sigma, na.rm = TRUE) / theta
      Esigma_k <- mean(tib_seg$sigma, na.rm = TRUE)
      sigma_k <- tib_seg$sigma
      theta = theta
      
      # Calcul de nb_k soit le nombre de rupture
      nb_k <- K-1
      
      # Calcul de la longueur des segments  
      segment_sizes <- as.numeric(table(tib$segment))
      min_segment <- which.min(segment_sizes)
      min_segment_length <- segment_sizes[min_segment]
      max_segment <- which.max(segment_sizes)
      max_segment_length <- segment_sizes[max_segment]
      
      
      return(list(data=tib,breakpoints=k))
    }
    
    result <- simulate(N = input$N, K = input$K, ratio = input$Ratio, theta = input$Theta,
                       delta = input$Delta, heterosced = input$heterosced)
  })
  
  
  personal_data = reactive({  # reactive pour obtenir les Import de données
    
    req(input$file1)
    read.csv(input$file1$datapath,
             header = input$header,
             sep = input$sep,
             dec = input$dec,
             na.strings = input$na.strings)
  })
  
  
  output$colonneUI <- renderUI({
    req(personal_data())  # Attendre que les données soient chargées
    colnames_data <- colnames(personal_data())  # Récupérer les noms de colonnes
    selectInput("selectedCol", "Choisissez une colonne de série temporelle", choices = colnames_data)
  })
  
  output$dataPreview <- renderTable({
    req(personal_data())  # Attendre que les données soient chargées
    head(personal_data())
  })
  
  
  get_values = reactive({  # recuperation les métriques + measure
    
    if (input$Data == "Exemple") {  # pour les Exemple
      data <- real_data() # recupere les données depuis la reactive
      values <- data[[input$metrique]]  # prends les valeurs de la metrique 
      measure <- data$measure/1000 # prends les distances de l'axe x
      geom <- data$geom # prends la geometrie pour afficher l'axe
      breakpoints <- NULL
    }
    
    else if (input$Data == "Simulations"){
      result <- simulate_data()
      serie_simulee <- result$data
      values <- serie_simulee$x
      breakpoints <- result$breakpoints # recupere les points de ruptures de référence
      measure <- seq(1, input$N)
      geom <- Data[Data$toponyme == "l'Isère", ]  # ne sert a rien mais obligatoire pour le code
      geom <- geom$geom  # ne sert a rien mais obligatoire pour le code
    }
    
    else if (input$Data == "Import de données"){
      # observe({
      #   req(personal_data())
      #   updateSelectInput(session, "selectedCol", 
      #                     choices = colnames(personal_data()))
      # })
      # 
      # selected_col <- reactiveVal(NULL)      # Valeur réactive pour stocker la colonne sélectionnée
      # 
      # observe({
      #   selected_col(input$selectedCol)
      # })
      # 
      # req(selected_col())
      # data <- personal_data()
      # values <- data[[input$selectedCol]]
      req(input$selectedCol)  # Attendre que l'utilisateur ait sélectionné une colonne
      data <- personal_data()  # Récupérer les données importées
      values <- data[[input$selectedCol]]  # Extraire les valeurs de la colonne sélectionnée
      values <- na.omit(values)  # Supprimer les valeurs manquantes
      
      measure <- 1:length(values)
      breakpoints <- NULL
      geom <- Data[Data$toponyme == "l'Isère", ]
      geom <- geom$geom
    }
    
    if (input$log) {
      values <- log(values+1)    #log transforme les donnes si condition
    }
    
    return(list(values = values, breakpoints = breakpoints, measure = measure, geom = geom))
    
  })
  
  
  get_penalty = reactive({
    
    values <- get_values()$values # récuperation des donné
    
    pen <- NULL
    if (input$penalty == "Manual") { 
      if (input$penalty_manual == "default") {
        pen <- var(values) * log(length(values))
      }
      else if (input$penalty_manual == "log (n)") {
        pen <- log(length(values)) * input$custom_penalty_log
      }
      else if (input$penalty_manual == "other") {
        pen <- input$custom_penalty_other
      }
      
      
    } else if (input$penalty == "Asymptotic") { 
      pen <- input$penalty_asymptotic
    } else if (input$penalty == "CROPS") { 
      pen <- c(input$penalty_crops)
    }
    
    if (is.null(pen)) {
      stop("La pénalité n'a pas été correctement définie.")
    }
    
    return(pen)
    
  })
  
  
  get_segments=reactive({   #reactive pour faire la segmentation selon les méthodes
    
    values <- get_values()$values
    
    if (input$methode_segmentation == "Changepoint"){  
      if (input$type_changement == "mean") {
        seg <- cpt.mean(values, method = input$methode, penalty = input$penalty,Q=50, pen.value = get_penalty(),minseglen = 4)
      }
      else if (input$type_changement == "var") {
        seg <- cpt.var(values, method = input$methode, penalty = input$penalty,Q=50, pen.value = get_penalty(),minseglen = 4)
      }
      else if (input$type_changement == "meanvar") {
        seg <- cpt.meanvar(values, method = input$methode, penalty = input$penalty, Q=50, pen.value = get_penalty(),minseglen = 4)
      }
      seg <- seg
      cpt_points <- seg@cpts
      cpt_points <- cpt_points[-length(cpt_points)]
    }
    
    else if (input$methode_segmentation == "Hubert"){  
      seg = Hubert_segmentation(values,alpha = input$alpha)
      cpt_points <- seg$locations
      cpt_points <- cpt_points[-length(cpt_points)]
      cpt_points <- cpt_points[-1]
    }
    
    else if (input$methode_segmentation == "E.divisive"){  
      seg = e.divisive(matrix(values), sig.lvl = 0.05 , R = input$R, k = NULL, min.size = input$min.size, alpha = input$alpha)
      cpt_points <- seg$estimates
      cpt_points <- cpt_points[-length(cpt_points)]
      cpt_points <- cpt_points[-1]
    }
    
    else if (input$methode_segmentation == "CPM"){  
      seg = processStream(values,cpmType = input$cpmType,ARL0 = input$ARL0, startup = input$startup)
      cpt_points <- seg$changePoints
    } 
    
    else if (input$methode_segmentation == "Rbeast"){
      seg = beast(values,season = "none")
      ncp_mode <- seg$trend$ncp_mode
      cp <- seg$trend$cp
      cpt_points <- cp[1:ncp_mode]
      cpt_points <- sort(cpt_points)
    }
    
    else if (input$methode_segmentation == "Cumseg"){
      seg <- jumpoints(values, output = "2")
      cpt_points <- seg$psi
    }
    
    else if (input$methode_segmentation == "E.agglo"){
      seg <- e.agglo(matrix(values),alpha = input$alpha_agglo, penalty = function(cps){0} )
      cpt_points <- seg$estimates
      cpt_points <- cpt_points[-length(cpt_points)]
      cpt_points <- cpt_points[-1]
    }
    
    else if (input$methode_segmentation == "EnvCpt"){
      seg <- envcpt(values,minseglen = 4)
      cpt_points <- seg[["trendcpt"]]@cpts
      cpt_points <- cpt_points[-length(cpt_points)]
    }
    
    else if (input$methode_segmentation == "KcpRS"){
      seg <- kcpRS(values, RS_fun = runMean, RS_name = "mean", wsize = input$wsize, nperm = input$nperm, Kmax = 30, alpha = 0.05 )
      cpt_points <- seg$changePoints
    }
    
    else if (input$methode_segmentation == "Pettitt"){
      seg <- pettitt_test(values,index_offset = 0, input$alpha_pettitt, input$min_seg_pettitt)
      cpt_points <- sort(unique(seg$change_points))
    }
    
    else if (input$methode_segmentation == "Buishands U test"){
      seg <- buishand_U_test(values,index_offset = 0, input$alpha_buishand_U, input$min_seg_buishand_U )
      cpt_points <- sort(unique(seg$change_points))
    }
    
    else if (input$methode_segmentation == "Buishands Range test"){
      seg <- buishand_range_test(values,index_offset = 0, input$alpha_buishand_range, input$min_seg_buishand_range )
      cpt_points <- sort(unique(seg$change_points))
    }
    
    else if (input$methode_segmentation == "SNHT"){
      seg <- snht_test(values,index_offset = 0, input$alpha_snht, input$min_seg_snht )
      cpt_points <- sort(unique(seg$change_points))
    }
    
    else if (input$methode_segmentation == "Bfast"){
      ts_values <- ts(values)
      seg <- bfast(ts_values, input$h_bfast , season = "none")
      cpt_points <- seg$output[[1]]$bp.Vt$breakpoints
    }
    
    else if (input$methode_segmentation == "Ruptures (python)"){
      
      values <- as.matrix(values)
      values <- r_to_py(values)
      
      if (input$algo_ruptures == "Pelt") {
        seg <- ruptures$Pelt(model = input$model_ruptures)$fit(values)
      }
      else if (input$algo_ruptures == "Binary Segmentation") {
        seg <- ruptures$Binseg(model = input$model_ruptures)$fit(values)
      }
      else if (input$algo_ruptures == "Bottom-Up") {
        seg <- ruptures$BottomUp(model = input$model_ruptures)$fit(values)
      }
      else if (input$algo_ruptures == "Window-Based") {
        seg <- ruptures$Window(model = input$model_ruptures)$fit(values)
      }
      
      cpt_points <- seg$predict(pen = input$pen_ruptures)
      cpt_points <- cpt_points[-length(cpt_points)]
      
    }
    
    return(list(seg = seg, cpt_points = cpt_points))
  })
  
  
  # reactive pour afficher le code de segmentation
  code <- reactive({
    
    methode <- input$methode
    pen <- input$penalty
    penalty_value <- get_penalty()
    alpha_hubert <- input$alpha
    alpha_agglo <- input$alpha_agglo
    alpha_pettitt <- input$alpha_pettitt
    R_value <- input$R
    min_size <- input$min.size
    wsize <- input$wsize
    nperm <- input$nperm
    cpmType <- input$cpmType
    ARL0 <- input$ARL0
    startup <- input$startup
    h_bfast <- input$h_bfast
    
    # Ajoute le code en fonction des méthodes
    if (input$methode_segmentation == "Changepoint"){  
      if (input$type_changement == "mean") {
        paste0( 
          "    cpt.mean(values,\n",
          "      method = \"", input$methode, "\",\n",  # Utilisation de guillemets doubles
          "      penalty = \"", input$penalty, "\",\n",     # Utilisation de guillemets doubles
          "      Q = 50,\n",
          "      pen.value = ", get_penalty(), ",\n",
          "      minseglen = 4)\n")
        
      } else if (input$type_changement == "var") {
        paste0( 
          "    cpt.var(values,\n",
          "      method = \"", input$methode, "\",\n",  # Utilisation de guillemets doubles
          "      penalty = \"", input$penalty, "\",\n",     # Utilisation de guillemets doubles
          "      Q = 50,\n",
          "      pen.value = ", get_penalty(), ",\n",
          "      minseglen = 4)\n")
        
      } else if (input$type_changement == "meanvar") {
        paste0( 
          "    cpt.meanvar(values,\n",
          "      method = \"", input$methode, "\",\n",  # Utilisation de guillemets doubles
          "      penalty = \"", input$penalty, "\",\n",     # Utilisation de guillemets doubles
          "      Q = 50,\n",
          "      pen.value = ", get_penalty(), ",\n",
          "      minseglen = 4)\n")
      }
      
    } else if (input$methode_segmentation == "Hubert"){  
      paste0(
        "    Hubert_segmentation(values,\n",
        "    alpha = ", input$alpha, ")\n")
      
    } else if (methode == "E.divisive") {  
      paste0(
        "    e.divisive(matrix(values),\n",
        "      sig.lvl = 0.05,\n",
        "      R = ", input$R, ",\n",
        "      k = NULL,\n",
        "      min.size = ", input$min.size, ",\n",
        "      alpha = ", input$alpha, ")\n")
      
    } else if (methode == "CPM") {  
      paste0(
        "    processStream(values,\n",
        "      cpmType = \"",input$cpmType, "\",\n",
        "      ARL0 = ", input$ARL0, ",\n",
        "      startup = \"", input$startup, "\")\n")
      
    } else if (methode == "Rbeast") {  
      paste0(
        "    beast(values,\n",
        "      season = \"none\")\n")
      
    } else if (methode == "Cumseg") {  
      paste0(
        "    jumpoints(values,\n",
        "      output = \"2\")\n")
      
    } else if (methode == "E.agglo") {  
      paste0(
        "    e.agglo(matrix(values),\n",
        "      alpha = ", input$alpha_agglo, ",\n",
        "      penalty = function(cps) {0})\n")
      
    } else if (methode == "EnvCpt") {  
      paste0(
        "    envcpt(values,\n",
        "      minseglen = 4)\n")
      
    } else if (methode == "KcpRS") {  
      paste0(
        "    kcpRS(values,\n",
        "      RS_fun = runMean,\n",
        "      RS_name = \"mean\",\n",
        "      wsize = ", wsize, ",\n",
        "      nperm = ", nperm, ",\n",
        "      Kmax = 30,\n",
        "      alpha = 0.05)\n")
      
    } else if (methode == "Pettitt") {  
      paste0(
        "    pettitt_test(values,\n",
        "      index_offset = 0,\n",
        "      alpha = ", input$alpha_pettitt, ",\n",
        "      min_seg_pettitt = ", input$min_seg_pettitt, ")\n")
      
    } else if (methode == "Buishands U test") {  
      paste0(
        "    buishand_U_test(values,\n",
        "      index_offset = 0,\n",
        "      alpha = ", input$alpha_buishand_U, ",\n",
        "      min_seg_buishand_U = ", input$min_seg_buishand_U, ")\n")
      
    } else if (methode == "Buishands Range test") {  
      paste0(
        "    buishand_range_test(values,\n",
        "      index_offset = 0,\n",
        "      alpha = ", input$alpha_buishand_range, ",\n",
        "      min_seg_buishand_range = ", input$min_seg_buishand_range, ")\n")
      
    } else if (methode == "SNHT") {  
      paste0(
        "    snht_test(values,\n",
        "      index_offset = 0,\n",
        "      alpha = ", input$alpha_snht, ",\n",
        "      min_seg_snht = ", input$min_seg_snht, ")\n")
      
    } else if (methode == "Bfast") {  
      paste0(
        "    bfast(ts(values),\n",
        "      h = ", input$h_bfast, ",\n",
        "      season = \"none\")\n")
    }
  })
  
  
  
  # Affichage dynamique dans la fenêtre
  output$code_display <- renderText({
    code()
  })
  
  
  # code pour copier le code afficher
  observeEvent(input$copy_button, {
    showNotification("Code copié dans le presse-papiers !", type = "message")
    session$sendCustomMessage("copyText", code())
  })
  
  
  
  
  
  
  get_objectives=reactive({  ## reactive pour obtenir les objectives/ graph pen ~ nb changepoint
    
    values <- get_values()$values
    
    if (input$methode_segmentation == "Changepoint"){
      Objective <- seq(0, 30, by = 0.5)
      Number_of_Breakpoints <- numeric(length(Objective))
      for (i in seq_along(Objective)) {
        cpt <- cpt.mean(values, method = input$methode, penalty = input$penalty, pen.value=Objective[i])
        Number_of_Breakpoints[i] <- length(cpt@cpts)
      }
      objective = data.frame(Number_of_Breakpoints,Objective)
    }
    
    return(objective)
  })
  
  
  get_data_axis <- reactive({  ## reactive pour aobtenir l'axe x correctement/ les coords
    
    geom <- get_values()$geom
    # geom <- structure(geom,class = "WKB")
    # geom <- st_as_sfc(geom,EWKB = TRUE)
    geom <- st_coordinates(geom)
    df <- as.data.frame(geom)
    colnames(df) <- c("x", "y","id")
    
    x <- tapply(df$x, df$id, mean)
    y <- tapply(df$y, df$id, mean)
    mean_coord <- data.frame(x,y)
    
    filtered_df <- df %>%
      filter(id %in% get_segments()$cpt_points)
    
    breakpoint <- filtered_df %>%
      group_by(id) %>%
      slice(1)
    return(list(df = df, breakpoint = breakpoint,mean_coord = mean_coord))
    
  })
  
  
  output$plot <- renderPlot({  #plot de la segmentation homogenize 
    
    values <- get_values()$values
    breakpoints = get_values()$breakpoints
    breakpoints_signal <- c(1,breakpoints,length(values)+1)
    measure <- get_values()$measure
    cpt <- get_segments()$cpt_points
    cpt_points <- c(1,get_segments()$cpt_points,length(values)+1)
    seg <- get_segments()$seg
    
    data <- data.frame(measure,values)
    
    p <- ggplot(data, aes(x = measure, y = values)) +
      geom_line() +
      theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, color = "black")
      ) +
      xlab("Distance") +
      ylab("Valeur") 
    
    
    if (input$segment_button) {
      
      if (input$Data == "Exemple") {
        p <- ggplot(data, aes(x = measure, y = values)) +
          geom_line() +
          geom_line(aes(y = model_signal(values, cpt_points)), color = "red", size = 1) +
          theme(
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(size = 0.5, color = "black")
          ) +
          xlab("Distance depuis l'exutoire en km") +
          ylab(input$metrique) +
          scale_x_continuous(breaks = seq(min(measure), max(measure), by = 20), labels = function(x) round(x, 1))
        
        
        
        
      } else if (input$Data == "Simulations") {
        p <- ggplot(data, aes(x = measure, y = values)) +
          geom_line() +
          geom_line(aes(y = model_signal(values, cpt_points)), color = "red", size = 1) +
          geom_line(aes(y = model_signal(values, breakpoints_signal)), col = "blue", size = 1) +
          theme(
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(size = 0.5, color = "black")
          ) +
          xlab("Distance") +
          ylab("valeur")+
          scale_x_continuous(breaks = seq(min(measure), max(measure), by = 20), labels = function(x) round(x, 1))
        
        
        
      } else if (input$Data == "Import de données") {
        p <- ggplot(data, aes(x = measure, y = values)) +
          geom_line() +
          geom_line(aes(y = model_signal(values, cpt_points)), color = "red", size = 1) +
          theme(
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(size = 0.5, color = "black")
          ) +
          xlab("Distance") +
          ylab("value") +
          scale_x_continuous(breaks = seq(min(measure), max(measure), by = 20), labels = function(x) round(x, 1))
      }
      
      
      else if (input$methode_segmentation == "EnvCpt") {
        p <-       plot(seg,type='fit')
      }
      
      if (input$penalty == "CROPS") {
        p <- plot(seg, ncpts = input$ncpts, type = "l", main = paste("Détection de points de changement pour la rivière :", input$riviere), xlab = "Distance à l'exutoire (km)", ylab = input$metrique_reelle)
      }
      
      
      
      p <- p + annotate("text", x = min(measure), y = max(values),
                        label = paste("Nb de ruptures:", length(cpt)), # permet d'afficher sur le graph le nombre de rupture de la segmentation
                        hjust = 0, vjust = 1, size = 4, color = "black")
      
    }
    print(p)
    
  })
  
  
  # plot 
  # output$plot2 <- renderPlot({
  #   
  #   seg <- get_segments()$seg
  #   values <- get_segments()$values
  #   
  #   if (input$penalty == "CROPS"){
  #     plot(seg,diagnostic=TRUE)
  #   }
  #   
  # })
  
  output$map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
      setView(lng = 2.2137, lat = 46.2276, zoom = 6)
  })
  
  observe({
    leafletProxy("map") %>%
      clearShapes() %>%
      clearMarkers()
    
    if (input$show_geology) {
      leafletProxy("map") %>%
        clearShapes() %>%
        clearMarkers() %>%
        addWMSTiles(
          baseUrl = "https://geoservices.brgm.fr/geologie",
          layers = "SCAN_F_GEOL1M",
          options = WMSTileOptions(format = "image/png", transparent = TRUE, opacity = 0.2),
          attribution = "&copy; BRGM"
        ) %>%
        addPolylines(data = get_data_axis()$df, lng = ~x, lat = ~y)
      
      if (input$segment_button) {
        leafletProxy("map") %>%
          addCircleMarkers(data = get_data_axis()$breakpoint, lng = ~x, lat = ~y, radius = 5, color = "red")
      }
      
    } else {
      leafletProxy("map") %>%
        clearShapes() %>%
        clearMarkers() %>%
        addProviderTiles(providers$CartoDB.Positron) %>%
        addPolylines(data = get_data_axis()$df, lng = ~x, lat = ~y)
      
      if (input$segment_button) {
        leafletProxy("map") %>%
          addCircleMarkers(data = get_data_axis()$breakpoint, lng = ~x, lat = ~y, radius = 5, color = "red")
      }
    }
  })
  
  
  tableau_df1 <- reactive({
    df <- get_values()$breakpoints
    df <- as.data.frame(df)
    colnames(df)[1] <- "Segmentation de référence (segment bleu)"
    return(df) 
  })
  
  tableau_df2 <- reactive({
    df <- get_segments()$cpt_points
    df <- as.data.frame(df)
    colnames(df)[1] <- "Segmentation de l'algorithme (segment rouge)"
    return(df)
  })
  
  output$tableaux <- renderUI({
    # Vérification si la checkbox "segment_button" est cochée
    if (input$segment_button) {
      fluidRow(
        column(6, tableOutput("tableau1")),  
        column(6, tableOutput("tableau2"))   
      )
    } else {
      # Retourner un UI vide si la checkbox n'est pas cochée
      return(NULL)
    }
  })
  
  output$tableau1 <- renderTable({
    tableau_df1()
  })
  
  output$tableau2 <- renderTable({
    tableau_df2()
  })
  
}


# Lancer l'application Shiny
shinyApp(ui = ui, server = server)

