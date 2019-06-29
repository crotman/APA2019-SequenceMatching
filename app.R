#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)




# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Sequence alignment"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(width = 2,
            textAreaInput("seq1",
                        "Sequência 1",
                        value = "mean"
                      ),
            textAreaInput("seq2",
                      "Sequência 2",
                      value = "name"),
            numericInput("penal_igual",
                      "Penalização mesma classe (vogal/consoante)",
                      value = 1,
                      min = 0,
                      step = 1,
                      max = 10
                      ),
            numericInput("penal_diferente",
                         "Penalização classe diferente (vogal/consoante)",
                         value = 3,
                         min = 0,
                         step = 1,
                         max = 10
            ),
            numericInput("penal_gap",
                         "Penalização gap",
                         value = 2,
                         min = 0,
                         step = 1,
                         max = 10
            ),
            
            uiOutput("passo")
            
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("matches"),
           plotOutput("algoritmo")
           
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {


    add_aresta <- function(arestas_param, origem_x_param, origem_y_param, destino_x_param, destino_y_param)
    {

        add <- tibble(
            origem_x = origem_x_param+1,
            origem_y = origem_y_param+1,
            destino_x = destino_x_param+1,
            destino_y = destino_y_param+1
        )


        arestas_param %>%
            rbind(add)

    }




    troca <- function(x, y, penal_vogal_diff, penal_vogal_igual)
    {

        vogalx <- x %in% c("a","e","i","o","u","A","E","I","O", "U")
        vogaly <- y %in% c("a","e","i","o","u","A","E","I","O", "U")

        if (x != y)
        {
            if(vogalx == vogaly)
            {
                penal_vogal_igual
            }
            else
            {
                penal_vogal_diff
            }
        }
        else
        {
            0
        }
    }


    reconstroix <- function(i, j, caminhox, caminhoy, x)
    {
        print("reconstroix")
        print(paste0(i,"-",j))
        if(i == 0 | j == 0 )
        {
            paste0(str_flatten(x[0:i]), str_flatten(rep("_",j)))
        }
        else if (caminhox[i, j] == 1 & caminhoy[i, j] == 1)
        {
            paste0(reconstroix(i-1, j-1, caminhox, caminhoy, x),x[i])
        }
        else if(caminhox[i, j] == 1 & caminhoy[i, j] == 0)
        {
            paste0(reconstroix(i-1, j, caminhox, caminhoy, x),x[i])
        }
        else if(caminhox[i, j] == 0 & caminhoy[i, j] == 1)
        {
            paste0(reconstroix(i, j-1, caminhox, caminhoy, x),"_")
        }

    }


    reconstroiy <- function(i, j, caminhox, caminhoy, y)
    {

        if(i == 0 | j==0)
        {
            print("reconstroiy")
            print(paste0(i,"-",j))
            paste0(str_flatten(y[0:j]), str_flatten(rep("_",i)))
        }
        else if (caminhox[i, j] == 1 & caminhoy[i, j] == 1)
            paste0(reconstroiy(i-1, j-1, caminhox, caminhoy, y),y[j])
        else if(caminhox[i, j] == 1 & caminhoy[i, j] == 0)
            paste0(reconstroiy(i-1, j, caminhox, caminhoy, y),"_")
        else if(caminhox[i, j] == 0 & caminhoy[i, j] == 1)
            paste0(reconstroiy(i, j-1, caminhox, caminhoy, y),y[j])

    }



    reconstroi_arestas <- function(arestas, i, j, caminhox, caminhoy)
    {
        if(j==0 | i ==0)
        {
            print(paste0(i,"-",j))
            if (i > 0)
            {
                print(paste0("i>0:",i,"-",j))
                add_aresta(arestas,0:(i-1),0,1:i,0)
            }
            else if (j > 0)
            {
                add_aresta(arestas,0,0:(j-1),0,1:j)
            }
            else
            {
                arestas
            }

        }
        else if (caminhox[i, j] == 1 & caminhoy[i, j] == 1)
        {
            arestas <- add_aresta(arestas,i-1,j-1,i,j)
            reconstroi_arestas(arestas,i-1,j-1, caminhox, caminhoy)
        }
        else if(caminhox[i, j] == 1 & caminhoy[i, j] == 0)
        {
            arestas <- add_aresta(arestas,i-1,j,i,j)
            reconstroi_arestas(arestas,i-1,j, caminhox, caminhoy)
        }
        else if(caminhox[i, j] == 0 & caminhoy[i, j] == 1)
        {
            arestas <- add_aresta(arestas,i,j-1,i,j)
            reconstroi_arestas(arestas,i,j-1, caminhox, caminhoy)
        }

    }




    matches <- reactive({


        arestas <- tibble()

        x <- unlist(strsplit(input$seq1, ""))
        y <- unlist(strsplit(input$seq2, ""))

        n <- length(x)
        m <- length(y)


        tabela <- matrix(0,nrow = length(x)+1, length(y)+1 )

        caminhox <- matrix(0,nrow = length(x), length(y) )
        caminhoy <- matrix(0,nrow = length(x), length(y) )

        arestas <- tibble()

        for(i in 1:n+1)
        {
            tabela[i,1] <- (i-1)*input$penal_gap
        }

        for(j in 1:m+1)
        {
            tabela[1,j] <- (j-1)*input$penal_gap
        }


        for (i in 2:(n+1) )
        {
            for (j in 2:(m+1))
            {
                op1 <- troca(x[i-1],y[j-1], input$penal_igual, input$penal_diferente) + tabela[i-1,j-1]
                op2 <- input$penal_gap + tabela[i-1, j]
                op3 <- input$penal_gap + tabela[i, j-1]
                tabela[i,j] <- min(op1, op2, op3)
                print(op1)
                print(op2)
                print(op3)
                tabela[i,j]
                if (tabela[i,j] == op1)
                {
                    caminhox[i-1,j-1] <- 1
                    caminhoy[i-1,j-1] <- 1
                }
                else if (tabela[i,j] == op2)
                {
                    caminhox[i-1,j-1] <- 1
                    caminhoy[i-1,j-1] <- 0
                }
                else if (tabela[i,j] == op3)
                {
                    caminhox[i-1,j-1] <- 0
                    caminhoy[i-1,j-1] <- 1
                }

            }
        }



        arestas <- reconstroi_arestas(arestas,n,m, caminhox, caminhoy)



        list(reconstroix(n,m, caminhox, caminhoy, x), reconstroiy(n,m, caminhox, caminhoy, y), tabela, arestas)


    })

    
    tamanho_match = reactive({
        
        matchx <- unlist(matches()[1])
        matchx_array <- unlist(strsplit(unlist(matchx), ""))
        length(matchx_array)
        
    })
    

    output$algoritmo <- renderPlot({

        matchx <- unlist(matches()[1])
        matchy <- unlist(matches()[2])
        tabela <- unlist(matches()[3])
        
        

        
        x <- unlist(strsplit(input$seq1, ""))
        y <- unlist(strsplit(input$seq2, ""))

        tabela <- matrix(tabela,nrow = length(x)+1, length(y)+1 )
        
                
        n <- length(x)
        m <- length(y)
        
        print(paste0("matchx:", matchx))
        print(paste0("tipo:", typeof(matchx)))
        print(paste0("matchy:", matchy))
        
        print("1")
        
        matchx_array <- unlist(strsplit(unlist(matchx), ""))
        matchy_array <- unlist(strsplit(unlist(matchy), ""))

        print("2")
        
        matchx_tib <- tibble(seq = "Sequência 1", letra = matchx_array, pos = 1:length(matchx_array) )

        matchy_tib <- tibble(seq = "Sequência 2", letra = matchy_array, pos = 1:length(matchy_array) )

        
        arestas <- matches()[4][[1]] %>% 
            mutate(um = 1) %>% 
            mutate(passo = cumsum(um))
            
        
        print("3")
        
        matches_comp <- rbind(matchx_tib, matchy_tib)

        print("arestas")
        
        print(arestas)
        
        tabela_tib <- as.tibble(as.matrix(tabela)) %>%
            mutate(um = 1) %>%
            mutate(linha = as.integer(cumsum(um))) %>%
            select(-um) %>%
            gather(coluna, OPT, -linha) %>%
            mutate(coluna = as.integer( str_replace(coluna,"V",""))) %>% 
            identity()  



        seq1_tib <- tibble(letra = c("-",x), pos = 1:(n+1))
        seq2_tib <- tibble(letra = c("-",y), pos = 1:(m+1))



        ggplot(tabela_tib, aes(y = linha, x = coluna)) +
            geom_tile(aes(fill = OPT)) +
            geom_text(aes(label = OPT), face = "bold", size=4, vjust = -0.8, hjust = 0.8 ) +
            geom_point()+
            scale_fill_gradient(low = "white", high = "skyblue") +
            labs(x = "Seq 2", y = "Seq 1") +
            geom_segment(data = arestas %>% filter(passo <= input$sliderPassos),
                         aes(
                             x = origem_y,
                             y = origem_x,
                             xend = destino_y,
                             yend = destino_x),
                         arrow = arrow( angle = 10, length = unit(0.5,"cm") , type = "closed",),
                         alpha = 0.2,
                         size = 1,
                         color = "black"

            ) +
            geom_text(data = seq1_tib, aes( label = letra, x = -.5, y = pos ), size = 6 ) +
            geom_text(data = seq2_tib, aes( label = letra, y = -.5, x = pos ), size = 6 ) +
            theme_void() +
            theme(
                panel.grid.minor = element_line(color = "gray90", size = 0.002, linetype = "dotted"),
                panel.grid.major = element_line(color = "gray90", size = 0.002, linetype = "dotted"),
                panel.ontop = TRUE
            ) +
            scale_x_continuous(minor_breaks = seq(1,100,1) ) +
            scale_y_continuous(minor_breaks = seq(1,100,1) )

    })


    output$matches <- renderPlot({


        matchx <- unlist(matches()[1])
        matchy <- unlist(matches()[2])

        matchx_array <- unlist(strsplit(unlist(matchx), ""))
        matchy_array <- unlist(strsplit(unlist(matchy), ""))


        matchx_tib <- tibble(seq = "Sequência 1", letra = matchx_array, pos = 1:length(matchx_array), passo = length(matchx_array):1  )

        matchy_tib <- tibble(seq = "Sequência 2", letra = matchy_array, pos = 1:length(matchy_array), passo = length(matchx_array):1 )


        matches_comp <- rbind(matchx_tib, matchy_tib) %>% 
            filter(passo <= input$sliderPassos)

        ggplot(matches_comp) +
            geom_text(aes(label = letra, x = pos, y = seq), size = 8) +
            theme_void() +
            theme(
                panel.grid.minor.x = element_line(color = "gray90"),
                panel.grid.major.x = element_line(color = "gray90")
            ) +
            scale_x_continuous(minor_breaks = seq(1,100,1) )

    })

    output$passo <- renderUI(
        
        sliderInput("sliderPassos", label = "Passos:", step = 1, value = tamanho_match(), min= 1, max = tamanho_match())
    )

    
}

# Run the application 
shinyApp(ui = ui, server = server)
