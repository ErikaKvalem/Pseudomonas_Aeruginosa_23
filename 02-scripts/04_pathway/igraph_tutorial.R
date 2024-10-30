#install.packages("igraph")

##Load libraries
library(igraph)
library(tidygraph)
library(ggraph)
library(tidyverse)
set_graph_style() # This sets the default style to the graph style

media.edge<-read.csv("data/Dataset1-Media-Example-EDGES.csv")
media.node<-read.csv("data/Dataset1-Media-Example-NODES.csv")

media<-tbl_graph(media.node,media.edge)
media

media_m <- media %>% activate(edges) %>% filter(type=="mention")
media_m %>%
  ggraph(layout = "sugiyama") +
  geom_node_text(aes(label = media, color = type.label), size=3) +     geom_edge_diagonal(color = "gray", alpha = 0.4) +theme_graph()

set.seed(100) 
media_m %>% 
  ggraph(layout = 'graphopt') + 
  geom_edge_link(arrow = arrow(length = unit(2, 'mm')), 
                 end_cap = circle(3, 'mm')) + 
  geom_node_point(size = 3) +  geom_node_text(aes(label = media, color = type.label), size=3,repel = T) + theme_graph()
