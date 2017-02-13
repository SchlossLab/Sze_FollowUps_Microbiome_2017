### Common OTUs Network
### Create data to be used in Network
## Marc Sze


###Load needed Libraries and functions
source('code/functions.R')

loadLibs(c("dplyr", "tidyr", "reshape2", "scales", "igraph"))

### Load in needed data tables
common_pred_data <- read.csv("results/tables/common_pred_data.csv", header = T, stringsAsFactors = F)

### Create data table
lesion_pos <- filter(common_pred_data, lesion == "Yes") %>% select(-lesion)
lesion_neg <- filter(common_pred_data, lesion == "No") %>% select(-lesion)

comparisons <- colnames(select(common_pred_data, -lesion))

key <- c(1:length(comparisons))
names(key) <- comparisons

### Create co-occurance percentages
neg_occurance <- as.data.frame(get_co_occurance(lesion_neg, comparisons))
pos_occurance <- as.data.frame(get_co_occurance(lesion_pos, comparisons))

### Create list of connections above a given threshold
neg_co_list <- get_list_coocur(neg_occurance, cutoff = 0.95)
pos_co_list <- get_list_coocur(pos_occurance, cutoff = 0.95)

### Create a network (igraph) object

total_nodes <- 24
test <- make_empty_graph(n=total_nodes, directed = F)

listNames <- as.numeric(names(neg_co_list))

for(i in 1:length(listNames)){
  
  for(j in 1:length(neg_co_list[[i]]))
    
    test <- add.edges(test, c(listNames[i],neg_co_list[[i]][j]))
}

set.seed(3)
plot(test, vertex.label=names(key))

set.seed(3)
plot(test)



test <- graph_from_literal()
make_empty_graph(n = 5, directed = F)
add_edges(g, c(1,2, 2,3, 3,4, 4,5))



