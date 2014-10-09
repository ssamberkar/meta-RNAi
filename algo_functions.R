###########################################################
#  Functions for generic meta-analysis of 7 genome-wide   #
#  RNAi screens                                           #
#  Author: Sandeep Amberkar                               #
##########################################################
require(GOSemSim)
require(gdata)
require(igraph)
require(org.Hs.eg.db)
#Map screen hits onto clusters
map_hits <- function(data1,hits,background)
{
  matched_hits <- matrix(nrow=nrow(data1),ncol=length(hits))
  background_match <- list()
  ptn_binom <- list()
  hit_stat <- matrix(NA,nrow=dim(data1)[1],ncol=3)
  rownames(hit_stat) <- row.names(data1)
  cat(paste("Initiating Cluster Mapping...\n"))
  #Sys.sleep(3)
  for (i in 1:(dim(data1)[1]))
  {
    match <-  which(hits%in%data1[i,])
    bg_match <- which(background%in%data1[i,])
    background_match[[i]] <- background[bg_match]
    le2 <- length(background_match[[i]])
    cat(paste("Processing Cluster No...", i,"\n"))
    if (length(match) != 0)
    {
      matched_hits[i,1:length(hits[match])] <- hits[match]
      tmp   <-  apply(as.matrix(unique(unlist(data1[i,]))),1,nchar)
      le <- 0
      if(any(tmp==0)){
        id <- which(tmp==0)
        le <- length(tmp[-id])
      }
      else le <- length(tmp)
      hit_stat[i,] <- c(length(match),le,le2)
      
    }
    #else matched_hits[[i]] <- NA
    
  }
  cat(paste("Mapped",length(hits),"hits to",dim(data1)[1],"clusters.\n",sep=" "))
  
  hit_stat[1,2]=length(data1[1,])
  rownames(hit_stat) <- row.names(data1)
  rownames(matched_hits) <- row.names(data1)
  hs_list <- list()
  hs_list[[1]] <- hit_stat
  hs_list[[2]] <- matched_hits
  return(hs_list);
}
#Function to separate singleton clusters
no_singleton_clusters <- function(hit_stat,data1)
{
  #
  cat(paste("\nRemoving Clusters with no hits...\n"))
  #tmp <- which(hit_stat[[1]][,1]!="NA")
  t1 <- data.frame(na.omit(hit_stat[[1]]))
  tmp1 <- which(t1[,1]==1|t1[,1]==2)
  #sing <- t1[tmp1,]
  t2 <- t1[-tmp1,]
  #   if(length((which(t2[,1]>=t2[,3])))==integer(0)){
  #     t2 <- t2[-which(t1[,1]>=t1[,3]),]
  #     cat(paste("\nSeperating Clusters with single hit...\n"))
  #     colnames(t2) <- c("Number of matched hits","Number of proteins in Cluster","Number of screened proteins")
  #   }
  #   
  #tmp2 <- which(t1[,1]>=t1[,3])
  cat(paste("\nSeperating Clusters with single hit...\n"))
  #colnames(sing) <- c("Number of matched hits","Number of proteins in Cluster","Number of screened proteins")
  colnames(t2) <- c("Number of matched hits","Number of proteins in Cluster","Number of screened proteins")
  cat(paste("Removed",length(tmp1),"singleton clusters.\n",sep=" "))
  return(t2);
}


# Determine hit enrichment in clusters
find_sig_clusters <- function(no_singletons,hits,bg_network,lib_ptns,data1){
  x1.vec=k1.vec=pvals=fil_sig_clust=list()
  clust_index=rownames(no_singletons)
  for(i in 1:length(clust_index)){
    cat(paste("Calculating hit overlap from cluster...",clust_index[i],"\n",sep=""))
    x1.vec[[i]]=length(intersect(unlist(hits),unlist(data1[clust_index[i],which(data1[clust_index[i],]!="")])))
  }
  for(j in 1:length(clust_index)){
    cat(paste("Calculating length of cluster...",clust_index[j],"\n",sep=""))
    #length(which(unlist(data1[clust_index[1],which(data1[clust_index[1],]!="")],use.names=F)%in%bg_ptns)
    k1.vec[[j]]=length(which(unlist(data1[clust_index[j],which(data1[clust_index[j],]!="")],use.names=F)%in%lib_ptns))  
  }
  #Calculating significance of clusters
  for(k in 1:length(clust_index)){
    cat(paste("Calculating significance of cluster...",clust_index[k],"\n",sep=""))
    pvals[[k]]=1-phyper(x1.vec[[k]]-1,length(unlist(hits,use.names=F)),length(lib_ptns)-length(unlist(hits,use.names=F)),k1.vec[[k]])
  }
  sig_clust=cbind(no_singletons,unlist(pvals))
  colnames(sig_clust)=c(colnames(no_singletons),"p-value")
  
  fil_sig_clust[[1]]=sig_clust[which(sig_clust[,4]<=0.05),]
  fil_sig_clust[[2]]=sig_clust[-which(sig_clust[,4]<=0.05),]
  return(fil_sig_clust);
  cat(paste("Done!!!",sep="\n"))
}

#Build protein subnetworks from clusters
build_sig_subgraphs=function(sig_clusts,bg_network,clusters){
  sig_nets_index=rownames(sig_clusts)
  seeds=subgraphs=list()
  for(i in 1:length(sig_nets_index)){
    cat(paste("Isolating seeds from ",sig_nets_index[i],"...\n",sep=""))
    seeds[[i]]=unlist(clusters[sig_nets_index[i],which(clusters[sig_nets_index[i],]!="")],use.names=F)
  }
  for(j in 1:length(seeds)){
    cat(paste("Creating subgraphs from ",seeds[j],"...\n",sep=""))
    subgraphs[[j]]=induced.subgraph(bg_network,na.omit(seeds[[j]]))
  }
  return(subgraphs)
  
}

#Paint screen hits within the subnetwork
paint_sig_subgraphs=function(subnet,screen_hits,lib_genes){
  cat(paste("Painting default colors..\n"))
  V(subnet)$color="grey" # Default node color
  cat(paste("Painting screen proteins..\n"))
  V(subnet)$color[which(V(subnet)$name%in%lib_genes)]="white" # Painting lib proteins as white
  cat(paste("Painting hit proteins..\n"))
  V(subnet)$color[which(V(subnet)$name%in%unlist(screen_hits[which(names(screen_hits)=="HCV")]))]=rep(paste(255,192,203,sep=","),length(which(V(subnet)$name%in%unlist(screen_hits[which(names(screen_hits)=="HCV")]))))#Painting HCV hits  
  V(subnet)$color[which(V(subnet)$name%in%unlist(screen_hits[which(names(screen_hits)=="HIV")]))]=rep(paste(255,165,0,sep=","),length(which(V(subnet)$name%in%unlist(screen_hits[which(names(screen_hits)=="HIV")]))))#Painting HIV hits
  V(subnet)$color[which(V(subnet)$name%in%unlist(screen_hits[which(names(screen_hits)=="WNV")]))]=rep(paste(255,255,0,sep=","),length(which(V(subnet)$name%in%unlist(screen_hits[which(names(screen_hits)=="WNV")]))))#Painting WNV hits
  
  cat(paste("Done!!\n"))
  return(subnet)
}

#Compute network centralities for subnetworks and rank them
calc_clust_ranks=function(sig_subgraphs_centralities,fil_sig_clust){
  uniprot_entrez_map=toTable(org.Hs.egUNIPROT)
  if(length(sig_subgraphs_centralities)==1){
    mean_deg_sig_subgraphs=mean(degree(sig_subgraphs_centralities[[1]]))
    cat(paste("Calculating Average path length of significant subgraph of length",length(sig_subgraphs_centralities),"...\n"))  
    mean_pathlength_sig_subgraphs=average.path.length(sig_subgraphs_centralities[[1]])
    cat(paste("Calculating Average clustering coefficient of significant subgraph of length",length(sig_subgraphs_centralities),"...\n"))
    mean_cc_sig_subgraphs=mean(transitivity(sig_subgraphs_centralities[[1]]))
    cat(paste("Calculating Average Betweenness of significant subgraph of length",length(sig_subgraphs_centralities),"...\n"))
    mean_bwn_sig_subgraphs=mean(betweenness(sig_subgraphs_centralities[[1]]))
    cat(paste("Calculating Average Closeness of significant subgraph of length",length(sig_subgraphs_centralities),"...\n"))
    mean_cls_sig_subgraphs=mean(closeness(sig_subgraphs_centralities[[1]]))
    cat(paste("Calculating Average Page.Rank of significant subgraph of length",length(sig_subgraphs_centralities),"...\n"))
    mean_pgr_sig_subgraphs=mean(page.rank(sig_subgraphs_centralities[[1]])$vector)
    clust_genes_entrez2=uniprot_entrez_map$gene_id[which(uniprot_entrez_map$uniprot_id%in%V(sig_subgraphs_centralities_data1[[1]])$name)]
    cat(paste("Calculating Average Dice.Similarity of significant subgraph of length",length(sig_subgraphs_centralities),"...\n"))
    mean_dice.sim_sig_subgraphs=mean(similarity.dice(graph=sig_subgraphs_centralities[[1]],vids=V(sig_subgraphs_centralities[[1]]),mode="all"))
    cat(paste("Calculating Average Wang_GO.BP Similarity of significant subgraph of length",length(sig_subgraphs_centralities),"...\n"))
    mean_wang.index_GO.BP=mean(mgeneSim(clust_genes_entrez2,ont="BP",organism="human",measure="Wang",combine="avg"))
    cat(paste("Calculating Average Wang_GO.CC Similarity of significant subgraph of length",length(sig_subgraphs_centralities),"...\n"))
    mean_wang.index_GO.CC=mean(mgeneSim(clust_genes_entrez2,ont="CC",organism="human",measure="Wang",combine="avg"))
    cat(paste("Calculating Average Wang_GO.MF Similarity of significant subgraph of length",length(sig_subgraphs_centralities),"...\n"))
    mean_wang.index_GO.MF=mean(mgeneSim(clust_genes_entrez2,ont="MF",organism="human",measure="Wang",combine="avg"))
    ranked_fil_sig_clust=cbind(fil_sig_clust,mean_bwn_sig_subgraphs,
                               mean_cls_sig_subgraphs,
                               mean_cc_sig_subgraphs,
                               mean_pgr_sig_subgraphs,
                               mean_deg_sig_subgraphs,
                               mean_pathlength_sig_subgraphs,
                               mean_dice.sim_sig_subgraphs,
                               mean_wang.index_GO.BP,
                               mean_wang.index_GO.CC,
                               mean_wang.index_GO.MF)
    colnames(ranked_fil_sig_clust)=c(colnames(fil_sig_clust),"Mean_Betweenness",
                                     "Mean_Closeness",
                                     "Mean_Clust_Coefficient",
                                     "Mean_PGR","Mean_Degree",
                                     "Mean_PathLength",
                                     "Mean_DiceSim",
                                     "Mean_WangSim_GO.BP",
                                     "Mean_WangSim_GO.CC",
                                     "Mean_WangSim_GO.MF")
    
    return(ranked_fil_sig_clust)
  }
  else
    cat(paste("Calculating Average Degree of significant subgraph...\n"))
    mean_deg_sig_subgraphs=unlist(lapply(lapply(sig_subgraphs_centralities,degree),mean))
    cat(paste("Calculating Average path length of significant subgraph...\n"))  
    mean_pathlength_sig_subgraphs=unlist(lapply(lapply(sig_subgraphs_centralities,average.path.length),mean))
    cat(paste("Calculating Average Clustering-coefficient of significant subgraph...\n"))
    mean_cc_sig_subgraphs=unlist(lapply(lapply(sig_subgraphs_centralities,transitivity),mean))
    cat(paste("Calculating Average Betweenness of significant subgraph...\n"))
    mean_bwn_sig_subgraphs=unlist(lapply(lapply(sig_subgraphs_centralities,betweenness),mean))
    cat(paste("Calculating Average Closeness of significant subgraph...\n"))
    mean_cls_sig_subgraphs=unlist(lapply(lapply(sig_subgraphs_centralities,closeness),mean))
    cat(paste("Calculating Average PageRank of significant subgraph...\n"))
    mean_pgr_sig_subgraphs=unlist(lapply(sapply(sig_subgraphs_centralities,page.rank)[,1:length(sig_subgraphs_centralities)][1,],mean))
  
  clust_genes_entrez=list()
  for(cle in 1:length(sig_subgraphs_centralities)){
    clust_genes_entrez[[cle]]=uniprot_entrez_map$gene_id[which(uniprot_entrez_map$uniprot_id%in%V(sig_subgraphs_centralities[[cle]])$name)]
  }
  cat(paste("Calculating Dice semantic similarity of genes in each cluster...\n"))
  mean_dice.sim_sig_subgraphs=unlist(lapply(lapply(sig_subgraphs_centralities,similarity.dice),mean))
  cat(paste("Calculating semantic similarity of genes for GO.BP ontologies...\n"))
  mean_wang.index_GO.BP=unlist(lapply(mapply(clust_genes_entrez,FUN=mgeneSim,MoreArgs=list(ont="BP",measure="Wang",organism="human",combine="avg")),mean))
  cat(paste("Calculating semantic similarity of genes for GO.CC ontologies in each cluster...\n"))
  mean_wang.index_GO.CC=unlist(lapply(mapply(clust_genes_entrez,FUN=mgeneSim,MoreArgs=list(ont="CC",measure="Wang",organism="human",combine="avg")),mean))
  cat(paste("Calculating semantic similarity of genes for GO.MF ontologies in each cluster...\n"))
  mean_wang.index_GO.MF=unlist(lapply(mapply(clust_genes_entrez,FUN=mgeneSim,MoreArgs=list(ont="MF",measure="Wang",organism="human",combine="avg")),mean))
  if(is.list(mean_wang.index_GO.BP)==FALSE){
    mean_wang.index_GO.BP=mean(mean_wang.index_GO.BP)
  }
  if(is.list(mean_wang.index_GO.CC)==FALSE){
    mean_wang.index_GO.CC=mean(mean_wang.index_GO.CC)
  }
  if(is.list(mean_wang.index_GO.MF)==FALSE){
    mean_wang.index_GO.MF=mean(mean_wang.index_GO.MF)
  }
  #----Please change after modified code run----#ranked_fil_sig_clust=cbind(fil_sig_clust,mean_bwn_sig_subgraphs,mean_cls_sig_subgraphs,mean_cc_sig_subgraphs,mean_pgr_sig_subgraphs,mean_sim.dice_sig_subgraphs,mean_wang.index_GO.CC,mean_wang.index_GO.MF)
  #----Please change after modified code run----colnames(ranked_fil_sig_clust)=c(colnames(fil_sig_clust),"Mean_Betweenness","Mean_Closeness","Mean_Clust_Coefficient","Mean_PGR","Mean_DiceSim","Mean_WangSim_GO.CC","Mean_WangSim_GO.MF")
  ranked_fil_sig_clust=cbind(fil_sig_clust,
                             unlist(mean_bwn_sig_subgraphs),
                             unlist(mean_cls_sig_subgraphs),
                             unlist(mean_cc_sig_subgraphs),
                             unlist(mean_pgr_sig_subgraphs),
                             mean_deg_sig_subgraphs,
                             mean_pathlength_sig_subgraphs,
                             mean_dice.sim_sig_subgraphs,
                             mean_wang.index_GO.BP,
                             mean_wang.index_GO.CC,
                             mean_wang.index_GO.MF)
  colnames(ranked_fil_sig_clust)=c(colnames(fil_sig_clust),
                                   "Mean_Betweenness",
                                   "Mean_Closeness",
                                   "Mean_Clust_Coefficient",
                                   "Mean_PGR",
                                   "Mean_Degree",
                                   "Mean_PathLength",
                                   "Mean_DiceSim",
                                   "Mean_WangSim_GO.BP",
                                   "Mean_WangSim_GO.CC",
                                   "Mean_WangSim_GO.MF")
  
  return(ranked_fil_sig_clust)
}
#Compute Jaccard coefficients between protein complex and subnetwork
jaccard=function(complex,clust){
  #Submit both complex and clusters ONLY in list format!
  jcd=matrix(NA,nrow=length(complex),ncol=length(clust))
  for(a in 1:length(complex)){
    for(b in 1:length(clust)){
      cat("Calculating Jaccard Score for Complex",a,"& Cluster",b,"\n")
      jcd[a,b]=length(intersect(complex[[a]],clust[[b]]))/length(union(complex[[a]],clust[[b]]))
    }
  }
  return(jcd)
}

#Assign node-wise centrality values to subnetworks in igraph format
assign_subgraph_centrality <- function(sig_subgraphs){
  for(sg in 1:length(sig_subgraphs)){
    #subgraph_index=which(hsa_UI3_cent$ProteinName%in%V(sig_subgraphs[[sg]])$name)
    cat(paste("Assigning Degree centrality of subnetwork",sg,"..\n"))
    V(sig_subgraphs[[sg]])$deg=degree(sig_subgraphs[[sg]],v=V(sig_subgraphs[[sg]])$name,mode="all",loops=F)
    cat(paste("Assigning Betweenness centrality of subnetwork",sg,"..\n"))
    V(sig_subgraphs[[sg]])$bwn=betweenness(sig_subgraphs[[sg]],v=V(sig_subgraphs[[sg]])$name,directed=F)
    cat(paste("Assigning Closeness centrality of",sg,"\n"))
    V(sig_subgraphs[[sg]])$cls=closeness(sig_subgraphs[[sg]],v=V(sig_subgraphs[[sg]])$name)
    cat(paste("Assigning PageRank centrality of",sg,"\n"))
    V(sig_subgraphs[[sg]])$pgr=page.rank(sig_subgraphs[[sg]],v=V(sig_subgraphs[[sg]])$name,directed=F)$vector
    
  }
  return(sig_subgraphs)
}
  





