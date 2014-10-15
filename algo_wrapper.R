#!/usr/local/bin/Rscript
#### This is an example PBS script for network-based meta-analysis of RNAi screens in R.
#### To define a name for the job (will be displayed in qstat, pbstop output):
#PBS -N Cl1-Random-WNV
### To change output files
#PBS -e Cl1-Random-WNV-error.log
#PBS -o Cl1-Random-WNV-output.txt
### Ressources requested: Memory and Time - for example 12 hours and 1 gb memory
#PBS -l cput=30:00:00,mem=8000mb
### Following are the R-commands to execute.
library(GOSemSim)
library(gdata)
library(igraph)
library(org.Hs.eg.db)
# Load analysis functions from algo_functions.R
source("/home/amberkar/Hades_Work/Work/Data/RNAi_Screens/NewAlgo/algo_functions.R")
# Read clusters obtained from HPIN
cl1_files=list.files("/home/amberkar/Hades_Work/Work/Data/RNAi_Screens/NewAlgo/Clust1_Output/Final_Clusters",pattern="*.csv",full.names=T)
# Read HPIN collated from diff. sources as tsv into a igraph object
iref9_graph <- graph.data.frame(read.table("/home/amberkar/Hades_Work/Work/Data/Protein_Int_Data/Unified_Interactome3/HSA_UnifiedInteractome2012_900.txt",sep="\t",header=T,as.is=T),directed=F)
file_prefix=painted_sig_subgraphs_data1=list()
# Read in hit lists from various RNAi screens
hiv_brass <- scan("/home/amberkar/Hades_Work/Work/Data/RNAi_Screens/Combined-3viruses/HIV-Brass-Hits-Uniprot.txt",sep="\n",what="char")
hiv_zhou <- scan("/home/amberkar/Hades_Work/Work/Data/RNAi_Screens/Combined-3viruses/HIV-Zhou-Hits-Uniprot.txt",sep="\n",what="char")
hiv_koenig <- scan("/home/amberkar/Hades_Work/Work/Data/RNAi_Screens/Combined-3viruses/HIV-Koenig-Hits-Uniprot.txt",sep="\n",what="char")
#hcv_mar <- scan("/home/amberkar/Hades_Work/Work/Data/RNAi_Screens/HCV_Drug_Genome/HCV_DG_Hits_Uniprot.txt",sep="\n",what="char")
hcv_brass <- scan("/home/amberkar/Hades_Work/Work/Data/RNAi_Screens/Combined-3viruses/HCV-Brass-Hits-Uniprot.txt",sep="\n",what="char")
hcv_tai <- scan("/home/amberkar/Hades_Work/Work/Data/RNAi_Screens/Combined-3viruses/HCV-Tai-Hits-Uniprot.txt",sep="\n",what="char")
hcv_lup <- scan("/home/amberkar/Hades_Work/Work/Data/RNAi_Screens/Combined-3viruses/Hits_Lupberger_Uniprot.txt",sep="\n",what="char")
#hcv_mar_val <- scan("../HCV-Marion-Val-Hits-Proteins.txt",sep="\n",what="char")
# Build virus-specific hits in a R lists object
hiv_combined <- c(hiv_brass,hiv_koenig,hiv_zhou)
hcv_combined <- c(hcv_brass,hcv_tai,hcv_lup)
wnv <- scan("/home/amberkar/Hades_Work/Work/Data/RNAi_Screens/WestNileVirus Screen/WNV-Hits-Uniprot.txt",sep="\n",what="char")
# Read Dharmacon siGENOME library for using as 'Background'
bg_ptns = scan("/home/amberkar/Hades_Work/Work/Data/RNAi_Screens/Dharmacon_Human_Genome_Library/G-005000_gene_uniprot.txt",sep="\n",what="char")
# Read centralities for individual proteins
hsa_UI3_cent=read.table("/home/amberkar/Hades_Work/Work/Data/Protein_Int_Data/Unified_Interactome3/HSA_UnifiedInteractome2012_900_Centralities.txt",sep="\t",header=T,as.is=T)
screens=list()
# Input virus hit lists into R list object
screens[[1]]=hcv_mar
screens[[2]]=hcv_brass
screens[[3]]=hcv_tai
screens[[4]]=hcv_lup
#hcv_combined_random=sample(V(iref9_random)$name,size=length(hcv_combined),replace=T)
#wnv_combined_random=sample(V(iref9_random)$name,size=length(wnv),replace=T)

viruses_folders=c("HCV_Marion","HCV_Brass","HCV_Tai","HCV_Lup")
names(screens)=viruses_folders
for(i in 1:16){
  file_prefix[[i]]=strsplit(strsplit(cl1_files[i],"Final_Clusters/")[[1]][2],".csv")
}
clust_size_range=seq(25,100,by=5) # Defining size range of clusters

setwd("/home/amberkar/Hades_Work/Work/Data/RNAi_Screens/NewAlgo/Cl1_Wrapper_Results_HCV_Meta/")
for(v in 1:length(viruses_folders)){ # Running loop over all viruses
  #if(length(clust_size_range)==length(cl1_files)){ # Check for defined clus
  #cat(paste("Cluster size range and associated files match! Now processing clusters...\n"))
  #hits=unique(unlist(screens))
  for(fl in length(cl1_files):1){ # Running analysis over all cluser files
    painted_sig_subgraphs_data1=sig_subgraphs_centralities_data1=painted_non_sig_subgraphs_data1=exact_clust_size=list() # Initialise various lists
    data1=read.table(cl1_files[fl],sep="\t",header=F,as.is=T,fill=T,row.names=1) # Reading a cluster file
    mp_hits=map_hits(data1,unlist(screens[v],use.names=F),bg_ptns) # Mapping hits onto clusters
    ns_data1=no_singleton_clusters(mp_hits) # Removing singleton clusters
    colnames(mp_hits[[1]])=colnames(ns_data1) 
    if(dim(ns_data1)>0){
      sig_clust_data1=find_sig_clusters(ns_data1,unlist(screens[v],use.names=F),iref9_graph,bg_ptns,data1) #Calculating significance of every cluster
      if(dim(sig_clust_data1[[1]])[1]>0){
        sig_subgraphs_data1=build_sig_subgraphs(sig_clust_data1[[1]],iref9_graph,data1) # Building subnetworks from clusters
        non_sig_subgraphs_data1=build_sig_subgraphs(sig_clust_data1[[2]],iref9_graph,data1) # Building subnetworks from non-sig clusters,
        for(k1 in 1:length(sig_subgraphs_data1)){
          cat(paste("Painting sig clusters",k1,"...\n"))
          painted_sig_subgraphs_data1[[k1]]=paint_sig_subgraphs(sig_subgraphs_data1[[k1]],unlist(screens[v],use.names=F),bg_ptns) #Painting hits in subnetwork
        }
        
        for(k2 in 1:length(non_sig_subgraphs_data1)){
          cat(paste("Painting non-sig clusters",k2,"...\n"))
          painted_non_sig_subgraphs_data1[[k2]]=paint_sig_subgraphs(non_sig_subgraphs_data1[[k2]],unlist(screens[v],use.names=F),bg_ptns) #Painting hits in subnetwork
        }
        cat(paste("Calculating centralities for sig clusters...\n"))
        sig_subgraphs_centralities_data1=assign_subgraph_centrality(painted_sig_subgraphs_data1) #Assigning node centralities to every protein in the subnetwork
        cat(paste("Calculating centralities for non-sig clusters...\n"))
        non_sig_subgraphs_centralities_data1=assign_subgraph_centrality(painted_non_sig_subgraphs_data1) #Assigning node centralities to every protein in the subnetwork
        sig_clusts_ranks_data1=calc_clust_ranks(sig_subgraphs_centralities_data1,sig_clust_data1[[1]])
        non_sig_clusts_ranks_data1=calc_clust_ranks(non_sig_subgraphs_centralities_data1,sig_clust_data1[[2]])
        rownames(sig_clust_data1[[1]])=rownames(sig_clusts_ranks_data1)#Table consisting of clusters with mean values of network centralities
        rownames(sig_clust_data1[[2]])=rownames(non_sig_clusts_ranks_data1)#Table consisting of clusters with mean values of network centralities
        names(painted_sig_subgraphs_data1)=rownames(sig_clust_data1[[1]])
        names(painted_non_sig_subgraphs_data1)=rownames(sig_clust_data1[[2]])
        cat(paste("Writing results to files...\n"))
        write.table(mp_hits[[1]],paste(viruses_folders[v],unlist(file_prefix)[fl],"SemSim_MappedHits.txt",sep="_"),sep="\t",col.names=T,row.names=T,quote=F)
        write.table(sig_clust_data1[[1]],paste(viruses_folders[v],unlist(file_prefix)[fl],"SemSim_SigClusters.txt",sep="_"),sep="\t",col.names=T,row.names=T,quote=F)
        write.table(sig_clust_data1[[2]],paste(viruses_folders[v],unlist(file_prefix)[fl],"SemSim_NonSigClusters.txt",sep="_"),sep="\t",col.names=T,row.names=T,quote=F)
        write.table(sig_clusts_ranks_data1,paste(viruses_folders[v],unlist(file_prefix)[fl],"SemSim_SigClustersRank.txt",sep="_"),sep="\t",col.names=T,row.names=T,quote=F)
        write.table(non_sig_clusts_ranks_data1,paste(viruses_folders[v],unlist(file_prefix)[fl],"SemSim_NonSigClustersRank.txt",sep="_"),sep="\t",col.names=T,row.names=T,quote=F)
       
  
        names(sig_subgraphs_centralities_data1)=rownames(sig_clust_data1[[1]])
        names(non_sig_subgraphs_centralities_data1)=rownames(sig_clust_data1[[2]])
        for(pn1 in 1:length(sig_subgraphs_centralities_data1)){
          write.graph(sig_subgraphs_centralities_data1[[pn1]],paste(viruses_folders[v],unlist(file_prefix)[fl],names(sig_subgraphs_centralities_data1)[pn1],"SemSim_PaintedSubgraph.gml",sep="_"),format="gml")  
        }
        
        
        
      }
      
      else  next;
      cat(paste("Skipped for cluster size",clust_size_range[fl]," for virus!..\n"))
    }
      
    }

}

