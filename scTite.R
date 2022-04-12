{
  library(data.table)
  library(expm)
  library(ggm)
  library(lpSolve)
  library(princurve)
  library(mclust)
  library(umap)
  library(igraph)
  library(SCENT)
}


setwd("path/scTite/example")

##example
#hsmm_0.15_MAGIC_computer_SR,cell*gene,processed dataset, gene names transformed according to SR
#hsmm_information, contain cell Classes and labels
#hsmm_SR,SR value of cells
#k, number of clusters
#transfer_paramter,Transition Cell Threshold
#startCluster,start cluster
#isNormalized , Sort after normalized pseudo-time
#Improve_efficiency, Whether to reduce the time cost, maybe reduce the accuracy


data<-"hsmm_0.15_MAGIC_computer_SR.txt"
data_filter_MAGIC<-read.table(data)
SR_name<-"hsmm_SR.txt"
SR_entropy<-read.table(SR_name)
k<-4
transfer_paramter<-0.2
startCluster<-1
isNormalized<-TRUE
Improve_efficiency<-FALSE

#Trajectory, contains trajectories and pseudo-time order
Trajectory<-sctite(data_filter_MAGIC,k,SR_entropy,transfer_paramter,startCluster,isNormalized)


sctite<-function(data_filter_MAGIC,k,SR_entropy,transfer_paramter,startCluster,isNormalized){
  
  
  #Clustering
  {
    #seed1 <- sample(1:10000, 1)
    seed1<-123#hsmm
    set.seed(seed1)
    #Rows are cells, columns are genes
    X<-umap(data_filter_MAGIC,config=umap.defaults,method=c("naive"))#X<-umap(X,min_dist=0.2)
    X<-as.matrix(X[["layout"]])
    use_reduc<-"Umap_mclust"
    red_mclust<-Mclust(X, G = k)#,modelNames = "VVV")
    mclust_label<-red_mclust[["classification"]]
    mclust_center<-t(red_mclust[["parameters"]][["mean"]])
    flag_stop<-0
    for(i in seq_len(k)){
      number_cluster_i<-length(which(mclust_label==i))
      if(number_cluster_i<10){
        flag_stop<-1
      }
    }
    if(length(unique(mclust_label))!=k){
      flag_stop<-1
    }
  }
  if(flag_stop==1){
    print("The clustering effect may be poor, please select an appropriate number of clusters or re-run the experiment!")
  }else{
    print("Clustering step completed!")
  }
  #Calculate transition entropy
  {
    probs<-t(red_mclust[["z"]])
    for(i in seq_len(nrow(probs))){
      for(j in seq_len(ncol(probs)))
        if(probs[i,j]!=0)
          probs[i,j]<-probs[i,j]*log(probs[i,j])/log(2)
    }
    cluster_entropy <-cbind(as.matrix( -apply(probs,2,sum)),c(1:nrow(data_filter_MAGIC)))
    temp_1<-cluster_entropy[order(cluster_entropy[,1],decreasing = TRUE),]
    #After obtaining the transition entropy, select the first 5%-20% of cells as transition cells
    Threshold<- round(transfer_paramter*nrow(data_filter_MAGIC))
    transfer_cell<-as.matrix(temp_1[1:Threshold,2])
    transfer_cell_coordinate<-cbind(X[transfer_cell,],transfer_cell,mclust_label[as.numeric(transfer_cell)])
    colnames(transfer_cell_coordinate)<-c("x","y","transfer_cell","label")
  }
  
  #Calculate the potency state and get the wasserstein distance between clusters
  {
    #Infer potency state
    pot.o <- InferPotencyStates(potest.v=SR_entropy,type=c("SR"), pheno.v =mclust_label)
    k_sr<-max(pot.o[["potS"]])
    X_label_SR<-cbind(X,mclust_label,pot.o[["potS"]])
    #Calculate the potency state distribution in each cluster
    dis_entropy<-matrix(0,nrow=k,ncol=k_sr)
    for(i in seq_len(k)){
      x_temp<-X_label_SR[X_label_SR[,3]==i,]
      for(j in seq_len(nrow(x_temp)))
        dis_entropy[i,x_temp[j,4]]<-dis_entropy[i,x_temp[j,4]]+1
      sum_i<-sum(dis_entropy[i,])
      dis_entropy[i,]<-round(dis_entropy[i,]/sum_i,digit=6)
    }
    Expand_parameters<-1e+6
    dis_entropy<-dis_entropy*Expand_parameters
    X_label<-rbind(mclust_center,X)
    rownames(X_label)<-c(1:nrow(X_label))
    #Calculate the wasserstein distance and build the inter-cluster distance matrix
    {
      cluster_sr_coordinate<-matrix(0,ncol=2)
      for(i in seq_len(k)){
        for(j in seq_len(k_sr)){
          x_temp<-matrix(0,ncol=4)
          x_temp<-as.matrix(rbind(x_temp,X_label_SR[X_label_SR[,3]==i&X_label_SR[,4]==j,]))
          if(nrow(x_temp)==1){
            cluster_sr_coordinate<-rbind(cluster_sr_coordinate,c(0,0))
          }else{
            nrow_temp<-nrow(x_temp)-1
            x_coordinate<-sum(x_temp[,1])/nrow_temp
            y_coordinate<-sum(x_temp[,2])/nrow_temp
            cluster_sr_coordinate<-rbind(cluster_sr_coordinate,c(x_coordinate,y_coordinate))
          }
        }
      }
      cluster_sr_coordinate<-cluster_sr_coordinate[-1,]
      entropy_dist<-matrix(0,nrow=k,ncol=k)
      for(i in seq_len(k)) {
        start<-(i-1)*k_sr+1
        cluster_i<-cluster_sr_coordinate[start:(start+k_sr-1),]
        cluster_i<-matrix(cluster_i[which(rowSums(cluster_i==0)==0),],ncol=2)
        for(j in seq_len(k)){
          if(i==j){
            next
          }
          start<-(j-1)*k_sr+1
          cluster_j<-cluster_sr_coordinate[start:(start+k_sr-1),]
          cluster_j<-matrix(cluster_j[which(rowSums(cluster_j==0)==0),],ncol=2)
          euclid_dist <- function(cluster_i, cluster_j) {
            return(sqrt(sum((cluster_i - cluster_j)^2)))
          }
          w1<-dis_entropy[i,]
          w1<-w1[which(w1!=0)]
          w2<-dis_entropy[j,]
          w2<-w2[which(w2!=0)]
          n1 = length(cluster_i[,1])
          n2 = length(cluster_j[,1])
          
          dist_entropy_temp = matrix(0, n1, n2)
          for (i1 in 1:n1) {
            for (j1 in 1:n2) 
              dist_entropy_temp[i1, j1] = euclid_dist(cluster_i[i1,], cluster_j[j1,])
          }
          entropy_dist[i,j]<-emd(dist_entropy_temp, w1, w2)
        }
      }
      for(i in seq_len(k)){
        for(j in seq_len(k)){
          if(entropy_dist[i,j]=="NaN"){
            entropy_dist[i,j]<-entropy_dist[j,i]
          }
        }
      }
    }
    distance_finary<-entropy_dist
    print(paste0("Matrix symmetry:",isSymmetric(distance_finary)))
    print("Get wasserstein between clusters!")
  }
  
  #lineages_temp save temporary lineage structure
  {
    # Minimum spanning tree constructed for identifing lineages
    colnames(distance_finary)<-c(1:k)
    rownames(distance_finary)<-c(1:k)
    mstree <- ape::mst(distance_finary)
    forest <- mstree
    # identify sub-trees 
    subtrees <- subtrees.update <- forest
    diag(subtrees) <- 1
    while(sum(subtrees.update) > 0){
      subtrees.new <- apply(subtrees,2,function(col){
        rowSums(subtrees[,as.logical(col), drop=FALSE]) > 0
      })
      subtrees.update <- subtrees.new - subtrees
      subtrees <- subtrees.new
    }
    subtrees <- unique(subtrees)
    trees <- lapply(seq_len(nrow(subtrees)),function(ri){
      colnames(forest)[subtrees[ri,]]
    })
    trees <- trees[order(vapply(trees,length,0),decreasing = TRUE)]
    ntree <- length(trees)
    
    lineages_temp<-list()
    start.clus <- as.character(startCluster)
    for(tree in trees){
      if(length(tree) == 1){
        lineages_temp[[length(lineages_temp)+1]] <- tree
        next
      }
      tree.ind <- rownames(forest) %in% tree
      tree.graph <- forest[tree.ind, tree.ind, drop = FALSE]
      degree <- rowSums(tree.graph)
      g <- graph.adjacency(tree.graph, mode="undirected")
      # if you have starting cluster(s) in this tree, draw lineages
      # to each leaf
      if(sum(start.clus %in% tree) > 0){
        starts <- start.clus[start.clus %in% tree]
        ends <- rownames(tree.graph)[
          degree == 1 & ! rownames(tree.graph) %in% starts]
        for(st in starts){
          paths <- shortest_paths(g, from = st, to = ends, 
                                  mode = 'out', 
                                  output = 'vpath')$vpath
          for(p in paths){
            lineages_temp[[length(lineages_temp)+1]] <- as.numeric(names(p))
          }
        }
      }
    }
  }
  
  #transfer_cell_coordinate_temp find the projected edge of the cell
  {
    lineages_all<-lineages_temp
    transfer_cell_coordinate_temp<-transfer_cell_coordinate
    nrow_transer<-nrow(transfer_cell_coordinate)
    lambad_side<-matrix(0,nrow=nrow_transer)
    edge_all<-matrix(0,ncol=2)
    for(i in seq_len(length(lineages_all))){
      for(j in seq_len(length(lineages_all[[i]])-1)){
        edge_all<-rbind(edge_all,c(as.numeric(lineages_all[[i]][[j]]),as.numeric(lineages_all[[i]][[j+1]])))
      }
    }
    edge_all<-edge_all[-1,]
    edge_all<-edge_all[!duplicated( edge_all),]
    for(ii in seq_len(nrow(edge_all))){
      vertix<-matrix(c(1:2),ncol=2)
      for(jj in seq_len(ncol(edge_all))){
        h<-edge_all[ii,jj]
        vertix<-rbind(vertix,X_label[h,])
      }
      vertix<-vertix[-1,]
      cluster_number<-c(edge_all[ii,])
      projection<-matrix(transfer_cell_coordinate[transfer_cell_coordinate[,4]==cluster_number[1]|transfer_cell_coordinate[,4]==cluster_number[2],],ncol=4)
      lambad_temp<-matrix(Inf,nrow=nrow_transer,ncol=1)
      if(nrow(projection)!=0){
        s <-as.matrix(vertix)
        s1<-matrix(projection[,1:2],ncol=2)
        fit1 <- project_to_curve(s1, s)
        b<-c(vertix[2,1]-vertix[1,1],vertix[2,2]-vertix[1,2])
        dis_temp<-sqrt(fit1[["dist_ind"]])
        jj<-1
        for(iii in seq_len(nrow_transer)){
          a<-c(transfer_cell_coordinate[iii,1]-vertix[1,1],
               transfer_cell_coordinate[iii,2]-vertix[1,2])
          label_temp<-transfer_cell_coordinate[iii,3]
          if(label_temp %in% projection[,3]&sum(a * b )/(sqrt(sum(a * a))* sqrt(sum(b * b)))>=0){#可以
            lambad_temp[iii,1]<-dis_temp[jj]
            jj<-jj+1
          }
        }
      }
      lambad_side<-cbind(lambad_side,lambad_temp)
    }
    lambad_side<-as.matrix(lambad_side[,-1])
    lambad_side_2<-matrix(0,nrow=nrow_transer)
    for(i5 in seq_len(nrow_transer)){
      lambad_side_2[i5,1]<- which.min(lambad_side[i5,] )
    }
    transfer_cell_coordinate_temp<-cbind(transfer_cell_coordinate_temp,lambad_side_2)
    colnames(transfer_cell_coordinate_temp)<-c("x","y","transfer_cell","label","projected_edge")
  }
  
  #transition_path save all transition paths，transition_path_final and lineages is the resulting lineage structure
  {
    transition_path<-list()
    list_length<-1
    data_t<-t(data_filter_MAGIC)
    for(i2 in seq_len(nrow(edge_all))){
      #print(paste("edge_now:",i2))
      temp<-matrix(1:3,ncol=3)
      transfer_cell_edge_i<-matrix(transfer_cell_coordinate_temp[transfer_cell_coordinate_temp[,5]==i2,],ncol=5)
      if(nrow(transfer_cell_edge_i)==0){
        transition_path[list_length]<-list(edge_all[i2,])
        list_length<-list_length+1
      }else if(nrow(transfer_cell_edge_i)==1){
        transition_path[list_length]<-list(c(edge_all[i2,1],transfer_cell_edge_i[1,3]+k,edge_all[i2,2]))
        list_length<-list_length+1
      }else{
        #Sort by SR
        transfer_cell_edge_i<-cbind(transfer_cell_edge_i,SR_entropy[transfer_cell_edge_i[,3],])
        colnames(transfer_cell_edge_i)<-c("x","y","number","label","projected edge","SR")
        transfer_cell_edge_i<-transfer_cell_edge_i[order(transfer_cell_edge_i[,6],decreasing = TRUE),]
        n_transfer_cell<-nrow(transfer_cell_edge_i)
        pcor_transfer_cell<-matrix(0,nrow=n_transfer_cell,ncol=n_transfer_cell)
        if(Improve_efficiency==TRUE){
          computer_cell<-50
          #Reduce time cost, only calculate partial correlation coefficients of two vectors and four indirect vectors
          for(i3 in seq(2,n_transfer_cell)){
            if((i3-1)<computer_cell){
              i_3_start<-1
            }else{
              i_3_start<-i3-computer_cell
            }
            
            for(i4 in seq(i_3_start,i3-1)){
              if((i3-1-i4)< 4){
                if((i3-1)==i4){
                  cov_save<-cov(data_t[,transfer_cell_edge_i[c(i4,i3),3]])
                  pcor_transfer_cell[i3,i4]<-pcor(c(1,2),cov_save)
                }else{
                  c1<-c(i4:(i3-1))
                  cov_save<-cov(data_t[,transfer_cell_edge_i[c(c1,i3),3]])
                  c2<-c(1:(ncol(cov_save)-1))
                  pcor_transfer_cell[i3,i4]<-pcor(c(1,ncol(cov_save),c2[-1]),cov_save)
                }
              }else{
                c1<-c(i4:(i4+4))
                cov_save<-cov(data_t[,transfer_cell_edge_i[c(c1,i3),3]])
                c2<-c(1:(ncol(cov_save)-1))
                pcor_transfer_cell[i3,i4]<-pcor(c(1,ncol(cov_save),c2[-1]),cov_save)
                
              }
              
            }
          }
        }else{
          cov_save<-cov(data_t[,transfer_cell_edge_i[,3]])
          for(i3 in seq_len(n_transfer_cell)){
            c1<-c(1:(i3-1))
            for(i4 in seq_len(i3-1)){
              if(i4==0||i4==i3){
                next
              }
              pcor_transfer_cell[i3,i4]<-pcor(c(i4,i3,c1[-i4]),cov_save)
            }
          }
        }
        for(i_transfer in seq_len(n_transfer_cell)){
          for(j_transfer in seq_len(n_transfer_cell))
            if(pcor_transfer_cell[j_transfer,i_transfer]>0)
              temp<-rbind(temp,
                          c(transfer_cell_edge_i[i_transfer,3],
                            transfer_cell_edge_i[j_transfer,3],
                            1/pcor_transfer_cell[j_transfer,i_transfer]))
        }
        temp = as.data.frame(temp)
        temp<-temp[-1,]
        names(temp) = c("start_id","end_id","newcost")
        g3<-graph.data.frame(temp, directed=FALSE)
        temp_path<-get.shortest.paths(g3, from=V(g3)[1], to=V(g3)[n_transfer_cell],weights=temp[,3])
        for(h in seq_len(length(temp_path[["vpath"]]))){
          path_save<-c()
          if(length(temp_path[["vpath"]][[h]])>1){
            for(h1 in seq_len(length(temp_path[["vpath"]][[h]]))){
              number<-as.integer( temp_path[["vpath"]][[h]][h1])
              path_save<-c(path_save,transfer_cell_edge_i[number,3])
            }
          }
          path_save<-c(edge_all[i2,1],path_save,edge_all[i2,2])
          
          for(temp_i in 2:(length(path_save)-1))
            path_save[temp_i]<-path_save[temp_i]+k
          
          transition_path[list_length]<-list(path_save)
          list_length<-list_length+1
        }
      }
    }
    
    
    
    
    
    
    transition_path_final<-list()
    final_length<-1
    for(i in seq_len(length(lineages_temp))){
      fina_temp<-c()
      temp_vertex<-lineages_temp[i]
      temp1<-c()
      for(j in seq_len(length(temp_vertex[[1]])))
        temp1<-c(temp1,as.numeric(temp_vertex[[1]][[j]]))
      #提出路径
      h<-1
      while(h<length(temp1)){
        for(j in seq_len(length(transition_path))){
          t<-length(transition_path[[j]])
          #找到路径
          if(as.numeric(transition_path[[j]][[1]])==temp1[h]&as.numeric(transition_path[[j]][[t]])==temp1[h+1]){
            for(h1 in seq_len(t))
              fina_temp<-c(fina_temp,as.numeric(transition_path[[j]][[h1]]))
          }
        }
        h<-h+1
      }
      fina_temp<-fina_temp[!duplicated(fina_temp)]
      transition_path_final[final_length]<-list(fina_temp)
      final_length<-final_length+1
    }
    
    #print(transition_path_final)
    lineages<-transition_path_final
  }
  print("Get transition path!")
  #get pseudo-time order
  {
    nTraj<-1
    Trajectory <- list()
    lanbda_temp<-matrix(0,nrow=nrow(X))
    #lineages_temp save paths between cluster centroid
    for(i1 in seq_len(length(lineages))){
      temp_router<-as.numeric(lineages[[i1]])
      vertix<-matrix(c(1:2),ncol=2)
      for(j in seq_len(length(temp_router))){
        h<-temp_router[j]
        vertix<-rbind(vertix,X_label[h,])
      }
      vertix<-vertix[-1,]
      projection_cluster<-as.numeric(lineages_temp[[i1]])
      vertix_1<-matrix(c(1:2),ncol=2)
      X_label_SR_index<-cbind(X_label_SR[,1:3],c(1:nrow(X_label_SR)))
      lambda_save_1<-c()
      for(j in seq_len(length(projection_cluster))){
        h<-projection_cluster[j]
        vertix_1<-rbind(vertix_1,X_label_SR[X_label_SR[,3]==h,1:2])
        lambda_save_1<-c(lambda_save_1,X_label_SR_index[X_label_SR[,3]==h,4])
      }
      vertix_1<-vertix_1[-1,]
      vertix_1<-as.matrix(vertix_1)
      if(nrow(vertix)<4){
        vertix<-as.matrix(vertix)
        s<-vertix
        fit1 <- project_to_curve(vertix_1, s)
      }else{
        fit <- principal_curve(vertix,start=vertix)
        s = fit$s[fit$ord, ,drop = FALSE]
        fit1 <- project_to_curve(vertix_1, s)
      }
      lambda_save_2<-as.matrix(fit1$lambda)
      lambda_save_2<-cbind(lambda_save_2,as.matrix(lambda_save_1))
      lambda_fit1<-matrix(0,nrow=nrow(X))
      for(i in seq_len(nrow(lambda_save_2))){
        lambda_fit1[lambda_save_2[i,2],1]<-lambda_save_2[i,1]
      }
      #lambda_fit1<-as.matrix(fit1$lambda)
      if(isNormalized==TRUE){
        max_lambda<-max(lambda_fit1)
        min_lambda<-min(lambda_fit1)
        lambda_fit1<-(lambda_fit1-min_lambda)/(max_lambda-min_lambda)
      }
      
      
      lanbda_temp<-cbind(lanbda_temp,lambda_fit1)
      
      Trajectory[[nTraj]] <- fit1
      names(Trajectory)[nTraj]<-paste0("trajectory",nTraj,sep='')
      nTraj<-nTraj+1
    }
    
    
    lanbda_temp<-as.matrix(lanbda_temp[,-1])
    
    #Determine the order of each point based on the projection index
    temp_order<-matrix(0,nrow=nrow(X))
    for(i11 in seq_len(nrow(temp_order))){
      temp_order[i11,1]<-max(lanbda_temp[i11,])
    }
    order<-temp_order
    Trajectory[nTraj]<-list(order)
    names(Trajectory)[nTraj]<-"pseudotime"
    nTraj<-nTraj+1
    
    Trajectory[nTraj]<-list(X)
    names(Trajectory)[nTraj]<-"X_two"
    nTraj<-nTraj+1
    
    Trajectory[nTraj]<-list(mclust_label)
    names(Trajectory)[nTraj]<-"mclust_label"
    nTraj<-nTraj+1
    
    Trajectory[nTraj]<-list(mclust_label)
    names(Trajectory)[nTraj]<-"mclust_label"
    nTraj<-nTraj+1
    
    Trajectory[nTraj]<-list(mclust_center)
    names(Trajectory)[nTraj]<-"mclust_center"
    nTraj<-nTraj+1
    
    Trajectory[nTraj]<-list(startCluster)
    names(Trajectory)[nTraj]<-"startCluster"
    nTraj<-nTraj+1
    
    
    Trajectory[nTraj]<-list(transfer_cell_coordinate)
    names(Trajectory)[nTraj]<-"transfer_cell_coordinate"
    nTraj<-nTraj+1
    
    Trajectory[nTraj]<-list(length(lineages))
    names(Trajectory)[nTraj]<-"number_trajectory"
    nTraj<-nTraj+1
    
    Trajectory[nTraj]<-list(seed1)
    names(Trajectory)[nTraj]<-"seed"
    nTraj<-nTraj+1
  }
  
  plot_tranfercell<-function(){
    cex=0.5
    plot(X[,1], X[,2],pch=16,cex=cex,col =mclust_label, yaxt="n",xaxt="n",xlab ="", ylab ="",)
    points(mclust_center,col = k+2, pch = 16, cex=cex+0.5)
    text(mclust_center,labels=1:nrow(mclust_center),pos="1",offset=-1,col="blue",cex=1.5)
    points(transfer_cell_coordinate[,1:2],col = "#8B4513", pch = 16, cex=cex)
    dims = seq_len(2)
    for(i1 in seq_len(length(lineages_temp))){
      for(j1 in seq(2,length(lineages_temp[[i1]]))){
        ii<-lineages_temp[[i1]][j1-1]
        h<-lineages_temp[[i1]][j1]
        lines(mclust_center[c(ii,h), dims], lwd = 4, col = "black")
      }
    }
  }
  plot_lineage<-function(){
    cex=0.5
    plot(X[,1], X[,2],cex=cex,col =mclust_label,bg=mclust_label,yaxt="n",xaxt="n", xlab ="", ylab ="")
    points(transfer_cell_coordinate[,1:2],pch=16,col = "#8B4513",cex=cex+0.5)
    text(mclust_center[startCluster,1],
         mclust_center[startCluster,2],
         labels="start",
         pos="1",
         offset=1,
         col="blue",cex=2)
    nTraj<-Trajectory[["number_trajectory"]]
    for(i1 in seq_len(nTraj)){
      lines(Trajectory[[i1]], lwd = 4, col = "black")
    }
  }
  plot_tranfercell()
  plot_lineage()
  return(Trajectory)
}

emd <- function(dist, w1, w2) {
  costs <- dist
  row.signs <- rep("<", length(w1))
  row.rhs <- w1
  col.signs <- rep(">", length(w2))
  col.rhs <- w2
  t <- lp.transport(costs, "min", row.signs, row.rhs, col.signs, col.rhs)
  flow <- t$solution
  work <- sum(flow * dist)
  e <- work / sum(flow)
  return(e)
}#EMD distance




